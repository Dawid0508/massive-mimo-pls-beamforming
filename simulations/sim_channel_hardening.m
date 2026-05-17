% =========================================================================
% SCENARIO: Channel hardening in Massive MIMO and its impact on PLS
% (Zaktualizowano do oficjalnego modelu nrCDLChannel z 5G Toolbox)
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
Nt_vec       = [4 8 16 32 64 128 256];   % Zmniejszono max do 256 dla optymalizacji
Nt_hist      = [4, 256];                 % values shown as histograms

% UWAGA: Zmniejszamy numIter z 4000 na 400. Model nrCDLChannel z 5G Toolbox 
% jest znacznie cięższy obliczeniowo niż zwykłe randn().
numIter      = 400;                      

fc           = p.fc_sub6;
dist_b       = 50;                       % Bob na 50 m
dist_e       = 40;                       % Ewa ukrywa się bliżej (40 m)
SNR_tx_dB    = 90;                       % Moc stacji bazowej [dB]
P_tx         = 10^(SNR_tx_dB / 10);
noise_var    = 1;
cdl_tag      = p.cdl_sub6;

fprintf('\n--- Channel Hardening (nrCDLChannel) @ %.0f GHz ---\n', fc/1e9);
fprintf('  Transmit SNR: %d dB\n', SNR_tx_dB);
fprintf('  Bob (d = %g m), Eve (d = %g m)\n', dist_b, dist_e);

[PL_lin_b, PL_dB_b] = compute_fspl(dist_b, fc);
[PL_lin_e, PL_dB_e] = compute_fspl(dist_e, fc);

var_norm_h   = zeros(size(Nt_vec));      % Var(||h||^2 / Nt)
mean_norm_h  = zeros(size(Nt_vec));
std_R_sec    = zeros(size(Nt_vec));      % Std(R_secrecy) under MRT
mean_R_sec   = zeros(size(Nt_vec));
hist_data    = cell(length(Nt_hist), 1);

% --- Sweep ---------------------------------------------------------------
for n_idx = 1:length(Nt_vec)
    Nt = Nt_vec(n_idx);
    fprintf('Simulating Nt = %d...\n', Nt);
    
    g_samples = zeros(numIter, 1);
    R_samples = zeros(numIter, 1);
    
    % Inicjalizujemy kanały dla danej liczby anten
    cdl_b = setup_matlab_cdl_stat(Nt, fc, -30, cdl_tag);
    cdl_e = setup_matlab_cdl_stat(Nt, fc, 30, cdl_tag);
    
    for it = 1:numIter
        % Generujemy nowe realizacje kanału
        release(cdl_b); cdl_b.Seed = randi([0 2^31-1]);
        [pg_b, ~] = cdl_b();
        
        release(cdl_e); cdl_e.Seed = randi([0 2^31-1]);
        [pg_e, ~] = cdl_e();
        
        % Bezpieczne spłaszczanie wyjścia do wektora kolumnowego
        h_raw_b = sum(pg_b, 2); h_raw_b = h_raw_b(:); 
        h_raw_e = sum(pg_e, 2); h_raw_e = h_raw_e(:);
        
        % --- KLUCZOWA ZMIANA: Brak norm(h_raw) w mianowniku ---
        % W modelu CDL, moc ścieżek sumuje się nominalnie do 1. 
        % Aplikujemy tylko fizyczne tłumienie.
        h_eff_b = sqrt(1/PL_lin_b) * h_raw_b;
        h_eff_e = sqrt(1/PL_lin_e) * h_raw_e;
        
        % Channel hardening opiera się na analizie surowego wzmocnienia małoskalowego.
        g_samples(it) = (h_raw_b' * h_raw_b) / Nt;          
        
        % Prekoder MRT (stacja kieruje wiązkę na Boba)
        w = h_eff_b / norm(h_eff_b);
        
        % Przepustowości z poprawną fizyką Mocy Transmisyjnej
        R_b = log2(1 + P_tx * abs(h_eff_b' * w)^2 / noise_var);
        R_e = log2(1 + P_tx * abs(h_eff_e' * w)^2 / noise_var);
        
        R_samples(it) = max(0, R_b - R_e);
    end
    
    mean_norm_h(n_idx) = mean(g_samples);
    var_norm_h(n_idx)  = var(g_samples);
    mean_R_sec(n_idx)  = mean(R_samples);
    std_R_sec(n_idx)   = std(R_samples);
    
    idx = find(Nt_hist == Nt, 1);
    if ~isempty(idx)
        hist_data{idx} = g_samples;
    end
end

% Theoretical reference: Var = 1/Nt dla czystego kanału i.i.d. Rayleigha
var_theory = 1 ./ Nt_vec;

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

% Top-left: hardening histograms
subplot(2, 2, 1);
colors = lines(length(Nt_hist));
edges  = 0:0.05:3.5;
for i = 1:length(Nt_hist)
    histogram(hist_data{i}, edges, 'Normalization', 'pdf', ...
        'FaceColor', colors(i,:), 'FaceAlpha', 0.55, ...
        'DisplayName', sprintf('N_t = %d', Nt_hist(i))); hold on;
end
c = pls_colors();
xline(1, '--', 'E[||h||^2/N_t] = 1', 'Color', c.ref, 'LineWidth', 1.2);
pls_axis_prefs(gca, 'refLabelV', 'top');
grid on; box on;
xlabel('||h||^2 / N_t'); ylabel('PDF');
title('Distribution of normalised gain (3GPP CDL)');
legend(arrayfun(@(n) sprintf('N_t = %d', n), Nt_hist, 'UniformOutput', false), ...
    'Location', 'NorthEast');

% Top-right: Var(||h||^2/Nt) vs Nt with theory line
subplot(2, 2, 2);
loglog(Nt_vec, var_norm_h, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b', ...
    'DisplayName', 'CDL (Skorelowany)'); hold on;
loglog(Nt_vec, var_theory, '--', 'Color', c.ref, 'LineWidth', 1.5, ...
    'DisplayName', 'Idealny i.i.d. (1/N_t)');
grid on; box on;
xlabel('Number of antennas N_t'); ylabel('Var(||h||^2 / N_t)');
title('Hardening rate: CDL vs Ideal i.i.d.');
legend('Location', 'NorthEast');

% Bottom-left: mean Secrecy Rate vs Nt
subplot(2, 2, 3);
semilogx(Nt_vec, mean_R_sec, '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
grid on; box on;
xlabel('Number of antennas N_t'); ylabel('E[R_{sec}] (bits/s/Hz)');
title(sprintf('Mean Secrecy Rate (Transmit SNR = %d dB)', SNR_tx_dB));

% Bottom-right: std of Secrecy Rate vs Nt - operational hardening
subplot(2, 2, 4);
semilogx(Nt_vec, std_R_sec, '-mo', 'LineWidth', 2, 'MarkerFaceColor', 'm');
grid on; box on;
xlabel('Number of antennas N_t'); ylabel('Std(R_{sec}) (bits/s/Hz)');
title('Outage sensitivity collapses with N_t');

sgtitle('Massive MIMO Channel Hardening (nrCDLChannel / P_tx Physics)');
save_figure(fig, 'fig_channel_hardening_cdl');

% =========================================================================
% FUNKCJA POMOCNICZA: Konfiguracja kanału statycznego
% =========================================================================
function cdl = setup_matlab_cdl_stat(Nt, fc, theta, band_tag)
    cdl = nrCDLChannel;
    cdl.DelayProfile = 'CDL-A';
    
    if fc < 10e9  
        cdl.DelaySpread = 30e-9;  
    else          
        cdl.DelaySpread = 10e-9;  
    end
    
    cdl.CarrierFrequency = fc;
    cdl.MaximumDopplerShift = 0; % Kanał statyczny do analizy pojemności
    
    cdl.TransmitAntennaArray.Size = [1 Nt 1 1 1]; 
    cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1]; 
    cdl.TransmitArrayOrientation = [-theta; 0; 0];
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
    
    cdl.NumTimeSamples = 1;
    cdl.ChannelFiltering = false; 
end