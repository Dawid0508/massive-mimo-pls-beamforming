% =========================================================================
% SCENARIO: Channel hardening in Massive MIMO and its impact on PLS
% -------------------------------------------------------------------------
% Zaktualizowano: Tłumienie 3GPP (compute_nr_pathloss), wspólny seed (CRN)
% oraz dynamiczne wyliczanie zysku Rx SNR wynikającego z rosnącego Nt.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
Nt_vec       = [4 8 16 32 64 128 256];   
Nt_hist      = [4, 256];                 
numIter      = 400;                      
fc           = p.fc_sub6;
dist_b       = 50;                       % Bob na 50 m
dist_e       = 40;                       % Ewa ukrywa się bliżej (40 m)
SNR_tx_dB    = 105;                       % Transmit SNR (rho_tx) [dB]
P_tx         = 10^(SNR_tx_dB / 10);
noise_var    = 1;
cdl_tag      = p.cdl_sub6;

% Zmiana z FSPL na model 3GPP
[PL_lin_b, PL_dB_b] = compute_nr_pathloss(dist_b, fc);
[PL_lin_e, PL_dB_e] = compute_nr_pathloss(dist_e, fc);

Rx_SNR_base_dB = 10 * log10(P_tx / PL_lin_b);

fprintf('\n--- Channel Hardening (nrCDLChannel) @ %.0f GHz ---\n', fc/1e9);
fprintf('  Transmit SNR (rho_tx): %d dB\n', SNR_tx_dB);
fprintf('  Bob (d = %g m, 3GPP PL = %.2f dB)\n', dist_b, PL_dB_b);
fprintf('  Eve (d = %g m, 3GPP PL = %.2f dB)\n', dist_e, PL_dB_e);

% Prezentacja zysku energetycznego Array Gain
fprintf('\n  [!] Rx SNR Scaling due to Array Gain (10*log10(Nt)):\n');
fprintf('      - Dla Nt = %3d: Rx SNR = %5.1f dB\n', Nt_vec(1), 10*log10(P_tx * Nt_vec(1) / PL_lin_b));
fprintf('      - Dla Nt = %3d: Rx SNR = %5.1f dB\n', Nt_vec(end), 10*log10(P_tx * Nt_vec(end) / PL_lin_b));
fprintf('-----------------------------------------------------------\n');

var_norm_h   = zeros(size(Nt_vec));      
mean_norm_h  = zeros(size(Nt_vec));
std_R_sec    = zeros(size(Nt_vec));      
mean_R_sec   = zeros(size(Nt_vec));
hist_data    = cell(length(Nt_hist), 1);

master_seeds = randi([0 2^31-1], numIter, 1); 

% --- Sweep ---------------------------------------------------------------
for n_idx = 1:length(Nt_vec)
    Nt = Nt_vec(n_idx);
    fprintf('Simulating Nt = %d...\n', Nt);
    
    g_samples = zeros(numIter, 1);
    R_samples = zeros(numIter, 1);
    
    cdl_b = setup_matlab_cdl_stat(Nt, fc, -30, cdl_tag);
    cdl_e = setup_matlab_cdl_stat(Nt, fc, 30, cdl_tag);
    
    for it = 1:numIter
        % Używamy stałego ziarna dla iteracji 'it' niezależnie od N_t
        common_seed = master_seeds(it); 
        
        release(cdl_b); cdl_b.Seed = common_seed;
        [pg_b, ~] = cdl_b();
        
        release(cdl_e); cdl_e.Seed = common_seed;
        [pg_e, ~] = cdl_e();
        
        h_raw_b = sum(pg_b, 2); h_raw_b = h_raw_b(:); 
        h_raw_e = sum(pg_e, 2); h_raw_e = h_raw_e(:);
        
        % Kanał efektywny z tłumieniem
        h_eff_b = sqrt(1/PL_lin_b) * h_raw_b;
        h_eff_e = sqrt(1/PL_lin_e) * h_raw_e;
        
        % Normalised gain do histogramu hardeningu
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

var_theory = 1 ./ Nt_vec;
plot_hardening_topology(dist_b, -30, dist_e, 30);

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
legend(arrayfun(@(n) sprintf('N_t = %d', n), Nt_hist, 'UniformOutput', false), 'Location', 'NorthEast');

% Top-right: Var(||h||^2/Nt) vs Nt with theory line
subplot(2, 2, 2);
loglog(Nt_vec, var_norm_h, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'CDL (Skorelowany)'); hold on;
loglog(Nt_vec, var_theory, '--', 'Color', c.ref, 'LineWidth', 1.5, 'DisplayName', 'Idealny i.i.d. (1/N_t)');
grid on; box on;
xlabel('Number of antennas N_t'); ylabel('Var(||h||^2 / N_t)');
title('Hardening rate: CDL vs Ideal i.i.d.');
legend('Location', 'NorthEast');

% Bottom-left: mean Secrecy Rate vs Nt
subplot(2, 2, 3);
semilogx(Nt_vec, mean_R_sec, '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
grid on; box on;
xlabel('Number of antennas N_t'); ylabel('E[R_{sec}] (bits/s/Hz)');
title(sprintf('Mean Secrecy Rate'));

% Bottom-right: std of Secrecy Rate vs Nt - operational hardening
subplot(2, 2, 4);
semilogx(Nt_vec, std_R_sec, '-mo', 'LineWidth', 2, 'MarkerFaceColor', 'm');
grid on; box on;
xlabel('Number of antennas N_t'); ylabel('Std(R_{sec}) (bits/s/Hz)');
title('Outage sensitivity collapses with N_t');

% Ujednolicony tytuł (Scenariusz 6)
sgtitle(sprintf('Scenario 6: Massive MIMO Channel Hardening (Base Rx SNR \\approx %.1f dB)', Rx_SNR_base_dB));
save_figure(fig, 'fig_channel_hardening');

% =========================================================================
% FUNKCJA POMOCNICZA: Konfiguracja kanału statycznego
% =========================================================================
function cdl = setup_matlab_cdl_stat(Nt, fc, theta, band_tag)
    cdl = nrCDLChannel;
    cdl.DelayProfile = 'CDL-A';
    
    if fc < 10e9, cdl.DelaySpread = 30e-9; else, cdl.DelaySpread = 10e-9; end
    
    cdl.CarrierFrequency = fc;
    cdl.MaximumDopplerShift = 0; 
    
    cdl.TransmitAntennaArray.Size = [1 Nt 1 1 1]; 
    cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1]; 
    cdl.TransmitArrayOrientation = [-theta; 0; 0];
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
    
    cdl.NumTimeSamples = 1;
    cdl.ChannelFiltering = false; 
end

% =========================================================================
% FUNKCJA POMOCNICZA: Generowanie topologii scenariusza (Channel Hardening)
% =========================================================================
function plot_hardening_topology(dist_b, theta_b, dist_e, theta_e)
    fig_top = figure('Color', 'w', 'Position', [150 150 700 700]);
    hold on; grid on; box on;
    
    x_bs = 0; y_bs = 0;
    max_d = max(dist_b, dist_e) + 15;
    
    plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    text(x_bs, y_bs - 4, 'BS (0,0)', 'HorizontalAlignment', 'center', 'Color', 'k');
    
    x_e = dist_e * sind(theta_e);
    y_e = dist_e * cosd(theta_e);
    plot(x_e, y_e, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Eve');
    text(x_e + 2, y_e, sprintf('Eve\n(%gm, %g\\circ)', dist_e, theta_e), 'Color', 'r', 'FontSize', 9);
    
    x_b = dist_b * sind(theta_b);
    y_b = dist_b * cosd(theta_b);
    plot(x_b, y_b, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Bob');
    text(x_b - 2, y_b - 2, sprintf('Bob\n(%gm, %g\\circ)', dist_b, theta_b), ...
         'Color', 'b', 'FontSize', 9, 'HorizontalAlignment', 'right');
    
    axis equal; xlim([-max_d, max_d]); ylim([-10, max_d]);
    xlabel('X [m]'); ylabel('Y [m]');
    title('Scenario 6: Channel Hardening');
    legend('Location', 'NorthWest');
    
    try save_figure(fig_top, '../topology/topology_channel_hardening'); catch; end
end