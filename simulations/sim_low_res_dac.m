% =========================================================================
% SCENARIO: Low-Resolution DAC quantisation in Massive MIMO PLS
% -------------------------------------------------------------------------
% Zaktualizowano: Zastosowano "Common Random Numbers" (Wspólne zmienne
% losowe). Pętla Monte Carlo (kanały) jest na zewnątrz, a pętla sprzętowa 
% (bity DAC) wewnątrz. Gwarantuje to idealnie uczciwe porównanie algorytmów.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- KONFIGURACJA SYSTEMU I FIZYKI ---------------------------------------
b_vec       = [1 2 3 4 5 Inf];               % Rozdzielczość DAC [bity]
b_show      = [1, 2, 4, Inf];                % Wybrane bity dla Sweepu B
SNR_tx_vec  = 80:10:140;                     % Transmit SNR sweep [dB]
SNR_fixed   = 90;                           % Stały SNR_tx dla Sweepu A
Nt          = 32;                            % Liczba anten (Massive MIMO)
K           = 4;                             % Liczba użytkowników
numIter     = 40;                            % Iteracje kanału
numSym      = 1000;                          % ZWIĘKSZONO! Lepsza statystyka
noise_var   = 1;                             % Znormalizowany szum tła

fc          = p.fc_sub6;                     % Pasmo 6 GHz
dist_b      = 50;                            % Dystans do Bobów [m]
dist_e      = 40;                            % Dystans do Ewy [m]
theta_e     = 25;                            % Kąt podsłuchiwacza

[PL_lin_b, PL_dB_b] = compute_fspl(dist_b, fc);
[PL_lin_e, PL_dB_e] = compute_fspl(dist_e, fc);

fprintf('\n--- Low-Res DACs @ 6 GHz (Nt=%d) ---\n', Nt);
fprintf('  Dystans Bob: %g m (FSPL: %.2f dB)\n', dist_b, PL_dB_b);
fprintf('  Dystans Eve: %g m (FSPL: %.2f dB)\n', dist_e, PL_dB_e);

% Przygotowanie obiektów kanałowych
theta_bobs = linspace(-60, 60, K);
cdl_b = cell(K, 1);
for k = 1:K, cdl_b{k} = setup_matlab_cdl(nrCDLChannel, Nt, fc, theta_bobs(k)); end
cdl_e = setup_matlab_cdl(nrCDLChannel, Nt, fc, theta_e);

% Akumulatory na wyniki (wypełniane wewnątrz pętli)
R_bob_acc = zeros(1, length(b_vec));
R_eve_acc = zeros(1, length(b_vec));
R_sec_acc = zeros(1, length(b_vec));
a_acc_vec = zeros(1, length(b_vec));
a_n_vec   = zeros(1, length(b_vec));

R_sec_snr_acc = zeros(length(b_show), length(SNR_tx_vec));

% =========================================================================
% --- SWEEP A: Wpływ rozdzielczości DAC (b) przy stałym Transmit SNR
% =========================================================================
fprintf('Rozpoczynam Sweep A (Rozdzielczość DAC)...\n');
P_tx = 10^(SNR_fixed / 10);

% GŁÓWNA PĘTLA ZEWNĘTRZNA (Losowanie "świata" raz na iterację)
for it = 1:numIter
    
    % 1. ŚRODOWISKO RADIOWE: Losujemy kanały tylko RAZ
    H_eff = zeros(Nt, K);
    for k = 1:K
        release(cdl_b{k}); cdl_b{k}.Seed = randi([0 2^31-1]);
        [pg_b, ~] = cdl_b{k}(); hb = squeeze(sum(pg_b, 2)); hb = hb(:);
        H_eff(:, k) = sqrt((1 / PL_lin_b) * Nt) * (hb / norm(hb));
    end
    release(cdl_e); cdl_e.Seed = randi([0 2^31-1]);
    [pg_e, ~] = cdl_e(); he = squeeze(sum(pg_e, 2)); he = he(:);
    h_eff_e = sqrt((1 / PL_lin_e) * Nt) * (he / norm(he));
    
    % 2. MATEMATYKA STACJI: Idealny Prekoder
    W_raw = H_eff * pinv(H_eff' * H_eff);
    W = W_raw / norm(W_raw, 'fro') * sqrt(P_tx);
    
    % 3. SYGNAŁY I SZUM: Identyczne dla każdego DAC!
    s = (randn(K, numSym) + 1j*randn(K, numSym)) / sqrt(2);
    x = W * s;
    n_th_b = sqrt(noise_var/2) * (randn(K, numSym) + 1j*randn(K, numSym));
    n_th_e = sqrt(noise_var/2) * (randn(1, numSym) + 1j*randn(1, numSym));
    
    % PĘTLA WEWNĘTRZNA: Testujemy sprzęt w tych samych warunkach
    for bi = 1:length(b_vec)
        b = b_vec(bi);
        
        % Kwantyzacja (Tylko ten element różnicuje nam symulację)
        x_q = local_uniform_quantize(x, b);       
        
        a_acc_vec(bi) = a_acc_vec(bi) + real(x(:)' * x_q(:)) / (x(:)' * x(:) + eps);
        a_n_vec(bi)   = a_n_vec(bi) + 1;
        
        % Propagacja w powietrzu (Używamy wspólnego szumu termicznego!)
        y_b = H_eff' * x_q + n_th_b;
        y_e = h_eff_e' * x_q + n_th_e;
        
        sum_R_b = 0; sum_R_e = 0; sum_R_s = 0;
        for k = 1:K
            useful_b = mean( real( y_b(k, :) .* conj(s(k, :)) ) );
            tot_b    = mean( abs(y_b(k, :)).^2 );
            sig_b    = useful_b^2;
            inr_b    = max(tot_b - sig_b, eps);
            R_b_k    = log2(1 + sig_b / inr_b);
            
            useful_e = mean( real( y_e .* conj(s(k, :)) ) );
            tot_e    = mean( abs(y_e).^2 );
            sig_e    = useful_e^2;
            inr_e    = max(tot_e - sig_e, eps);
            R_e_k    = log2(1 + sig_e / inr_e);
            
            sum_R_b = sum_R_b + R_b_k;
            sum_R_e = sum_R_e + R_e_k;
            sum_R_s = sum_R_s + max(0, R_b_k - R_e_k);
        end
        % Zrzut do akumulatorów globalnych
        R_bob_acc(bi) = R_bob_acc(bi) + sum_R_b;
        R_eve_acc(bi) = R_eve_acc(bi) + sum_R_e;
        R_sec_acc(bi) = R_sec_acc(bi) + sum_R_s;
    end
end
% Uśrednianie wyników Sweep A
R_bob = R_bob_acc / numIter;
R_eve = R_eve_acc / numIter;
R_sec = R_sec_acc / numIter;
a_emp = a_acc_vec ./ max(a_n_vec, 1);

% =========================================================================
% --- SWEEP B: Secrecy Rate vs Transmit SNR (Również wspólne kanały!)
% =========================================================================
fprintf('Rozpoczynam Sweep B (Transmit SNR)...\n');
for s_idx = 1:length(SNR_tx_vec)
    P_tx_s = 10^(SNR_tx_vec(s_idx)/10);
    
    for it = 1:numIter
        % 1. Wylosowanie środowiska
        H_eff = zeros(Nt, K);
        for k = 1:K
            release(cdl_b{k}); cdl_b{k}.Seed = randi([0 2^31-1]);
            [pg_b, ~] = cdl_b{k}(); hb = squeeze(sum(pg_b, 2)); hb = hb(:);
            H_eff(:, k) = sqrt((1 / PL_lin_b) * Nt) * (hb / norm(hb));
        end
        release(cdl_e); cdl_e.Seed = randi([0 2^31-1]);
        [pg_e, ~] = cdl_e(); he = squeeze(sum(pg_e, 2)); he = he(:);
        h_eff_e = sqrt((1 / PL_lin_e) * Nt) * (he / norm(he));
        
        W_raw = H_eff * pinv(H_eff' * H_eff);
        W     = W_raw / norm(W_raw, 'fro') * sqrt(P_tx_s);
        
        s = (randn(K, numSym) + 1j*randn(K, numSym)) / sqrt(2);
        x = W * s;
        n_th_b = sqrt(noise_var/2) * (randn(K, numSym) + 1j*randn(K, numSym));
        n_th_e = sqrt(noise_var/2) * (randn(1, numSym) + 1j*randn(1, numSym));
        
        % Pętla po rozdzielczościach DAC na wspólnym środowisku
        for ii = 1:length(b_show)
            b = b_show(ii);
            x_q = local_uniform_quantize(x, b);
            
            y_b = H_eff' * x_q + n_th_b;
            y_e = h_eff_e' * x_q + n_th_e;
            
            sum_R_s = 0;
            for k = 1:K
                useful_b = mean(real(y_b(k,:) .* conj(s(k,:))));
                tot_b    = mean(abs(y_b(k,:)).^2);
                R_b_k    = log2(1 + useful_b^2 / max(tot_b - useful_b^2, eps));
                
                useful_e = mean(real(y_e .* conj(s(k,:))));
                tot_e    = mean(abs(y_e).^2);
                R_e_k    = log2(1 + useful_e^2 / max(tot_e - useful_e^2, eps));
                
                sum_R_s  = sum_R_s + max(0, R_b_k - R_e_k);
            end
            R_sec_snr_acc(ii, s_idx) = R_sec_snr_acc(ii, s_idx) + sum_R_s;
        end
    end
end
R_sec_snr = R_sec_snr_acc / numIter;

plot_dac_topology(dist_b, theta_bobs, dist_e, theta_e);

% =========================================================================
% --- WIZUALIZACJA --------------------------------------------------------
% =========================================================================
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);
bar_x = 1:length(b_vec);
bar_labels = arrayfun(@(b) ternary(isinf(b), 'Inf (Ideal)', sprintf('%d bit', b)), b_vec, 'UniformOutput', false);

% Górny Lewy: Pojemność Boba i Ewy vs Rozdzielczość
subplot(2, 2, 1);
bar(bar_x, [R_bob(:) R_eve(:)], 1.0); grid on; box on;
set(gca, 'XTick', bar_x, 'XTickLabel', bar_labels);
xlabel('DAC Resolution (Bits per branch)');
ylabel('Sum-Rate (bits/s/Hz)');
title(sprintf('Bob vs Eve Capacity (Tx SNR = %d dB)', SNR_fixed));
legend('Bob (Target)', 'Eve (Eavesdropper)', 'Location', 'NorthWest');

% Górny Prawy: Secrecy Rate vs Rozdzielczość
subplot(2, 2, 2);
bar(bar_x, R_sec, 0.6, 'FaceColor', [0.2 0.6 0.3]);
grid on; box on;
set(gca, 'XTick', bar_x, 'XTickLabel', bar_labels);
xlabel('DAC Resolution (Bits per branch)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Overall Network Security (Secrecy Sum-Rate)');

% Dolny Lewy: Empiryczny czynnik Bussganga
subplot(2, 2, 3);
b_finite = b_vec(~isinf(b_vec));
a_finite = a_emp(~isinf(b_vec));
bar(1:length(b_finite), a_finite, 0.6, 'FaceColor', [0.6 0.4 0.8]); hold on;
c = pls_colors();
yline(1, '--', 'Ideal Hardware (a = 1.0)', 'Color', c.perfect, 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
grid on; box on;
set(gca, 'XTick', 1:length(b_finite), 'XTickLabel', arrayfun(@num2str, b_finite, 'Uniform', false));
xlabel('DAC bits (b)'); ylabel('Bussgang Gain (a)');
title('Signal Survival Rate (Bussgang Theorem)');
ylim([0 1.1]);

% Dolny Prawy: Sweep SNR
subplot(2, 2, 4);
markers = {'-rs', '-ms', '-go', '-bo'};
hold on;
for ii = 1:length(b_show)
    plot(SNR_tx_vec, R_sec_snr(ii,:), markers{ii}, 'LineWidth', 2, 'MarkerFaceColor', markers{ii}(2));
end
grid on; box on;
xlabel('Transmit SNR at Base Station (dB)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Secrecy Capacity vs Transmit Power');
labels = arrayfun(@(b) ternary(isinf(b), 'Ideal DACs (Inf)', sprintf('%d-bit DACs', b)), b_show, 'UniformOutput', false);
legend(labels, 'Location', 'NorthWest');

sgtitle(sprintf('Massive MIMO Hardware Limits: Low-Resolution DACs (Nt = %d, K = %d)', Nt, K));
save_figure(fig, 'fig_low_res_dac');

% =========================================================================
% FUNKCJE POMOCNICZE
% =========================================================================
function xq = local_uniform_quantize(x, b)
    if isinf(b)
        xq = x; return;
    end
    sigma = sqrt(mean(abs(x(:)).^2));
    clip_val = 3 * sigma / sqrt(2); 
    step = 2 * clip_val / (2^b);

    xr = real(x); xi = imag(x);
    xr(xr > clip_val) = clip_val; xr(xr < -clip_val) = -clip_val;
    xi(xi > clip_val) = clip_val; xi(xi < -clip_val) = -clip_val;
    
    xqr = floor(xr / step) * step + step/2;
    xqi = floor(xi / step) * step + step/2;
    xq = xqr + 1j * xqi;
    
    gain = sqrt(mean(abs(x(:)).^2) / mean(abs(xq(:)).^2));
    xq = xq * gain;
end

function cdl = setup_matlab_cdl(cdl, Nt, fc, theta)
    release(cdl);
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

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

% =========================================================================
% FUNKCJA POMOCNICZA: Generowanie topologii scenariusza (DAC Quantisation)
% =========================================================================
function plot_dac_topology(dist_b, theta_bobs, dist_e, theta_e)
    fig_top = figure('Color', 'w', 'Position', [150 150 700 700]);
    hold on; grid on; box on;
    
    % Konwersja na współrzędne kartezjańskie (BS w 0,0)
    x_bs = 0; y_bs = 0;
    max_d = max(dist_b, dist_e) + 15;
    
    % Rysowanie BS
    p_bs = plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    text(x_bs, y_bs - 4, 'BS (0,0)', 'HorizontalAlignment', 'center', 'Color', 'k');
    
    % Rysowanie Ewy
    x_e = dist_e * sind(theta_e);
    y_e = dist_e * cosd(theta_e);
    p_e = plot(x_e, y_e, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Eve');
    text(x_e + 2, y_e, sprintf('Eve\n(%gm, %g\\circ)', dist_e, theta_e), 'Color', 'r', 'FontSize', 9);
    
    % Rysowanie Bobów
    p_b = [];
    for k = 1:length(theta_bobs)
        x_b = dist_b * sind(theta_bobs(k));
        y_b = dist_b * cosd(theta_bobs(k));
        
        p_b = plot(x_b, y_b, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        text(x_b - 2.5, y_b - 2.5, sprintf('B_{%d}\n(%gm, %g\\circ)', k, dist_b, theta_bobs(k)), ...
             'Color', 'b', 'FontSize', 8, 'HorizontalAlignment', 'right');
    end
    if ~isempty(p_b)
        set(p_b, 'DisplayName', sprintf('Bobs (K=%d)', length(theta_bobs)));
    end
    
    % Ustawienia osi
    axis equal;
    xlim([-max_d, max_d]);
    ylim([-10, max_d]);
    xlabel('X [m]'); ylabel('Y [m]');
    title(sprintf('Scenario 5: Hardware Limits (K=%d Bobs, 1 Eve)', length(theta_bobs)));
    legend([p_bs, p_b, p_e], 'Location', 'NorthWest');
    
    % Zapis do pliku
    try
        save_figure(fig_top, '../topology/topology_low_res_dac');
    catch
        warning('Funkcja save_figure nie jest dostępna.');
    end
end