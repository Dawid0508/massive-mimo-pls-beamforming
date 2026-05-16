% =========================================================================
% SCENARIO: Vector vs Matrix ZF normalization - Sum-Rate vs Fairness
% Oparty na modelu 3GPP nrCDLChannel i fizycznym Path Loss (FSPL)
% ZAWARTY FIX: Dynamiczny i prawidłowy rozkład przestrzenny dla każdego K!
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- PARAMETRY SYSTEMU ---
Nt          = 32;                 
fc          = p.fc_sub6;          
K_fixed     = 8;                  
SNR_fixed   = 80;                 % Dostosowano do fizycznego FSPL (bezwzględnego)
SNR_tx_vec  = 60:5:120;           % Oś X: Transmit SNR (Moc stacji bazowej w dB)
K_vec       = 2:2:16;             
numIter     = 50;                 

dist_e  = 40;
theta_e = 30;

% Kanał podsłuchiwacza Ewy (jedna, stała pozycja)
cdl_e = setup_matlab_cdl(Nt, fc, theta_e, 'CDL-A');

SR_vs_SNR = zeros(2, length(SNR_tx_vec));
J_vs_SNR  = zeros(2, length(SNR_tx_vec));
SR_vs_K   = zeros(2, length(K_vec));
J_vs_K    = zeros(2, length(K_vec));

% =========================================================================
% --- 1. SWEEP TRANSMIT SNR (Dla stałego K = 8) ---
% =========================================================================
fprintf('Rozpoczynam sweep Transmit SNR dla K = %d...\n', K_fixed);

% Generujemy dedykowaną geometrię, która idealnie pokrywa CAŁY sektor
dist_snr_sweep  = linspace(20, 150, K_fixed);
theta_snr_sweep = linspace(-60, 60, K_fixed);

cdl_b_snr = cell(K_fixed, 1);
for k = 1:K_fixed
    cdl_b_snr{k} = setup_matlab_cdl(Nt, fc, theta_snr_sweep(k), 'CDL-A');
end

for s_idx = 1:length(SNR_tx_vec)
    P_tx_lin = 10^(SNR_tx_vec(s_idx) / 10);
    [SR_vs_SNR(:, s_idx), J_vs_SNR(:, s_idx)] = run_cdl_sweep( ...
        cdl_b_snr, cdl_e, dist_snr_sweep, dist_e, fc, Nt, K_fixed, P_tx_lin, numIter);
end

% =========================================================================
% --- 2. SWEEP K (Dla stałego Transmit SNR) ---
% =========================================================================
P_tx_fixed_lin = 10^(SNR_fixed / 10);
fprintf('Rozpoczynam sweep K dla Transmit SNR = %d dB...\n', SNR_fixed);

for k_idx = 1:length(K_vec)
    K_current = K_vec(k_idx);
    
    % KLUCZOWA POPRAWKA: Dynamiczne skalowanie siatki geometrycznej dla danego K!
    dist_k_sweep  = linspace(20, 150, K_current);
    theta_k_sweep = linspace(-60, 60, K_current);
    
    cdl_b_k = cell(K_current, 1);
    for k = 1:K_current
        cdl_b_k{k} = setup_matlab_cdl(Nt, fc, theta_k_sweep(k), 'CDL-A');
    end
    
    [SR_vs_K(:, k_idx), J_vs_K(:, k_idx)] = run_cdl_sweep( ...
        cdl_b_k, cdl_e, dist_k_sweep, dist_e, fc, Nt, K_current, P_tx_fixed_lin, numIter);
end

% =========================================================================
% --- WIZUALIZACJA ---
% =========================================================================
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

subplot(2, 2, 1);
plot(SNR_tx_vec, SR_vs_SNR(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(SNR_tx_vec, SR_vs_SNR(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Transmit SNR at Base Station (dB)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Sum-Rate vs SNR (K = %d)', K_fixed));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'NorthWest');

subplot(2, 2, 2);
plot(SNR_tx_vec, J_vs_SNR(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(SNR_tx_vec, J_vs_SNR(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Transmit SNR at Base Station (dB)'); ylabel("Jain's index (Fairness)");
title(sprintf('Fairness vs SNR (K = %d)', K_fixed));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthWest');

subplot(2, 2, 3);
plot(K_vec, SR_vs_K(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(K_vec, SR_vs_K(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Number of users (K)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Sum-Rate vs K (Transmit SNR = %d dB)', SNR_fixed));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthWest');

subplot(2, 2, 4);
plot(K_vec, J_vs_K(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(K_vec, J_vs_K(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Number of users (K)'); ylabel("Jain's index (Fairness)");
title(sprintf('Fairness vs K (Transmit SNR = %d dB)', SNR_fixed));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthWest');

sgtitle(sprintf('ZF Normalization Trade-off: 6 GHz CDL-A (Nt = %d, Full Sector Spread)', Nt));
save_figure(fig, 'fig_fairness_normalization_cdl');

% =========================================================================
% FUNKCJA WYKONAWCZA SWEEPÓW
% =========================================================================
function [SR_out, J_out] = run_cdl_sweep(cdl_b, cdl_e, dist_b, dist_e, fc, Nt, K, P_tx, numIter)
    SR_acc = zeros(2, 1);
    J_acc  = zeros(2, 1);
    
    for it = 1:numIter
        H_eff = zeros(Nt, K);
        
        % Kanały użytkowników (Bobowie)
        for k = 1:K
            release(cdl_b{k});
            cdl_b{k}.Seed = randi([0 2^31-1]);
            [pg_b, ~] = cdl_b{k}();
            h_b = squeeze(sum(pg_b, 2)); h_b = h_b(:);
            
            % Wyliczanie fizycznego tłumienia przestrzennego dla każdego Boba
            [PL_lin_b, ~] = compute_fspl(dist_b(k), fc);
            
            % Kanał efektywny z twardą fizyką tłumienia amplitudy sygnału
            H_eff(:, k) = sqrt((1 / PL_lin_b) * Nt) * (h_b / norm(h_b));
        end
        
        % Kanał podsłuchiwacza (Ewa)
        release(cdl_e);
        cdl_e.Seed = randi([0 2^31-1]);
        [pg_e, ~] = cdl_e();
        h_e = squeeze(sum(pg_e, 2)); h_e = h_e(:);
        
        % Wyliczanie fizycznego tłumienia przestrzennego dla Ewy
        [PL_lin_e, ~] = compute_fspl(dist_e, fc);
        h_e_eff = sqrt((1 / PL_lin_e) * Nt) * (h_e / norm(h_e));
        
        % Obliczanie surowego prekodera pseudoodwrotności (ZF)
        W_raw = H_eff * pinv(H_eff' * H_eff);
        
        % 1. Normalizacja Macierzowa (Frobeniusa)
        W_mat = W_raw / norm(W_raw, 'fro');
        
        % 2. Normalizacja Wektorowa (per-user)
        W_vec = zeros(Nt, K);
        for k = 1:K
            if norm(W_raw(:, k)) > 1e-9
                W_vec(:, k) = W_raw(:, k) / norm(W_raw(:, k)) * sqrt(1/K);
            end
        end
        
        % Ewaluacja metryk PLS
        [sr_m, j_m] = compute_pls_metrics(H_eff, h_e_eff, W_mat, P_tx, K);
        [sr_v, j_v] = compute_pls_metrics(H_eff, h_e_eff, W_vec, P_tx, K);
        
        SR_acc(1) = SR_acc(1) + sr_m;
        SR_acc(2) = SR_acc(2) + sr_v;
        J_acc(1)  = J_acc(1)  + j_m;
        J_acc(2)  = J_acc(2)  + j_v;
    end
    
    SR_out = SR_acc / numIter;
    J_out  = J_acc  / numIter;
end

% =========================================================================
% FUNKCJA POMOCNICZA: Metryki Sum-Rate oraz Jain's Fairness Index
% =========================================================================
function [sum_secrecy, jains_index] = compute_pls_metrics(H_eff, h_e_eff, W, P_tx, K)
    R_b = zeros(1, K);
    R_e = zeros(1, K);
    
    for k = 1:K
        S_b = P_tx * abs(H_eff(:, k)' * W(:, k))^2;
        I_b = 0;
        for j = 1:K
            if j ~= k
                I_b = I_b + P_tx * abs(H_eff(:, k)' * W(:, j))^2;
            end
        end
        R_b(k) = log2(1 + S_b / (I_b + 1)); % Szum tła znormalizowany do 1
        
        S_e = P_tx * abs(h_e_eff' * W(:, k))^2;
        I_e = 0;
        for j = 1:K
            if j ~= k
                I_e = I_e + P_tx * abs(h_e_eff' * W(:, j))^2;
            end
        end
        R_e(k) = log2(1 + S_e / (I_e + 1));
    end
    
    sum_secrecy = sum(max(0, R_b - R_e));
    
    if sum(R_b.^2) == 0
        jains_index = 0;
    else
        jains_index = (sum(R_b))^2 / (K * sum(R_b.^2));
    end
end

% =========================================================================
% FUNKCJA POMOCNICZA: Konfiguracja kanału nrCDLChannel
% =========================================================================
function cdl = setup_matlab_cdl(Nt, fc, theta, band_tag)
    cdl = nrCDLChannel;
    cdl.DelayProfile = 'CDL-A';      
    cdl.DelaySpread = 30e-9; 
    cdl.CarrierFrequency = fc;
    cdl.MaximumDopplerShift = 0;         
    cdl.TransmitAntennaArray.Size = [1 Nt 1 1 1]; 
    cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1]; 
    cdl.TransmitArrayOrientation = [-theta; 0; 0];
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
    cdl.NumTimeSamples = 1;
    cdl.ChannelFiltering = false; 
end