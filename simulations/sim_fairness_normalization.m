% =========================================================================
% SCENARIO: Vector vs Matrix ZF normalization - Sum-Rate vs Fairness
% Oparty na modelu 3GPP nrCDLChannel i fizycznym Path Loss (3GPP UMi)
% ZAWARTY FIX: Dynamiczny rozkład przestrzenny oraz CRN (Common Random Numbers)!
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- PARAMETRY SYSTEMU ---
Nt          = 32;                 
fc          = p.fc_sub6;          
K_fixed     = 8;                  
SNR_fixed   = 110;                 % Bazowy Transmit SNR do fizycznego path loss
SNR_tx_vec  = 80:5:140;           % Transmit SNR (Moc stacji bazowej w dB)
K_vec       = 2:2:16;             
numIter     = 100;                 
dist_e  = 40;
theta_e = 30;

% Kanał podsłuchiwacza Ewy (jedna, stała pozycja)
cdl_e = setup_matlab_cdl(Nt, fc, theta_e, 'CDL-A');

SR_vs_SNR  = zeros(2, length(SNR_tx_vec));
J_vs_SNR   = zeros(2, length(SNR_tx_vec));
SNR_rx_vec = zeros(1, length(SNR_tx_vec)); % Wektor do zbierania uśrednionego Rx SNR
SR_vs_K    = zeros(2, length(K_vec));
J_vs_K     = zeros(2, length(K_vec));
Rx_SNR_K_vec = zeros(1, length(K_vec));    % Średni Rx SNR przy zmianie K

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

% [CRN]: Wektor stałych seedów dla Sweepu SNR
seeds_snr = randi([0 2^31-1], numIter, 1);

for s_idx = 1:length(SNR_tx_vec)
    P_tx_lin = 10^(SNR_tx_vec(s_idx) / 10);
    
    % Przekazujemy seeds_snr jako dodatkowy argument do funkcji sweepu
    [SR_vs_SNR(:, s_idx), J_vs_SNR(:, s_idx), SNR_rx_vec(s_idx)] = run_cdl_sweep( ...
        cdl_b_snr, cdl_e, dist_snr_sweep, dist_e, fc, Nt, K_fixed, P_tx_lin, numIter, seeds_snr);
end

% =========================================================================
% --- 2. SWEEP K (Dla stałego Transmit SNR) ---
% =========================================================================
P_tx_fixed_lin = 10^(SNR_fixed / 10);
fprintf('Rozpoczynam sweep K dla Transmit SNR = %d dB...\n', SNR_fixed);

% Przywracamy bazowy stan generatora przed losowaniem seedów dla Sweepu K (powtarzalność)
rng(p.rng_seed); 
seeds_k = randi([0 2^31-1], numIter, 1);

for k_idx = 1:length(K_vec)
    K_current = K_vec(k_idx);
    
    % Dynamiczne skalowanie siatki geometrycznej dla danego K
    dist_k_sweep  = linspace(20, 150, K_current);
    theta_k_sweep = linspace(-60, 60, K_current);
    
    cdl_b_k = cell(K_current, 1);
    for k = 1:K_current
        cdl_b_k{k} = setup_matlab_cdl(Nt, fc, theta_k_sweep(k), 'CDL-A');
    end
    
    % Przekazujemy seeds_k jako dodatkowy argument do funkcji sweepu
    [SR_vs_K(:, k_idx), J_vs_K(:, k_idx), Rx_SNR_K_vec(k_idx)] = run_cdl_sweep( ...
        cdl_b_k, cdl_e, dist_k_sweep, dist_e, fc, Nt, K_current, P_tx_fixed_lin, numIter, seeds_k);
end

% =========================================================================
% --- WIZUALIZACJA ---
% =========================================================================
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);
% 1. Sum-Rate vs Received SNR
subplot(2, 2, 1);
plot(SNR_rx_vec, SR_vs_SNR(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(SNR_rx_vec, SR_vs_SNR(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlim([min(SNR_rx_vec) max(SNR_rx_vec)]); % <-- DOPASOWANIE OSI DO DANYCH
xlabel('Average Received SNR at Clients (dB)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Sum-Rate vs Received SNR (K = %d)', K_fixed));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'NorthWest');
% 2. Fairness vs Received SNR
subplot(2, 2, 2);
plot(SNR_rx_vec, J_vs_SNR(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(SNR_rx_vec, J_vs_SNR(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlim([min(SNR_rx_vec) max(SNR_rx_vec)]); % <-- DOPASOWANIE OSI DO DANYCH
xlabel('Average Received SNR at Clients (dB)'); ylabel("Jain's index (Fairness)");
title(sprintf('Fairness vs Received SNR (K = %d)', K_fixed));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthWest');
% 3. Sum-Rate vs K
subplot(2, 2, 3);
plot(K_vec, SR_vs_K(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(K_vec, SR_vs_K(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Number of users (K)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Sum-Rate vs K (Mean Rx SNR \\approx %.1f dB)', mean(Rx_SNR_K_vec)));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthWest');
% 4. Fairness vs K
subplot(2, 2, 4);
plot(K_vec, J_vs_K(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(K_vec, J_vs_K(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Number of users (K)'); ylabel("Jain's index (Fairness)");
title(sprintf('Fairness vs K (Mean Rx SNR \\approx %.1f dB)', mean(Rx_SNR_K_vec)));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthWest');
sgtitle(sprintf('Scenario 2: ZF Normalization Trade-off'));
save_figure(fig, 'fig_fairness_normalization');
plot_multi_user_topology(dist_k_sweep, theta_k_sweep, dist_e, theta_e);

% =========================================================================
% FUNKCJA WYKONAWCZA SWEEPÓW (Nagłówek rozszerzony o argument seeds_vec)
% =========================================================================
function [SR_out, J_out, Rx_SNR_dB_avg] = run_cdl_sweep(cdl_b, cdl_e, dist_b, dist_e, fc, Nt, K, P_tx, numIter, seeds_vec)
    SR_acc = zeros(2, 1);
    J_acc  = zeros(2, 1);
    
    % Wyliczanie tłumienia 3GPP wektorowo poza główną pętlą Monte Carlo (Szybkość!)
    [PL_lin_b_vec, ~] = compute_nr_pathloss(dist_b, fc);
    [PL_lin_e, ~]     = compute_nr_pathloss(dist_e, fc);
    
    % Wyznaczenie nominalnego uśrednionego Received SNR u klientów (szum tła = 1)
    % Uwzględnia moc nadawczą P_tx oraz zysk z wieloantenowości Nt (zgodnie z normowaniem kanału)
    Rx_SNR_dB_avg = 10 * log10(mean(P_tx * Nt ./ PL_lin_b_vec));
    
    for it = 1:numIter
        H_eff = zeros(Nt, K);
        
        % [CRN]: Pobranie unikalnego, ale stałego seeda przypisanego do danej iteracji 'it'
        iter_seed = seeds_vec(it);
        
        % Kanały użytkowników (Bobowie)
        for k = 1:K
            release(cdl_b{k});
            cdl_b{k}.Seed = iter_seed;
            [pg_b, ~] = cdl_b{k}();
            h_b = squeeze(sum(pg_b, 2)); h_b = h_b(:);
            
            % Kanał efektywny z tłumieniem 3GPP NR UMi
            H_eff(:, k) = sqrt(1 / PL_lin_b_vec(k)) * h_b;
        end
        
        % Kanał podsłuchiwacza (Ewa)
        release(cdl_e);
        cdl_e.Seed = iter_seed;
        [pg_e, ~] = cdl_e();
        h_e = squeeze(sum(pg_e, 2)); h_e = h_e(:);
        
        % Kanał efektywny Ewy z tłumieniem 3GPP NR UMi
        h_e_eff = sqrt(1 / PL_lin_e) * h_e;
        
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
    cdl.DelaySpread = 92e-9; 
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
% FUNKCJA POMOCNICZA: Generowanie topologii scenariusza (Multi-User)
% =========================================================================
function plot_multi_user_topology(dist_b_vec, theta_b_vec, dist_e, theta_e)
    fig_top = figure('Color', 'w', 'Position', [150 150 700 700]);
    hold on; grid on; box on;
    
    x_bs = 0; y_bs = 0;
    max_d = max([dist_b_vec, dist_e]) + 20;
    
    % Rysowanie Ewy
    x_e = dist_e * sind(theta_e);
    y_e = dist_e * cosd(theta_e);
    p_e = plot(x_e, y_e, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Eve');
    text(x_e + 3, y_e, sprintf('Eve\n(%gm, %g\\circ)', dist_e, theta_e), 'Color', 'r', 'FontSize', 9);
    
    % Rysowanie Bobów
    p_b = [];
    for k = 1:length(dist_b_vec)
        x_b = dist_b_vec(k) * sind(theta_b_vec(k));
        y_b = dist_b_vec(k) * cosd(theta_b_vec(k));
        
        p_b = plot(x_b, y_b, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        
        if k == 1 || k == length(dist_b_vec)
            text(x_b + 3, y_b, sprintf('B_{%d}\n(%gm, %g\\circ)', k, dist_b_vec(k), theta_b_vec(k)), 'Color', 'b', 'FontSize', 8);
        end
    end
    set(p_b, 'DisplayName', sprintf('Bobs (K=%d)', length(dist_b_vec)));
    
    % Rysowanie stacji bazowej
    p_bs = plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    text(x_bs, y_bs - 5, 'BS (0,0)', 'HorizontalAlignment', 'center', 'Color', 'k');
    
    axis equal;
    xlim([-max_d, max_d]);
    ylim([-20, max_d]);
    xlabel('X [m]'); ylabel('Y [m]');
    title(sprintf('Scenario 2: Normalization fairness (K=%d)', length(dist_b_vec)));
    legend([p_bs, p_b, p_e], 'Location', 'NorthWest');
    
    try
        save_figure(fig_top, '../topology/topology_fairness');
    catch
        warning('Funkcja save_figure nie jest dostępna. Wykres nie został zapisany automatycznie.');
    end
end