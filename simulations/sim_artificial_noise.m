% =========================================================================
% SCENARIO: Artificial Noise (AN) injection in the null space of H'
% -------------------------------------------------------------------------
% PEŁNY FIX: Podejście "Zamrożonego Kanału" i "Rosnącej Koalicji" 
% zastosowane dla OBU SWEEPÓW. Maksymalna szybkość i idealna gładkość.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- PARAMETRY FIZYCZNE I SYSTEMOWE --------------------------------------
dist        = 50;                     
L_values    = 1:2:15;                 
phi_values  = 0.1:0.1:1.0;            
L_fixed     = 7;                      
phi_fixed   = 0.7;                    
K           = 4;                      
numIter     = 100;                     
SNR_tx_dB   = 110;                    
P_tx        = 10^(SNR_tx_dB / 10);
noise_var   = 1;                      
fc          = p.fc_sub6;              
Nt          = p.Nt_sub6;              

[PL_lin, PL_dB] = compute_nr_pathloss(dist, fc);
Rx_SNR_dB = 10 * log10(P_tx * Nt / PL_lin); 

fprintf('\n--- Sztuczny Szum (AN) @ 6 GHz ---\n');
fprintf('  Dystans: %g m (FSPL: %.2f dB)\n', dist, PL_dB);
fprintf('  Transmit SNR (rho_tx): %d dB\n', SNR_tx_dB);
fprintf('  Average Rx SNR: %.1f dB\n', Rx_SNR_dB);

% Wstępna alokacja obiektów 3GPP za pomocą oryginalnej funkcji
max_L = max([max(L_values), L_fixed]);
cdl_b = cell(K, 1);
for k = 1:K, cdl_b{k} = setup_matlab_cdl(Nt, fc); end
cdl_e = cell(max_L, 1);
for e = 1:max_L, cdl_e{e} = setup_matlab_cdl(Nt, fc); end

% =========================================================================
% --- SWEEP A: Wpływ liczby Ew (L) przy stałym podziale mocy (phi = 0.7)
% =========================================================================
fprintf('Rozpoczynam Sweep A (Liczba Ew - Wektorowo!)...\n');
SR_no_AN_A_acc   = zeros(1, length(L_values));
SR_with_AN_A_acc = zeros(1, length(L_values));

for it = 1:numIter
    % Przekazujemy cały wektor L_values do jednej iteracji
    [SR_no_vec, SR_an_vec] = run_AN_sweep_L(cdl_b, cdl_e, PL_lin, Nt, K, ...
                                            L_values, P_tx, phi_fixed, noise_var, max_L);
    SR_no_AN_A_acc   = SR_no_AN_A_acc + SR_no_vec;
    SR_with_AN_A_acc = SR_with_AN_A_acc + SR_an_vec;
end
SR_no_AN_A   = SR_no_AN_A_acc / numIter;
SR_with_AN_A = SR_with_AN_A_acc / numIter;

% =========================================================================
% --- SWEEP B: Wpływ podziału mocy (phi) przy stałej liczbie Ew (L = 7)
% =========================================================================
fprintf('Rozpoczynam Sweep B (Podział mocy phi - Wektorowo!)...\n');
SR_with_AN_B_acc = zeros(1, length(phi_values));

for it = 1:numIter
    % Przekazujemy cały wektor phi_values do jednej iteracji
    [~, SR_an_vec] = run_AN_sweep_phi(cdl_b, cdl_e, PL_lin, Nt, K, ...
                                      L_fixed, P_tx, phi_values, noise_var);
    SR_with_AN_B_acc = SR_with_AN_B_acc + SR_an_vec;
end
SR_with_AN_B = SR_with_AN_B_acc / numIter;

% Znalezienie optymalnego phi
[max_SR, opt_idx] = max(SR_with_AN_B);
opt_phi = phi_values(opt_idx);
fprintf('=> Znaleziono optymalne phi: %.1f (Secrecy Sum-Rate: %.2f bits/s/Hz)\n', opt_phi, max_SR);

% =========================================================================
% --- WIZUALIZACJA --------------------------------------------------------
% =========================================================================
fig = figure('Color', 'w', 'Position', [100 100 1100 450]);

subplot(1, 2, 1);
plot(L_values, SR_no_AN_A,   '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, SR_with_AN_A, '-g^', 'LineWidth', 2, 'MarkerFaceColor', 'g');
grid on; box on;
xlabel('Number of colluding eavesdroppers (L)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Robustness vs L  (\\phi = %.1f)', phi_fixed));
legend('Plain ZF', 'ZF + Artificial Noise', 'Location', 'SouthWest');

subplot(1, 2, 2);
plot(phi_values, SR_with_AN_B, '-mh', 'LineWidth', 2, 'MarkerFaceColor', 'm');
grid on; box on;
xlabel('Data-power fraction \phi (1.0 = No AN)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Power-allocation trade-off  (L = %d Eves)', L_fixed));

xline(phi_fixed, 'k:', sprintf('Referencja (\\phi = %.1f)', phi_fixed), 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(opt_phi, 'r--', sprintf('Optimum (\\phi = %.1f)', opt_phi), 'LineWidth', 2, 'Color', 'r', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center');

sgtitle(sprintf('Scenario 4: Artificial Noise (Rx SNR \\approx %.1f dB)', Rx_SNR_dB));
try save_figure(fig, 'fig_artificial_noise'); catch; end
plot_an_topology(dist, K, 8);

% =========================================================================
% FUNKCJA WYKONAWCZA DLA SWEEPU A: Wektorowa ocena liczby Ew (L)
% =========================================================================
function [SR_no_AN_vec, SR_with_AN_vec] = run_AN_sweep_L(cdl_b, cdl_e, PL_lin, Nt, K, L_vec, P_tx, phi, noise_var, max_L)
    common_seed = randi([0 2^31-1]);
    
    % Generujemy geometrię od razu dla maksymalnej koalicji Ew
    theta_bobs = -60 + 120*rand(1, K);
    target_bob_idx = randi(K, 1, max_L);
    angular_error  = -6 + 12 * rand(1, max_L);
    theta_eves     = theta_bobs(target_bob_idx) + angular_error;
    
    H_eff = zeros(Nt, K);
    for k = 1:K
        release(cdl_b{k});
        cdl_b{k}.TransmitArrayOrientation = [-theta_bobs(k); 0; 0];
        cdl_b{k}.Seed = common_seed; 
        [pg_b, ~] = cdl_b{k}();
        H_eff(:, k) = sqrt(1 / PL_lin) * squeeze(sum(pg_b, 2)); 
    end
    
    G_max = zeros(Nt, max_L);
    for e = 1:max_L
        release(cdl_e{e});
        cdl_e{e}.TransmitArrayOrientation = [-theta_eves(e); 0; 0];
        cdl_e{e}.Seed = common_seed; 
        [pg_e, ~] = cdl_e{e}();
        G_max(:, e) = sqrt(1 / PL_lin) * squeeze(sum(pg_e, 2));
    end
    
    W_raw = H_eff * pinv(H_eff' * H_eff + 1e-9 * eye(K));
    W = W_raw / norm(W_raw, 'fro');
    Z = null(H_eff'); 
    
    % Wstępne wyliczenie bazy dla Boba (stałe niezależnie od L)
    S_b_no = zeros(K, 1); I_b_no = zeros(K, 1);
    S_b_an = zeros(K, 1); I_b_an = zeros(K, 1);
    for k = 1:K
        S_b_no(k) = P_tx * abs(H_eff(:, k)' * W(:, k))^2;
        for j = 1:K, if j ~= k, I_b_no(k) = I_b_no(k) + P_tx * abs(H_eff(:, k)' * W(:, j))^2; end; end
        S_b_an(k) = phi * S_b_no(k);
        I_b_an(k) = phi * I_b_no(k);
    end
    
    SR_no_AN_vec   = zeros(1, length(L_vec));
    SR_with_AN_vec = zeros(1, length(L_vec));
    
    % Pętla wewnętrzna krojąca "zamrożoną" macierz Ew
    for l_idx = 1:length(L_vec)
        L = L_vec(l_idx);
        G_eff = G_max(:, 1:L); % ROSNĄCA KOALICJA - wycinamy podmacierz
        
        R_no_eve_inv = eye(L);
        P_AN = (1 - phi) * P_tx;
        Q_AN = (P_AN / (Nt - K)) * (Z * Z');
        R_an_eve_inv = inv(G_eff' * Q_AN * G_eff + noise_var * eye(L));
        
        R_b_no = log2(1 + S_b_no ./ (I_b_no + noise_var));
        R_b_an = log2(1 + S_b_an ./ (I_b_an + noise_var));
        R_e_no = zeros(K, 1); R_e_an = zeros(K, 1);
        
        for k = 1:K
            h_ek_no = sqrt(P_tx) * G_eff' * W(:, k);
            R_e_no(k) = log2(1 + real(h_ek_no' * R_no_eve_inv * h_ek_no));
            
            h_ek_an = sqrt(phi * P_tx) * G_eff' * W(:, k);
            R_e_an(k) = log2(1 + real(h_ek_an' * R_an_eve_inv * h_ek_an));
        end
        SR_no_AN_vec(l_idx)   = sum(max(0, R_b_no - R_e_no));
        SR_with_AN_vec(l_idx) = sum(max(0, R_b_an - R_e_an));
    end
end

% =========================================================================
% FUNKCJA WYKONAWCZA DLA SWEEPU B: Wektorowa ocena podziału mocy (phi)
% =========================================================================
function [SR_no_AN, SR_with_AN_vec] = run_AN_sweep_phi(cdl_b, cdl_e, PL_lin, Nt, K, num_eve, P_tx, phi_vec, noise_var)
    common_seed = randi([0 2^31-1]);
    
    theta_bobs = -60 + 120*rand(1, K);
    target_bob_idx = randi(K, 1, num_eve);
    angular_error  = -6 + 12 * rand(1, num_eve);
    theta_eves     = theta_bobs(target_bob_idx) + angular_error;
    
    H_eff = zeros(Nt, K);
    for k = 1:K
        release(cdl_b{k});
        cdl_b{k}.TransmitArrayOrientation = [-theta_bobs(k); 0; 0];
        cdl_b{k}.Seed = common_seed; 
        [pg_b, ~] = cdl_b{k}();
        H_eff(:, k) = sqrt(1 / PL_lin) * squeeze(sum(pg_b, 2)); 
    end
    
    G_eff = zeros(Nt, num_eve);
    for e = 1:num_eve
        release(cdl_e{e});
        cdl_e{e}.TransmitArrayOrientation = [-theta_eves(e); 0; 0];
        cdl_e{e}.Seed = common_seed; 
        [pg_e, ~] = cdl_e{e}();
        G_eff(:, e) = sqrt(1 / PL_lin) * squeeze(sum(pg_e, 2));
    end
    
    W_raw = H_eff * pinv(H_eff' * H_eff + 1e-9 * eye(K));
    W = W_raw / norm(W_raw, 'fro');
    Z = null(H_eff'); 
    
    R_no_eve_inv = eye(num_eve); 
    S_b_base = zeros(K, 1); I_b_base = zeros(K, 1);
    R_b_no = zeros(K, 1); R_e_no = zeros(K, 1);
    
    for k = 1:K
        S_b_base(k) = P_tx * abs(H_eff(:, k)' * W(:, k))^2;
        for j = 1:K, if j ~= k, I_b_base(k) = I_b_base(k) + P_tx * abs(H_eff(:, k)' * W(:, j))^2; end; end
        R_b_no(k) = log2(1 + S_b_base(k) / (I_b_base(k) + noise_var));
        h_ek_no = sqrt(P_tx) * G_eff' * W(:, k);
        R_e_no(k) = log2(1 + real(h_ek_no' * R_no_eve_inv * h_ek_no));
    end
    SR_no_AN = sum(max(0, R_b_no - R_e_no));
    
    SR_with_AN_vec = zeros(1, length(phi_vec));
    for ph_idx = 1:length(phi_vec)
        phi = phi_vec(ph_idx);
        P_AN = (1 - phi) * P_tx; 
        Q_AN = (P_AN / (Nt - K)) * (Z * Z'); 
        R_an_eve_inv = inv(G_eff' * Q_AN * G_eff + noise_var * eye(num_eve));
        
        R_b_an = zeros(K, 1); R_e_an = zeros(K, 1);
        for k = 1:K
            R_b_an(k) = log2(1 + (phi * S_b_base(k)) / (phi * I_b_base(k) + noise_var)); 
            h_ek_an = sqrt(phi * P_tx) * G_eff' * W(:, k);
            R_e_an(k) = log2(1 + real(h_ek_an' * R_an_eve_inv * h_ek_an));
        end
        SR_with_AN_vec(ph_idx) = sum(max(0, R_b_an - R_e_an));
    end
end

% =========================================================================
% FUNKCJA POMOCNICZA: Inicjalizacja bazowego modelu CDL
% =========================================================================
function cdl = setup_matlab_cdl(Nt, fc)
    cdl = nrCDLChannel;
    cdl.DelayProfile = 'CDL-A';
    if fc < 10e9, cdl.DelaySpread = 92e-9; else, cdl.DelaySpread = 30e-9; end
    cdl.CarrierFrequency = fc;
    cdl.MaximumDopplerShift = 0;         
    cdl.TransmitAntennaArray.Size = [1 Nt 1 1 1]; 
    cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1]; 
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
    cdl.NumTimeSamples = 1;
    cdl.ChannelFiltering = false; 
end

% =========================================================================
% FUNKCJA POMOCNICZA: Generowanie topologii scenariusza
% =========================================================================
function plot_an_topology(dist, K, num_eve)
    fig_top = figure('Color', 'w', 'Position', [150 150 700 700]);
    hold on; grid on; box on;
    
    theta_bobs = -60 + 120 * rand(1, K);
    target_bob_idx = randi(K, 1, num_eve);
    angular_error  = -6 + 12 * rand(1, num_eve);
    theta_eves     = theta_bobs(target_bob_idx) + angular_error;
    
    x_bs = 0; y_bs = 0;
    max_d = dist + 15;
    
    p_bs = plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    text(x_bs, y_bs - 3, 'BS (0,0)', 'HorizontalAlignment', 'center', 'Color', 'k');
    
    x_eves = dist * sind(theta_eves); y_eves = dist * cosd(theta_eves);
    p_e = [];
    for e = 1:num_eve
        p_e = plot(x_eves(e), y_eves(e), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        text(x_eves(e) + 1.5, y_eves(e) + 1.5, sprintf('E_{%d}', e), 'Color', 'r', 'FontSize', 9);
    end
    if ~isempty(p_e), set(p_e, 'DisplayName', sprintf('Eves (L=%d)', num_eve)); end
    
    x_b = dist * sind(theta_bobs); y_b = dist * cosd(theta_bobs);
    p_b = [];
    for k = 1:K
        p_b = plot(x_b(k), y_b(k), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        text(x_b(k) - 1.5, y_b(k) - 1.5, sprintf('B_{%d}', k), 'Color', 'b', 'FontSize', 9, 'HorizontalAlignment', 'right');
    end
    if ~isempty(p_b), set(p_b, 'DisplayName', sprintf('Bobs (K=%d)', K)); end
    
    axis equal; xlim([-max_d, max_d]); ylim([-10, max_d]);
    xlabel('X [m]'); ylabel('Y [m]');
    title(sprintf('Scenario 4: Artificial Noise'));
    legend([p_bs, p_b, p_e], 'Location', 'NorthWest');
    try save_figure(fig_top, '../topology/topology_artificial_noise'); catch; end
end