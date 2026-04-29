% =========================================================================
% SCENARIO: 6 GHz vs 28 GHz under colluding eavesdroppers (worst-case MRC)
% -------------------------------------------------------------------------
% L geographically scattered eavesdroppers cooperate via a back-haul
% link and apply Maximum-Ratio Combining (MRC) — the worst case for the
% transmitter. The combined receive vector for stream k is
%       g_eq(k) = sum_e g_e * (g_e^H w_k)^*       (matched to leakage)
% which gives  SNR_eve_k = sum_e |g_e^H w_k|^2 / noise_var.
% This is the canonical MRC bound used in PLS literature.
%
% We compare a sub-6 GHz Massive-MIMO array (Nt = 32) to a 28 GHz
% Ultra-Massive-MIMO array (Nt = 512). The latter wins because much
% narrower beams reduce per-Eve leakage despite the higher path loss.
% =========================================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
dist        = 30;                  % link distance [m]
L_values    = 1:2:15;              % # colluding Eves
K           = 4;                   % legitimate users
numIter     = 80;
SNR_rx_dB   = 30;                  % desired received SNR after path loss
noise_var   = p.noise_var;

bands = struct( ...
    'name', {'6 GHz (Massive MIMO)', '28 GHz (Ultra-Massive MIMO)'}, ...
    'fc',   {p.fc_sub6,              p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6,              p.Nt_mmwave});

results_SR    = zeros(2, length(L_values));
results_Fair  = zeros(2, length(L_values));

% --- Sweep ---------------------------------------------------------------
for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;
    PL_lin = compute_fspl(dist, fc);
    P_tx   = 10^(SNR_rx_dB/10) * PL_lin;     % transmit power for SNR_rx_dB at user
    [~, sv] = setup_ula(Nt, fc);

    for l_idx = 1:length(L_values)
        num_eve = L_values(l_idx);
        SR_acc = 0; F_acc = 0;
        for it = 1:numIter
            theta_bobs = -60 + 120*rand(1, K);
            theta_eves = -90 + 180*rand(1, num_eve);

            H = step(sv, fc, theta_bobs);   % Nt x K
            G = step(sv, fc, theta_eves);   % Nt x num_eve

            W = H * pinv(H' * H);
            W = W / norm(W, 'fro');         % Frobenius normalisation
            P_eff = P_tx / PL_lin;          % received power scale per stream

            R_b = zeros(K, 1);
            R_e = zeros(K, 1);
            for k = 1:K
                R_b(k) = log2(1 + P_eff * abs(H(:,k)' * W(:,k))^2 / noise_var);

                % MRC bound: sum of leakage powers across all Eves
                snr_eve_k = P_eff * sum(abs(G' * W(:,k)).^2);
                R_e(k)    = log2(1 + snr_eve_k / noise_var);
            end

            R_s = secrecy_rate(R_b, R_e);
            SR_acc = SR_acc + sum(R_s);
            F_acc  = F_acc  + jains_fairness(R_s);
        end
        results_SR(b, l_idx)   = SR_acc / numIter;
        results_Fair(b, l_idx) = F_acc  / numIter;
    end
end

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1100 420]);

subplot(1, 2, 1);
plot(L_values, results_SR(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, results_SR(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Number of colluding eavesdroppers (L)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Security under MRC attack');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

subplot(1, 2, 2);
plot(L_values, results_Fair(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, results_Fair(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
ylim([0 1.05]);
xlabel('Number of colluding eavesdroppers (L)');
ylabel("Jain's fairness index");
title('Fairness across legitimate users');
legend(bands(1).name, bands(2).name, 'Location', 'SouthWest');

sgtitle(sprintf('6 GHz vs 28 GHz under colluding Eves  (K = %d, d = %d m, SNR = %d dB)', ...
    K, dist, SNR_rx_dB));

save_figure(fig, 'fig_colluding_eavesdroppers');
