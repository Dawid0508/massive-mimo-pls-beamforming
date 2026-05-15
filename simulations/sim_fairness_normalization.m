% =========================================================================
% SCENARIO: Vector vs Matrix ZF normalization - Sum-Rate vs Fairness
% -------------------------------------------------------------------------
% Matrix: ||W||_F^2 = P_rx  (more power to weak ZF columns)
% Vector: ||W(:,k)||^2 = P_rx/K  (equal power per stream)
%
% Per-user gain spread is applied so the two constraints do not collapse
% to the same precoder under symmetric i.i.d. Rayleigh channels.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

Nt          = 64;
K_fixed     = 8;
SNR_fixed   = 20;
SNR_rx_vec  = 0:5:30;
K_vec       = 2:2:16;
gain_spread_dB = 12;              % total spread across users [dB]
numIter     = 200;
noise_var   = p.noise_var;

print_scenario_snr('title', 'Fairness / ZF normalization', ...
    'SNR_rx_dB', SNR_fixed, ...
    'actors', {sprintf('Bob (K=%d)', K_fixed), 'Eve'}, ...
    'notes', sprintf('SNR sweep %d:%d:%d dB; K sweep %d:%d:%d', ...
        SNR_rx_vec(1), SNR_rx_vec(2)-SNR_rx_vec(1), SNR_rx_vec(end), ...
        K_vec(1), K_vec(2)-K_vec(1), K_vec(end)));

SR_vs_SNR   = zeros(2, length(SNR_rx_vec));
J_vs_SNR    = zeros(2, length(SNR_rx_vec));
SR_vs_K     = zeros(2, length(K_vec));
J_vs_K      = zeros(2, length(K_vec));

for s_idx = 1:length(SNR_rx_vec)
    P_rx = rx_snr_power('linear', SNR_rx_vec(s_idx));
    [SR_vs_SNR(:, s_idx), J_vs_SNR(:, s_idx)] = ...
        run_sweep_point(Nt, K_fixed, P_rx, noise_var, numIter, gain_spread_dB, p.eve_attn_dB);
end

P_rx_fixed = rx_snr_power('linear', SNR_fixed);
for k_idx = 1:length(K_vec)
    [SR_vs_K(:, k_idx), J_vs_K(:, k_idx)] = ...
        run_sweep_point(Nt, K_vec(k_idx), P_rx_fixed, noise_var, numIter, gain_spread_dB, p.eve_attn_dB);
end

fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

subplot(2, 2, 1);
plot(SNR_rx_vec, SR_vs_SNR(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(SNR_rx_vec, SR_vs_SNR(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Received SNR (dB, norm. Rayleigh)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Sum-Rate vs SNR  (K = %d)', K_fixed));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthEast');

subplot(2, 2, 2);
plot(SNR_rx_vec, J_vs_SNR(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(SNR_rx_vec, J_vs_SNR(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Received SNR (dB, norm. Rayleigh)'); ylabel("Jain's index (Bob rates)");
title(sprintf('Fairness vs SNR  (K = %d)', K_fixed));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthWest');

subplot(2, 2, 3);
plot(K_vec, SR_vs_K(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(K_vec, SR_vs_K(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Number of users K'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Sum-Rate vs K  (SNR = %d dB)', SNR_fixed));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthEast');

subplot(2, 2, 4);
plot(K_vec, J_vs_K(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(K_vec, J_vs_K(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Number of users K'); ylabel("Jain's index (Bob rates)");
title(sprintf('Fairness vs K  (SNR = %d dB)', SNR_fixed));
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthWest');

sgtitle(sprintf(['ZF normalization trade-off  (Nt = %d, user gain spread = %d dB)'], ...
    Nt, gain_spread_dB));

save_figure(fig, 'fig_fairness_normalization');


function [SR, J] = run_sweep_point(Nt, K, P_rx, noise_var, numIter, gain_spread_dB, eve_attn_dB)
    SR_acc = zeros(2, 1);
    J_acc  = zeros(2, 1);

    for it = 1:numIter
        H = (randn(Nt, K) + 1j*randn(Nt, K)) / sqrt(2);
        H = apply_user_gain_spread(H, gain_spread_dB);
        h_eve = (randn(Nt, 1) + 1j*randn(Nt, 1)) / sqrt(2);
        h_eve = attenuate_eve_channel(h_eve, eve_attn_dB);

        W_raw = H / (H' * H);

        W_mat = zf_precoder_normalize(W_raw, P_rx, 'matrix');
        W_vec = zf_precoder_normalize(W_raw, P_rx, 'vector');

        [sr_m, j_m] = compute_zf_secrecy_metrics(H, h_eve, W_mat, noise_var);
        [sr_v, j_v] = compute_zf_secrecy_metrics(H, h_eve, W_vec, noise_var);
        SR_acc(1) = SR_acc(1) + sr_m;
        SR_acc(2) = SR_acc(2) + sr_v;
        J_acc(1)  = J_acc(1)  + j_m;
        J_acc(2)  = J_acc(2)  + j_v;
    end

    SR = SR_acc / numIter;
    J  = J_acc  / numIter;
end
