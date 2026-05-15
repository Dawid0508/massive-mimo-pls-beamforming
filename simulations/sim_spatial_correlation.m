% =========================================================================
% SCENARIO 1: Spatial correlation impact on Secrecy Rate and Fairness
% -------------------------------------------------------------------------
% Compares matrix (Frobenius) vs vector (per-user) ZF normalization under
% exponential spatial correlation rho. Per-user gain spread prevents the
% two normalizations from coinciding when channels are nearly symmetric.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

Nt          = 64;
K           = 8;
rho_values  = 0:0.1:0.9;
gain_spread_dB = 12;
numIter     = 200;
SNR_rx_dB   = 20;
P_rx        = rx_snr_power('linear', SNR_rx_dB);
noise_var   = p.noise_var;

print_scenario_snr('title', 'Spatial correlation', ...
    'SNR_rx_dB', SNR_rx_dB, ...
    'actors', {sprintf('Bob (K=%d)', K), 'Eve'});

SR_sum      = zeros(2, length(rho_values));
fairness    = zeros(2, length(rho_values));

for r_idx = 1:length(rho_values)
    rho = rho_values(r_idx);
    R = toeplitz(rho.^(0:Nt-1));
    R_sqrt = sqrtm(R);

    SR_sum_acc   = zeros(2, 1);
    fairness_acc = zeros(2, 1);

    for mc = 1:numIter
        H_iid = (randn(Nt, K) + 1j*randn(Nt, K)) / sqrt(2);
        H     = R_sqrt * H_iid;
        H     = apply_user_gain_spread(H, gain_spread_dB);
        h_eve = (randn(Nt, 1) + 1j*randn(Nt, 1)) / sqrt(2);
        h_eve = attenuate_eve_channel(h_eve, p.eve_attn_dB);

        W_raw = H / (H' * H);
        W_mat = zf_precoder_normalize(W_raw, P_rx, 'matrix');
        W_vec = zf_precoder_normalize(W_raw, P_rx, 'vector');

        [sr_m, j_m] = compute_zf_secrecy_metrics(H, h_eve, W_mat, noise_var);
        [sr_v, j_v] = compute_zf_secrecy_metrics(H, h_eve, W_vec, noise_var);
        SR_sum_acc(1)   = SR_sum_acc(1)   + sr_m;
        SR_sum_acc(2)   = SR_sum_acc(2)   + sr_v;
        fairness_acc(1) = fairness_acc(1) + j_m;
        fairness_acc(2) = fairness_acc(2) + j_v;
    end

    SR_sum(:, r_idx)   = SR_sum_acc   / numIter;
    fairness(:, r_idx) = fairness_acc / numIter;
end

fig = figure('Color', 'w', 'Position', [100 100 1100 420]);

subplot(1, 2, 1);
plot(rho_values, SR_sum(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(rho_values, SR_sum(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Spatial correlation \rho'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Secrecy Sum-Rate vs \rho');
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthWest');

subplot(1, 2, 2);
plot(rho_values, fairness(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(rho_values, fairness(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
ylim([0 1.05]);
xlabel('Spatial correlation \rho'); ylabel("Jain's index (Bob rates)");
title("Fairness vs \rho");
legend('Matrix (Frobenius)', 'Vector (per-user)', 'Location', 'SouthWest');

sgtitle(sprintf(['Spatial correlation  (Nt=%d, K=%d, SNR=%d dB, gain spread %d dB)'], ...
    Nt, K, SNR_rx_dB, gain_spread_dB));

save_figure(fig, 'fig_spatial_correlation');
