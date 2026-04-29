% =========================================================================
% SCENARIO: Vector vs Matrix ZF normalization - Sum-Rate vs Fairness
% -------------------------------------------------------------------------
% Reference: Scenariusze_pomysly.docx (Scenario 3 - Jain's Fairness Index);
%            PLS in Massive MIMO Challenges (ResearchGate 399962830).
%
% Two ZF normalization strategies are compared in a clean Rayleigh
% multi-user setting (no correlation) so the difference is *only* the
% per-user power constraint:
%
%   Matrix normalization : W = ZF / ||ZF||_F * sqrt(P_tot)
%       - joint power constraint, exploits weakest user as a sink
%       - typically higher Secrecy Sum-Rate but uneven per-user rates
%
%   Vector normalization : each column scaled to sqrt(P_tot/K)
%       - equal per-stream power, "fair" allocation
%       - usually lower Sum-Rate but a Jain index close to 1
%
% The script sweeps SNR at fixed K (left column) and K at fixed SNR
% (right column) and reports both Secrecy Sum-Rate and Jain's index.
% =========================================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
Nt          = 64;                       % BS antennas
K_fixed     = 8;                        % users for the SNR sweep
SNR_fixed   = 20;                       % dB, used in the K sweep
SNR_dB_vec  = 0:5:30;
K_vec       = 2:2:16;                   % must satisfy K <= Nt
numIter     = 200;
noise_var   = p.noise_var;

% Result storage:  rows = {Matrix, Vector}, cols = sweep index
SR_vs_SNR   = zeros(2, length(SNR_dB_vec));
J_vs_SNR    = zeros(2, length(SNR_dB_vec));
SR_vs_K     = zeros(2, length(K_vec));
J_vs_K      = zeros(2, length(K_vec));

% --- Sweep A: SNR at K = K_fixed ----------------------------------------
for s_idx = 1:length(SNR_dB_vec)
    P_tot = 10^(SNR_dB_vec(s_idx)/10);
    [SR_vs_SNR(:, s_idx), J_vs_SNR(:, s_idx)] = ...
        run_sweep_point(Nt, K_fixed, P_tot, noise_var, numIter);
end

% --- Sweep B: K at SNR = SNR_fixed --------------------------------------
P_tot_fixed = 10^(SNR_fixed/10);
for k_idx = 1:length(K_vec)
    [SR_vs_K(:, k_idx), J_vs_K(:, k_idx)] = ...
        run_sweep_point(Nt, K_vec(k_idx), P_tot_fixed, noise_var, numIter);
end

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

subplot(2, 2, 1);
plot(SNR_dB_vec, SR_vs_SNR(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(SNR_dB_vec, SR_vs_SNR(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Transmit SNR (dB)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Sum-Rate vs SNR  (K = %d)', K_fixed));
legend('Matrix', 'Vector', 'Location', 'NorthWest');

subplot(2, 2, 2);
plot(SNR_dB_vec, J_vs_SNR(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(SNR_dB_vec, J_vs_SNR(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Transmit SNR (dB)'); ylabel("Jain's fairness index");
title(sprintf('Fairness vs SNR  (K = %d)', K_fixed));
legend('Matrix', 'Vector', 'Location', 'SouthEast');

subplot(2, 2, 3);
plot(K_vec, SR_vs_K(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(K_vec, SR_vs_K(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Number of users K'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Sum-Rate vs K  (SNR = %d dB)', SNR_fixed));
legend('Matrix', 'Vector', 'Location', 'NorthWest');

subplot(2, 2, 4);
plot(K_vec, J_vs_K(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(K_vec, J_vs_K(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Number of users K'); ylabel("Jain's fairness index");
title(sprintf('Fairness vs K  (SNR = %d dB)', SNR_fixed));
legend('Matrix', 'Vector', 'Location', 'SouthEast');

sgtitle(sprintf('ZF normalization trade-off: throughput vs fairness  (Nt = %d)', Nt));

save_figure(fig, 'fig_fairness_normalization');


% =========================================================================
%                          Local helper
% =========================================================================
function [SR, J] = run_sweep_point(Nt, K, P_tot, noise_var, numIter)
% Returns SR=[SR_matrix; SR_vector], J=[J_matrix; J_vector] averaged over
% numIter Monte-Carlo realisations.

    SR = zeros(2, 1);
    J  = zeros(2, 1);

    for it = 1:numIter
        % Bobs: i.i.d. Rayleigh.   Eve: i.i.d. Rayleigh, single antenna.
        H     = (randn(Nt, K) + 1j*randn(Nt, K)) / sqrt(2);
        h_eve = (randn(Nt, 1) + 1j*randn(Nt, 1)) / sqrt(2);

        W_raw = H * pinv(H' * H);

        W_mat = W_raw / norm(W_raw, 'fro') * sqrt(P_tot);
        W_vec = zeros(Nt, K);
        for k = 1:K
            col = W_raw(:, k);
            if norm(col) > 1e-9
                W_vec(:, k) = col / norm(col) * sqrt(P_tot / K);
            end
        end

        for n_idx = 1:2
            if n_idx == 1, W = W_mat; else, W = W_vec; end

            R_b = zeros(K, 1); R_e = zeros(K, 1);
            for k = 1:K
                sig    = abs(H(:,k)' * W(:,k))^2;
                intf   = sum(abs(H(:,k)' * W).^2) - sig;
                R_b(k) = log2(1 + sig / (intf + noise_var));

                sig_e  = abs(h_eve' * W(:,k))^2;
                intf_e = sum(abs(h_eve' * W).^2) - sig_e;
                R_e(k) = log2(1 + sig_e / (intf_e + noise_var));
            end

            R_s = secrecy_rate(R_b, R_e);
            SR(n_idx) = SR(n_idx) + sum(R_s);
            J(n_idx)  = J(n_idx)  + jains_fairness(R_s);
        end
    end

    SR = SR / numIter;
    J  = J  / numIter;
end
