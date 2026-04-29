% =========================================================================
% SCENARIO 1: Spatial correlation impact on Secrecy Rate and Fairness
% -------------------------------------------------------------------------
% Reference: PLS in Massive MIMO Challenges (ResearchGate 399962830);
%            Scenariusze_pomysły.docx, Scenario 1.
%
% This script compares two ZF normalization strategies under increasing
% spatial correlation rho in {0, 0.1, ..., 0.9}:
%   * Matrix normalization  : W = ZF / ||ZF||_F  (joint power constraint)
%   * Vector normalization  : each column of W normalised to unit norm
%
% A passive eavesdropper (Eve) observes the same downlink with an
% independent Rayleigh channel. We report:
%   * Secrecy Sum-Rate         (sum_k max(0, R_b(k) - R_e(k)))
%   * Jain's fairness index    over per-user secrecy rates
% =========================================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- System parameters ---------------------------------------------------
Nt          = 64;                  % BS antennas
K           = 8;                   % legitimate users (Bobs)
rho_values  = 0:0.1:0.9;           % spatial-correlation sweep
numIter     = 200;                 % Monte-Carlo iterations
SNR_dB      = 20;                  % transmit SNR
P_tot       = 10^(SNR_dB/10);      % total transmit power (linear)
noise_var   = p.noise_var;

% Result storage:  rows = {Matrix, Vector} normalization
SR_sum      = zeros(2, length(rho_values));
fairness    = zeros(2, length(rho_values));

% --- Main sweep ----------------------------------------------------------
for r_idx = 1:length(rho_values)
    rho = rho_values(r_idx);

    % Exponential correlation matrix (Toeplitz, rho^|i-j|)
    R = toeplitz(rho.^(0:Nt-1));
    R_sqrt = sqrtm(R);

    SR_sum_acc   = zeros(2, 1);
    fairness_acc = zeros(2, 1);

    for mc = 1:numIter
        % --- Channels -------------------------------------------------
        % Bob: correlated Rayleigh.   Eve: independent Rayleigh.
        H_iid = (randn(Nt, K) + 1j*randn(Nt, K)) / sqrt(2);
        H     = R_sqrt * H_iid;                                    % Nt x K

        h_eve = (randn(Nt, 1) + 1j*randn(Nt, 1)) / sqrt(2);

        % --- ZF precoder (regularised pseudo-inverse) ----------------
        % W_raw = H * (H'*H)^(-1)  ;  pinv handles near-singular H'*H.
        Gram   = H' * H;
        W_raw  = H * pinv(Gram);

        % (1) Matrix normalization  -> joint Frobenius constraint
        W_mat = W_raw / norm(W_raw, 'fro') * sqrt(P_tot);

        % (2) Vector normalization -> equal per-user power P_tot/K
        W_vec = zeros(Nt, K);
        for k = 1:K
            col = W_raw(:, k);
            if norm(col) > 1e-9
                W_vec(:, k) = col / norm(col) * sqrt(P_tot / K);
            end
        end

        for n_idx = 1:2
            if n_idx == 1, W = W_mat; else, W = W_vec; end

            R_b = zeros(K, 1);
            R_e = zeros(K, 1);
            for k = 1:K
                % Bob k: his own signal + intra-cell interference
                sig   = abs(H(:,k)' * W(:,k))^2;
                intf  = sum(abs(H(:,k)' * W).^2) - sig;
                R_b(k) = log2(1 + sig / (intf + noise_var));

                % Eve attempts to decode user k via her own channel
                sig_e  = abs(h_eve' * W(:,k))^2;
                intf_e = sum(abs(h_eve' * W).^2) - sig_e;
                R_e(k) = log2(1 + sig_e / (intf_e + noise_var));
            end

            R_s = secrecy_rate(R_b, R_e);
            SR_sum_acc(n_idx)   = SR_sum_acc(n_idx)   + sum(R_s);
            fairness_acc(n_idx) = fairness_acc(n_idx) + jains_fairness(R_s);
        end
    end

    SR_sum(:, r_idx)   = SR_sum_acc   / numIter;
    fairness(:, r_idx) = fairness_acc / numIter;
end

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1100 420]);

subplot(1, 2, 1);
plot(rho_values, SR_sum(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(rho_values, SR_sum(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Spatial correlation \rho'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Secrecy Sum-Rate vs \rho');
legend('Matrix normalization', 'Vector normalization', 'Location', 'SouthWest');

subplot(1, 2, 2);
plot(rho_values, fairness(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(rho_values, fairness(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
ylim([0 1.05]);
xlabel('Spatial correlation \rho'); ylabel("Jain's fairness index");
title("Fairness vs \rho");
legend('Matrix normalization', 'Vector normalization', 'Location', 'SouthWest');

sgtitle(sprintf('Spatial correlation in Massive MIMO PLS  (Nt=%d, K=%d, SNR=%d dB)', ...
    Nt, K, SNR_dB));

save_figure(fig, 'fig_spatial_correlation');
