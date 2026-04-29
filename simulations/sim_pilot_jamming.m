% =========================================================================
% SCENARIO: Pilot Jamming attack (denial-of-service on the training phase)
% -------------------------------------------------------------------------
% Reference: PLS in Massive MIMO Challenges (ResearchGate 399962830,
%            "Pilot signal jamming"); contrast with sim_pilot_contamination.
%
% Unlike pilot *contamination* (Eve replays Bob's pilot to hijack the
% beam), here Eve emits a *random* Gaussian noise burst during the
% training slot. The BS receives:
%     y_train = sqrt(P_p) * h_b * p + j + n_train
% where ||p|| = sqrt(tau) and j ~ CN(0, P_j * I_Nt) is the jammer.
% LS channel estimation gives
%     h_hat = h_b + (j + n_train) / (sqrt(P_p) * conj(p))
%           ~ h_b + CN(0, sigma_eps^2 * I)
% with sigma_eps^2 = (P_j + sigma^2) / (P_p * tau).
%
% Effect on PLS:
%   * Bob's beam degrades (the precoder points at h_hat, not h_b)
%   * Eve gains nothing directly - she merely creates pointing error
% So Pilot Jamming is a *denial-of-service* attack, while Pilot
% Contamination is a *steering* attack. Both are needed for a full
% threat picture.
%
% Two sweeps:
%   (A) Jammer-to-Pilot Ratio (JPR) on the x-axis at fixed K, Nt
%   (B) JPR vs Nt heat-map showing how Massive MIMO mitigates jamming
%       through training-sequence correlation gain (longer tau).
% =========================================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
Nt          = 64;
K           = 4;
tau         = K;                      % minimum orthogonal pilot length
JPR_dB_vec  = -10:5:30;               % jammer-to-pilot ratio [dB]
SNR_dB      = 20;
P_pilot     = 10^(SNR_dB/10);
P_data      = P_pilot;                % data SNR = pilot SNR (typical)
noise_var   = p.noise_var;
numIter     = 200;

precoders   = {'MRT', 'ZF'};

% Result storage at fixed Nt
R_b   = zeros(2, length(JPR_dB_vec));   % rows = {MRT, ZF}
R_e   = zeros(2, length(JPR_dB_vec));
R_s   = zeros(2, length(JPR_dB_vec));
R_s_perfect = zeros(2, 1);              % perfect-CSI baselines

% Sweep B: JPR vs Nt heat-map (ZF only)
Nt_vec      = [16 32 64 128 256];
JPR_grid_dB = -10:5:30;
SR_grid     = zeros(length(Nt_vec), length(JPR_grid_dB));

% --- Sweep A: JPR sweep at Nt = 64 --------------------------------------
for j_idx = 1:length(JPR_dB_vec)
    P_jam = 10^(JPR_dB_vec(j_idx)/10) * P_pilot;
    sigma_eps2 = (P_jam + noise_var) / (P_pilot * tau);

    acc = zeros(2, 3);  % rows = {MRT,ZF}, cols = {Bob, Eve, Sec}
    for it = 1:numIter
        H     = (randn(Nt, K) + 1j*randn(Nt, K)) / sqrt(2);
        h_eve = (randn(Nt, 1) + 1j*randn(Nt, 1)) / sqrt(2);

        % LS channel estimate per user (independent jamming per pilot)
        E    = sqrt(sigma_eps2/2) * (randn(Nt, K) + 1j*randn(Nt, K));
        Hhat = H + E;

        % --- Build precoders from h_hat ---
        % MRT: matched-filter on the *estimate*
        W_mrt = Hhat ./ vecnorm(Hhat) * sqrt(P_data / K);

        % ZF: pseudo-inverse on the estimate, Frobenius normalised
        W_raw = Hhat * pinv(Hhat' * Hhat);
        W_zf  = W_raw / norm(W_raw, 'fro') * sqrt(P_data);

        for pidx = 1:2
            if pidx == 1, W = W_mrt; else, W = W_zf; end
            R_b_k = zeros(K, 1); R_e_k = zeros(K, 1);
            for k = 1:K
                sig    = abs(H(:,k)' * W(:,k))^2;
                intf   = sum(abs(H(:,k)' * W).^2) - sig;
                R_b_k(k) = log2(1 + sig / (intf + noise_var));

                sig_e  = abs(h_eve' * W(:,k))^2;
                intf_e = sum(abs(h_eve' * W).^2) - sig_e;
                R_e_k(k) = log2(1 + sig_e / (intf_e + noise_var));
            end
            R_s_k = secrecy_rate(R_b_k, R_e_k);
            acc(pidx, 1) = acc(pidx, 1) + sum(R_b_k);
            acc(pidx, 2) = acc(pidx, 2) + sum(R_e_k);
            acc(pidx, 3) = acc(pidx, 3) + sum(R_s_k);
        end
    end
    R_b(:, j_idx) = acc(:, 1) / numIter;
    R_e(:, j_idx) = acc(:, 2) / numIter;
    R_s(:, j_idx) = acc(:, 3) / numIter;
end

% --- Perfect-CSI baselines (no jamming) ---------------------------------
acc = zeros(2, 1);
for it = 1:numIter
    H = (randn(Nt, K) + 1j*randn(Nt, K)) / sqrt(2);
    h_eve = (randn(Nt, 1) + 1j*randn(Nt, 1)) / sqrt(2);

    W_mrt = H ./ vecnorm(H) * sqrt(P_data / K);
    W_raw = H * pinv(H' * H);
    W_zf  = W_raw / norm(W_raw, 'fro') * sqrt(P_data);

    for pidx = 1:2
        if pidx == 1, W = W_mrt; else, W = W_zf; end
        Rsk = 0;
        for k = 1:K
            sig    = abs(H(:,k)' * W(:,k))^2;
            intf   = sum(abs(H(:,k)' * W).^2) - sig;
            R_b_k  = log2(1 + sig / (intf + noise_var));

            sig_e  = abs(h_eve' * W(:,k))^2;
            intf_e = sum(abs(h_eve' * W).^2) - sig_e;
            R_e_k  = log2(1 + sig_e / (intf_e + noise_var));

            Rsk = Rsk + secrecy_rate(R_b_k, R_e_k);
        end
        acc(pidx) = acc(pidx) + Rsk;
    end
end
R_s_perfect = acc / numIter;

% --- Sweep B: heat-map over (Nt, JPR), ZF only -------------------------
for n_idx = 1:length(Nt_vec)
    Nt_b = Nt_vec(n_idx);
    K_b  = K;  tau_b = K_b;
    for j_idx = 1:length(JPR_grid_dB)
        P_jam_b = 10^(JPR_grid_dB(j_idx)/10) * P_pilot;
        sigma_eps2 = (P_jam_b + noise_var) / (P_pilot * tau_b);
        acc_s = 0;
        for it = 1:60               % lighter MC for the heat-map
            H = (randn(Nt_b, K_b) + 1j*randn(Nt_b, K_b)) / sqrt(2);
            h_eve = (randn(Nt_b, 1) + 1j*randn(Nt_b, 1)) / sqrt(2);
            Hhat = H + sqrt(sigma_eps2/2) * (randn(Nt_b, K_b) + 1j*randn(Nt_b, K_b));
            W_raw = Hhat * pinv(Hhat' * Hhat);
            W = W_raw / norm(W_raw, 'fro') * sqrt(P_data);
            Rsk = 0;
            for k = 1:K_b
                sig    = abs(H(:,k)' * W(:,k))^2;
                intf   = sum(abs(H(:,k)' * W).^2) - sig;
                R_b_k  = log2(1 + sig / (intf + noise_var));
                sig_e  = abs(h_eve' * W(:,k))^2;
                intf_e = sum(abs(h_eve' * W).^2) - sig_e;
                R_e_k  = log2(1 + sig_e / (intf_e + noise_var));
                Rsk = Rsk + secrecy_rate(R_b_k, R_e_k);
            end
            acc_s = acc_s + Rsk;
        end
        SR_grid(n_idx, j_idx) = acc_s / 60;
    end
end

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

% Top-left: Bob, Eve, Secrecy under jamming, ZF
subplot(2, 2, 1);
plot(JPR_dB_vec, R_b(2,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(JPR_dB_vec, R_e(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
plot(JPR_dB_vec, R_s(2,:), '-g^', 'LineWidth', 2, 'MarkerFaceColor', 'g');
yline(R_s_perfect(2), 'k:', 'Perfect CSI', 'LineWidth', 1.5);
grid on; box on;
xlabel('Jammer-to-Pilot Ratio (dB)'); ylabel('Rate (bits/s/Hz)');
title('ZF under pilot jamming');
legend('Bob Sum-Rate', 'Eve Sum-Rate', 'Secrecy Sum-Rate', 'Location', 'NorthEast');

% Top-right: MRT vs ZF Secrecy under jamming
subplot(2, 2, 2);
plot(JPR_dB_vec, R_s(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(JPR_dB_vec, R_s(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
yline(R_s_perfect(1), 'b:', sprintf('MRT perfect: %.2f', R_s_perfect(1)), 'LineWidth', 1.2);
yline(R_s_perfect(2), 'r:', sprintf('ZF perfect: %.2f',  R_s_perfect(2)), 'LineWidth', 1.2);
grid on; box on;
xlabel('Jammer-to-Pilot Ratio (dB)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('MRT vs ZF robustness to jamming');
legend('MRT', 'ZF', 'Location', 'NorthEast');

% Bottom-left: heat-map Nt vs JPR
subplot(2, 2, 3);
imagesc(JPR_grid_dB, Nt_vec, SR_grid);
set(gca, 'YDir', 'normal'); colormap(parula); colorbar;
xlabel('Jammer-to-Pilot Ratio (dB)'); ylabel('Antennas N_t');
title('Secrecy Sum-Rate (ZF) vs N_t and JPR');

% Bottom-right: how much does adding antennas help at high JPR?
subplot(2, 2, 4);
plot(Nt_vec, SR_grid(:, end), '-mo', 'LineWidth', 2, 'MarkerFaceColor', 'm'); hold on;
plot(Nt_vec, SR_grid(:, ceil(end/2)), '-co', 'LineWidth', 2, 'MarkerFaceColor', 'c');
plot(Nt_vec, SR_grid(:, 1),   '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
grid on; box on;
xlabel('Antennas N_t'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Massive MIMO as defence against pilot jamming');
legend(sprintf('JPR = %d dB (heavy)', JPR_grid_dB(end)), ...
       sprintf('JPR = %d dB (medium)', JPR_grid_dB(ceil(end/2))), ...
       sprintf('JPR = %d dB (mild)', JPR_grid_dB(1)), ...
       'Location', 'NorthWest');

sgtitle(sprintf('Pilot jamming (DoS attack on training)  (K = %d, \\tau = %d, SNR = %d dB)', ...
    K, tau, SNR_dB));

save_figure(fig, 'fig_pilot_jamming');
