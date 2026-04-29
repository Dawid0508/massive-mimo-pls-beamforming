% =========================================================================
% SCENARIO: CSI aging (Doppler) and L-tap Wiener prediction as defence
% -------------------------------------------------------------------------
% Reference: J. Zhu, R. Schober, V. Bhargava,
%            "Secure Massive MIMO Under Imperfect CSI:
%             Performance Analysis and Channel Prediction",
%            IEEE TWC 2018 (doc/8543651);
%            Scenariusze_pomysly.docx (CSI aging scenario).
%
% Asymmetry argument (paper Sec. III): outdated CSI hurts Bob (the BS
% beam misses) but Eve estimates her own channel locally and is
% unaffected by the BS's stale view. The paper's remedy is an L-tap
% Wiener predictor of h[tau] from past samples h[0], h[-1], ...,
% h[-L+1] using the Jakes correlation J0(2*pi*f_d*dt).
%
% This script compares three strategies under increasing UE velocity:
%   (1) No prediction      - precoder built on h[0], used at t = tau
%   (2) Wiener predictor   - precoder built on h_hat[tau] from L taps
%   (3) Perfect CSI        - oracle baseline (precoder built on h[tau])
% =========================================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
fc          = p.fc_sub6;              % use 6 GHz for tractable Nt = 32
Nt          = p.Nt_sub6;
K           = 4;
v_kmh_vec   = 0:10:120;               % UE velocity sweep [km/h]
v_fixed     = 60;                     % km/h, used for the L sweep
L_vec       = 1:6;                    % Wiener predictor order
L_fixed     = 4;
delta_t     = 1e-3;                   % time slot duration [s]
tau         = 1;                      % prediction horizon (slots ahead)
SNR_dB      = 20;
P_tot       = 10^(SNR_dB/10);
noise_var   = p.noise_var;
numIter     = 250;

% Result storage
SR_no       = zeros(size(v_kmh_vec));
SR_wiener   = zeros(size(v_kmh_vec));
SR_perfect  = zeros(size(v_kmh_vec));
rho_curve   = zeros(size(v_kmh_vec));      % single-lag Jakes correlation

SR_vs_L_no      = zeros(size(L_vec));
SR_vs_L_wiener  = zeros(size(L_vec));
SR_vs_L_perfect = zeros(size(L_vec));

% --- Sweep A: velocity at L = L_fixed -----------------------------------
for v_idx = 1:length(v_kmh_vec)
    v = v_kmh_vec(v_idx);
    rho_curve(v_idx) = jakes_correlation(v, fc, delta_t);
    [SR_no(v_idx), SR_wiener(v_idx), SR_perfect(v_idx)] = ...
        run_aging_trial(Nt, K, P_tot, noise_var, numIter, ...
                        v, fc, delta_t, tau, L_fixed);
end

% --- Sweep B: predictor order L at v = v_fixed --------------------------
for l_idx = 1:length(L_vec)
    [SR_vs_L_no(l_idx), SR_vs_L_wiener(l_idx), SR_vs_L_perfect(l_idx)] = ...
        run_aging_trial(Nt, K, P_tot, noise_var, numIter, ...
                        v_fixed, fc, delta_t, tau, L_vec(l_idx));
end

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

% Top-left: Jakes correlation vs velocity
subplot(2, 2, 1);
plot(v_kmh_vec, rho_curve, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b');
grid on; box on;
xlabel('UE velocity (km/h)'); ylabel('J_0(2\pi f_d \Delta t)');
title(sprintf('Jakes correlation @ %.0f GHz, \\Delta t = %.0f ms', fc/1e9, delta_t*1e3));

% Top-right: Secrecy vs velocity for the three strategies
subplot(2, 2, 2);
plot(v_kmh_vec, SR_no,      '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r'); hold on;
plot(v_kmh_vec, SR_wiener,  '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
plot(v_kmh_vec, SR_perfect, '-k^', 'LineWidth', 2, 'MarkerFaceColor', 'k');
grid on; box on;
xlabel('UE velocity (km/h)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Three strategies under CSI aging');
legend(sprintf('No prediction  (L=%d-tap Wiener disabled)', L_fixed), ...
       sprintf('Wiener (L = %d)', L_fixed), 'Perfect CSI (oracle)', ...
       'Location', 'NorthEast');

% Bottom-left: Secrecy vs predictor order L at v_fixed
subplot(2, 2, 3);
plot(L_vec, SR_vs_L_no,      '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r'); hold on;
plot(L_vec, SR_vs_L_wiener,  '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
plot(L_vec, SR_vs_L_perfect, '-k^', 'LineWidth', 2, 'MarkerFaceColor', 'k');
grid on; box on;
xlabel('Wiener predictor order L'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Diminishing returns of predictor order  (v = %d km/h)', v_fixed));
legend('No prediction', 'Wiener', 'Perfect CSI', 'Location', 'East');

% Bottom-right: relative gain of Wiener over No-prediction
subplot(2, 2, 4);
gain = (SR_wiener - SR_no) ./ max(SR_perfect - SR_no, eps);
plot(v_kmh_vec, gain*100, '-mo', 'LineWidth', 2, 'MarkerFaceColor', 'm');
grid on; box on; ylim([0 110]);
xlabel('UE velocity (km/h)'); ylabel('Recovered fraction (%)');
title('Wiener as a "shield": fraction of perfect-CSI gap closed');

sgtitle(sprintf('CSI aging at %.0f GHz: Wiener prediction vs Doppler', fc/1e9));

save_figure(fig, 'fig_csi_aging');


% =========================================================================
%                          Local helper
% =========================================================================
function [SR_no, SR_w, SR_p] = run_aging_trial(Nt, K, P_tot, noise_var, ...
                                               numIter, v, fc, dt, tau, L)
% Runs numIter Monte-Carlo trials and returns the three Secrecy Sum-Rates.

    SR_no = 0; SR_w = 0; SR_p = 0;

    % Pre-compute the Wiener weights once per call (geometry-independent
    % under the Jakes model: covariance is the same for every antenna).
    rr = zeros(L, 1);          % cross-corr  E[h[t+tau] h*[t-l]] / E[|h|^2]
    R  = zeros(L, L);          % auto-corr   E[h[t-i] h*[t-j]] / E[|h|^2]
    for i = 1:L
        for j = 1:L
            R(i, j) = jakes_correlation(v, fc, (i-j) * dt);
        end
        % Past samples are at lags 0, -1, ..., -(L-1)*dt
        rr(i) = jakes_correlation(v, fc, (tau + (i-1)) * dt);
    end
    R = R + 1e-9 * eye(L);     % regularisation for v -> 0
    w_pred = R \ rr;           % L x 1 Wiener filter

    for it = 1:numIter
        % --- Generate a length-(L+tau) Jakes-correlated process per
        %     antenna and per user. We synthesise it via the Cholesky
        %     factor of the (L+tau) x (L+tau) Jakes covariance.
        T = L + tau;
        Rt = zeros(T, T);
        for i = 1:T
            for j = 1:T
                Rt(i, j) = jakes_correlation(v, fc, (i-j)*dt);
            end
        end
        Rt = Rt + 1e-9 * eye(T);
        Lt = chol(Rt, 'lower');

        % H_time(:, :, t) is the K-user channel at slot t  (t = 1..T).
        H_time = zeros(Nt, K, T);
        for k = 1:K
            for n = 1:Nt
                z = (randn(T, 1) + 1j*randn(T, 1)) / sqrt(2);
                H_time(n, k, :) = Lt * z;
            end
        end
        h_eve_time = zeros(Nt, T);
        for n = 1:Nt
            z = (randn(T, 1) + 1j*randn(T, 1)) / sqrt(2);
            h_eve_time(n, :) = (Lt * z).';
        end

        % Past observations at t = 1..L,  prediction horizon at t = L+tau
        H_past = H_time(:, :, 1:L);     % Nt x K x L
        H_now  = squeeze(H_time(:, :, L+tau));  % Nt x K (the *true* channel)
        h_eve  = h_eve_time(:, L+tau);

        % --- Strategy 1: no prediction (use the most recent sample) ----
        H_used_no = squeeze(H_past(:, :, L));   % h[0]

        % --- Strategy 2: L-tap Wiener prediction ---------------------
        H_used_w = zeros(Nt, K);
        for k = 1:K
            % Past samples for user k stacked oldest -> newest
            samples = squeeze(H_past(:, k, end:-1:1));   % Nt x L
            H_used_w(:, k) = samples * w_pred;
        end

        % --- Strategy 3: perfect CSI ---------------------------------
        H_used_p = H_now;

        SR_no = SR_no + secrecy_with_estimate(H_used_no, H_now, h_eve, P_tot, noise_var);
        SR_w  = SR_w  + secrecy_with_estimate(H_used_w,  H_now, h_eve, P_tot, noise_var);
        SR_p  = SR_p  + secrecy_with_estimate(H_used_p,  H_now, h_eve, P_tot, noise_var);
    end

    SR_no = SR_no / numIter;
    SR_w  = SR_w  / numIter;
    SR_p  = SR_p  / numIter;
end


function SR = secrecy_with_estimate(H_est, H_true, h_eve, P_tot, noise_var)
% Builds a ZF precoder from H_est, evaluates it on the true channel
% H_true (Bob) and the eavesdropper h_eve, and returns Secrecy Sum-Rate.

    K = size(H_est, 2);
    W_raw = H_est * pinv(H_est' * H_est + 1e-9*eye(K));
    W = W_raw / norm(W_raw, 'fro') * sqrt(P_tot);

    R_b = zeros(K, 1); R_e = zeros(K, 1);
    for k = 1:K
        sig    = abs(H_true(:,k)' * W(:,k))^2;
        intf   = sum(abs(H_true(:,k)' * W).^2) - sig;
        R_b(k) = log2(1 + sig / (intf + noise_var));

        sig_e  = abs(h_eve' * W(:,k))^2;
        intf_e = sum(abs(h_eve' * W).^2) - sig_e;
        R_e(k) = log2(1 + sig_e / (intf_e + noise_var));
    end
    SR = sum(secrecy_rate(R_b, R_e));
end
