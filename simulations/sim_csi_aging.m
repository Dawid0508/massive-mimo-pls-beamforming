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
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

fc          = p.fc_sub6;
Nt          = p.Nt_sub6;
K           = 4;
v_kmh_vec   = 0:10:120;
v_fixed     = 60;
L_vec       = 1:6;
L_fixed     = 4;
delta_t     = 1e-3;
tau         = 1;
SNR_rx_dB   = 20;
P_rx        = rx_snr_power('linear', SNR_rx_dB);
noise_var   = p.noise_var;
numIter     = p.numIter;

print_scenario_snr('title', 'CSI aging @ 6 GHz', ...
    'SNR_rx_dB', SNR_rx_dB, ...
    'actors', {sprintf('Bob (K=%d)', K), 'Eve'}, ...
    'notes', 'Abstract Jakes fading; Eve attenuated (eve_attn_dB)');

SR_no       = zeros(size(v_kmh_vec));
SR_wiener   = zeros(size(v_kmh_vec));
SR_perfect  = zeros(size(v_kmh_vec));
rho_curve   = zeros(size(v_kmh_vec));

SR_vs_L_no      = zeros(size(L_vec));
SR_vs_L_wiener  = zeros(size(L_vec));
SR_vs_L_perfect = zeros(size(L_vec));

for v_idx = 1:length(v_kmh_vec)
    v = v_kmh_vec(v_idx);
    rho_curve(v_idx) = jakes_correlation(v, fc, delta_t);
    [SR_no(v_idx), SR_wiener(v_idx), SR_perfect(v_idx)] = ...
        run_aging_trial(Nt, K, P_rx, noise_var, numIter, ...
                        v, fc, delta_t, tau, L_fixed);
end

for l_idx = 1:length(L_vec)
    [SR_vs_L_no(l_idx), SR_vs_L_wiener(l_idx), SR_vs_L_perfect(l_idx)] = ...
        run_aging_trial(Nt, K, P_rx, noise_var, numIter, ...
                        v_fixed, fc, delta_t, tau, L_vec(l_idx));
end

c = pls_colors();
fig = figure('Color', c.bg, 'Position', [100 100 1200 760]);

subplot(2, 2, 1);
plot(v_kmh_vec, rho_curve, '-o', 'LineWidth', 2, 'Color', c.sub6, ...
    'MarkerFaceColor', c.sub6);
grid on; box on;
xlabel('UE velocity (km/h)'); ylabel('J_0(2\pi f_d \Delta t)');
title(sprintf('Jakes correlation @ %.0f GHz, \\Delta t = %.0f ms', fc/1e9, delta_t*1e3));

subplot(2, 2, 2);
plot(v_kmh_vec, SR_no, '-s', 'LineWidth', 2, 'Color', c.eve, ...
    'MarkerFaceColor', c.eve); hold on;
plot(v_kmh_vec, SR_wiener, '-o', 'LineWidth', 2, 'Color', c.bob, ...
    'MarkerFaceColor', c.bob);
plot(v_kmh_vec, SR_perfect, '-^', 'Color', c.perfect, 'LineWidth', 2, ...
    'MarkerFaceColor', c.perfect, 'DisplayName', 'Perfect CSI');
grid on; box on;
xlabel('UE velocity (km/h)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Three strategies under CSI aging');
legend('No prediction', sprintf('Wiener (L=%d)', L_fixed), 'Perfect CSI', ...
       'Location', 'East');

subplot(2, 2, 3);
plot(L_vec, SR_vs_L_no, '-s', 'LineWidth', 2, 'Color', c.eve, ...
    'MarkerFaceColor', c.eve); hold on;
plot(L_vec, SR_vs_L_wiener, '-o', 'LineWidth', 2, 'Color', c.bob, ...
    'MarkerFaceColor', c.bob);
plot(L_vec, SR_vs_L_perfect, '-^', 'Color', c.perfect, 'LineWidth', 2, ...
    'MarkerFaceColor', c.perfect, 'DisplayName', 'Perfect CSI');
grid on; box on;
xlabel('Wiener predictor order L'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Diminishing returns of predictor order  (v = %d km/h)', v_fixed));
legend('No prediction', 'Wiener', 'Perfect CSI', 'Location', 'West');

subplot(2, 2, 4);
gain = (SR_wiener - SR_no) ./ max(SR_perfect - SR_no, eps);
plot(v_kmh_vec, gain*100, '-o', 'LineWidth', 2, 'Color', c.beam, ...
    'MarkerFaceColor', c.beam);
grid on; box on; ylim([0 110]);
xlabel('UE velocity (km/h)'); ylabel('Recovered fraction (%)');
title('Wiener as a "shield": fraction of perfect-CSI gap closed');

sgtitle(sprintf('CSI aging at %.0f GHz: Wiener prediction vs Doppler', fc/1e9));

save_figure(fig, 'fig_csi_aging');


% =========================================================================
%                          Local helpers
% =========================================================================
function [SR_no, SR_w, SR_p] = run_aging_trial(Nt, K, P_rx, noise_var, ...
                                               numIter, v, fc, dt, tau, L)
    SR_no = 0; SR_w = 0; SR_p = 0;

    rr = zeros(L, 1);
    R  = zeros(L, L);
    for i = 1:L
        for j = 1:L
            R(i, j) = jakes_correlation(v, fc, (i-j) * dt);
        end
        rr(i) = jakes_correlation(v, fc, (tau + (i-1)) * dt);
    end
    R = R + 1e-9 * eye(L);
    w_pred = R \ rr;

    T = L + tau;
    Rt = zeros(T, T);
    for i = 1:T
        for j = 1:T
            Rt(i, j) = jakes_correlation(v, fc, (i-j)*dt);
        end
    end
    Rt = Rt + 1e-9 * eye(T);
    Lt = chol(Rt, 'lower');

    for it = 1:numIter
        H_time = zeros(Nt, K, T);
        for k = 1:K
            Z = (randn(T, Nt) + 1j*randn(T, Nt)) / sqrt(2);
            H_time(:, k, :) = reshape((Lt * Z).', [Nt, 1, T]);
        end
        Z_e = (randn(T, Nt) + 1j*randn(T, Nt)) / sqrt(2);
        h_eve_time = (Lt * Z_e).';

        H_past = H_time(:, :, 1:L);
        H_now  = H_time(:, :, L+tau);
        h_eve  = attenuate_eve_channel(h_eve_time(:, L+tau));

        H_used_no = H_time(:, :, L);
        H_used_w = zeros(Nt, K);
        for k = 1:K
            samples = squeeze(H_past(:, k, end:-1:1));
            H_used_w(:, k) = samples * w_pred;
        end
        H_used_p = H_now;

        SR_no = SR_no + secrecy_with_estimate(H_used_no, H_now, h_eve, P_rx, noise_var);
        SR_w  = SR_w  + secrecy_with_estimate(H_used_w,  H_now, h_eve, P_rx, noise_var);
        SR_p  = SR_p  + secrecy_with_estimate(H_used_p,  H_now, h_eve, P_rx, noise_var);
    end

    SR_no = SR_no / numIter;
    SR_w  = SR_w  / numIter;
    SR_p  = SR_p  / numIter;
end


function SR = secrecy_with_estimate(H_est, H_true, h_eve, P_rx, noise_var)
    K = size(H_est, 2);
    W_raw = H_est * pinv(H_est' * H_est + 1e-9*eye(K));
    W = zf_precoder_normalize(W_raw, P_rx, 'matrix');
    SR = compute_zf_secrecy_metrics(H_true, h_eve, W, noise_var);
end
