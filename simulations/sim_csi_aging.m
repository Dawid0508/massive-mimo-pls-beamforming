% =========================================================================
% SCENARIO: CSI aging (Doppler) and L-tap Wiener prediction as defence
% -------------------------------------------------------------------------
% Zaktualizowano: nrCDLChannel, FSPL, Transmit SNR (wzorzec zespołu).
% Wiener weights from Jakes model; time-varying channels from nrCDL + Doppler.
%
% Three strategies under increasing UE velocity:
%   (1) No prediction  — ZF from latest CSI slot
%   (2) Wiener (L taps)— ZF from predicted CSI
%   (3) Perfect CSI    — ZF from true CSI at prediction horizon
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

fc            = p.fc_sub6;
Nt            = p.Nt_sub6;
K             = 4;
dist_b        = 50;
dist_e        = 75;                      
theta_b_ref   = -10;
theta_e       =  20;
v_kmh_vec     = 0:10:120;
v_fixed       = 60;
L_vec         = 1:6;
L_fixed       = 4;
delta_t       = 1e-3;
tau           = 1;
sample_rate   = 1 / delta_t;

SNR_tx_dB     = 100;
P_tx          = 10^(SNR_tx_dB / 10);
noise_var     = 1;
numIter       = p.numIter;

[PL_lin_b, PL_dB_b] = compute_fspl(dist_b, fc);
[PL_lin_e, PL_dB_e] = compute_fspl(dist_e, fc);

fprintf('\n--- CSI aging @ %.0f GHz (nrCDLChannel) ---\n', fc / 1e9);
fprintf('  Transmit SNR: %d dB\n', SNR_tx_dB);
fprintf('  Bob: %g m (FSPL: %.2f dB)\n', dist_b, PL_dB_b);
fprintf('  Eve: %g m (FSPL: %.2f dB)\n', dist_e, PL_dB_e);

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
        run_aging_trial(Nt, K, P_tx, noise_var, numIter, ...
        v, fc, delta_t, tau, L_fixed, sample_rate, ...
        PL_lin_b, PL_lin_e, theta_b_ref, theta_e);
end

for l_idx = 1:length(L_vec)
    [SR_vs_L_no(l_idx), SR_vs_L_wiener(l_idx), SR_vs_L_perfect(l_idx)] = ...
        run_aging_trial(Nt, K, P_tx, noise_var, numIter, ...
        v_fixed, fc, delta_t, tau, L_vec(l_idx), sample_rate, ...
        PL_lin_b, PL_lin_e, theta_b_ref, theta_e);
end

fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

subplot(2, 2, 1);
plot(v_kmh_vec, rho_curve, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b');
grid on; box on;
xlabel('UE velocity (km/h)'); ylabel('J_0(2\pi f_d \Delta t)');
title(sprintf('Jakes correlation @ %.0f GHz, \\Delta t = %.0f ms', fc/1e9, delta_t*1e3));

subplot(2, 2, 2);
plot(v_kmh_vec, SR_no, '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r'); hold on;
plot(v_kmh_vec, SR_wiener, '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
plot(v_kmh_vec, SR_perfect, '-^', 'LineWidth', 2, 'MarkerFaceColor', [0.85 0.75 0.2]);
grid on; box on;
xlabel('UE velocity (km/h)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Three strategies under CSI aging');
legend('No prediction', sprintf('Wiener (L=%d)', L_fixed), 'Perfect CSI', ...
       'Location', 'East');

subplot(2, 2, 3);
plot(L_vec, SR_vs_L_no, '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r'); hold on;
plot(L_vec, SR_vs_L_wiener, '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
plot(L_vec, SR_vs_L_perfect, '-^', 'LineWidth', 2, 'MarkerFaceColor', [0.85 0.75 0.2]);
grid on; box on;
xlabel('Wiener predictor order L'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Diminishing returns of predictor order  (v = %d km/h)', v_fixed));
legend('No prediction', 'Wiener', 'Perfect CSI', 'Location', 'West');

subplot(2, 2, 4);
gain = (SR_wiener - SR_no) ./ max(SR_perfect - SR_no, eps);
plot(v_kmh_vec, gain*100, '-mo', 'LineWidth', 2, 'MarkerFaceColor', 'm');
grid on; box on; ylim([0 110]);
xlabel('UE velocity (km/h)'); ylabel('Recovered fraction (%)');
title('Wiener as a "shield": fraction of perfect-CSI gap closed');

sgtitle(sprintf('CSI aging @ %.0f GHz (nrCDL, SNR_{tx} = %d dB)', fc/1e9, SNR_tx_dB));

save_figure(fig, 'fig_csi_aging');


% =========================================================================
%                          Local helpers
% =========================================================================
function [SR_no, SR_w, SR_p] = run_aging_trial(Nt, K, P_tx, noise_var, numIter, ...
    v_kmh, fc, dt, tau, L, sample_rate, PL_b, PL_e, theta_b_ref, theta_e)

    SR_no = 0; SR_w = 0; SR_p = 0;
    T = L + tau;

    rr = zeros(L, 1);
    R  = zeros(L, L);
    for i = 1:L
        for j = 1:L
            R(i, j) = jakes_correlation(v_kmh, fc, (i-j) * dt);
        end
        rr(i) = jakes_correlation(v_kmh, fc, (tau + (i-1)) * dt);
    end
    w_pred = (R + 1e-9 * eye(L)) \ rr;

    for it = 1:numIter
        theta_bobs = -60 + 120 * rand(1, K);
        theta_bobs(1) = theta_b_ref;

        H_time = zeros(Nt, K, T);
        for k = 1:K
            cdl_k = setup_matlab_cdl_doppler(nrCDLChannel, Nt, fc, ...
                theta_bobs(k), v_kmh, sample_rate, T);
            release(cdl_k);
            cdl_k.Seed = randi([0 2^31-1]);
            [pg_b, ~] = cdl_k();
            for t = 1:T
                ti = min(t, size(pg_b, 1));
                H_time(:, k, t) = cdl_to_heff(pg_b(ti, :, :, :), PL_b, Nt);
            end
        end

        cdl_e = setup_matlab_cdl_doppler(nrCDLChannel, Nt, fc, theta_e, ...
            0, sample_rate, T);
        release(cdl_e);
        cdl_e.Seed = randi([0 2^31-1]);
        [pg_e, ~] = cdl_e();
        h_eve_time = zeros(Nt, T);
        for t = 1:T
            ti = min(t, size(pg_e, 1));
            h_eve_time(:, t) = cdl_to_heff(pg_e(ti, :, :, :), PL_e, Nt);
        end

        H_past    = H_time(:, :, 1:L);
        H_now     = H_time(:, :, L + tau);
        h_eve     = attenuate_eve_channel(h_eve_time(:, L + tau));

        H_used_no = H_time(:, :, L);
        H_used_w  = zeros(Nt, K);
        for k = 1:K
            samples = squeeze(H_past(:, k, end:-1:1));
            H_used_w(:, k) = samples * w_pred;
        end
        H_used_p = H_now;

        SR_no = SR_no + secrecy_zf_tx(H_used_no, H_now, h_eve, P_tx, noise_var);
        SR_w  = SR_w  + secrecy_zf_tx(H_used_w,  H_now, h_eve, P_tx, noise_var);
        SR_p  = SR_p  + secrecy_zf_tx(H_used_p,  H_now, h_eve, P_tx, noise_var);
    end

    SR_no = SR_no / numIter;
    SR_w  = SR_w  / numIter;
    SR_p  = SR_p  / numIter;
end


function h_eff = cdl_to_heff(pg_slice, PL_lin, Nt)
    h = squeeze(sum(pg_slice, 2));
    h = h(:);
    nh = norm(h);
    if nh < 1e-12
        nh = 1;
    end
    h_eff = sqrt((1 / PL_lin) * Nt) * (h / nh);
end


function SR = secrecy_zf_tx(H_est, H_true, h_eve, P_tx, noise_var)
% ZF precoder + secrecy sum-rate (same rate loop as sim_phase_noise.m).
    K = size(H_est, 2);
    W_raw = H_est * pinv(H_est' * H_est + 1e-9 * eye(K));
    W = W_raw / norm(W_raw, 'fro');

    SR = 0;
    for k = 1:K
        sig_b  = P_tx * abs(H_true(:, k)' * W(:, k))^2;
        intf_b = 0;
        for j = 1:K
            if j ~= k
                intf_b = intf_b + P_tx * abs(H_true(:, k)' * W(:, j))^2;
            end
        end
        R_b = log2(1 + sig_b / (intf_b + noise_var));

        sig_e  = P_tx * abs(h_eve' * W(:, k))^2;
        intf_e = 0;
        for j = 1:K
            if j ~= k
                intf_e = intf_e + P_tx * abs(h_eve' * W(:, j))^2;
            end
        end
        R_e = log2(1 + sig_e / (intf_e + noise_var));

        SR = SR + max(0, R_b - R_e);
    end
end


function cdl = setup_matlab_cdl_doppler(cdl, Nt, fc, theta, v_kmh, sample_rate, num_samples)
    release(cdl);
    cdl.DelayProfile = 'CDL-A';
    if fc < 10e9
        cdl.DelaySpread = 30e-9;
    else
        cdl.DelaySpread = 10e-9;
    end
    cdl.CarrierFrequency = fc;
    c = physconst('LightSpeed');
    v_ms = v_kmh / 3.6;
    cdl.MaximumDopplerShift = (v_ms * fc) / c;
    cdl.SampleRate = sample_rate;
    cdl.NumTimeSamples = num_samples;
    cdl.TransmitAntennaArray.Size = [1 Nt 1 1 1];
    cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1];
    cdl.TransmitArrayOrientation = [-theta; 0; 0];
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
    cdl.ChannelFiltering = false;
end
