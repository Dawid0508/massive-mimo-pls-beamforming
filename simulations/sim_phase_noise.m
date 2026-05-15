% =========================================================================
% SCENARIO: Hardware phase noise impact on Massive MIMO PLS
% -------------------------------------------------------------------------
% Reference: Scenariusze_pomysly.docx (Scenario 4 - Phase Noise).
%
% Real RF chains add per-antenna phase jitter. We model it as a
% multiplicative perturbation on the precoder columns:
%
%       W_err = W .* exp(1j * sigma_phi * randn(size(W)))
%
% where sigma_phi is the std-dev of the per-element phase error in
% radians. For mmWave Ultra-Massive MIMO this is critical because the
% beam is so narrow that even degree-scale jitter can spray energy
% towards the eavesdropper. We compare 6 GHz (Nt = 32) and 28 GHz
% (Nt = 512) and contrast Matrix vs Vector ZF normalization.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
sigma_deg_vec = 0:1:15;                     % phase-error std [deg]
K             = 4;                          % users
numIter       = 80;
dist          = p.link_dist_m;
SNR_rx_dB     = 20;
noise_var     = p.noise_var;

bands = struct( ...
    'name',  {'6 GHz (Nt=32)', '28 GHz (Nt=512)'}, ...
    'fc',    {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',    {p.Nt_sub6, p.Nt_mmwave}, ...
    'cdl',   {p.cdl_sub6, p.cdl_mmwave});

SR_results = zeros(2, 2, length(sigma_deg_vec));   % {6/28 GHz} x {Matrix,Vector} x sigma

% Beam-pattern snapshot at sigma in {0, 5, 10} deg (band 2 only - mmWave)
sigma_snap_deg = [0, 5, 10];
angles    = -90:0.25:90;
bp_snap   = zeros(length(sigma_snap_deg), length(angles));
theta_b   = -10;
theta_e   =  15;

% --- Sweep ---------------------------------------------------------------
for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;  cdl = bands(b).cdl;
    P_rx = rx_snr_power('linear', SNR_rx_dB);
    print_scenario_snr('title', sprintf('Phase noise @ %s', bands(b).name), ...
        'SNR_rx_dB', SNR_rx_dB, 'dist_m', dist, 'fc_Hz', fc, ...
        'actors', {sprintf('Bob (K=%d)', K), 'Eve'});
    [~, sv] = setup_ula(Nt, fc);

    for s_idx = 1:length(sigma_deg_vec)
        sigma_phi = deg2rad(sigma_deg_vec(s_idx));
        SR_acc = zeros(2, 1);

        for it = 1:numIter
            theta_bobs = -60 + 120*rand(1, K);
            H = zeros(Nt, K);
            for k = 1:K
                H(:, k) = channel_3gpp_ula(sv, fc, theta_bobs(k), cdl);
            end
            h_eve = channel_3gpp_ula(sv, fc, -90 + 180*rand, cdl);

            W_raw = H * pinv(H' * H);
            W_mat = W_raw / norm(W_raw, 'fro') * sqrt(P_rx);
            W_vec = zeros(Nt, K);
            for k = 1:K
                col = W_raw(:, k);
                if norm(col) > 1e-9
                    W_vec(:, k) = col / norm(col) * sqrt(P_rx / K);
                end
            end

            % Inject independent per-antenna phase error per stream
            err_mat = exp(1j * sigma_phi * randn(Nt, K));
            W_mat_e = W_mat .* err_mat;
            W_vec_e = W_vec .* err_mat;

            for n_idx = 1:2
                if n_idx == 1, W = W_mat_e; else, W = W_vec_e; end
                R_b = zeros(K, 1); R_e = zeros(K, 1);
                for k = 1:K
                    sig    = P_rx * abs(H(:,k)' * W(:,k))^2;
                    intf   = P_rx * (sum(abs(H(:,k)' * W).^2) - abs(H(:,k)' * W(:,k))^2);
                    R_b(k) = log2(1 + sig / (intf + noise_var));

                    sig_e  = P_rx * abs(h_eve' * W(:,k))^2;
                    intf_e = P_rx * (sum(abs(h_eve' * W).^2) - abs(h_eve' * W(:,k))^2);
                    R_e(k) = log2(1 + sig_e / (intf_e + noise_var));
                end
                SR_acc(n_idx) = SR_acc(n_idx) + sum(secrecy_rate(R_b, R_e));
            end
        end
        SR_results(b, :, s_idx) = SR_acc / numIter;
    end

    % Beam-pattern snapshots (mmWave only)
    if b == 2
        h_b = channel_3gpp_ula(sv, fc, theta_b, cdl);
        h_e = channel_3gpp_ula(sv, fc, theta_e, cdl);
        % Single-stream MRT-style design just to visualise distortion
        w0 = h_b / norm(h_b);
        a_sweep = step(sv, fc, angles);
        for ss = 1:length(sigma_snap_deg)
            sigma_phi = deg2rad(sigma_snap_deg(ss));
            bp_acc = zeros(length(angles), 1);
            for it = 1:200
                w = w0 .* exp(1j * sigma_phi * randn(Nt, 1));
                bp_acc = bp_acc + abs(a_sweep' * w).^2;
            end
            bp_snap(ss, :) = 10*log10(bp_acc.' / 200);
        end
    end
end

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

for b = 1:2
    subplot(2, 2, b);
    plot(sigma_deg_vec, squeeze(SR_results(b,1,:)), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
    plot(sigma_deg_vec, squeeze(SR_results(b,2,:)), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    grid on; box on;
    xlabel('Phase-error std \sigma_\phi (deg)');
    ylabel('Secrecy Sum-Rate (bits/s/Hz)');
    title(['Robustness to phase noise: ', bands(b).name]);
    legend('Matrix normalization', 'Vector normalization', 'Location', 'NorthEast');
end

% Bottom-left: relative degradation (normalised to sigma = 0, Matrix norm.)
subplot(2, 2, 3);
deg6  = squeeze(SR_results(1, 1, :)) / max(SR_results(1, 1, 1), eps);
deg28 = squeeze(SR_results(2, 1, :)) / max(SR_results(2, 1, 1), eps);
plot(sigma_deg_vec, deg6,  '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(sigma_deg_vec, deg28, '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Phase-error std \sigma_\phi (deg)');
ylabel('Normalised Secrecy Sum-Rate');
title('Relative degradation (Matrix norm.)');
legend('6 GHz (Nt=32)', '28 GHz (Nt=512)', 'Location', 'SouthWest');

% Bottom-right: beam-pattern degradation at mmWave
subplot(2, 2, 4);
colors = lines(length(sigma_snap_deg));
for ss = 1:length(sigma_snap_deg)
    plot(angles, bp_snap(ss,:) - max(bp_snap(ss,:)), 'LineWidth', 2, 'Color', colors(ss,:)); hold on;
end
mark_eve(theta_e, sprintf('Eve (+%d^{\\circ})', theta_e), 'left');
mark_bob(theta_b, sprintf('Bob (%d^{\\circ})', theta_b), 'right');
pls_axis_prefs(gca, 'refLabelV', 'top');
grid on; box on;
xlim([-90 90]); ylim([-40 5]);
xlabel('Angle (deg)'); ylabel('Normalised gain (dB)');
title('mmWave beam smearing under phase noise');
legend(arrayfun(@(s) sprintf('\\sigma_\\phi = %d^{\\circ}', s), sigma_snap_deg, 'UniformOutput', false), ...
       'Location', 'SouthWest');

sgtitle(sprintf('Phase noise (3GPP CDL, K = %d, received SNR = %d dB, d = %d m)', ...
    K, SNR_rx_dB, dist));

save_figure(fig, 'fig_phase_noise');
