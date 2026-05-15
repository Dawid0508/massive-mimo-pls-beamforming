% =========================================================================
% SCENARIO: Location-Error Amplification (the "Narrow-Beam" paradox)
% -------------------------------------------------------------------------
% 3GPP CDL channels; received SNR at Bob after FSPL.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

sigma_loc_vec = 0:0.25:5;
theta_b       = -10;
theta_e       =  15;
dist          = p.link_dist_m;
SNR_rx_dB     = 20;
P_rx          = rx_snr_power('linear', SNR_rx_dB);
noise_var     = p.noise_var;
numIter       = 250;

bands = struct( ...
    'name', {'6 GHz (Nt=32)', '28 GHz (Nt=512)'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave}, ...
    'cdl',  {p.cdl_sub6, p.cdl_mmwave});

SR        = zeros(2, length(sigma_loc_vec));
P_on_tgt  = zeros(2, length(sigma_loc_vec));

sigma_snap = [0, 0.5, 2];
angles     = -30:0.05:30;
bp_snap    = zeros(length(sigma_snap), length(angles));
bw_3dB     = 102 ./ [bands.Nt];

for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;  cdl = bands(b).cdl;
    PL_lin = compute_fspl(dist, fc);
    print_scenario_snr('title', sprintf('Location error @ %s', bands(b).name), ...
        'SNR_rx_dB', SNR_rx_dB, 'dist_m', dist, 'fc_Hz', fc, ...
        'actors', {'Bob', 'Eve'});
    [~, sv] = setup_ula(Nt, fc);
    a_sweep = step(sv, fc, angles);

    for s_idx = 1:length(sigma_loc_vec)
        sigma_loc = sigma_loc_vec(s_idx);
        SR_acc = 0; on_acc = 0;
        for it = 1:numIter
            theta_hat = theta_b + sigma_loc * randn;
            h_b = channel_3gpp_ula(sv, fc, theta_b, cdl);
            h_e = channel_3gpp_ula(sv, fc, theta_e, cdl);
            h_hat = channel_3gpp_ula(sv, fc, theta_hat, cdl);

            w = h_hat / norm(h_hat) * sqrt(P_rx);

            R_b = log2(1 + abs(h_b' * w)^2 / noise_var);
            R_e = log2(1 + abs(h_e' * w)^2 / noise_var);
            SR_acc = SR_acc + secrecy_rate(R_b, R_e);

            pat = abs(a_sweep' * w).^2;
            [pmax, ~] = max(pat);
            [~, b_bin] = min(abs(angles - theta_b));
            on_acc = on_acc + (pat(b_bin) >= pmax/2);
        end
        SR(b, s_idx)       = SR_acc / numIter;
        P_on_tgt(b, s_idx) = on_acc / numIter;
    end

    if b == 2
        for ss = 1:length(sigma_snap)
            sigma_loc = sigma_snap(ss);
            bp_acc = zeros(length(angles), 1);
            for it = 1:300
                theta_hat = theta_b + sigma_loc * randn;
                h_hat = channel_3gpp_ula(sv, fc, theta_hat, cdl);
                w = h_hat / norm(h_hat);
                bp_acc = bp_acc + abs(a_sweep' * w).^2;
            end
            bp_snap(ss, :) = 10*log10(bp_acc.' / 300);
        end
    end
end

fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

subplot(2, 2, 1);
plot(sigma_loc_vec, SR(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(sigma_loc_vec, SR(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Location std-dev \sigma_{loc} (deg)');
ylabel('Secrecy Rate (bits/s/Hz)');
title('Secrecy Rate vs pointing error');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

subplot(2, 2, 2);
plot(sigma_loc_vec, P_on_tgt(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(sigma_loc_vec, P_on_tgt(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hbw6 = xline(bw_3dB(1)/2, 'b:', sprintf('BW_{3dB}/2 (6 GHz) \\approx %.2f^{\\circ}', bw_3dB(1)/2), ...
    'Color', [0.35 0.65 1.0], 'LineWidth', 1.2);
hbw28 = xline(bw_3dB(2)/2, 'r:', sprintf('BW_{3dB}/2 (28 GHz) \\approx %.2f^{\\circ}', bw_3dB(2)/2), ...
    'Color', [1.0 0.42 0.35], 'LineWidth', 1.2);
setappdata(hbw6, 'plsConstLabelSide', 'left');
setappdata(hbw28, 'plsConstLabelSide', 'left');
pls_axis_prefs(gca, 'refLabelV', 'bottom', 'refLabelOrient', 'aligned', 'staggerRef', true);
grid on; box on; ylim([0 1.05]);
xlabel('Location std-dev \sigma_{loc} (deg)');
ylabel('P(main lobe covers Bob)');
title('Beam-on-target probability (3 dB criterion)');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

subplot(2, 2, 3);
SR_rel = [SR(1,:)/max(SR(1,1), eps); SR(2,:)/max(SR(2,1), eps)];
plot(sigma_loc_vec, SR_rel(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(sigma_loc_vec, SR_rel(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Location std-dev \sigma_{loc} (deg)');
ylabel('Normalised Secrecy Rate');
title('Relative collapse: narrow-beam paradox');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

subplot(2, 2, 4);
colors = lines(length(sigma_snap));
for ss = 1:length(sigma_snap)
    plot(angles, bp_snap(ss,:) - max(bp_snap(ss,:)), 'LineWidth', 2, 'Color', colors(ss,:)); hold on;
end
mark_bob(theta_b, sprintf('Bob (%d^{\\circ})', theta_b));
mark_eve(theta_e, sprintf('Eve (+%d^{\\circ})', theta_e));
pls_axis_prefs(gca, 'refLabelV', 'top', 'staggerRef', true);
grid on; box on; ylim([-40 5]);
xlabel('Angle (deg)'); ylabel('Normalised gain (dB)');
title('mmWave beam under location error');
legend(arrayfun(@(s) sprintf('\\sigma_{loc} = %.2g^{\\circ}', s), sigma_snap, 'UniformOutput', false), ...
       'Location', 'NorthWest');

sgtitle(sprintf('Location-error (3GPP CDL, received SNR = %d dB, d = %d m)', ...
    SNR_rx_dB, dist));

save_figure(fig, 'fig_location_error');
