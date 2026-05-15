% =========================================================================
% SCENARIO: Active pilot contamination ("beam hijacking")
% -------------------------------------------------------------------------
% 3GPP CDL channels at 6 GHz vs 28 GHz; received SNR at Bob after FSPL.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

beta_values = 0:0.05:1;
theta_bob   = -30;
theta_eve   =  20;
dist        = p.link_dist_m;
numIter     = 120;
SNR_rx_dB   = 20;
noise_var   = p.noise_var;
sigma_est   = 0.05;

bands = struct( ...
    'name',  {'6 GHz', '28 GHz'}, ...
    'fc',    {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',    {p.Nt_sub6, p.Nt_mmwave}, ...
    'cdl',   {p.cdl_sub6, p.cdl_mmwave});

angles = -90:0.5:90;
SR_curves   = zeros(2, length(beta_values));
beam_clean  = zeros(2, length(angles));
beam_hijack = zeros(2, length(angles));

for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;  cdl = bands(b).cdl;
    PL_lin = compute_fspl(dist, fc);
    P_rx   = rx_snr_power('linear', SNR_rx_dB);
    print_scenario_snr('title', sprintf('Pilot contamination @ %s', bands(b).name), ...
        'SNR_rx_dB', SNR_rx_dB, 'dist_m', dist, 'fc_Hz', fc, ...
        'actors', {'Bob', 'Eve'});
    [~, sv] = setup_ula(Nt, fc);
    a_sweep = step(sv, fc, angles);

    for be_idx = 1:length(beta_values)
        beta = beta_values(be_idx);
        SR_acc = 0;
        bp_acc = zeros(length(angles), 1);
        for it = 1:numIter
            h_b = channel_3gpp_ula(sv, fc, theta_bob, cdl);
            h_e = channel_3gpp_ula(sv, fc, theta_eve, cdl);

            n_est = sigma_est * (randn(Nt,1) + 1j*randn(Nt,1))/sqrt(2);
            h_est = h_b + sqrt(beta) * h_e + n_est;
            w = h_est / norm(h_est);

            R_b = log2(1 + P_rx * abs(h_b' * w)^2 / noise_var);
            R_e = log2(1 + P_rx * abs(h_e' * w)^2 / noise_var);
            SR_acc = SR_acc + secrecy_rate(R_b, R_e);
            bp_acc = bp_acc + abs(a_sweep' * w).^2;
        end
        SR_curves(b, be_idx) = SR_acc / numIter;
        bp_avg = bp_acc / numIter;

        if abs(beta) < 1e-6
            beam_clean(b, :) = 10*log10(bp_avg).';
        elseif abs(beta - 1) < 1e-6
            beam_hijack(b, :) = 10*log10(bp_avg).';
        end
    end
end

fig = figure('Color', 'w', 'Position', [100 100 1100 760]);

subplot(2, 2, 1);
plot(beta_values*100, SR_curves(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(beta_values*100, SR_curves(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel("Eve's pilot power \beta (% of Bob)");
ylabel('Secrecy Rate (bits/s/Hz)');
title('Secrecy collapse under pilot contamination');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

subplot(2, 2, 2);
ratio = SR_curves ./ (SR_curves(:,1) + eps);
plot(beta_values*100, ratio(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(beta_values*100, ratio(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel("Eve's pilot power \beta (% of Bob)");
ylabel('Normalised Secrecy Rate');
title('Relative collapse');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

for b = 1:2
    subplot(2, 2, 2 + b);
    plot(angles, beam_clean(b,:)  - max(beam_clean(b,:)),  'b-',  'LineWidth', 1.7); hold on;
    plot(angles, beam_hijack(b,:) - max(beam_hijack(b,:)), 'r--', 'LineWidth', 2.0);
    mark_bob(theta_bob, sprintf('Bob (%d^{\\circ})', theta_bob), 'right');
    mark_eve(theta_eve, sprintf('Eve (+%d^{\\circ})', theta_eve), 'left');
    pls_axis_prefs(gca, 'refLabelV', 'top');
    grid on; box on;
    xlim([-90 90]); ylim([-40 5]);
    xlabel('Angle (deg)'); ylabel('Normalised gain (dB)');
    title(['Beam pattern @ ', bands(b).name]);
    legend('Clean (\beta=0)', 'Hijacked (\beta=1)', 'Location', 'South');
end

sgtitle(sprintf('Pilot contamination (3GPP CDL, received SNR = %d dB, d = %d m)', ...
    SNR_rx_dB, dist));

save_figure(fig, 'fig_pilot_contamination');
