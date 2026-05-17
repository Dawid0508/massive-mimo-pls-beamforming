% =========================================================================
% SCENARIO: Location-Error Amplification (the "Narrow-Beam" paradox)
% -------------------------------------------------------------------------
% Zaktualizowano: nrCDLChannel, fizyczny FSPL, Transmit SNR (wzorzec zespołu).
% BS forms MRT from estimate at mis-pointed angle
%   theta_hat = theta_b + N(0, sigma_loc^2).
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

sigma_loc_vec = 0:0.25:5;
theta_b       = -10;
theta_e       =  15;
dist          = p.link_dist_m;
SNR_tx_dB     = 100;
P_tx          = 10^(SNR_tx_dB / 10);
noise_var     = 1;
numIter       = p.numIter;

bands = struct( ...
    'name', {'6 GHz (Nt=32)', '28 GHz (Nt=512)'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave});

SR        = zeros(2, length(sigma_loc_vec));
P_on_tgt  = zeros(2, length(sigma_loc_vec));

sigma_snap = [0, 0.5, 2];
angles     = -30:0.05:30;
bp_snap    = zeros(length(sigma_snap), length(angles));
bw_3dB     = 102 ./ [bands.Nt];

cdl_b   = nrCDLChannel;
cdl_e   = nrCDLChannel;
cdl_hat = nrCDLChannel;

for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;
    [PL_lin, PL_dB] = compute_fspl(dist, fc);

    fprintf('\n--- Location error @ %s ---\n', bands(b).name);
    fprintf('  Transmit SNR: %d dB\n', SNR_tx_dB);
    fprintf('  d = %g m (FSPL: %.2f dB)\n', dist, PL_dB);

    cdl_b   = setup_matlab_cdl(cdl_b, Nt, fc, theta_b);
    cdl_e   = setup_matlab_cdl(cdl_e, Nt, fc, theta_e);
    [~, sv] = setup_ula(Nt, fc);
    a_sweep = step(sv, fc, angles);

    for s_idx = 1:length(sigma_loc_vec)
        sigma_loc = sigma_loc_vec(s_idx);
        SR_acc = 0; on_acc = 0;
        for it = 1:numIter
            theta_hat = theta_b + sigma_loc * randn;

            h_eff_b = draw_cdl_eff(cdl_b, PL_lin, Nt);
            h_eff_e = draw_cdl_eff(cdl_e, PL_lin, Nt);
            cdl_hat = setup_matlab_cdl(cdl_hat, Nt, fc, theta_hat);
            h_eff_hat = draw_cdl_eff(cdl_hat, PL_lin, Nt);

            w = h_eff_hat / norm(h_eff_hat);

            R_b = log2(1 + P_tx * abs(h_eff_b' * w)^2 / noise_var);
            R_e = log2(1 + P_tx * abs(h_eff_e' * w)^2 / noise_var);
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
                cdl_hat = setup_matlab_cdl(cdl_hat, Nt, fc, theta_hat);
                h_eff_hat = draw_cdl_eff(cdl_hat, PL_lin, Nt);
                w = h_eff_hat / norm(h_eff_hat);
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
    'LineWidth', 1.2);
hbw28 = xline(bw_3dB(2)/2, 'r:', sprintf('BW_{3dB}/2 (28 GHz) \\approx %.2f^{\\circ}', bw_3dB(2)/2), ...
    'LineWidth', 1.2);
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

sgtitle(sprintf('Location-error (nrCDL + FSPL, SNR_{tx} = %d dB, d = %d m)', ...
    SNR_tx_dB, dist));

save_figure(fig, 'fig_location_error');


% =========================================================================
%                          Local helpers
% =========================================================================
function h_eff = draw_cdl_eff(cdl, PL_lin, Nt)
    release(cdl);
    cdl.Seed = randi([0 2^31-1]);
    [pg, ~] = cdl();
    h = squeeze(sum(pg, 2)); h = h(:);
    h_eff = sqrt((1 / PL_lin) * Nt) * (h / norm(h));
end

function cdl = setup_matlab_cdl(cdl, Nt, fc, theta)
    release(cdl);
    cdl.DelayProfile = 'CDL-A';
    if fc < 10e9
        cdl.DelaySpread = 30e-9;
    else
        cdl.DelaySpread = 10e-9;
    end
    cdl.CarrierFrequency = fc;
    cdl.MaximumDopplerShift = 0;
    cdl.TransmitAntennaArray.Size = [1 Nt 1 1 1];
    cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1];
    cdl.TransmitArrayOrientation = [-theta; 0; 0];
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
    cdl.ChannelFiltering = false;
end
