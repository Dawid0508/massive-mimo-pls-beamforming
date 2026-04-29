% =========================================================================
% SCENARIO: Location-Error Amplification (the "Narrow-Beam" paradox)
% -------------------------------------------------------------------------
% Reference: PLS in Massive MIMO Challenges (ResearchGate 399962830,
%            "Location-Error Amplification Effect");
%            Scenariusze_pomysly.docx (narrow-beam paradox).
%
% In LoS-dominant Massive MIMO the BS often steers a beam from a
% reported user location (e.g. GNSS / radar fix). A location estimate
% with std-dev sigma_loc translates directly into a pointing error.
% The half-power beam-width of a uniform-illumination ULA is
%       BW_3dB ~ 102 / Nt   degrees
% so a 32-antenna array (BW ~ 3.2 deg) tolerates several degrees of
% pointing error, while a 512-antenna mmWave array (BW ~ 0.2 deg) is
% catastrophically sensitive: even 1 deg of localisation error sends
% the main lobe completely off Bob.
%
% This script demonstrates the paradox by sweeping sigma_loc on both
% bands and reporting Secrecy Rate, beam-on-target probability, and
% snapshot beam patterns.
% =========================================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
sigma_loc_vec = 0:0.25:5;             % localisation std-dev [deg]
theta_b       = -10;                  % Bob's true bearing [deg]
theta_e       =  15;                  % Eve's true bearing [deg]
SNR_dB        = 20;
P_tx          = 10^(SNR_dB/10);
noise_var     = p.noise_var;
numIter       = 250;

bands = struct( ...
    'name', {'6 GHz (Nt=32)', '28 GHz (Nt=512)'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave});

SR        = zeros(2, length(sigma_loc_vec));
P_on_tgt  = zeros(2, length(sigma_loc_vec));   % main-lobe-on-Bob probability

% Beam patterns at sigma in {0, 0.5, 2} deg (mmWave)
sigma_snap = [0, 0.5, 2];
angles     = -30:0.05:30;
bp_snap    = zeros(length(sigma_snap), length(angles));

% Theoretical 3 dB beam-width per band (informational)
bw_3dB = 102 ./ [bands.Nt];

% --- Sweep ---------------------------------------------------------------
for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;
    [~, sv] = setup_ula(Nt, fc);

    h_b = step(sv, fc, theta_b);
    h_e = step(sv, fc, theta_e);

    % Reference beam-pattern across angles (used for "on-target" check)
    a_sweep = step(sv, fc, angles);

    for s_idx = 1:length(sigma_loc_vec)
        sigma_loc = sigma_loc_vec(s_idx);
        SR_acc = 0; on_acc = 0;
        for it = 1:numIter
            theta_hat = theta_b + sigma_loc * randn;     % noisy estimate
            h_hat     = step(sv, fc, theta_hat);

            % Single-stream MRT precoder built on the *estimate*
            w = h_hat / norm(h_hat) * sqrt(P_tx);

            R_b = log2(1 + abs(h_b' * w)^2 / noise_var);
            R_e = log2(1 + abs(h_e' * w)^2 / noise_var);
            SR_acc = SR_acc + secrecy_rate(R_b, R_e);

            % "Beam-on-target": is Bob within 3 dB of the peak?
            pat = abs(a_sweep' * w).^2;
            [pmax, ~] = max(pat);
            % Find Bob's angular bin
            [~, b_bin] = min(abs(angles - theta_b));
            on_acc = on_acc + (pat(b_bin) >= pmax/2);
        end
        SR(b, s_idx)       = SR_acc / numIter;
        P_on_tgt(b, s_idx) = on_acc / numIter;
    end

    % --- mmWave beam-pattern snapshots -------------------------------
    if b == 2
        for ss = 1:length(sigma_snap)
            sigma_loc = sigma_snap(ss);
            bp_acc = zeros(length(angles), 1);
            for it = 1:300
                theta_hat = theta_b + sigma_loc * randn;
                h_hat = step(sv, fc, theta_hat);
                w = h_hat / norm(h_hat);
                bp_acc = bp_acc + abs(a_sweep' * w).^2;
            end
            bp_snap(ss, :) = 10*log10(bp_acc.' / 300);
        end
    end
end

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

% Top-left: Secrecy Rate vs sigma_loc on both bands
subplot(2, 2, 1);
plot(sigma_loc_vec, SR(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(sigma_loc_vec, SR(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Location std-dev \sigma_{loc} (deg)');
ylabel('Secrecy Rate (bits/s/Hz)');
title('Secrecy Rate vs pointing error');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

% Top-right: probability the main lobe still covers Bob (3 dB criterion)
subplot(2, 2, 2);
plot(sigma_loc_vec, P_on_tgt(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(sigma_loc_vec, P_on_tgt(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xline(bw_3dB(1)/2, 'b:', sprintf('BW_{3dB}/2 (6 GHz) \\approx %.2f^{\\circ}',  bw_3dB(1)/2));
xline(bw_3dB(2)/2, 'r:', sprintf('BW_{3dB}/2 (28 GHz) \\approx %.3f^{\\circ}', bw_3dB(2)/2));
grid on; box on; ylim([0 1.05]);
xlabel('Location std-dev \sigma_{loc} (deg)');
ylabel('P(main lobe covers Bob)');
title('Beam-on-target probability (3 dB criterion)');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

% Bottom-left: relative degradation, normalised to sigma_loc = 0
subplot(2, 2, 3);
SR_rel = [SR(1,:)/max(SR(1,1), eps); SR(2,:)/max(SR(2,1), eps)];
plot(sigma_loc_vec, SR_rel(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(sigma_loc_vec, SR_rel(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Location std-dev \sigma_{loc} (deg)');
ylabel('Normalised Secrecy Rate');
title('Relative collapse: narrow-beam paradox');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

% Bottom-right: mmWave beam patterns at three sigma values
subplot(2, 2, 4);
colors = lines(length(sigma_snap));
for ss = 1:length(sigma_snap)
    plot(angles, bp_snap(ss,:) - max(bp_snap(ss,:)), 'LineWidth', 2, 'Color', colors(ss,:)); hold on;
end
xline(theta_b, 'g:', sprintf('Bob (%d^{\\circ})', theta_b), 'LineWidth', 1.5);
xline(theta_e, 'k:', sprintf('Eve (+%d^{\\circ})', theta_e), 'LineWidth', 1.5);
grid on; box on; ylim([-40 5]);
xlabel('Angle (deg)'); ylabel('Normalised gain (dB)');
title('mmWave beam under location error');
legend(arrayfun(@(s) sprintf('\\sigma_{loc} = %.2g^{\\circ}', s), sigma_snap, 'UniformOutput', false), ...
       'Location', 'South');

sgtitle(sprintf('Location-error amplification  (\\theta_B = %d^{\\circ}, \\theta_E = +%d^{\\circ})', ...
    theta_b, theta_e));

save_figure(fig, 'fig_location_error');
