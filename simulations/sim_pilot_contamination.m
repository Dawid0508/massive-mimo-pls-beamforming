% =========================================================================
% SCENARIO: Active pilot contamination ("beam hijacking")
% -------------------------------------------------------------------------
% During the uplink training phase Eve transmits the same pilot as Bob,
% scaled by sqrt(beta). The BS estimates the *combined* channel and
% steers a Matched-Filter beam towards a phantom user located between
% Bob and Eve. As beta -> 1 the beam swings towards Eve and the
% Secrecy Rate collapses to zero — the canonical PLS failure mode of
% Massive MIMO with TDD reciprocity.
%
% We run the attack at both 6 GHz (Nt = 32) and 28 GHz (Nt = 512). The
% mmWave array's pencil-thin beams are *more* sensitive to a small
% pointing error, so the secrecy collapse there is more abrupt.
% =========================================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
beta_values = 0:0.05:1;            % relative power of Eve's fake pilot
theta_bob   = -30;                 % Bob's bearing  [deg]
theta_eve   =  20;                 % Eve's bearing  [deg]
numIter     = 120;
SNR_dB      = 20;
P_tx        = 10^(SNR_dB/10);
noise_var   = p.noise_var;
sigma_est   = 0.05;                % training noise std

bands = struct( ...
    'name',  {'6 GHz', '28 GHz'}, ...
    'fc',    {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',    {p.Nt_sub6, p.Nt_mmwave});

angles = -90:0.5:90;
SR_curves   = zeros(2, length(beta_values));
beam_clean  = zeros(2, length(angles));
beam_hijack = zeros(2, length(angles));

% --- Sweep ---------------------------------------------------------------
for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;
    [~, sv] = setup_ula(Nt, fc);

    a_sweep = step(sv, fc, angles);    % Nt x #angles, for beam patterns

    for be_idx = 1:length(beta_values)
        beta = beta_values(be_idx);
        SR_acc = 0;
        bp_acc = zeros(length(angles), 1);
        for it = 1:numIter
            % True channels with mild scattering
            h_b = step(sv, fc, theta_bob);
            h_e = step(sv, fc, theta_eve);
            h_b_c = 0.8*h_b + 0.2*(randn(Nt,1) + 1j*randn(Nt,1))/sqrt(2);
            h_e_c = 0.8*h_e + 0.2*(randn(Nt,1) + 1j*randn(Nt,1))/sqrt(2);

            % Pilot phase: BS estimates Bob's channel but receives
            % Eve's contamination on top of it.
            n_est = sigma_est * (randn(Nt,1) + 1j*randn(Nt,1))/sqrt(2);
            h_est = h_b_c + sqrt(beta) * h_e_c + n_est;

            % MF precoder built on the *contaminated* estimate
            w = h_est / norm(h_est);

            R_b = log2(1 + P_tx * abs(h_b_c' * w)^2 / noise_var);
            R_e = log2(1 + P_tx * abs(h_e_c' * w)^2 / noise_var);
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

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1100 760]);

% Top-left: secrecy collapse vs beta
subplot(2, 2, 1);
plot(beta_values*100, SR_curves(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(beta_values*100, SR_curves(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel("Eve's pilot power \beta (% of Bob)");
ylabel('Secrecy Rate (bits/s/Hz)');
title('Secrecy collapse under pilot contamination');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

% Top-right: secrecy ratio (collapse rate) - useful contrast
subplot(2, 2, 2);
ratio = SR_curves ./ (SR_curves(:,1) + eps);
plot(beta_values*100, ratio(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(beta_values*100, ratio(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel("Eve's pilot power \beta (% of Bob)");
ylabel('Normalised Secrecy Rate');
title('Relative collapse');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

% Bottom: beam patterns clean vs hijacked
for b = 1:2
    subplot(2, 2, 2 + b);
    plot(angles, beam_clean(b,:)  - max(beam_clean(b,:)),  'b-',  'LineWidth', 1.7); hold on;
    plot(angles, beam_hijack(b,:) - max(beam_hijack(b,:)), 'r--', 'LineWidth', 2.0);
    xline(theta_bob, 'g:', sprintf('Bob (%d^{\\circ})', theta_bob), 'LineWidth', 1.5);
    xline(theta_eve, 'k:', sprintf('Eve (+%d^{\\circ})', theta_eve), 'LineWidth', 1.5);
    grid on; box on;
    xlim([-90 90]); ylim([-40 5]);
    xlabel('Angle (deg)'); ylabel('Normalised gain (dB)');
    title(['Beam pattern @ ', bands(b).name]);
    legend('Clean (\beta=0)', 'Hijacked (\beta=1)', 'Location', 'South');
end

sgtitle('Pilot contamination: 6 GHz vs 28 GHz');

save_figure(fig, 'fig_pilot_contamination');
