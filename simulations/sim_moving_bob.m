% =========================================================================
% SCENARIO: Fixed beam, moving Bob (6 GHz vs 28 GHz, 3GPP CDL)
% -------------------------------------------------------------------------
% The BS steers a beam to a fixed bearing (expected Bob location). Bob's
% true angle is swept; Eve stays at a fixed off-boresight angle. Secrecy
% collapses when Bob drifts out of the narrow mmWave main lobe.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

dist           = p.link_dist_m;
theta_beam     = 0;                    % fixed beam boresight [deg]
theta_eve      = 22;                   % Eve bearing [deg]
theta_bob_vec  = -50:2:50;             % Bob true angle sweep [deg]
SNR_rx_dB      = 25;
P_rx           = rx_snr_power('linear', SNR_rx_dB);
noise_var      = p.noise_var;
numIter        = 100;

bands = struct( ...
    'name', {'6 GHz', '28 GHz'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave}, ...
    'cdl',  {p.cdl_sub6, p.cdl_mmwave});

SR_mrt = zeros(2, length(theta_bob_vec));
SR_zf  = zeros(2, length(theta_bob_vec));

for b = 1:2
    fc = bands(b).fc; Nt = bands(b).Nt; cdl = bands(b).cdl;
    print_scenario_snr('title', sprintf('Moving Bob @ %s', bands(b).name), ...
        'SNR_rx_dB', SNR_rx_dB, 'dist_m', dist, 'fc_Hz', fc, ...
        'actors', {'Bob', 'Eve'}, ...
        'notes', sprintf('Fixed beam @ %d deg, Eve @ %d deg', theta_beam, theta_eve));
    [~, sv] = setup_ula(Nt, fc);

    % Fixed beam from estimated geometry at theta_beam
    h_beam = channel_3gpp_ula(sv, fc, theta_beam, cdl);
    w_fixed = h_beam / norm(h_beam);

    for t_idx = 1:length(theta_bob_vec)
        theta_b = theta_bob_vec(t_idx);
        acc_mrt = 0; acc_zf = 0;
        for it = 1:numIter
            h_b = channel_3gpp_ula(sv, fc, theta_b, cdl);
            h_e = channel_3gpp_ula(sv, fc, theta_eve, cdl);

            % MRT-style fixed beam (no re-steering when Bob moves)
            w_mrt = w_fixed;
            R_b = log2(1 + P_rx * abs(h_b' * w_mrt)^2 / noise_var);
            R_e = log2(1 + P_rx * abs(h_e' * w_mrt)^2 / noise_var);
            acc_mrt = acc_mrt + secrecy_rate(R_b, R_e);

            % ZF null toward Eve, but still anchored to beam direction
            P_null = eye(Nt) - (h_e * (h_e' / (h_e' * h_e)));
            w_zf = P_null * h_beam;
            if norm(w_zf) > 1e-9
                w_zf = w_zf / norm(w_zf);
            else
                w_zf = w_fixed;
            end
            R_b = log2(1 + P_rx * abs(h_b' * w_zf)^2 / noise_var);
            R_e = log2(1 + P_rx * abs(h_e' * w_zf)^2 / noise_var);
            acc_zf = acc_zf + secrecy_rate(R_b, R_e);
        end
        SR_mrt(b, t_idx) = acc_mrt / numIter;
        SR_zf(b, t_idx)  = acc_zf  / numIter;
    end
end

% Beam-pattern snapshot (28 GHz) with Bob markers at three angles
angles = -60:0.1:60;
col = pls_colors();
fig = figure('Color', col.bg, 'Position', [100 100 1200 520]);

subplot(1, 2, 1);
plot(theta_bob_vec, SR_mrt(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(theta_bob_vec, SR_mrt(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
plot(theta_bob_vec, SR_zf(1,:),  '--b', 'LineWidth', 1.5);
plot(theta_bob_vec, SR_zf(2,:),  '--r', 'LineWidth', 1.5);
hb = xline(theta_beam, '--', 'Fixed beam', 'Color', col.beam, 'LineWidth', 1.2);
setappdata(hb, 'plsConstLabelSide', 'right');
pls_axis_prefs(gca, 'refLabelV', 'top');
grid on; box on;
xlabel('Bob true angle (deg)');
ylabel('Secrecy Rate (bits/s/Hz)');
title('Secrecy vs Bob position (fixed beam)');
legend('6 GHz MRT', '28 GHz MRT', '6 GHz ZF', '28 GHz ZF', 'Location', 'best');

fc2 = bands(2).fc; Nt2 = bands(2).Nt; cdl2 = bands(2).cdl;
[~, sv2] = setup_ula(Nt2, fc2);
h_beam2 = channel_3gpp_ula(sv2, fc2, theta_beam, cdl2);
w_fix2 = h_beam2 / norm(h_beam2);
a_sweep = step(sv2, fc2, angles);
pat = 10*log10(abs(w_fix2' * a_sweep).^2);
pat = pat - max(pat);

subplot(1, 2, 2);
plot(angles, pat, '-', 'Color', col.fg, 'LineWidth', 2); hold on;
bob_mark = [-20, 0, 20];
for bm = bob_mark
    mark_bob(bm, sprintf('Bob @ %d^\\circ', bm));
end
xline(theta_beam, '--', 'Beam', 'Color', col.beam, 'LineWidth', 1.5);
mark_eve(theta_eve, 'Eve');
grid on; box on;
xlabel('Angle (deg)'); ylabel('Gain (dB)');
title('28 GHz fixed beam (snapshot)');
ylim([-40 5]);

sgtitle(sprintf('Moving Bob: beam @ %d^\\circ, Eve @ %d^\\circ (received SNR = %d dB)', ...
    theta_beam, theta_eve, SNR_rx_dB));

save_figure(fig, 'fig_moving_bob');
