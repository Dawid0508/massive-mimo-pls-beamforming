% =========================================================================
% BASELINE: 6 GHz Massive MIMO vs 28 GHz Ultra-Massive MIMO with FSPL
% -------------------------------------------------------------------------
%   * sub-6 GHz : Nt = 32,  TR 38.901 CDL-A-like (rich NLOS)
%   * mmWave   : Nt = 512, TR 38.901 CDL-D-like (strong LOS)
%   * Received SNR swept at Bob after path loss
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

dist       = 1;                          % link distance [m]
SNR_rx_dB  = 40:2:80;                    % received SNR at Bob [dB]
SNR_rx_lin = rx_snr_power('linear', SNR_rx_dB);
numIter    = 200;

bands = struct( ...
    'name', {'6 GHz', '28 GHz'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave}, ...
    'cdl',  {p.cdl_sub6, p.cdl_mmwave});

for b = 1:2
    [bands(b).PL_lin, bands(b).PL_dB] = compute_fspl(dist, bands(b).fc);
    [~, bands(b).sv] = setup_ula(bands(b).Nt, bands(b).fc);
    print_scenario_snr('title', sprintf('Baseline @ %s', bands(b).name), ...
        'SNR_rx_dB', SNR_rx_dB(1), 'dist_m', dist, 'fc_Hz', bands(b).fc, ...
        'actors', {'Bob', 'Eve'}, ...
        'notes', sprintf('SNR sweep (received at Bob): %d:%d:%d dB', ...
            SNR_rx_dB(1), SNR_rx_dB(2)-SNR_rx_dB(1), SNR_rx_dB(end)));
end
fprintf('mmWave path-loss penalty vs 6 GHz: %.2f dB\n', ...
    bands(2).PL_dB - bands(1).PL_dB);

SecrecyCap = zeros(2, 2, length(SNR_rx_dB));    % bands x {MRT,ZF} x SNR

theta_b = -30; theta_e = -40;

for f_idx = 1:2
    fc = bands(f_idx).fc;  Nt = bands(f_idx).Nt;
    sv = bands(f_idx).sv;  L  = bands(f_idx).PL_lin;
    cdl = bands(f_idx).cdl;

    for it = 1:numIter
        hb = channel_3gpp_ula(sv, fc, theta_b, cdl);
        he = channel_3gpp_ula(sv, fc, theta_e, cdl);

        w_mrt = hb / norm(hb);
        P_null = eye(Nt) - (he * (he' / (he' * he)));
        w_zf   = P_null * hb;
        if norm(w_zf) > 1e-9
            w_zf = w_zf / norm(w_zf);
        else
            w_zf = zeros(Nt, 1);
        end

        for s = 1:length(SNR_rx_lin)
            P_rx = SNR_rx_lin(s);   % received signal power scale at Bob

            R_b = log2(1 + P_rx * abs(hb' * w_mrt)^2);
            R_e = log2(1 + P_rx * abs(he' * w_mrt)^2);
            SecrecyCap(f_idx,1,s) = SecrecyCap(f_idx,1,s) + max(0, R_b - R_e);

            R_b = log2(1 + P_rx * abs(hb' * w_zf)^2);
            R_e = log2(1 + P_rx * abs(he' * w_zf)^2);
            SecrecyCap(f_idx,2,s) = SecrecyCap(f_idx,2,s) + max(0, R_b - R_e);
        end
    end
end
SecrecyCap = SecrecyCap / numIter;

snap_theta_b = -30; snap_theta_e = 20;
angles = -90:0.05:90;
fig = figure('Color', 'w', 'Position', [100 100 1100 800]);

for f_idx = 1:2
    fc = bands(f_idx).fc; Nt = bands(f_idx).Nt; sv = bands(f_idx).sv;
    cdl = bands(f_idx).cdl;
    bandStr = [bands(f_idx).name, ' (CDL, PL included)'];

    hb_s = channel_3gpp_ula(sv, fc, snap_theta_b, cdl);
    he_s = channel_3gpp_ula(sv, fc, snap_theta_e, cdl);
    w_mrt_s = hb_s / norm(hb_s);
    P_null_s = eye(Nt) - (he_s * (he_s' / (he_s' * he_s)));
    w_zf_s = P_null_s * hb_s; w_zf_s = w_zf_s / norm(w_zf_s);

    a_sweep = step(sv, fc, angles);
    pat_mrt = 10*log10(abs(w_mrt_s' * a_sweep).^2);
    pat_zf  = 10*log10(abs(w_zf_s'  * a_sweep).^2);

    subplot(2, 2, f_idx);
    plot(SNR_rx_dB, squeeze(SecrecyCap(f_idx,1,:)), 'b-o', 'LineWidth', 1.5); hold on;
    plot(SNR_rx_dB, squeeze(SecrecyCap(f_idx,2,:)), 'r--s', 'LineWidth', 1.5);
    grid on; box on;
    title(['Secrecy: ', bandStr]);
    xlabel('Received SNR at Bob (dB)'); ylabel('bits/s/Hz');
    legend('MRT', 'ZF', 'Location', 'NorthWest');

    subplot(2, 2, f_idx + 2);
    plot(angles, pat_mrt - max(pat_mrt), 'b',  'LineWidth', 1.7); hold on;
    plot(angles, pat_zf  - max(pat_zf),  'r--','LineWidth', 1.5);
    mark_bob(snap_theta_b, 'Bob');
    mark_eve(snap_theta_e, 'Eve');
    pls_axis_prefs(gca, 'refLabelV', 'top', 'staggerRef', true);
    grid on; box on;
    title(['Normalised beam pattern: ', bandStr]);
    xlabel('Angle (deg)'); ylabel('Gain (dB)');
    ylim([-40 5]);
end
sgtitle(sprintf('6G PLS baseline (3GPP CDL) at d = %g m', dist));

save_figure(fig, 'fig_baseline_6GHz_vs_28GHz');
