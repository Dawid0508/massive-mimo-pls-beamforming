% =========================================================================
% BASELINE: 6 GHz Massive MIMO vs 28 GHz Ultra-Massive MIMO with FSPL
% -------------------------------------------------------------------------
% Reference setup that fixes the physical-layer parameters used by the
% rest of the project:
%   * sub-6 GHz : Nt = 32   antennas, FSPL @ 1 m baseline
%   * mmWave   : Nt = 512  antennas, FSPL @ 1 m baseline
%
% Two precoders are evaluated per band:
%   - MRT (matched-filter, max signal towards Bob)
%   - ZF  (null-projection towards Eve, full energy on Bob)
% =========================================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

dist     = 1;                          % link distance [m]
SNR_dB   = 40:2:80;                    % transmit SNR range [dB]
SNR_lin  = 10.^(SNR_dB/10);
numIter  = 200;

bands = struct( ...
    'name', {'6 GHz', '28 GHz'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave});

% Path loss + ULA per band
for b = 1:2
    [bands(b).PL_lin, bands(b).PL_dB] = compute_fspl(dist, bands(b).fc);
    [~, bands(b).sv] = setup_ula(bands(b).Nt, bands(b).fc);
end

fprintf('Path loss @  6 GHz : %.2f dB\n', bands(1).PL_dB);
fprintf('Path loss @ 28 GHz : %.2f dB\n', bands(2).PL_dB);
fprintf('mmWave penalty     : %.2f dB\n', bands(2).PL_dB - bands(1).PL_dB);

SecrecyCap = zeros(2, 2, length(SNR_dB));    % bands x {MRT,ZF} x SNR

%% Monte-Carlo
for f_idx = 1:2
    fc = bands(f_idx).fc;  Nt = bands(f_idx).Nt;
    sv = bands(f_idx).sv;  L  = bands(f_idx).PL_lin;

    for it = 1:numIter
        theta_b = -30; theta_e = -40;
        hb = step(sv, fc, theta_b)';   % 1 x Nt
        he = step(sv, fc, theta_e)';

        % Sub-6 GHz: add scattering. mmWave: keep LoS-dominant.
        if f_idx == 1
            hb = 0.8*hb + 0.2*(randn(1,Nt) + 1j*randn(1,Nt))/sqrt(2);
            he = 0.8*he + 0.2*(randn(1,Nt) + 1j*randn(1,Nt))/sqrt(2);
        end

        % --- Precoders ---
        w_mrt = hb' / norm(hb);
        P_null = eye(Nt) - (he' * pinv(he*he') * he);
        w_zf   = P_null * hb';
        if norm(w_zf) > 1e-9
            w_zf = w_zf / norm(w_zf);
        else
            w_zf = zeros(Nt, 1);
        end

        for s = 1:length(SNR_lin)
            P_eff = SNR_lin(s) / L;

            R_b = log2(1 + P_eff * abs(hb * w_mrt)^2);
            R_e = log2(1 + P_eff * abs(he * w_mrt)^2);
            SecrecyCap(f_idx,1,s) = SecrecyCap(f_idx,1,s) + max(0, R_b - R_e);

            R_b = log2(1 + P_eff * abs(hb * w_zf)^2);
            R_e = log2(1 + P_eff * abs(he * w_zf)^2);
            SecrecyCap(f_idx,2,s) = SecrecyCap(f_idx,2,s) + max(0, R_b - R_e);
        end
    end
end
SecrecyCap = SecrecyCap / numIter;

%% Visualisation
snap_theta_b = -30; snap_theta_e = 20;
angles = -90:0.05:90;
fig = figure('Color', 'w', 'Position', [100 100 1100 800]);

for f_idx = 1:2
    fc = bands(f_idx).fc; Nt = bands(f_idx).Nt; sv = bands(f_idx).sv;
    bandStr = [bands(f_idx).name, ' (PL included)'];

    hb_s = step(sv, fc, snap_theta_b)';
    he_s = step(sv, fc, snap_theta_e)';
    w_mrt_s = hb_s' / norm(hb_s);
    P_null_s = eye(Nt) - (he_s' * pinv(he_s*he_s') * he_s);
    w_zf_s = P_null_s * hb_s'; w_zf_s = w_zf_s / norm(w_zf_s);

    a_sweep = step(sv, fc, angles);
    pat_mrt = 10*log10(abs(w_mrt_s' * a_sweep).^2);
    pat_zf  = 10*log10(abs(w_zf_s'  * a_sweep).^2);

    subplot(2, 2, f_idx);
    plot(SNR_dB, squeeze(SecrecyCap(f_idx,1,:)), 'b-o', 'LineWidth', 1.5); hold on;
    plot(SNR_dB, squeeze(SecrecyCap(f_idx,2,:)), 'r--s', 'LineWidth', 1.5);
    grid on; box on;
    title(['Secrecy: ', bandStr]);
    xlabel('SNR (dB)'); ylabel('bits/s/Hz');
    legend('MRT', 'ZF', 'Location', 'NorthWest');

    subplot(2, 2, f_idx + 2);
    plot(angles, pat_mrt - max(pat_mrt), 'b',  'LineWidth', 1.7); hold on;
    plot(angles, pat_zf  - max(pat_zf),  'r--','LineWidth', 1.5);
    xline(snap_theta_b, 'g:', 'Bob');
    xline(snap_theta_e, 'k:', 'Eve');
    grid on; box on;
    title(['Normalised beam pattern: ', bandStr]);
    xlabel('Angle (deg)'); ylabel('Gain (dB)');
    ylim([-40 5]);
end
sgtitle(sprintf('6G PLS baseline at d = %g m', dist));

save_figure(fig, 'fig_baseline_6GHz_vs_28GHz');
