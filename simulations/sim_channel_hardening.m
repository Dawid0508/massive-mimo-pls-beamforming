% =========================================================================
% SCENARIO: Channel hardening in Massive MIMO and its impact on PLS
% -------------------------------------------------------------------------
% For h ~ CN(0, I_Nt) the law of large numbers gives
%       ||h||^2 / Nt  ->  1     (a.s. as Nt -> inf)
%       Var(||h||^2 / Nt) = 1 / Nt
%
% In other words a Massive MIMO channel "hardens": the random gain
% becomes deterministic. Two PLS-relevant consequences:
%   1) Outage probability of the legitimate user collapses, lower
%      transmit-power margins are needed.
%   2) Bob's rate variance shrinks but Eve's does not - this is what
%      gives Massive MIMO its native PLS advantage.
%
% This script:
%   - shows histograms of the normalised channel gain at small vs
%     large Nt (visual "hardening proof")
%   - sweeps Nt and reports both Var(||h||^2/Nt) (fundamental) and
%     Std(R_secrecy) (operational consequence)
% =========================================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
Nt_vec       = [4 8 16 32 64 128 256 512];
Nt_hist      = [4, 256];                 % values shown as histograms
numIter      = 4000;                     % heavy MC for clean variance curves
SNR_dB       = 20;
P_tot        = 10^(SNR_dB/10);
noise_var    = p.noise_var;

var_norm_h   = zeros(size(Nt_vec));      % Var(||h||^2 / Nt)
mean_norm_h  = zeros(size(Nt_vec));
std_R_sec    = zeros(size(Nt_vec));      % Std(R_secrecy) under MRT
mean_R_sec   = zeros(size(Nt_vec));

hist_data    = cell(length(Nt_hist), 1);

% --- Sweep ---------------------------------------------------------------
for n_idx = 1:length(Nt_vec)
    Nt = Nt_vec(n_idx);
    g_samples = zeros(numIter, 1);
    R_samples = zeros(numIter, 1);
    for it = 1:numIter
        h_b = (randn(Nt,1) + 1j*randn(Nt,1)) / sqrt(2);
        h_e = (randn(Nt,1) + 1j*randn(Nt,1)) / sqrt(2);

        g_samples(it) = (h_b' * h_b) / Nt;          % normalised gain

        % MRT towards Bob (single-user illustrative case)
        w   = h_b / norm(h_b);
        R_b = log2(1 + P_tot * abs(h_b' * w)^2 / noise_var);
        R_e = log2(1 + P_tot * abs(h_e' * w)^2 / noise_var);
        R_samples(it) = secrecy_rate(R_b, R_e);
    end
    mean_norm_h(n_idx) = mean(g_samples);
    var_norm_h(n_idx)  = var(g_samples);
    mean_R_sec(n_idx)  = mean(R_samples);
    std_R_sec(n_idx)   = std(R_samples);

    idx = find(Nt_hist == Nt, 1);
    if ~isempty(idx)
        hist_data{idx} = g_samples;
    end
end

% Theoretical reference: Var = 1/Nt
var_theory = 1 ./ Nt_vec;

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

% Top-left: hardening histograms
subplot(2, 2, 1);
colors = lines(length(Nt_hist));
edges  = 0:0.05:3.5;
for i = 1:length(Nt_hist)
    histogram(hist_data{i}, edges, 'Normalization', 'pdf', ...
        'FaceColor', colors(i,:), 'FaceAlpha', 0.55, ...
        'DisplayName', sprintf('Nt = %d', Nt_hist(i))); hold on;
end
xline(1, 'k--', 'E[||h||^2/N_t] = 1', 'LineWidth', 1.2);
grid on; box on;
xlabel('||h||^2 / N_t'); ylabel('PDF');
title('Channel hardening: distribution of normalised gain');
legend('Location', 'NorthEast');

% Top-right: Var(||h||^2/Nt) vs Nt with theory line
subplot(2, 2, 2);
loglog(Nt_vec, var_norm_h, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
loglog(Nt_vec, var_theory, 'k--', 'LineWidth', 1.5);
grid on; box on;
xlabel('Number of antennas N_t'); ylabel('Var(||h||^2 / N_t)');
title('Hardening rate: Var \propto 1/N_t');
legend('Monte-Carlo', 'Theory: 1/N_t', 'Location', 'NorthEast');

% Bottom-left: mean Secrecy Rate vs Nt
subplot(2, 2, 3);
semilogx(Nt_vec, mean_R_sec, '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
grid on; box on;
xlabel('Number of antennas N_t'); ylabel('E[R_{sec}] (bits/s/Hz)');
title(sprintf('Mean Secrecy Rate vs N_t  (SNR = %d dB)', SNR_dB));

% Bottom-right: std of Secrecy Rate vs Nt - operational hardening
subplot(2, 2, 4);
semilogx(Nt_vec, std_R_sec, '-mo', 'LineWidth', 2, 'MarkerFaceColor', 'm');
grid on; box on;
xlabel('Number of antennas N_t'); ylabel('Std(R_{sec}) (bits/s/Hz)');
title('Outage sensitivity collapses with N_t');

sgtitle('Massive MIMO channel hardening and its PLS consequence');

save_figure(fig, 'fig_channel_hardening');
