% =========================================================================
% SCENARIO: Low-Resolution DAC quantisation in Massive MIMO PLS
% -------------------------------------------------------------------------
% Reference: Y. Xu et al., "Secure Massive MIMO Communication With
%            Low-Resolution DACs", IEEE TWC 2018 (doc/8626548).
%            Scenariusze_pomysly.docx (Low-Resolution DACs scenario).
%
% Massive arrays burn power on per-antenna RF chains, so 6G hardware
% trends towards 1-3 bit DACs per branch. Quantisation distortion can
% be modelled via Bussgang's theorem:
%
%       x_q = a * x + n_q ,    Cov(n_q) = (1 - a^2) * P_x * I
%
% with n_q approximately uncorrelated with x (Gaussian input). The
% quantisation noise is *isotropic*, so it leaks energy in every
% direction, including towards Eve. The interesting question: does
% the relative degradation hurt Bob or Eve more?
%
% This script uses the *direct* uniform quantiser (utils/uniform_quantize)
% and reports Bob, Eve, and Secrecy Sum-Rates for b = {1, 2, 3, 4, Inf}.
% =========================================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
b_vec       = [1 2 3 4 5 Inf];               % DAC bits per I/Q branch
SNR_dB_vec  = 0:5:30;                        % SNR sweep for the second panel
SNR_fixed   = 20;                            % dB, used for the b sweep
Nt          = 64;
K           = 4;
numIter     = 80;
numSym      = 64;                            % symbols per Monte-Carlo run
noise_var   = p.noise_var;

% Result storage at SNR_fixed
R_bob   = zeros(1, length(b_vec));
R_eve   = zeros(1, length(b_vec));
R_sec   = zeros(1, length(b_vec));

% Result storage for SNR sweep at b in {1, 2, Inf}
b_show  = [1, 2, Inf];
R_sec_snr = zeros(length(b_show), length(SNR_dB_vec));

% Empirical Bussgang factor a(b) - useful auxiliary curve
a_emp   = zeros(size(b_vec));

% --- Sweep A: Secrecy decomposition vs DAC bits --------------------------
P_tot = 10^(SNR_fixed/10);
for bi = 1:length(b_vec)
    b = b_vec(bi);
    acc_b = 0; acc_e = 0; acc_s = 0;
    a_acc = 0; a_n = 0;

    for it = 1:numIter
        H     = (randn(Nt, K) + 1j*randn(Nt, K)) / sqrt(2);
        h_eve = (randn(Nt, 1) + 1j*randn(Nt, 1)) / sqrt(2);

        W_raw = H * pinv(H' * H);
        W = W_raw / norm(W_raw, 'fro') * sqrt(P_tot);

        % i.i.d. Gaussian symbols
        s   = (randn(K, numSym) + 1j*randn(K, numSym)) / sqrt(2);
        x   = W * s;                              % Nt x numSym before DAC
        x_q = uniform_quantize(x, b);             % per-antenna quantisation

        % Empirical Bussgang factor: a = E[<x,x_q>]/E[|x|^2]
        a_acc = a_acc + real(x(:)' * x_q(:)) / (x(:)' * x(:) + eps);
        a_n   = a_n + 1;

        % Received samples (Bob and Eve)
        y_b = H' * x_q + sqrt(noise_var/2) * (randn(K, numSym) + 1j*randn(K, numSym));
        y_e = h_eve' * x_q + sqrt(noise_var/2) * (randn(1, numSym) + 1j*randn(1, numSym));

        % Per-user empirical SINR via signal/(residual) decomposition
        % Bob k:   useful = h_k^H W e_k s_k   ;   anything else = noise+leakage
        for k = 1:K
            % Useful and aggregate signals through the *quantised* x
            % estimate via least-squares fit on the symbol sequence
            % (works because s is i.i.d. and unit-variance):
            useful_b = mean( real( y_b(k, :) .* conj(s(k, :)) ) );
            % Total received power minus the estimated useful component
            tot_b    = mean( abs(y_b(k, :)).^2 );
            sig_b    = useful_b^2;
            inr_b    = max(tot_b - sig_b, eps);
            R_b_k    = log2(1 + sig_b / inr_b);

            % Eve treats user-k symbol as the desired one (worst case)
            useful_e = mean( real( y_e .* conj(s(k, :)) ) );
            tot_e    = mean( abs(y_e).^2 );
            sig_e    = useful_e^2;
            inr_e    = max(tot_e - sig_e, eps);
            R_e_k    = log2(1 + sig_e / inr_e);

            acc_b = acc_b + R_b_k;
            acc_e = acc_e + R_e_k;
            acc_s = acc_s + secrecy_rate(R_b_k, R_e_k);
        end
    end

    R_bob(bi) = acc_b / numIter;
    R_eve(bi) = acc_e / numIter;
    R_sec(bi) = acc_s / numIter;
    a_emp(bi) = a_acc / max(a_n, 1);
end

% --- Sweep B: SNR vs Secrecy at three DAC resolutions --------------------
for ii = 1:length(b_show)
    b = b_show(ii);
    for s_idx = 1:length(SNR_dB_vec)
        P_tot = 10^(SNR_dB_vec(s_idx)/10);
        acc_s = 0;
        for it = 1:numIter
            H     = (randn(Nt, K) + 1j*randn(Nt, K)) / sqrt(2);
            h_eve = (randn(Nt, 1) + 1j*randn(Nt, 1)) / sqrt(2);
            W_raw = H * pinv(H' * H);
            W     = W_raw / norm(W_raw, 'fro') * sqrt(P_tot);

            s   = (randn(K, numSym) + 1j*randn(K, numSym)) / sqrt(2);
            x   = W * s;
            x_q = uniform_quantize(x, b);
            y_b = H' * x_q + sqrt(noise_var/2) * (randn(K, numSym) + 1j*randn(K, numSym));
            y_e = h_eve' * x_q + sqrt(noise_var/2) * (randn(1, numSym) + 1j*randn(1, numSym));

            for k = 1:K
                useful_b = mean(real(y_b(k,:) .* conj(s(k,:))));
                tot_b    = mean(abs(y_b(k,:)).^2);
                R_b_k    = log2(1 + useful_b^2 / max(tot_b - useful_b^2, eps));
                useful_e = mean(real(y_e .* conj(s(k,:))));
                tot_e    = mean(abs(y_e).^2);
                R_e_k    = log2(1 + useful_e^2 / max(tot_e - useful_e^2, eps));
                acc_s    = acc_s + secrecy_rate(R_b_k, R_e_k);
            end
        end
        R_sec_snr(ii, s_idx) = acc_s / numIter;
    end
end

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

bar_x = 1:length(b_vec);
bar_labels = arrayfun(@(b) ternary(isinf(b), 'Inf', sprintf('%d', b)), b_vec, 'UniformOutput', false);

% Top-left: Bob vs Eve sum-rate per b
subplot(2, 2, 1);
bar(bar_x, [R_bob(:) R_eve(:)], 1.0); grid on; box on;
set(gca, 'XTick', bar_x, 'XTickLabel', bar_labels);
xlabel('DAC resolution b (bits / I-or-Q branch)');
ylabel('Sum-Rate (bits/s/Hz)');
title(sprintf('Bob vs Eve under quantisation  (SNR = %d dB)', SNR_fixed));
legend('Bob (legitimate)', 'Eve (eavesdropper)', 'Location', 'NorthWest');

% Top-right: Secrecy Sum-Rate vs b
subplot(2, 2, 2);
bar(bar_x, R_sec, 0.6, 'FaceColor', [0.2 0.6 0.3]);
grid on; box on;
set(gca, 'XTick', bar_x, 'XTickLabel', bar_labels);
xlabel('DAC resolution b (bits / branch)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Secrecy Sum-Rate vs DAC bits');

% Bottom-left: Empirical Bussgang factor a(b)
subplot(2, 2, 3);
b_finite = b_vec(~isinf(b_vec));
a_finite = a_emp(~isinf(b_vec));
bar(b_finite, a_finite, 0.6, 'FaceColor', [0.6 0.4 0.8]); hold on;
yline(1, 'k--', 'Ideal DAC: a = 1', 'LineWidth', 1.2);
grid on; box on;
xlabel('DAC bits b'); ylabel('Bussgang gain a');
title('Empirical Bussgang factor');
ylim([0 1.05]);

% Bottom-right: Secrecy Sum-Rate vs SNR for selected b
subplot(2, 2, 4);
markers = {'-bo', '-ms', '-k^'};
hold on;
for ii = 1:length(b_show)
    plot(SNR_dB_vec, R_sec_snr(ii,:), markers{ii}, 'LineWidth', 2, 'MarkerFaceColor', markers{ii}(2));
end
grid on; box on;
xlabel('Transmit SNR (dB)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Secrecy Sum-Rate vs SNR per DAC resolution');
labels = arrayfun(@(b) ternary(isinf(b), 'b = Inf (ideal)', sprintf('b = %d bit', b)), b_show, 'UniformOutput', false);
legend(labels, 'Location', 'NorthWest');

sgtitle(sprintf('Low-resolution DAC quantisation in Massive MIMO PLS  (Nt = %d, K = %d)', Nt, K));

save_figure(fig, 'fig_low_res_dac');


% --- tiny inline ternary helper for clean labels -------------------------
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
