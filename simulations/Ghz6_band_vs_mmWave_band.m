% =========================================================================
% SCENARIO: Artificial Noise (AN) injection in the null space of H'
% -------------------------------------------------------------------------
% Goh & Hong (2011) showed that splitting BS power between data
% (fraction phi) and isotropic noise in the null space of H' (fraction
% 1-phi) leaks zero AN to the legitimate users while jamming any
% off-subspace eavesdropper. We replicate the effect at 6 GHz under a
% colluding multi-Eve attack and (a) sweep L = #Eves at fixed phi, (b)
% sweep phi at fixed L to expose the data/noise trade-off optimum.
% =========================================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
dist        = 30;                     % Alice <-> users distance [m]
L_values    = 1:2:15;                 % # colluding eavesdroppers (sweep A)
phi_values  = 0.1:0.1:1.0;            % data-power fraction (sweep B)
L_fixed     = 7;                      % # Eves used in sweep B
phi_fixed   = 0.7;                    % phi used in sweep A
K           = 4;                      % legitimate users
numIter     = 80;
noise_var   = p.noise_var;

fc = p.fc_sub6;  Nt = p.Nt_sub6;
PL_lin = compute_fspl(dist, fc);

% Operating point: 30 dB *received* SNR after path loss
P_tx = 10^(30/10) * PL_lin;            % transmit power so that P_tx/PL = 30 dB

[ula, sv] = setup_ula(Nt, fc); %#ok<ASGLU>

% --- Sweep A: number of Eves at phi = phi_fixed -------------------------
SR_no_AN_A   = zeros(1, length(L_values));
SR_with_AN_A = zeros(1, length(L_values));

for l_idx = 1:length(L_values)
    num_eve = L_values(l_idx);
    acc_no = 0; acc_an = 0;
    for it = 1:numIter
        [SR_no, SR_an] = run_AN_trial(sv, fc, Nt, K, num_eve, ...
                                      P_tx, PL_lin, phi_fixed, noise_var);
        acc_no = acc_no + SR_no;
        acc_an = acc_an + SR_an;
    end
    SR_no_AN_A(l_idx)   = acc_no / numIter;
    SR_with_AN_A(l_idx) = acc_an / numIter;
end

% --- Sweep B: phi sweep at L = L_fixed ----------------------------------
SR_with_AN_B = zeros(1, length(phi_values));
for ph_idx = 1:length(phi_values)
    acc = 0;
    for it = 1:numIter
        [~, SR_an] = run_AN_trial(sv, fc, Nt, K, L_fixed, ...
                                  P_tx, PL_lin, phi_values(ph_idx), noise_var);
        acc = acc + SR_an;
    end
    SR_with_AN_B(ph_idx) = acc / numIter;
end

% --- Visualisation -------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1100 420]);

subplot(1, 2, 1);
plot(L_values, SR_no_AN_A,   '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, SR_with_AN_A, '-g^', 'LineWidth', 2, 'MarkerFaceColor', 'g');
grid on; box on;
xlabel('Number of colluding eavesdroppers (L)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Robustness vs L  (\\phi = %.1f)', phi_fixed));
legend('Plain ZF', 'ZF + AN', 'Location', 'NorthEast');

subplot(1, 2, 2);
plot(phi_values, SR_with_AN_B, '-mh', 'LineWidth', 2, 'MarkerFaceColor', 'm');
grid on; box on;
xlabel('Data-power fraction \phi');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Power-allocation trade-off  (L = %d)', L_fixed));
xline(phi_fixed, 'k:', sprintf('\\phi = %.1f', phi_fixed));

sgtitle(sprintf('Artificial Noise at 6 GHz (Nt = %d, K = %d, d = %d m)', ...
    Nt, K, dist));

save_figure(fig, 'fig_artificial_noise');


% =========================================================================
%                          Local helper
% =========================================================================
function [SR_no_AN, SR_with_AN] = run_AN_trial(sv, fc, Nt, K, num_eve, ...
                                               P_tx, PL_lin, phi, noise_var)
    % Random user / eavesdropper geometries (uniform in azimuth)
    theta_bobs = -60 + 120*rand(1, K);
    theta_eves = -90 + 180*rand(1, num_eve);

    H = step(sv, fc, theta_bobs);              % Nt x K
    G = step(sv, fc, theta_eves);              % Nt x num_eve

    % Standard ZF precoder (Frobenius-normalised)
    W = H * pinv(H' * H);
    W = W / norm(W, 'fro');

    % Null space of H'   (any vector z s.t. H' z = 0  =>  Bob hears nothing)
    Z = null(H');                              % Nt x m
    m = size(Z, 2);

    % Per-stream effective received power (after path loss)
    P_eff = P_tx / PL_lin;

    R_b_no  = zeros(K, 1); R_b_an = zeros(K, 1);
    R_e_no  = zeros(K, 1); R_e_an = zeros(K, 1);
    for k = 1:K
        % Bob (no AN, full power on data)
        sig   = P_eff * abs(H(:,k)' * W(:,k))^2;
        R_b_no(k) = log2(1 + sig / noise_var);

        % Bob (AN active: data power = phi*P_eff). AN is in null(H'),
        % therefore H(:,k)'*Z = 0  =>  no leakage onto Bob.
        sig_an = phi * P_eff * abs(H(:,k)' * W(:,k))^2;
        R_b_an(k) = log2(1 + sig_an / noise_var);

        % Eve(s): co-located colluding receivers, energy combining
        leak_no = 0; leak_an = 0; an_at_eve = 0;
        for e = 1:num_eve
            leak_no = leak_no + P_eff * abs(G(:,e)' * W(:,k))^2;
            leak_an = leak_an + (phi * P_eff) * abs(G(:,e)' * W(:,k))^2;
            an_at_eve = an_at_eve + ((1-phi) * P_eff / m) * sum(abs(G(:,e)' * Z).^2);
        end
        R_e_no(k) = log2(1 + leak_no / noise_var);
        R_e_an(k) = log2(1 + leak_an / (an_at_eve + noise_var));
    end

    SR_no_AN   = sum(secrecy_rate(R_b_no, R_e_no));
    SR_with_AN = sum(secrecy_rate(R_b_an, R_e_an));
end
