% =========================================================================
% SCENARIO: 6 GHz vs 28 GHz under colluding eavesdroppers (worst-case MRC)
% -------------------------------------------------------------------------
% L cooperating Eves, 3GPP CDL channels, received SNR at Bob after FSPL.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

dist        = p.link_dist_m;
L_values    = 1:2:15;
K           = 4;
numIter     = 80;
SNR_rx_dB   = p.SNR_rx_dB;
noise_var   = p.noise_var;

bands = struct( ...
    'name', {'6 GHz (Massive MIMO)', '28 GHz (Ultra-Massive MIMO)'}, ...
    'fc',   {p.fc_sub6,              p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6,              p.Nt_mmwave}, ...
    'cdl',  {p.cdl_sub6,             p.cdl_mmwave});

results_SR    = zeros(2, length(L_values));
results_Fair  = zeros(2, length(L_values));

for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;  cdl = bands(b).cdl;
    PL_lin = compute_fspl(dist, fc);
    P_rx   = rx_snr_power('linear', SNR_rx_dB);
    print_scenario_snr('title', sprintf('Colluding Eves @ %s', bands(b).name), ...
        'SNR_rx_dB', SNR_rx_dB, 'dist_m', dist, 'fc_Hz', fc, ...
        'actors', {sprintf('Bob (K=%d)', K), sprintf('Eve (L eaves, colluding)')});
    [~, sv] = setup_ula(Nt, fc);

    for l_idx = 1:length(L_values)
        num_eve = L_values(l_idx);
        SR_acc = 0; F_acc = 0;
        for it = 1:numIter
            theta_bobs = -60 + 120*rand(1, K);
            theta_eves = -90 + 180*rand(1, num_eve);

            H = zeros(Nt, K);
            for k = 1:K
                H(:, k) = channel_3gpp_ula(sv, fc, theta_bobs(k), cdl);
            end
            G = zeros(Nt, num_eve);
            for e = 1:num_eve
                G(:, e) = channel_3gpp_ula(sv, fc, theta_eves(e), cdl);
            end

            W = H * pinv(H' * H);
            W = W / norm(W, 'fro');

            R_b = zeros(K, 1);
            R_e = zeros(K, 1);
            for k = 1:K
                R_b(k) = log2(1 + P_rx * abs(H(:,k)' * W(:,k))^2 / noise_var);
                snr_eve_k = P_rx * sum(abs(G' * W(:,k)).^2);
                R_e(k)    = log2(1 + snr_eve_k / noise_var);
            end

            R_s = secrecy_rate(R_b, R_e);
            SR_acc = SR_acc + sum(R_s);
            F_acc  = F_acc  + jains_fairness(R_s);
        end
        results_SR(b, l_idx)   = SR_acc / numIter;
        results_Fair(b, l_idx) = F_acc  / numIter;
    end
end

fig = figure('Color', 'w', 'Position', [100 100 1100 420]);

subplot(1, 2, 1);
plot(L_values, results_SR(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, results_SR(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Number of colluding eavesdroppers (L)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Security under MRC attack');
legend(bands(1).name, bands(2).name, 'Location', 'SouthWest');

subplot(1, 2, 2);
plot(L_values, results_Fair(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, results_Fair(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
ylim([0 1.05]);
xlabel('Number of colluding eavesdroppers (L)');
ylabel("Jain's fairness index");
title('Fairness across legitimate users');
legend(bands(1).name, bands(2).name, 'Location', 'SouthWest');

sgtitle(sprintf('6 GHz vs 28 GHz (3GPP CDL)  (K = %d, d = %d m, received SNR = %d dB)', ...
    K, dist, SNR_rx_dB));

save_figure(fig, 'fig_colluding_eavesdroppers');
