% =========================================================================
% SCENARIO: Colluding eavesdroppers (6 GHz vs 28 GHz)
% -------------------------------------------------------------------------
% nrCDLChannel + FSPL + Transmit SNR.
% BS: ZF across K Bobs. Eve coalition: cooperative MRC (R = I_L).
% Sweep: number of cooperating eavesdroppers L.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

dist      = 50;
L_values  = 1:2:15;
K         = 4;
numIter   = 10;
SNR_tx_dB = 100;
P_tx      = 10^(SNR_tx_dB / 10);
noise_var = 1;

bands = struct( ...
    'name', {'6 GHz', '28 GHz'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave});

results_SR   = zeros(2, length(L_values));
results_Fair = zeros(2, length(L_values));
max_L = max(L_values);

for b = 1:2
    fc = bands(b).fc;
    Nt = bands(b).Nt;
    [PL_lin, PL_dB] = compute_fspl(dist, fc);

    fprintf('\n--- Colluding eavesdroppers @ %s ---\n', bands(b).name);
    fprintf('  Transmit SNR: %d dB\n', SNR_tx_dB);
    fprintf('  d = %g m (FSPL: %.2f dB)\n', dist, PL_dB);

    cdl_b = cell(K, 1);
    for k = 1:K
        cdl_b{k} = nrCDLChannel;
    end
    cdl_e = cell(max_L, 1);
    for e = 1:max_L
        cdl_e{e} = nrCDLChannel;
    end

    for l_idx = 1:length(L_values)
        num_eve = L_values(l_idx);
        SR_acc = 0;
        F_acc  = 0;

        for it = 1:numIter
            theta_bobs = -60 + 120 * rand(1, K);
            theta_eves = -60 + 120 * rand(1, num_eve);

            H_eff = zeros(Nt, K);
            for k = 1:K
                cdl_b{k} = setup_matlab_cdl(cdl_b{k}, Nt, fc, theta_bobs(k));
                release(cdl_b{k});
                cdl_b{k}.Seed = randi([0 2^31 - 1]);
                [pg_b, ~] = cdl_b{k}();
                hb = squeeze(sum(pg_b, 2));
                hb = hb(:);
                H_eff(:, k) = sqrt((1 / PL_lin) * Nt) * (hb / norm(hb));
            end

            G_eff = zeros(Nt, num_eve);
            for e = 1:num_eve
                cdl_e{e} = setup_matlab_cdl(cdl_e{e}, Nt, fc, theta_eves(e));
                release(cdl_e{e});
                cdl_e{e}.Seed = randi([0 2^31 - 1]);
                [pg_e, ~] = cdl_e{e}();
                he = squeeze(sum(pg_e, 2));
                he = he(:);
                G_eff(:, e) = sqrt((1 / PL_lin) * Nt) * (he / norm(he));
            end

            W_raw = H_eff * pinv(H_eff' * H_eff + 1e-9 * eye(K));
            W = W_raw / norm(W_raw, 'fro');
            R_inv = eye(num_eve);

            R_b = zeros(K, 1);
            R_e = zeros(K, 1);
            for k = 1:K
                sig_b = P_tx * abs(H_eff(:, k)' * W(:, k))^2;
                intf_b = 0;
                for j = 1:K
                    if j ~= k
                        intf_b = intf_b + P_tx * abs(H_eff(:, k)' * W(:, j))^2;
                    end
                end
                R_b(k) = log2(1 + sig_b / (intf_b + noise_var));

                h_ek = sqrt(P_tx) * G_eff' * W(:, k);
                R_e(k) = log2(1 + real(h_ek' * R_inv * h_ek));
            end

            R_s = max(0, R_b - R_e);
            SR_acc = SR_acc + sum(R_s);
            rs_sum = sum(R_s);
            if rs_sum > 0
                F_acc = F_acc + (rs_sum^2) / (K * sum(R_s.^2));
            end
        end

        results_SR(b, l_idx)   = SR_acc / numIter;
        results_Fair(b, l_idx) = F_acc / numIter;
    end
end

fig = figure('Color', 'w', 'Position', [100 100 1100 450]);

subplot(1, 2, 1);
plot(L_values, results_SR(1, :), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;
plot(L_values, results_SR(2, :), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on;
box on;
xlabel('Number of colluding eavesdroppers (L)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Secrecy vs coalition size');
legend('6 GHz (N_t = 32)', '28 GHz (N_t = 512)', 'Location', 'NorthEast');

subplot(1, 2, 2);
plot(L_values, results_Fair(1, :), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;
plot(L_values, results_Fair(2, :), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on;
box on;
ylim([0 1.05]);
xlabel('Number of colluding eavesdroppers (L)');
ylabel("Jain's fairness index (secrecy rates)");
title('Fairness of per-user secrecy');
legend('6 GHz', '28 GHz', 'Location', 'SouthWest');

sgtitle(sprintf('Colluding eavesdroppers (nrCDL + FSPL, K = %d, SNR_{tx} = %d dB)', K, SNR_tx_dB));
save_figure(fig, 'fig_colluding_eavesdroppers');
plot_colluding_topology(dist, K, max_L);


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
    cdl.NumTimeSamples = 1;
    cdl.ChannelFiltering = false;
end

% =========================================================================
% FUNKCJA POMOCNICZA: Generowanie topologii scenariusza (Zmowa Ew) - UPROSZCZONA
% =========================================================================
function plot_colluding_topology(dist, K, L)
    fig_top = figure('Color', 'w', 'Position', [150 150 700 700]);
    hold on; grid on; box on;
    
    % Losowanie przykładowych kątów w sektorze [-60, 60] stopni
    theta_bobs = -60 + 120 * rand(1, K);
    theta_eves = -60 + 120 * rand(1, L);
    
    % Konwersja na kartezjańskie (BS w 0,0, oś Y to broadside 0 st.)
    x_bs = 0; y_bs = 0;
    max_d = dist + 15;
    
    % Rysowanie BS
    p_bs = plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    text(x_bs, y_bs - 3, 'BS (0,0)', 'HorizontalAlignment', 'center', 'Color', 'k');
    
    % Rysowanie zmowy Ew (Colluding Eves)
    x_eves = dist * sind(theta_eves);
    y_eves = dist * cosd(theta_eves);
    
    p_e = [];
    for e = 1:L
        p_e = plot(x_eves(e), y_eves(e), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        text(x_eves(e) + 1.5, y_eves(e) + 1.5, sprintf('E_{%d}', e), 'Color', 'r', 'FontSize', 9);
    end
    if ~isempty(p_e)
        set(p_e, 'DisplayName', sprintf('Colluding Eves (L=%d)', L));
    end
    
    % Rysowanie Bobów
    p_b = [];
    for k = 1:K
        x_b = dist * sind(theta_bobs(k));
        y_b = dist * cosd(theta_bobs(k));
        p_b = plot(x_b, y_b, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        text(x_b - 1.5, y_b - 1.5, sprintf('B_{%d}', k), 'Color', 'b', 'FontSize', 9, 'HorizontalAlignment', 'right');
    end
    if ~isempty(p_b)
        set(p_b, 'DisplayName', sprintf('Bobs (K=%d)', K));
    end
    
    % Ustawienia osi
    axis equal;
    xlim([-max_d, max_d]);
    ylim([-10, max_d]);
    xlabel('X [m]'); ylabel('Y [m]');
    title(sprintf('Scenario 3: Colluding eavesdroppers (K=%d, L=%d, d=%gm)', K, L, dist));
    legend([p_bs, p_b, p_e], 'Location', 'NorthWest');
    
    save_topology(fig_top, 'topology_colluding_eves');
end