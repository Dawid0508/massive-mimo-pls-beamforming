% =========================================================================
% SCENARIO: Colluding eavesdroppers (6 GHz vs 28 GHz)
% -------------------------------------------------------------------------
% ZAWARTY FIX: Podejście "Rosnącej Koalicji" (Wektorowa ocena L na tym samym
% seedzie i kanale), optymalizacja prędkości prekodera ZF oraz 3GPP Path Loss.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- PARAMETRY SYSTEMU ---
dist      = 50;
L_values  = 1:2:15;
K         = 4;
numIter   = 100;          % Dzięki wektoryzacji pętli 100 iteracji wykona się błyskawicznie
SNR_tx_dB = 100;          % Traktowane jako rho_tx = P_tx / N_0 (Transmit SNR)
P_tx      = 10^(SNR_tx_dB / 10);
noise_var = 1;            % Znormalizowany szum tła

bands = struct( ...
    'name', {'6 GHz', '28 GHz'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave});

results_SR   = zeros(2, length(L_values));
results_Fair = zeros(2, length(L_values));
Rx_SNR_dB    = zeros(2, 1); 
max_L = max(L_values);

for b = 1:2
    fc = bands(b).fc;
    Nt = bands(b).Nt;
    
    % Wyliczenie tłumienia z modelu 3GPP NR Path Loss
    [PL_lin, PL_dB] = compute_nr_pathloss(dist, fc); 
    
    % Wyliczenie średniego docelowego Received SNR
    Rx_SNR_dB(b) = 10 * log10(P_tx * Nt / PL_lin);
    
    fprintf('\n--- Colluding eavesdroppers @ %s ---\n', bands(b).name);
    fprintf('  Transmit SNR (rho_tx): %d dB\n', SNR_tx_dB);
    fprintf('  d = %g m (3GPP PL: %.2f dB)\n', dist, PL_dB);
    fprintf('  Average Rx SNR: %.1f dB\n', Rx_SNR_dB(b));
    
    % Alokacja obiektów 3GPP przed pętlami
    cdl_b = cell(K, 1);
    for k = 1:K, cdl_b{k} = create_base_cdl(Nt, fc); end
    cdl_e = cell(max_L, 1);
    for e = 1:max_L, cdl_e{e} = create_base_cdl(Nt, fc); end
    
    % Akumulatory dla całego sweepu w danej częstotliwości
    SR_sweep_acc = zeros(1, length(L_values));
    F_sweep_acc  = zeros(1, length(L_values));
    
    % PĘTLA GŁÓWNA MONTE CARLO (Generowanie "świata" raz na iterację)
    for it = 1:numIter
        iter_seed = randi([0 2^31 - 1]);
        
        % Geometria sektora
        theta_bobs = -60 + 120 * rand(1, K);
        theta_eves_max = -60 + 120 * rand(1, max_L); % Losujemy pozycje od razu dla max_L Ew
        
        % --- 1. KANAŁY LEGALNE (BOBOWIE) ---
        H_eff = zeros(Nt, K);
        for k = 1:K
            release(cdl_b{k});
            cdl_b{k}.TransmitArrayOrientation = [-theta_bobs(k); 0; 0];
            cdl_b{k}.Seed = iter_seed;
            [pg_b, ~] = cdl_b{k}();
            hb = squeeze(sum(pg_b, 2)); hb = hb(:);
            H_eff(:, k) = sqrt(1 / PL_lin) * hb; 
        end
        
        % --- 2. KANAŁY PODSŁUCHIWACZY (MAKSYMALNA KOALICJA EW) ---
        G_eff_max = zeros(Nt, max_L);
        for e = 1:max_L
            release(cdl_e{e});
            cdl_e{e}.TransmitArrayOrientation = [-theta_eves_max(e); 0; 0];
            cdl_e{e}.Seed = iter_seed;
            [pg_e, ~] = cdl_e{e}();
            he = squeeze(sum(pg_e, 2)); he = he(:);
            G_eff_max(:, e) = sqrt(1 / PL_lin) * he;
        end
        
        % --- 3. PREKODOWANIE ZF (Wyliczane RAZ na całą iterację!) ---
        W_raw = H_eff * pinv(H_eff' * H_eff + 1e-9 * eye(K));
        W = W_raw / norm(W_raw, 'fro'); 
        
        % --- 4. PRZEPUSTOWOŚĆ BOBÓW (Stała niezależnie od liczby Ew) ---
        R_b = zeros(K, 1);
        S_b_base = zeros(K, 1);
        I_b_base = zeros(K, 1);
        for k = 1:K
            S_b_base(k) = P_tx * abs(H_eff(:, k)' * W(:, k))^2;
            for j = 1:K
                if j ~= k, I_b_base(k) = I_b_base(k) + P_tx * abs(H_eff(:, k)' * W(:, j))^2; end
            end
            R_b(k) = log2(1 + S_b_base(k) / (I_b_base(k) + noise_var));
        end
        
        % --- 5. PĘTLA WEWNĘTRZNA: Testowanie rozmiaru koalicji Ew (L) ---
        for l_idx = 1:length(L_values)
            num_eve = L_values(l_idx);
            
            % ROSNĄCA KOALICJA: Wycinamy podmacierz kanałów dla aktualnego L
            G_eff = G_eff_max(:, 1:num_eve);
            R_inv = eye(num_eve);
            R_e = zeros(K, 1);
            
            for k = 1:K
                % Wyliczenie kooperatywnego odbioru MRC dla aktualnej liczby Ew
                h_ek = sqrt(P_tx) * G_eff' * W(:, k);
                R_e(k) = log2(1 + real(h_ek' * R_inv * h_ek));
            end
            
            % Obliczenie metryk PLS dla danego L na zamrożonym kanale
            R_s = max(0, R_b - R_e);
            rs_sum = sum(R_s);
            
            SR_sweep_acc(l_idx) = SR_sweep_acc(l_idx) + rs_sum;
            if sum(R_s.^2) > 0
                F_sweep_acc(l_idx) = F_sweep_acc(l_idx) + (rs_sum^2) / (K * sum(R_s.^2));
            end
        end
    end
    
    % Zapisanie uśrednionych wyników całej krzywej sweepu
    results_SR(b, :)   = SR_sweep_acc / numIter;
    results_Fair(b, :) = F_sweep_acc / numIter;
end

% =========================================================================
% --- WIZUALIZACJA WYNIKÓW ---
% =========================================================================
legend_6GHz  = sprintf('6 GHz (N_t = %d, Rx SNR \\approx %.1f dB)', bands(1).Nt, Rx_SNR_dB(1));
legend_28GHz = sprintf('28 GHz (N_t = %d, Rx SNR \\approx %.1f dB)', bands(2).Nt, Rx_SNR_dB(2));

fig = figure('Color', 'w', 'Position', [100 100 1100 450]);

subplot(1, 2, 1);
plot(L_values, results_SR(1, :), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, results_SR(2, :), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlim([min(L_values) max(L_values)]); % Usunięcie przerw po bokach
xlabel('Number of colluding eavesdroppers (L)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Secrecy vs coalition size');
legend(legend_6GHz, legend_28GHz, 'Location', 'NorthEast');

subplot(1, 2, 2);
plot(L_values, results_Fair(1, :), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, results_Fair(2, :), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlim([min(L_values) max(L_values)]); % Usunięcie przerw po bokach
xlabel('Number of colluding eavesdroppers (L)'); ylabel("Jain's fairness index");
title('Fairness of per-user secrecy');
legend(legend_6GHz, legend_28GHz, 'Location', 'SouthWest');

sgtitle('Scenario 3: Colluding Eavesdroppers Trade-off');
try save_figure(fig, 'fig_colluding_eavesdroppers'); catch; end
plot_colluding_topology(dist, K, max_L);

% =========================================================================
% FUNKCJE POMOCNICZE
% =========================================================================
function cdl = create_base_cdl(Nt, fc)
    cdl = nrCDLChannel;
    cdl.DelayProfile = 'CDL-A';
    if fc < 10e9, cdl.DelaySpread = 30e-9; else, cdl.DelaySpread = 10e-9; end
    cdl.CarrierFrequency = fc;
    cdl.MaximumDopplerShift = 0;
    cdl.TransmitAntennaArray.Size = [1 Nt 1 1 1];
    cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1];
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
    cdl.NumTimeSamples = 1;
    cdl.ChannelFiltering = false;
end

function plot_colluding_topology(dist, K, L)
    fig_top = figure('Color', 'w', 'Position', [150 150 700 700]);
    hold on; grid on; box on;
    theta_bobs = -60 + 120 * rand(1, K);
    theta_eves = -60 + 120 * rand(1, L);
    
    x_bs = 0; y_bs = 0;
    max_d = dist + 15;
    
    p_bs = plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    text(x_bs, y_bs - 3, 'BS (0,0)', 'HorizontalAlignment', 'center', 'Color', 'k');
    
    x_eves = dist * sind(theta_eves); y_eves = dist * cosd(theta_eves);
    p_e = [];
    for e = 1:L
        p_e = plot(x_eves(e), y_eves(e), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        text(x_eves(e) + 1.5, y_eves(e) + 1.5, sprintf('E_{%d}', e), 'Color', 'r', 'FontSize', 9);
    end
    if ~isempty(p_e), set(p_e, 'DisplayName', sprintf('Colluding Eves (L=%d)', L)); end
    
    x_b = dist * sind(theta_bobs); y_b = dist * cosd(theta_bobs);
    p_b = [];
    for k = 1:K
        p_b = plot(x_b(k), y_b(k), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        text(x_b(k) - 1.5, y_b(k) - 1.5, sprintf('B_{%d}', k), 'Color', 'b', 'FontSize', 9, 'HorizontalAlignment', 'right');
    end
    if ~isempty(p_b), set(p_b, 'DisplayName', sprintf('Bobs (K=%d)', K)); end
    
    axis equal; xlim([-max_d, max_d]); ylim([-10, max_d]);
    xlabel('X [m]'); ylabel('Y [m]');
    title(sprintf('Scenario 3: Colluding eavesdroppers'));
    legend([p_bs, p_b, p_e], 'Location', 'NorthWest');
    try save_figure(fig_top, '../topology/topology_colluding_eves'); catch; end
end