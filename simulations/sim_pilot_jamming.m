% =========================================================================
% SCENARIO: Pilot Jamming attack (denial-of-service on the training phase)
% -------------------------------------------------------------------------
% Zaktualizowano: Tłumienie trasy compute_nr_pathloss (model 3GPP),
% wyliczanie nominalnego Received SNR, udoskonalona wizualizacja mapy ciepła
% oraz wdrożenie Common Random Numbers (CRN) dla stabilności Monte Carlo.
%
% Ewa wysyła losowy szum podczas fazy treningowej (Jammer).
% Błąd estymacji kanału LS zależy od JPR (Jammer-to-Pilot Ratio) oraz
% naturalnie skaluje się z fizycznym tłumieniem trasy 3GPP.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Konfiguracja -------------------------------------------------------
Nt          = 32;                     % Bazowa liczba anten dla Sweep A
K           = 4;
tau         = K;                      % Długość sekwencji treningowej
JPR_dB_vec  = -10:5:30;               % Jammer-to-Pilot Ratio [dB] na odbiorniku stacji!
SNR_tx_dB   = 110;                    % Transmit SNR [dB]
P_tx        = 10^(SNR_tx_dB / 10);
noise_var   = 1;                      % Znormalizowany szum tła
numIter     = 40;                     % Iteracje Monte Carlo dla Sweep A
fc          = p.fc_sub6;              % Pasmo (np. 6 GHz)
dist_b      = 50;                     % Dystans do Bobów [m]
dist_e      = 40;                     % Dystans do Ewy [m] (podsłuch/jammer)
theta_e     = 30;                     % Kąt Ewy

% --- TŁUMIENIE 3GPP ---
[PL_lin_b, PL_dB_b] = compute_nr_pathloss(dist_b, fc);
[PL_lin_e, PL_dB_e] = compute_nr_pathloss(dist_e, fc);

% Wyliczenie nominalnego Received SNR u Bobów (dla stałego P_tx i bazowego Nt)
Rx_SNR_dB_b = 10 * log10(P_tx * Nt / PL_lin_b);

% Moce odbierane na stacji bazowej podczas fazy treningowej
P_rx_pilot = P_tx / PL_lin_b; 

fprintf('\n--- Pilot Jamming (Denial-of-Service) @ %.1f GHz ---\n', fc/1e9);
fprintf('  Dystans Bob: %g m (3GPP PL: %.2f dB)\n', dist_b, PL_dB_b);
fprintf('  Dystans Eve: %g m (3GPP PL: %.2f dB)\n', dist_e, PL_dB_e);
fprintf('  Transmit SNR (rho_tx): %d dB\n', SNR_tx_dB);
fprintf('  Average Rx SNR at Bobs (Nt=%d): %.1f dB\n', Nt, Rx_SNR_dB_b);
fprintf('  Moc odbierana pilota (P_rx_pilot): %g liniowo\n', P_rx_pilot);

% Macierze wynikowe
R_b   = zeros(2, length(JPR_dB_vec));   % wiersze = {MRT, ZF}
R_e   = zeros(2, length(JPR_dB_vec));
R_s   = zeros(2, length(JPR_dB_vec));
R_s_perfect = zeros(2, 1);              % Baseline idealnego CSI

% Parametry Sweep B
Nt_vec      = [16 32 64 128 256];
JPR_grid_dB = -10:5:30;
SR_grid     = zeros(length(Nt_vec), length(JPR_grid_dB));

% --- Wstępna alokacja dla Sweep A ---------------------------------------
theta_bobs = linspace(-60, 60, K);
cdl_b = cell(K, 1);
for k = 1:K
    cdl_b{k} = setup_matlab_cdl(nrCDLChannel, Nt, fc, theta_bobs(k)); 
end
cdl_e = setup_matlab_cdl(nrCDLChannel, Nt, fc, theta_e);

% =========================================================================
% --- SWEEP A: Zmiana JPR dla stałego Nt
% =========================================================================
fprintf('Rozpoczynam Sweep A (JPR przy Nt = %d)...\n', Nt);

% [CRN]: Wektor stałych seedów dla wszystkich iteracji w Sweep A
seeds_A = randi([0 2^31-1], numIter, 1); 

for j_idx = 1:length(JPR_dB_vec)
    % Moc Jammera odbierana przez Stację Bazową
    P_jam_rx = P_rx_pilot * 10^(JPR_dB_vec(j_idx)/10); 
    var_err = (1 / PL_lin_b) * (P_jam_rx + noise_var) / (P_rx_pilot * tau);
    
    acc = zeros(2, 3);  % {MRT,ZF} x {Bob, Eve, Sec}
    for it = 1:numIter
        H_eff = zeros(Nt, K);
        
        iter_seed = seeds_A(it); % Pobieramy seed dla danej iteracji
        
        for k = 1:K
            release(cdl_b{k});
            cdl_b{k}.Seed = iter_seed;
            [pg_b, ~] = cdl_b{k}();
            hb = squeeze(sum(pg_b, 2)); hb = hb(:);
            H_eff(:, k) = sqrt(1 / PL_lin_b) * hb;
        end
        release(cdl_e);
        cdl_e.Seed = iter_seed;
        [pg_e, ~] = cdl_e();
        he = squeeze(sum(pg_e, 2)); he = he(:);
        he_eff = sqrt(1 / PL_lin_e) * he;
        
        % [CRN]: Wymuszamy seed przed losowaniem szumu estymacji
        % Dzięki temu kierunek błędu jest stały dla danej iteracji, zmienia się tylko wariancja.
        rng(iter_seed); 
        E = sqrt(var_err / 2) * (randn(Nt, K) + 1j*randn(Nt, K));
        Hhat = H_eff + E;
        
        % Prekodery tworzone w oparciu o "zatrutą" estymatę
        W_mrt = Hhat / norm(Hhat, 'fro'); 
        W_raw = Hhat * pinv(Hhat' * Hhat);
        W_zf  = W_raw / norm(W_raw, 'fro');
        
        for pidx = 1:2
            if pidx == 1, W = W_mrt; else, W = W_zf; end
            R_b_k = zeros(K, 1); R_e_k = zeros(K, 1);
            for k = 1:K
                sig    = P_tx * abs(H_eff(:,k)' * W(:,k))^2;
                intf   = P_tx * (sum(abs(H_eff(:,k)' * W).^2) - abs(H_eff(:,k)' * W(:,k))^2);
                R_b_k(k) = log2(1 + sig / (intf + noise_var));
                
                sig_e  = P_tx * abs(he_eff' * W(:,k))^2;
                intf_e = P_tx * (sum(abs(he_eff' * W).^2) - abs(he_eff' * W(:,k))^2);
                R_e_k(k) = log2(1 + sig_e / (intf_e + noise_var));
            end
            R_s_k = max(0, R_b_k - R_e_k);
            acc(pidx, 1) = acc(pidx, 1) + sum(R_b_k);
            acc(pidx, 2) = acc(pidx, 2) + sum(R_e_k);
            acc(pidx, 3) = acc(pidx, 3) + sum(R_s_k);
        end
    end
    R_b(:, j_idx) = acc(:, 1) / numIter;
    R_e(:, j_idx) = acc(:, 2) / numIter;
    R_s(:, j_idx) = acc(:, 3) / numIter;
end

% --- Baseline: Idealne CSI (Brak Jammera) ---
acc_perf = zeros(2, 1);
% Odświeżamy stan generatora po ewentualnych zmianach z CRN
rng(p.rng_seed); 
seeds_perf = randi([0 2^31-1], numIter, 1);

for it = 1:numIter
    H_eff = zeros(Nt, K);
    iter_seed = seeds_perf(it);
    
    for k = 1:K
        release(cdl_b{k});
        cdl_b{k}.Seed = iter_seed;
        [pg_b, ~] = cdl_b{k}();
        hb = squeeze(sum(pg_b, 2)); hb = hb(:);
        H_eff(:, k) = sqrt(1 / PL_lin_b) * hb;
    end
    release(cdl_e);
    cdl_e.Seed = iter_seed;
    [pg_e, ~] = cdl_e();
    he = squeeze(sum(pg_e, 2)); he = he(:);
    he_eff = sqrt(1 / PL_lin_e) * he;
    
    W_mrt = H_eff / norm(H_eff, 'fro');
    W_raw = H_eff * pinv(H_eff' * H_eff);
    W_zf  = W_raw / norm(W_raw, 'fro');
    
    for pidx = 1:2
        if pidx == 1, W = W_mrt; else, W = W_zf; end
        Rsk = 0;
        for k = 1:K
            sig    = P_tx * abs(H_eff(:,k)' * W(:,k))^2;
            intf   = P_tx * (sum(abs(H_eff(:,k)' * W).^2) - abs(H_eff(:,k)' * W(:,k))^2);
            R_b_k  = log2(1 + sig / (intf + noise_var));
            
            sig_e  = P_tx * abs(he_eff' * W(:,k))^2;
            intf_e = P_tx * (sum(abs(he_eff' * W).^2) - abs(he_eff' * W(:,k))^2);
            R_e_k  = log2(1 + sig_e / (intf_e + noise_var));
            
            Rsk = Rsk + max(0, R_b_k - R_e_k);
        end
        acc_perf(pidx) = acc_perf(pidx) + Rsk;
    end
end
R_s_perfect = acc_perf / numIter;

% =========================================================================
% --- SWEEP B: Mapa ciepła (Nt vs JPR) - Tylko ZF
% =========================================================================
fprintf('Rozpoczynam Sweep B (Mapa Ciepła Nt vs JPR)...\n');

numIter_B = 40; % Zredukowane MC dla szybkości mapy ciepła
rng(p.rng_seed); 
seeds_B = randi([0 2^31-1], numIter_B, 1); % [CRN] dla mapy ciepła

for n_idx = 1:length(Nt_vec)
    Nt_b = Nt_vec(n_idx);
    
    % Re-konfiguracja modeli dla nowego Nt
    cdl_b_b = cell(K, 1);
    for k = 1:K, cdl_b_b{k} = setup_matlab_cdl(nrCDLChannel, Nt_b, fc, theta_bobs(k)); end
    cdl_e_b = setup_matlab_cdl(nrCDLChannel, Nt_b, fc, theta_e);
    
    for j_idx = 1:length(JPR_grid_dB)
        P_jam_rx = P_rx_pilot * 10^(JPR_grid_dB(j_idx)/10);
        var_err = (1 / PL_lin_b) * (P_jam_rx + noise_var) / (P_rx_pilot * tau);
        
        acc_s = 0;
        for it = 1:numIter_B  
            H_eff = zeros(Nt_b, K);
            iter_seed = seeds_B(it);
            
            for k = 1:K
                release(cdl_b_b{k});
                cdl_b_b{k}.Seed = iter_seed;
                [pg_b, ~] = cdl_b_b{k}();
                hb = squeeze(sum(pg_b, 2)); hb = hb(:);
                H_eff(:, k) = sqrt(1 / PL_lin_b) * hb;
            end
            release(cdl_e_b);
            cdl_e_b.Seed = iter_seed;
            [pg_e, ~] = cdl_e_b();
            he = squeeze(sum(pg_e, 2)); he = he(:);
            he_eff = sqrt(1 / PL_lin_e) * he;
            
            rng(iter_seed); % [CRN] Spójny szum przestrzenny
            E = sqrt(var_err / 2) * (randn(Nt_b, K) + 1j*randn(Nt_b, K));
            Hhat = H_eff + E;
            
            W_raw = Hhat * pinv(Hhat' * Hhat);
            W = W_raw / norm(W_raw, 'fro');
            
            Rsk = 0;
            for k = 1:K
                sig    = P_tx * abs(H_eff(:,k)' * W(:,k))^2;
                intf   = P_tx * (sum(abs(H_eff(:,k)' * W).^2) - abs(H_eff(:,k)' * W(:,k))^2);
                R_b_k  = log2(1 + sig / (intf + noise_var));
                
                sig_e  = P_tx * abs(he_eff' * W(:,k))^2;
                intf_e = P_tx * (sum(abs(he_eff' * W).^2) - abs(he_eff' * W(:,k))^2);
                R_e_k  = log2(1 + sig_e / (intf_e + noise_var));
                
                Rsk = Rsk + max(0, R_b_k - R_e_k);
            end
            acc_s = acc_s + Rsk;
        end
        SR_grid(n_idx, j_idx) = acc_s / numIter_B;
    end
end

plot_jamming_topology(dist_b, theta_bobs, dist_e, theta_e);

% =========================================================================
% --- WIZUALIZACJA --------------------------------------------------------
% =========================================================================
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

% Top-left: Bob, Eve, Secrecy under jamming (ZF)
subplot(2, 2, 1);
plot(JPR_dB_vec, R_b(2,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(JPR_dB_vec, R_e(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
plot(JPR_dB_vec, R_s(2,:), '-g^', 'LineWidth', 2, 'MarkerFaceColor', 'g');
c = pls_colors();
yline(R_s_perfect(2), ':', 'Perfect CSI', 'Color', c.perfect, 'LineWidth', 1.5);
pls_axis_prefs(gca, 'refLabelV', 'bottom', 'staggerRef', true);
grid on; box on;
xlabel('Jammer-to-Pilot Ratio (dB)'); ylabel('Rate (bits/s/Hz)');
title(sprintf('ZF Capacity under Pilot Jamming (Nt = %d)', Nt));
legend('Bob Sum-Rate', 'Eve Sum-Rate', 'Secrecy Sum-Rate', 'Location', 'SouthWest');

% Top-right: MRT vs ZF
subplot(2, 2, 2);
plot(JPR_dB_vec, R_s(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(JPR_dB_vec, R_s(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
yline(R_s_perfect(1), 'b:', sprintf('MRT perfect: %.2f', R_s_perfect(1)), 'LineWidth', 1.2);
yline(R_s_perfect(2), 'r:', sprintf('ZF perfect: %.2f',  R_s_perfect(2)), 'LineWidth', 1.2);
pls_axis_prefs(gca, 'refLabelV', 'bottom', 'staggerRef', true);
grid on; box on;
xlabel('Jammer-to-Pilot Ratio (dB)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('MRT vs ZF Robustness');
legend('MRT', 'ZF', 'Location', 'SouthWest');

% Bottom-left: Heat-map Nt vs JPR (UDOSKONALONA)
subplot(2, 2, 3);
imagesc(JPR_grid_dB, 1:length(Nt_vec), SR_grid);
set(gca, 'YDir', 'normal', 'YTick', 1:length(Nt_vec), 'YTickLabel', Nt_vec);
colormap(parula); 

% Pasek kolorów z wyraźnym podpisem jednostki obok
cb = colorbar;
cb.Label.String = 'Secrecy Sum-Rate (bits/s/Hz)';
cb.Label.FontSize = 10;
cb.Label.FontWeight = 'bold';

xlabel('Jammer-to-Pilot Ratio (dB)'); ylabel('Antennas (N_t)');
title('Secrecy Sum-Rate (ZF): N_t vs JPR');
% Bottom-right: Array Gain vs Jamming
subplot(2, 2, 4);
plot(Nt_vec, SR_grid(:, end), '-mo', 'LineWidth', 2, 'MarkerFaceColor', 'm'); hold on;
plot(Nt_vec, SR_grid(:, ceil(end/2)), '-co', 'LineWidth', 2, 'MarkerFaceColor', 'c');
plot(Nt_vec, SR_grid(:, 1),   '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
grid on; box on;
xlabel('Antennas N_t'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Massive MIMO as Defence Against Jamming');
legend(sprintf('JPR = %d dB (Heavy)', JPR_grid_dB(end)), ...
       sprintf('JPR = %d dB (Medium)', JPR_grid_dB(ceil(end/2))), ...
       sprintf('JPR = %d dB (Mild)', JPR_grid_dB(1)), ...
       'Location', 'NorthWest');

sgtitle(sprintf('Scenario 7: Pilot Jamming DoS Attack (Base Rx SNR \\approx %.1f dB)', Rx_SNR_dB_b));

save_figure(fig, 'fig_pilot_jamming');

% =========================================================================
% FUNKCJE POMOCNICZE
% =========================================================================
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

function plot_jamming_topology(dist_b, theta_bobs, dist_e, theta_e)
    fig_top = figure('Color', 'w', 'Position', [150 150 700 700]);
    hold on; grid on; box on;
    
    x_bs = 0; y_bs = 0;
    max_d = max(dist_b, dist_e) + 15;
    
    x_e = dist_e * sind(theta_e);
    y_e = dist_e * cosd(theta_e);
    p_e = plot(x_e, y_e, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Eve (Jammer)');
    text(x_e + 2, y_e, sprintf('Eve (Jammer)\n(%gm, %g\\circ)', dist_e, theta_e), 'Color', 'r', 'FontSize', 9);
    
    p_b = [];
    for k = 1:length(theta_bobs)
        x_b = dist_b * sind(theta_bobs(k));
        y_b = dist_b * cosd(theta_bobs(k));
        p_b = plot(x_b, y_b, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        
        if k == 1 || k == length(theta_bobs)
            text(x_b - 2.5, y_b - 2.5, sprintf('B_{%d}\n(%gm, %g\\circ)', k, dist_b, theta_bobs(k)), ...
                 'Color', 'b', 'FontSize', 8, 'HorizontalAlignment', 'right');
        end
    end
    if ~isempty(p_b)
        set(p_b, 'DisplayName', sprintf('Bobs (K=%d)', length(theta_bobs)));
    end
    
    p_bs = plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    text(x_bs, y_bs - 4, 'BS (0,0)', 'HorizontalAlignment', 'center', 'Color', 'k');
    
    axis equal;
    xlim([-max_d, max_d]);
    ylim([-10, max_d]);
    xlabel('X [m]'); ylabel('Y [m]');
    title('Scenario 7: Pilot Jamming Attack (Denial-of-Service)');
    legend([p_bs, p_b, p_e], 'Location', 'NorthWest');
    
    try
        save_figure(fig_top, '../topology/topology_pilot_jamming');
    catch
        warning('Funkcja save_figure nie jest dostępna.');
    end
end