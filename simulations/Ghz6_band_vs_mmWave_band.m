% =========================================================================
% BASELINE: 6 GHz Massive MIMO vs 28 GHz Ultra-Massive MIMO
% -------------------------------------------------------------------------
% ZAŁOŻENIA MODELU (Do obrony projektu):
% 1. Kanał: 3GPP TR 38.901 UMi NLOS (nrCDLChannel).
% 2. Zagrożenie: Active/Untrusted Eavesdropper. Ewa jest prawowitym 
%    użytkownikiem sieci, stąd BS posiada idealną estymatę jej kanału 
%    (Perfect CSI) wykorzystywaną w rzutowaniu ortogonalnym ZF.
% 3. Szum: Znormalizowany do N_0 = 1 (0 dB). 
% 4. Oś X: Reprezentuje uśredniony odebrany SNR u Boba (Received SNR).
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- FIZYKA DYSTANSU I MOCY ---
dist_b     = 50;                        % Bob jest na 50 m
dist_e     = 40;                        % Ewa ukrywa się na 40 m (bliżej stacji!)

% PARAMETR WEJŚCIOWY: Znormalizowana Moc Nadawania (P_tx / N_0)
P_tx_norm_dB  = 70:2:140;                  
P_tx_norm_lin = 10.^(P_tx_norm_dB / 10);
numIter       = 50;

bands = struct( ...
    'name', {'6 GHz', '28 GHz'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave}, ...
    'cdl',  {p.cdl_sub6, p.cdl_mmwave});

% Obliczamy tłumienia z wyprzedzeniem używając 3GPP Path Loss
for b = 1:2
    [bands(b).PL_lin_b, bands(b).PL_dB_b] = compute_nr_pathloss(dist_b, bands(b).fc);
    [bands(b).PL_lin_e, bands(b).PL_dB_e] = compute_nr_pathloss(dist_e, bands(b).fc);
    
    rx_snr_start = P_tx_norm_dB(1) + 10*log10(bands(b).Nt) - bands(b).PL_dB_b;
    rx_snr_end   = P_tx_norm_dB(end) + 10*log10(bands(b).Nt) - bands(b).PL_dB_b;
    
    fprintf('\n--- Baseline @ %s ---\n', bands(b).name);
    fprintf('  Normalized Transmit Power: %d:%d:%d dB\n', P_tx_norm_dB(1), P_tx_norm_dB(2)-P_tx_norm_dB(1), P_tx_norm_dB(end));
    fprintf('  Equivalent Bob Rx SNR: %.1f to %.1f dB\n', rx_snr_start, rx_snr_end);
    fprintf('  Bob (d = %g m, PL = %.2f dB)\n', dist_b, bands(b).PL_dB_b);
    fprintf('  Eve (d = %g m, PL = %.2f dB)\n', dist_e, bands(b).PL_dB_e);
end

fprintf('\nmmWave path-loss penalty vs 6 GHz (Bob): %.2f dB\n', ...
    bands(2).PL_dB_b - bands(1).PL_dB_b);

SecrecyCap = zeros(2, 2, length(P_tx_norm_dB)); 
theta_b = -30; 
theta_e = 40;

for f_idx = 1:2
    fc = bands(f_idx).fc;  
    Nt = bands(f_idx).Nt;
    cdl_tag = bands(f_idx).cdl;
    
    PL_b = bands(f_idx).PL_lin_b;
    PL_e = bands(f_idx).PL_lin_e;
    
    % --- INICJALIZACJA KANAŁÓW 3GPP (5G Toolbox) ---
    cdl_b = setup_matlab_cdl(Nt, fc, theta_b, cdl_tag);
    cdl_e = setup_matlab_cdl(Nt, fc, theta_e, cdl_tag);
    
    for it = 1:numIter
        release(cdl_b);
        release(cdl_e);
        iter_seed = randi([0 2^31-1]); 
        cdl_b.Seed = iter_seed;
        cdl_e.Seed = iter_seed;
        
        [pg_b, ~] = cdl_b();
        [pg_e, ~] = cdl_e();
        
        hb = squeeze(sum(pg_b, 2)); hb = hb(:);
        he = squeeze(sum(pg_e, 2)); he = he(:);
        
        % Kanały efektywne uwzględniające potężny Path Loss UMi NLOS
        hb_eff = sqrt(1/PL_b) * hb;
        he_eff = sqrt(1/PL_e) * he;
        
        w_mrt = hb_eff / norm(hb_eff);
        P_null = eye(Nt) - (he_eff * (he_eff' / (he_eff' * he_eff)));
        w_zf   = P_null * hb_eff;
        
        if norm(w_zf) > 1e-9
            w_zf = w_zf / norm(w_zf);
        else
            w_zf = zeros(Nt, 1);
        end
        
        for s = 1:length(P_tx_norm_lin)
            P_tx = P_tx_norm_lin(s); 
            
            % Obliczenia przepustowości przy N_0 = 1
            R_b = log2(1 + P_tx * abs(hb_eff' * w_mrt)^2);
            R_e = log2(1 + P_tx * abs(he_eff' * w_mrt)^2);
            SecrecyCap(f_idx,1,s) = SecrecyCap(f_idx,1,s) + max(0, R_b - R_e);
            
            R_b = log2(1 + P_tx * abs(hb_eff' * w_zf)^2);
            R_e = log2(1 + P_tx * abs(he_eff' * w_zf)^2);
            SecrecyCap(f_idx,2,s) = SecrecyCap(f_idx,2,s) + max(0, R_b - R_e);
        end
    end
end
SecrecyCap = SecrecyCap / numIter;

% =========================================================================
% TWORZENIE WYKRESÓW POJEMNOŚCI (PANORAMA 1x2)
% =========================================================================
% Panoramiczne proporcje figury dopasowane do układu dwóch wykresów obok siebie
fig = figure('Color', 'w', 'Position', [100 100 1200 480]);

for f_idx = 1:2
    Nt = bands(f_idx).Nt; 
    PL_dB_b = bands(f_idx).PL_dB_b;
    bandStr = [bands(f_idx).name, ' (nrCDLChannel)'];
    
    % Przeliczenie osi X na rzeczywisty Received SNR u Boba
    SNR_rx_dB_plot = P_tx_norm_dB + 10*log10(Nt) - PL_dB_b;
    
    % --- WYKRES POJEMNOŚCI (Układ 1x2) ---
    ax1 = subplot(1, 2, f_idx);
    plot(SNR_rx_dB_plot, squeeze(SecrecyCap(f_idx,1,:)), 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
    plot(SNR_rx_dB_plot, squeeze(SecrecyCap(f_idx,2,:)), 'r--s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    grid on; box on;
    
    title(['Secrecy: ', bandStr], 'Color', 'k');
    xlabel('Average Receiver SNR at Bob (dB)', 'Color', 'k'); 
    ylabel('Secrecy Rate (bits/s/Hz)', 'Color', 'k');
    
    lgd = legend('MRT', 'ZF', 'Location', 'NorthWest');
    set(lgd, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
    xlim([min(SNR_rx_dB_plot), max(SNR_rx_dB_plot)]);
    set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.15);
end

sgt = sgtitle(sprintf('Scenario1: 6GHz vs mmWave Comparison'));
set(sgt, 'Color', 'k', 'FontWeight', 'bold');
set(fig, 'InvertHardcopy', 'off');
set(fig, 'Color', 'w');

try
    save_figure(fig, 'fig_baseline_6GHz_vs_28GHz');
catch
    warning('Funkcja save_figure nie jest dostępna. Wykres nie został zapisany automatycznie.');
end

plot_topology(dist_b, theta_b, dist_e, theta_e);

% =========================================================================
% FUNKCJE POMOCNICZE
% =========================================================================
function cdl = setup_matlab_cdl(Nt, fc, theta, band_tag)
    cdl = nrCDLChannel;
    cdl.DelayProfile = 'CDL-A';
    if fc < 10e9  
        cdl.DelaySpread = 30e-9;  
    else          
        cdl.DelaySpread = 10e-10; % dostosowane do mmWave  
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

function plot_topology(dist_b, theta_b, dist_e, theta_e)
    fig_top = figure('Color', 'w', 'Position', [150 150 600 600]);
    x_bs = 0; y_bs = 0;
    x_b = dist_b * sind(theta_b); y_b = dist_b * cosd(theta_b);
    x_e = dist_e * sind(theta_e); y_e = dist_e * cosd(theta_e);
    
    ax = gca;
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.15);
    
    hold on; grid on; box on;
    
    % Rysowanie wyłącznie punktów obiektów fizycznych
    p_bs = plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    p_b  = plot(x_b, y_b, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Bob');
    p_e  = plot(x_e, y_e, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Eve');
    
    % Tekstowe etykiety współrzędnych i odległości
    text(x_b + 2, y_b, sprintf('Bob\n(%gm, %g\\circ)', dist_b, theta_b), 'Color', 'b', 'FontSize', 10);
    text(x_e + 2, y_e, sprintf('Eve\n(%gm, %g\\circ)', dist_e, theta_e), 'Color', 'r', 'FontSize', 10);
    text(x_bs, y_bs - 3, 'BS (0,0)', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k');
    
    axis equal;
    max_d = max(dist_b, dist_e) + 15;
    xlim([-max_d, max_d]); ylim([-max_d/2, max_d]);
    xlabel('X [m]', 'Color', 'k'); ylabel('Y [m]', 'Color', 'k');
    title('Scenario 1: 6GHz vs mmWave Topology', 'Color', 'k');
    
    lgd = legend([p_bs, p_b, p_e], 'Location', 'NorthEast');
    set(lgd, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
    
    set(fig_top, 'InvertHardcopy', 'off');
    set(fig_top, 'Color', 'w');
    try
        save_figure(fig_top, '../topology/topology_baseline');
    catch
    end
end