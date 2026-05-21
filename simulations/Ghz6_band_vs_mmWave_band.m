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
    [~, bands(b).sv] = setup_ula(bands(b).Nt, bands(b).fc);
    
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
% TWORZENIE WYKRESÓW I ZDJĘCIE SYTUACYJNE (SNAPSHOT)
% =========================================================================
snap_theta_b = theta_b; snap_theta_e = theta_e;
angles = -90:0.05:90;

% Tworzymy figurę i twardo ustawiamy białe tło środowiska
fig = figure('Color', 'w', 'Position', [100 100 1100 800]);

for f_idx = 1:2
    fc = bands(f_idx).fc; 
    Nt = bands(f_idx).Nt; 
    sv = bands(f_idx).sv;
    cdl_tag = bands(f_idx).cdl;
    PL_dB_b = bands(f_idx).PL_dB_b;
    bandStr = [bands(f_idx).name, ' (nrCDLChannel)'];
    
    % Przeliczenie osi X na rzeczywisty Received SNR u Boba
    SNR_rx_dB_plot = P_tx_norm_dB + 10*log10(Nt) - PL_dB_b;
    
    cdl_b_s = setup_matlab_cdl(Nt, fc, snap_theta_b, cdl_tag);
    cdl_e_s = setup_matlab_cdl(Nt, fc, snap_theta_e, cdl_tag);
    cdl_b_s.Seed = 102;
    cdl_e_s.Seed = 102; 
    
    [pg_b_s, ~] = cdl_b_s();
    [pg_e_s, ~] = cdl_e_s();
    
    hb_s = squeeze(sum(pg_b_s, 2)); hb_s = hb_s(:); 
    he_s = squeeze(sum(pg_e_s, 2)); he_s = he_s(:); 
    
    w_mrt_s = hb_s / norm(hb_s);
    P_null_s = eye(Nt) - (he_s * (he_s' / (he_s' * he_s)));
    w_zf_s = P_null_s * hb_s; w_zf_s = w_zf_s / norm(w_zf_s);
    
    a_sweep = step(sv, fc, angles);
    
    % --- WYKRES POJEMNOŚCI ---
    ax1 = subplot(2, 2, f_idx);
    plot(SNR_rx_dB_plot, squeeze(SecrecyCap(f_idx,1,:)), 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
    plot(SNR_rx_dB_plot, squeeze(SecrecyCap(f_idx,2,:)), 'r--s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    grid on; box on;
    
    title(['Secrecy: ', bandStr], 'Color', 'k');
    % OŚ X USTAWIENIE: Received SNR at Bob (dB)
    xlabel('Average Receiver SNR at Bob (dB)', 'Color', 'k'); 
    ylabel('bits/s/Hz', 'Color', 'k');
    
    lgd = legend('MRT', 'ZF', 'Location', 'NorthWest');
    set(lgd, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
    xlim([min(SNR_rx_dB_plot), max(SNR_rx_dB_plot)]);
    set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.15);
    
    % --- WYKRES WIĄZKI ABSOLUTNEJ (Antenna Gain, PSLL, HPBW) ---
    pat_mrt_abs = 10*log10(abs(w_mrt_s' * a_sweep).^2);
    pat_zf_abs  = 10*log10(abs(w_zf_s'  * a_sweep).^2);
    
    ax2 = subplot(2, 2, f_idx + 2);
    plot(angles, pat_mrt_abs, 'b',  'LineWidth', 1.7); hold on;
    plot(angles, pat_zf_abs,  'r--','LineWidth', 1.5);
    
    % Szukanie parametrów głównej wiązki (na przykładzie wiązki MRT)
    [peak_gain, peak_idx] = max(pat_mrt_abs);
    
    % --- OBLICZANIE HPBW (3dB Beamwidth) ---
    level_3db = peak_gain - 3.0;
    
    % Szukaj w lewo od piku
    idx_left = find(pat_mrt_abs(1:peak_idx) <= level_3db, 1, 'last');
    if isempty(idx_left), idx_left = 1; end
    
    % Szukaj w prawo od piku
    idx_right = find(pat_mrt_abs(peak_idx:end) <= level_3db, 1, 'first') + peak_idx - 1;
    if isempty(idx_right), idx_right = length(angles); end
    
    hpbw_deg = angles(idx_right) - angles(idx_left);
    
    % Zaznaczenie HPBW na wykresie
    plot([angles(idx_left), angles(idx_right)], [level_3db, level_3db], 'k|-', 'LineWidth', 1.5, 'MarkerSize', 6);
    
    % --- OBLICZANIE PSLL ---
    [pks, locs] = findpeaks(pat_mrt_abs);
    exclusion_zone_deg = 5; 
    angle_step = angles(2) - angles(1);
    exclusion_samples = round(exclusion_zone_deg / angle_step);
    
    valid_idx = abs(locs - peak_idx) > exclusion_samples;
    valid_pks = pks(valid_idx);
    valid_locs = locs(valid_idx);
    
    if ~isempty(valid_pks)
        [sidelobe_gain, max_sidelobe_idx_temp] = max(valid_pks);
        sidelobe_idx = valid_locs(max_sidelobe_idx_temp);
        psll = peak_gain - sidelobe_gain;
        
        plot(angles(sidelobe_idx), sidelobe_gain, 'kv', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
        plot([angles(sidelobe_idx), angles(peak_idx)], [sidelobe_gain, sidelobe_gain], 'k:', 'LineWidth', 1.5);
        
        txt = sprintf('Max Gain: %.1f dB\nHPBW (3dB): %.1f\\circ\nPSLL: %.1f dB', peak_gain, hpbw_deg, psll);
    else
        txt = sprintf('Max Gain: %.1f dB\nHPBW (3dB): %.1f\\circ\nPSLL: N/A', peak_gain, hpbw_deg);
    end
    
    text(0.96, 0.94, txt, 'Units', 'normalized', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
         'FontSize', 10, 'Color', 'k', 'BackgroundColor', 'w', ...
         'EdgeColor', 'k', 'Margin', 4);
         
    mark_bob(snap_theta_b, 'Bob');
    mark_eve(snap_theta_e, 'Eve');
    pls_axis_prefs(gca, 'refLabelV', 'top', 'staggerRef', true);
    
    grid on; box on;
    title(['Absolute beam pattern: ', bandStr], 'Color', 'k');
    xlabel('Angle (deg)', 'Color', 'k'); 
    ylabel('Gain (dB)', 'Color', 'k');
    
    ylim([-60, ceil(peak_gain/5)*5 + 5]);
    xlim([-90, 90]);
    set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.15);
end

sgt = sgtitle(sprintf('6G PLS Baseline (Bob = %gm, Eve = %gm)', dist_b, dist_e));
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

function plot_topology(dist_b, theta_b, dist_e, theta_e)
    fig_top = figure('Color', 'w', 'Position', [150 150 600 600]);
    x_bs = 0; y_bs = 0;
    x_b = dist_b * sind(theta_b); y_b = dist_b * cosd(theta_b);
    x_e = dist_e * sind(theta_e); y_e = dist_e * cosd(theta_e);
    
    ax = gca;
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.15);
    
    hold on; grid on; box on;
    t = linspace(0, 2*pi, 100);
    plot(dist_b * sin(t), dist_b * cos(t), 'b:', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(dist_e * sin(t), dist_e * cos(t), 'r:', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot([x_bs, x_b], [y_bs, y_b], 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot([x_bs, x_e], [y_bs, y_e], 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    plot(x_b, y_b, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Bob');
    plot(x_e, y_e, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Eve');
    
    text(x_b + 2, y_b, sprintf('Bob\n(%gm, %g\\circ)', dist_b, theta_b), 'Color', 'b', 'FontSize', 10);
    text(x_e + 2, y_e, sprintf('Eve\n(%gm, %g\\circ)', dist_e, theta_e), 'Color', 'r', 'FontSize', 10);
    text(x_bs, y_bs - 3, 'BS (0,0)', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k');
    
    axis equal;
    max_d = max(dist_b, dist_e) + 15;
    xlim([-max_d, max_d]); ylim([-max_d/2, max_d]);
    xlabel('X [m]', 'Color', 'k'); ylabel('Y [m]', 'Color', 'k');
    title('Scenario 1: 6GHz vs mmWave', 'Color', 'k');
    
    lgd = legend('Location', 'NorthEast');
    set(lgd, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
    
    set(fig_top, 'InvertHardcopy', 'off');
    set(fig_top, 'Color', 'w');
    try
        save_figure(fig_top, '../topology/topology_baseline');
    catch
    end
end