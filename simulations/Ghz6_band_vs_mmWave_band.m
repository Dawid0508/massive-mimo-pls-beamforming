% =========================================================================
% BASELINE: 6 GHz Massive MIMO vs 28 GHz Ultra-Massive MIMO with FSPL
% (Zaktualizowano do oficjalnego modelu nrCDLChannel z 5G Toolbox)
% Poprawiono logikę Asymetrii Dystansu i Transmit SNR (Brak "Magicznej Ewy")
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- FIZYKA DYSTANSU I MOCY ---
dist_b     = 50;                        % Bob jest na 50 m
dist_e     = 40;                        % Ewa ukrywa się na 40 m (bliżej stacji!)
SNR_tx_dB  = 60:2:130;                  % Oś X: Transmit SNR (Moc stacji bazowej)
SNR_tx_lin = 10.^(SNR_tx_dB / 10);
numIter    = 200;

bands = struct( ...
    'name', {'6 GHz', '28 GHz'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave}, ...
    'cdl',  {p.cdl_sub6, p.cdl_mmwave});

% Obliczamy tłumienia z wyprzedzeniem
for b = 1:2
    [bands(b).PL_lin_b, bands(b).PL_dB_b] = compute_fspl(dist_b, bands(b).fc);
    [bands(b).PL_lin_e, bands(b).PL_dB_e] = compute_fspl(dist_e, bands(b).fc);
    [~, bands(b).sv] = setup_ula(bands(b).Nt, bands(b).fc);
    
    % --- NOWY, POPRAWNY WYDRUK W KONSOLI ---
    fprintf('\n--- Baseline @ %s ---\n', bands(b).name);
    fprintf('  Transmit SNR sweep: %d:%d:%d dB\n', SNR_tx_dB(1), SNR_tx_dB(2)-SNR_tx_dB(1), SNR_tx_dB(end));
    fprintf('  Bob (d = %g m, PL = %.2f dB)\n', dist_b, bands(b).PL_dB_b);
    fprintf('  Eve (d = %g m, PL = %.2f dB)\n', dist_e, bands(b).PL_dB_e);
end
fprintf('\nmmWave path-loss penalty vs 6 GHz (Bob): %.2f dB\n', ...
    bands(2).PL_dB_b - bands(1).PL_dB_b);

SecrecyCap = zeros(2, 2, length(SNR_tx_dB));    % bands x {MRT,ZF} x SNR
theta_b = -30; 
theta_e = -40;

for f_idx = 1:2
    fc = bands(f_idx).fc;  
    Nt = bands(f_idx).Nt;
    cdl_tag = bands(f_idx).cdl;
    
    % Pobieramy wyliczone wczesniej fizyczne tlumienie przestrzenne
    PL_b = bands(f_idx).PL_lin_b;
    PL_e = bands(f_idx).PL_lin_e;
    
    % --- INICJALIZACJA KANAŁÓW 3GPP (5G Toolbox) ---
    cdl_b = setup_matlab_cdl(Nt, fc, theta_b, cdl_tag);
    cdl_e = setup_matlab_cdl(Nt, fc, theta_e, cdl_tag);
    
    for it = 1:numIter
        release(cdl_b);
        release(cdl_e);
        cdl_b.Seed = randi([0 2^31-1]);
        cdl_e.Seed = randi([0 2^31-1]);
        
        [pg_b, ~] = cdl_b();
        [pg_e, ~] = cdl_e();
        
        hb = squeeze(sum(pg_b, 2)); hb = hb(:);
        he = squeeze(sum(pg_e, 2)); he = he(:);
        
        % TWORZYMY KANAŁY EFEKTYWNE (Zysk anten + Tłumienie w powietrzu)
        % To tutaj zaszyta jest cała twarda fizyka dystansu (Near-Far effect)
        hb_eff = sqrt((1/PL_b) * Nt) * (hb / norm(hb));
        he_eff = sqrt((1/PL_e) * Nt) * (he / norm(he));
        
        % Prekodery wyliczane z kanału efektywnego
        w_mrt = hb_eff / norm(hb_eff);
        P_null = eye(Nt) - (he_eff * (he_eff' / (he_eff' * he_eff)));
        w_zf   = P_null * hb_eff;
        
        if norm(w_zf) > 1e-9
            w_zf = w_zf / norm(w_zf);
        else
            w_zf = zeros(Nt, 1);
        end
        
        % OBLICZENIA POJEMNOŚCI
        for s = 1:length(SNR_tx_lin)
            P_tx = SNR_tx_lin(s); % Stacja bazowa wypuszcza P_tx z anteny
            
            % Tłumienie trasy zjada sygnał naturalnie dzięki wektorom hb_eff/he_eff
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

% --- TWORZENIE WYKRESÓW I ZDJĘCIE SYTUACYJNE (SNAPSHOT) ---
snap_theta_b = theta_b; snap_theta_e = theta_e;
angles = -90:0.05:90;
fig = figure('Color', 'w', 'Position', [100 100 1100 800]);

for f_idx = 1:2
    fc = bands(f_idx).fc; Nt = bands(f_idx).Nt; sv = bands(f_idx).sv;
    cdl_tag = bands(f_idx).cdl;
    bandStr = [bands(f_idx).name, ' (nrCDLChannel)'];
    
    cdl_b_s = setup_matlab_cdl(Nt, fc, snap_theta_b, cdl_tag);
    cdl_e_s = setup_matlab_cdl(Nt, fc, snap_theta_e, cdl_tag);
    cdl_b_s.Seed = 101;
    cdl_e_s.Seed = 101; 
    
    [pg_b_s, ~] = cdl_b_s();
    [pg_e_s, ~] = cdl_e_s();
    
    hb_s = squeeze(sum(pg_b_s, 2)); hb_s = hb_s(:); 
    he_s = squeeze(sum(pg_e_s, 2)); he_s = he_s(:); 
    
    % Do narysowania migawki wiązki interesują nas tylko kierunki, więc tu
    % używamy zwykłej normalizacji, aby wykres kątowy był czysty i wyraźny.
    w_mrt_s = hb_s / norm(hb_s);
    P_null_s = eye(Nt) - (he_s * (he_s' / (he_s' * he_s)));
    w_zf_s = P_null_s * hb_s; w_zf_s = w_zf_s / norm(w_zf_s);
    
    a_sweep = step(sv, fc, angles);
    pat_mrt = 10*log10(abs(w_mrt_s' * a_sweep).^2);
    pat_zf  = 10*log10(abs(w_zf_s'  * a_sweep).^2);
    
    subplot(2, 2, f_idx);
    plot(SNR_tx_dB, squeeze(SecrecyCap(f_idx,1,:)), 'b-o', 'LineWidth', 1.5); hold on;
    plot(SNR_tx_dB, squeeze(SecrecyCap(f_idx,2,:)), 'r--s', 'LineWidth', 1.5);
    grid on; box on;
    title(['Secrecy: ', bandStr]);
    % Oś X jest teraz poprawnie podpisana jako moc stacji bazowej
    xlabel('Transmit SNR at Base Station (dB)'); 
    ylabel('bits/s/Hz');
    legend('MRT', 'ZF', 'Location', 'NorthWest');
    
    subplot(2, 2, f_idx + 2);
    plot(angles, pat_mrt - max(pat_mrt), 'b',  'LineWidth', 1.7); hold on;
    plot(angles, pat_zf  - max(pat_zf),  'r--','LineWidth', 1.5);
    mark_bob(snap_theta_b, 'Bob');
    mark_eve(snap_theta_e, 'Eve');
    pls_axis_prefs(gca, 'refLabelV', 'top', 'staggerRef', true);
    grid on; box on;
    title(['Normalised beam pattern: ', bandStr]);
    xlabel('Angle (deg)'); ylabel('Gain (dB)');
    ylim([-40 5]);
end
sgtitle(sprintf('6G PLS Baseline (Bob = %gm, Eve = %gm)', dist_b, dist_e));
save_figure(fig, 'fig_baseline_6GHz_vs_28GHz');
plot_topology(dist_b, theta_b, dist_e, theta_e);

% =========================================================================
% FUNKCJA POMOCNICZA: Konfiguracja kanału nrCDLChannel
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

% =========================================================================
% FUNKCJA POMOCNICZA: Generowanie topologii scenariusza
% =========================================================================
function plot_topology(dist_b, theta_b, dist_e, theta_e)
    fig_top = figure('Color', 'w', 'Position', [150 150 600 600]);
    
    % Przeliczenie współrzędnych biegunowych na kartezjańskie
    % Zakładamy, że stacja bazowa (BS) jest w punkcie (0,0)
    % Zgodnie z osią broadside ULA (0 stopni na osi Y, kąty rosną względem niej)
    x_bs = 0; y_bs = 0;
    x_b = dist_b * sind(theta_b); y_b = dist_b * cosd(theta_b);
    x_e = dist_e * sind(theta_e); y_e = dist_e * cosd(theta_e);
    
    hold on; grid on; box on;
    
    % Okręgi dystansu (promienie dla Boba i Ewy)
    t = linspace(0, 2*pi, 100);
    plot(dist_b * sin(t), dist_b * cos(t), 'b:', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(dist_e * sin(t), dist_e * cos(t), 'r:', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Linie kierunkowe (LOS path)
    plot([x_bs, x_b], [y_bs, y_b], 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot([x_bs, x_e], [y_bs, y_e], 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Rysowanie punktów
    plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    plot(x_b, y_b, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Bob');
    plot(x_e, y_e, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Eve');
    
    % Dodanie etykiet tekstowych
    text(x_b + 2, y_b, sprintf('Bob\n(%gm, %g\\circ)', dist_b, theta_b), 'Color', 'b', 'FontSize', 10);
    text(x_e + 2, y_e, sprintf('Eve\n(%gm, %g\\circ)', dist_e, theta_e), 'Color', 'r', 'FontSize', 10);
    text(x_bs, y_bs - 3, 'BS (0,0)', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k');
    
    % Formatowanie osi układu
    axis equal;
    max_d = max(dist_b, dist_e) + 15;
    xlim([-max_d, max_d]);
    ylim([-max_d/2, max_d]); % Dopasowane do widoku "z przodu" anteny
    xlabel('X [m]'); ylabel('Y [m]');
    title('Scenario 1: 6GHz vs mmWave');
    legend('Location', 'NorthEast');

    save_figure(fig_top, '../topology/topology_baseline');
end