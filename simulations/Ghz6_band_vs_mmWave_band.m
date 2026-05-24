% =========================================================================
% BASELINE: 6 GHz Massive MIMO vs 28 GHz Ultra-Massive MIMO
% -------------------------------------------------------------------------
% ZAŁOŻENIA MODELU (Passive Eve & K=4 Bobów):
% 1. Kanał: 3GPP TR 38.901 UMi NLOS (nrCDLChannel).
% 2. Zagrożenie: Passive Eavesdropper. BS NIE ZNA kanału Ewy.
%    Prekodowanie (ZF) eliminuje jedynie interferencje międzysystemowe 
%    (MUI) dla K=4 Bobów, ignorując Ewę.
% 3. Szum: Znormalizowany do N_0 = 1 (0 dB). 
% 4. Oś X: Reprezentuje uśredniony odebrany SNR u Boba (Received SNR).
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- FIZYKA DYSTANSU I MOCY ---
K          = 4;                         % Liczba prawowitych użytkowników (Bob)
dist_b     = 50;                        % Bobowie są na 50 m
theta_b    = -60 + 120 * rand(1, K);    % Losowe kąty Bobów od -60 do 60 stopni
dist_e     = 40;                        % Ewa ukrywa się na 40 m
theta_e    = 40;                        % Kąt Ewy

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
    fprintf('  Bob (d = %g m, PL = %.2f dB, K = %d)\n', dist_b, bands(b).PL_dB_b, K);
    fprintf('  Eve (d = %g m, PL = %.2f dB)\n', dist_e, bands(b).PL_dB_e);
end
fprintf('\nmmWave path-loss penalty vs 6 GHz (Bob): %.2f dB\n', ...
    bands(2).PL_dB_b - bands(1).PL_dB_b);

SecrecyCap = zeros(2, 2, length(P_tx_norm_dB)); 

for f_idx = 1:2
    fc = bands(f_idx).fc;  
    Nt = bands(f_idx).Nt;
    cdl_tag = bands(f_idx).cdl;
    
    PL_b = bands(f_idx).PL_lin_b;
    PL_e = bands(f_idx).PL_lin_e;
    
    % --- INICJALIZACJA KANAŁÓW 3GPP (5G Toolbox) ---
    cdl_b = cell(1, K);
    for k = 1:K
        cdl_b{k} = setup_matlab_cdl(Nt, fc, theta_b(k), cdl_tag);
    end
    cdl_e = setup_matlab_cdl(Nt, fc, theta_e, cdl_tag);
    
    for it = 1:numIter
        iter_seed = randi([0 2^31-1]); 
        theta_b_iter = -60 + 120 * rand(1, K);
        
        Hb_eff = zeros(Nt, K);
        for k = 1:K
            release(cdl_b{k});
            cdl_b{k}.TransmitArrayOrientation = [-theta_b_iter(k); 0; 0];
            cdl_b{k}.Seed = iter_seed;
            [pg_b, ~] = cdl_b{k}();
            hb = squeeze(sum(pg_b, 2)); hb = hb(:);
            Hb_eff(:, k) = sqrt(1/PL_b) * hb;
        end
        
        release(cdl_e);
        cdl_e.Seed = iter_seed;
        [pg_e, ~] = cdl_e();
        he = squeeze(sum(pg_e, 2)); he = he(:);
        he_eff = sqrt(1/PL_e) * he;
        
        % --- PREKODOWANIE (Tylko dla K Bobów, bez wiedzy o Ewie) ---
        % 1. MRT
        W_mrt = Hb_eff; 
        W_mrt = W_mrt / norm(W_mrt, 'fro'); % Normalizacja mocy sumarycznej
        
        % 2. ZF
        % Zamiast rzutować w przestrzeń zerową Ewy, odwracamy interferencje między Bobami
        W_zf = Hb_eff / (Hb_eff' * Hb_eff); 
        W_zf = W_zf / norm(W_zf, 'fro'); % Normalizacja mocy sumarycznej
        
        for s = 1:length(P_tx_norm_lin)
            P_tx = P_tx_norm_lin(s); 
            p_k = P_tx / K; % Równy podział mocy na K użytkowników
            
            sum_sec_mrt = 0;
            sum_sec_zf = 0;
            
            for k = 1:K
                % --- Analiza dla MRT ---
                sig_b_mrt = p_k * abs(Hb_eff(:,k)' * W_mrt(:,k))^2;
                inf_b_mrt = p_k * sum(abs(Hb_eff(:,k)' * W_mrt(:,[1:k-1, k+1:K])).^2);
                R_b_mrt   = log2(1 + sig_b_mrt / (1 + inf_b_mrt));
                
                sig_e_mrt = p_k * abs(he_eff' * W_mrt(:,k))^2;
                % Ewa traktuje strumienie innych Bobów jako szum (interferencje)
                inf_e_mrt = p_k * sum(abs(he_eff' * W_mrt(:,[1:k-1, k+1:K])).^2);
                R_e_mrt   = log2(1 + sig_e_mrt / (1 + inf_e_mrt));
                
                sum_sec_mrt = sum_sec_mrt + max(0, R_b_mrt - R_e_mrt);
                
                % --- Analiza dla ZF ---
                sig_b_zf = p_k * abs(Hb_eff(:,k)' * W_zf(:,k))^2;
                inf_b_zf = p_k * sum(abs(Hb_eff(:,k)' * W_zf(:,[1:k-1, k+1:K])).^2);
                R_b_zf   = log2(1 + sig_b_zf / (1 + inf_b_zf));
                
                sig_e_zf = p_k * abs(he_eff' * W_zf(:,k))^2;
                inf_e_zf = p_k * sum(abs(he_eff' * W_zf(:,[1:k-1, k+1:K])).^2);
                R_e_zf   = log2(1 + sig_e_zf / (1 + inf_e_zf));
                
                sum_sec_zf = sum_sec_zf + max(0, R_b_zf - R_e_zf);
            end
            
            SecrecyCap(f_idx,1,s) = SecrecyCap(f_idx,1,s) + sum_sec_mrt;
            SecrecyCap(f_idx,2,s) = SecrecyCap(f_idx,2,s) + sum_sec_zf;
        end
    end
end
SecrecyCap = SecrecyCap / numIter;

% =========================================================================
% TWORZENIE WYKRESÓW POJEMNOŚCI (PANORAMA 1x2)
% =========================================================================
fig = figure('Color', 'w', 'Position', [100 100 1200 480]);
for f_idx = 1:2
    Nt = bands(f_idx).Nt; 
    PL_dB_b = bands(f_idx).PL_dB_b;
    bandStr = [bands(f_idx).name, ' (nrCDLChannel)'];
    
    % Przeliczenie osi X na rzeczywisty Received SNR u Boba
    SNR_rx_dB_plot = P_tx_norm_dB + 10*log10(Nt) - PL_dB_b - 10*log10(K);
    
    ax1 = subplot(1, 2, f_idx);
    plot(SNR_rx_dB_plot, squeeze(SecrecyCap(f_idx,1,:)), 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
    plot(SNR_rx_dB_plot, squeeze(SecrecyCap(f_idx,2,:)), 'r--s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    grid on; box on;
    
    title(['Secrecy: ', bandStr], 'Color', 'k');
    xlabel('Average Receiver SNR per Bob (dB)', 'Color', 'k'); 
    ylabel('Secrecy Sum Rate (bits/s/Hz)', 'Color', 'k');
    
    lgd = legend('MRT', 'ZF', 'Location', 'NorthWest');
    set(lgd, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
    xlim([min(SNR_rx_dB_plot), max(SNR_rx_dB_plot)]);
    set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.15);
end
sgt = sgtitle(sprintf('Scenario 1: 6GHz vs mmWave Topology', K));
set(sgt, 'Color', 'k', 'FontWeight', 'bold');
set(fig, 'InvertHardcopy', 'off', 'Color', 'w');

try
    save_figure(fig, 'fig_baseline_6GHz_vs_28GHz');
catch
end
plot_topology(dist_b, theta_b, dist_e, theta_e);

% =========================================================================
% FUNKCJE POMOCNICZE
% =========================================================================
function cdl = setup_matlab_cdl(Nt, fc, theta, band_tag)
    cdl = nrCDLChannel;
    cdl.DelayProfile = 'CDL-A';
    if fc < 10e9  
        cdl.DelaySpread = 92e-9;  
    else          
        cdl.DelaySpread = 30e-9;
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

function plot_topology(dist_b, theta_b_vec, dist_e, theta_e)
    fig_top = figure('Color', 'w', 'Position', [150 150 600 600]);
    x_bs = 0; y_bs = 0;
    
    ax = gca;
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.15);
    hold on; grid on; box on;
    
    p_bs = plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    text(x_bs, y_bs - 4, 'BS (0,0)', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k');
    
    % Rysowanie Ewy
    x_e = dist_e * sind(theta_e); y_e = dist_e * cosd(theta_e);
    p_e = plot(x_e, y_e, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Passive Eve');
    text(x_e + 2, y_e, sprintf('Eve\n(%gm, %g\\circ)', dist_e, theta_e), 'Color', 'r', 'FontSize', 10);
    
    % Rysowanie Bobów
    p_b = [];
    for k = 1:length(theta_b_vec)
        x_b = dist_b * sind(theta_b_vec(k)); 
        y_b = dist_b * cosd(theta_b_vec(k));
        p_b = plot(x_b, y_b, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
        text(x_b + 2, y_b, sprintf('B_{%d} (%g\\circ)', k, round(theta_b_vec(k))), 'Color', 'b', 'FontSize', 9);
    end
    
    % Przypisanie do legendy tylko raz dla Bobów
    if ~isempty(p_b)
        set(p_b, 'DisplayName', sprintf('Bobs (K=%d)', length(theta_b_vec)));
    end
    
    axis equal;
    max_d = max(dist_b, dist_e) + 15;
    xlim([-max_d, max_d]); ylim([-max_d/2, max_d]);
    xlabel('X [m]', 'Color', 'k'); ylabel('Y [m]', 'Color', 'k');
    title('Scenario 1: 6GHz vs mmWave Comparison', 'Color', 'k');
    
    lgd = legend([p_bs, p_b, p_e], 'Location', 'NorthEast');
    set(lgd, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
    
    set(fig_top, 'InvertHardcopy', 'off', 'Color', 'w');
    try
        save_figure(fig_top, '../topology/topology_baseline');
    catch
    end
end