% =========================================================================
% SCENARIO: Spatial Correlation & MRT vs ZF Precoding
% -------------------------------------------------------------------------
% TOPOLOGIA: Ewa stoi na środku (0 stopni). 4 użytkowników początkowo 
% jest szeroko rozstawionych (niska korelacja), ale z czasem zbijają się 
% w ciasny tłum wokół Ewy (ekstremalna korelacja przestrzenna).
%
% OBSERWUJEMY:
% 1. Upadek algorytmu ZF (Noise Enhancement) przy wysokiej korelacji.
% 2. Skok pojemności Ewy, gdy wiązki Bobów nakładają się na jej pozycję.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(123); % Ustawiony stały seed dla idealnie powtarzalnych wykresów

% --- KONFIGURACJA SYSTEMU ---
K           = 4;                             % Liczba użytkowników
spread_vec  = 120:-10:10;                    % Szerokość tłumu [stopnie]
numIter     = 80;                            % Iteracje kanału (CDL-D)
SNR_tx_dB   = 100;                           
P_tx        = 10^(SNR_tx_dB / 10);
noise_var   = 1;

dist_b      = 50;                            % Bobowie są w równej odległości
dist_e      = 50;                            
theta_eve   = 0;                             % EWA CZEKA NA ŚRODKU!

bands = struct('name', {'6 GHz (Nt=32)', '28 GHz (Nt=256)'}, 'fc', {p.fc_sub6, p.fc_mmwave}, 'Nt', {32, 256});

% Macierze wynikowe: [Pasmo x Algorytm(1=MRT, 2=ZF) x Spread]
R_bob_sum = zeros(2, 2, length(spread_vec));
R_eve_sum = zeros(2, 2, length(spread_vec));
R_sec_sum = zeros(2, 2, length(spread_vec));

% Modele CDL (Line-of-Sight, żeby kąty miały znaczenie!)
cdl_b = cell(K, 1);
for k = 1:K, cdl_b{k} = setup_matlab_cdl(nrCDLChannel, 64, bands(1).fc, 0); end
cdl_e = setup_matlab_cdl(nrCDLChannel, 64, bands(1).fc, theta_eve);

% --- PĘTLA GŁÓWNA ---
for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;  
    [PL_lin_b, ~] = compute_fspl(dist_b, fc);
    [PL_lin_e, ~] = compute_fspl(dist_e, fc);
    
    for k = 1:K, cdl_b{k} = setup_matlab_cdl(cdl_b{k}, Nt, fc, 0); end
    cdl_e = setup_matlab_cdl(cdl_e, Nt, fc, theta_eve);
    
    fprintf('Obliczam pasmo %s...\n', bands(b).name);
    
    for s_idx = 1:length(spread_vec)
        spread = spread_vec(s_idx);
        theta_bobs = linspace(-spread/2, spread/2, K); % Rozstawienie symetryczne wokół 0
        
        acc_b = zeros(1, 2); acc_e = zeros(1, 2); acc_s = zeros(1, 2);
        
        for it = 1:numIter
            % 1. Generowanie fizycznych kanałów
            H_eff = zeros(Nt, K);
            for k = 1:K
                cdl_b{k}.TransmitArrayOrientation = [-theta_bobs(k); 0; 0];
                release(cdl_b{k}); cdl_b{k}.Seed = randi([0 2^31-1]);
                [pg_b, ~] = cdl_b{k}(); hb = squeeze(sum(pg_b, 2));
                H_eff(:, k) = sqrt((1 / PL_lin_b) * Nt) * (hb(:) / norm(hb(:)));
            end
            
            release(cdl_e); cdl_e.Seed = randi([0 2^31-1]);
            [pg_e, ~] = cdl_e(); he = squeeze(sum(pg_e, 2));
            h_eff_e = sqrt((1 / PL_lin_e) * Nt) * (he(:) / norm(he(:)));
            
            % 2. Prekodery: MRT vs ZF
            W_mrt_raw = H_eff;
            W_mrt = (W_mrt_raw / norm(W_mrt_raw, 'fro')) * sqrt(P_tx);
            
            W_zf_raw = H_eff * pinv(H_eff' * H_eff + 1e-9*eye(K)); % pinv chroni przed NaN
            W_zf = (W_zf_raw / norm(W_zf_raw, 'fro')) * sqrt(P_tx);
            
            W_all = {W_mrt, W_zf};
            
            % 3. Ewaluacja pojemności
            for alg = 1:2
                W = W_all{alg};
                sum_rb = 0; sum_re = 0; sum_rs = 0;
                
                for k = 1:K
                    % Pojemność dla Boba
                    sig_b  = abs(H_eff(:,k)' * W(:,k))^2;
                    intf_b = sum(abs(H_eff(:,k)' * W).^2) - sig_b;
                    R_b_k  = log2(1 + sig_b / (intf_b + noise_var));
                    
                    % Ewa "wysysa" dane z wiązki k
                    sig_e  = abs(h_eff_e' * W(:,k))^2;
                    intf_e = sum(abs(h_eff_e' * W).^2) - sig_e;
                    R_e_k  = log2(1 + sig_e / (intf_e + noise_var));
                    
                    sum_rb = sum_rb + R_b_k;
                    sum_re = sum_re + R_e_k;
                    sum_rs = sum_rs + max(0, R_b_k - R_e_k);
                end
                acc_b(alg) = acc_b(alg) + sum_rb;
                acc_e(alg) = acc_e(alg) + sum_re;
                acc_s(alg) = acc_s(alg) + sum_rs;
            end
        end
        R_bob_sum(b, :, s_idx) = acc_b / numIter;
        R_eve_sum(b, :, s_idx) = acc_e / numIter;
        R_sec_sum(b, :, s_idx) = acc_s / numIter;
    end
end

% --- WIZUALIZACJA ---
fig = figure('Color', 'w', 'Position', [100 50 1400 800]);

for b = 1:2
    % Wykres 1: Total Bob Capacity (Co potrafi stacja?)
    subplot(3, 2, b);
    plot(spread_vec, squeeze(R_bob_sum(b, 1, :)), '-k^', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
    plot(spread_vec, squeeze(R_bob_sum(b, 2, :)), '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
    grid on; box on; set(gca, 'XDir', 'reverse');
    ylabel('Bob Sum-Rate'); title(sprintf('%s: Network Capacity', bands(b).name));
    if b==1, legend('MRT (Max Signal)', 'ZF (Zero Interference)', 'Location', 'SouthWest'); end
    
    % Wykres 2: Eve Capacity (Kradzież danych)
    subplot(3, 2, b + 2);
    plot(spread_vec, squeeze(R_eve_sum(b, 1, :)), '--k^', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
    plot(spread_vec, squeeze(R_eve_sum(b, 2, :)), '--ro', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    grid on; box on; set(gca, 'XDir', 'reverse');
    ylabel('Eve Sum-Rate'); title('Information Leakage (Eavesdropping)');
    
    % Wykres 3: SECRECY RATE (Finał - Sieć chroniona)
    subplot(3, 2, b + 4);
    plot(spread_vec, squeeze(R_sec_sum(b, 1, :)), '-k^', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
    plot(spread_vec, squeeze(R_sec_sum(b, 2, :)), '-bo', 'LineWidth', 2.5, 'MarkerFaceColor', 'b');
    grid on; box on; set(gca, 'XDir', 'reverse');
    xlabel('Angular Spread (degrees) -> Users crowding'); ylabel('Secrecy Rate');
    title('Final Secure Capacity');
end
sgtitle('The Spatial Correlation Trap: MRT vs ZF Precoding (Eve hiding at 0^{\circ})');
save_figure(fig, 'fig_spatial_mrt_zf');

function cdl = setup_matlab_cdl(cdl, Nt, fc, theta)
    release(cdl);
    cdl.DelayProfile = 'CDL-D'; % Czysty Line-of-Sight
    if fc < 10e9, cdl.DelaySpread = 30e-9; else, cdl.DelaySpread = 10e-9; end
    cdl.CarrierFrequency = fc;
    cdl.MaximumDopplerShift = 0;         
    cdl.TransmitAntennaArray.Size = [1 Nt 1 1 1]; 
    cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1]; 
    cdl.TransmitArrayOrientation = [-theta; 0; 0];
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
    cdl.NumTimeSamples = 1;
    cdl.ChannelFiltering = false; 
end