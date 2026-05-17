% =========================================================================
% SCENARIO: 6 GHz vs 28 GHz under Targeted Attack with ARTIFICIAL NOISE
% -------------------------------------------------------------------------
% Porównanie odporności systemu z włączonym i wyłączonym Sztucznym Szumem.
% BS dzieli moc: 70% na dane, 30% na zagłuszanie (AN) w null-space Bobów.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- PARAMETRY FIZYCZNE I SYSTEMOWE ---
dist        = 50;                     
L_values    = 1:2:15;                 
K           = 4;                      
numIter     = 10;                     
SNR_tx_dB   = 100;                    
P_tx        = 10^(SNR_tx_dB / 10);

% Parametry Sztucznego Szumu
phi_AN      = 0.7; % 70% mocy na przesył danych, 30% na Sztuczny Szum

bands = struct( ...
    'name', {'6 GHz', '28 GHz'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave}, ...
    'cdl',  {p.cdl_sub6, p.cdl_mmwave});

% Macierze wynikowe: (Pasmo, Tryb_AN, L_values)
results_SR    = zeros(2, 2, length(L_values));
results_Fair  = zeros(2, 2, length(L_values));

max_L = max(L_values);

for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;  
    [PL_lin, PL_dB] = compute_fspl(dist, fc);
    
    fprintf('\n--- Symulacja Artificial Noise @ %s ---\n', bands(b).name);
    
    cdl_b = cell(K, 1);
    for k = 1:K, cdl_b{k} = nrCDLChannel; end
    cdl_e = cell(max_L, 1);
    for e = 1:max_L, cdl_e{e} = nrCDLChannel; end

    for l_idx = 1:length(L_values)
        num_eve = L_values(l_idx);
        % Akumulatory: [Bez AN, Z AN]
        SR_acc = zeros(1, 2); 
        F_acc  = zeros(1, 2);
        
        for it = 1:numIter
            theta_bobs = -60 + 120*rand(1, K);
            
            % Atak z bliska (Ewy skradają się pod kąty Bobów)
            target_bob_idx = randi(K, 1, num_eve); 
            angular_error  = -6 + 12 * rand(1, num_eve); 
            theta_eves     = theta_bobs(target_bob_idx) + angular_error;
            
            H_eff = zeros(Nt, K);
            for k = 1:K
                cdl_b{k} = setup_matlab_cdl(cdl_b{k}, Nt, fc, theta_bobs(k));
                cdl_b{k}.Seed = randi([0 2^31-1]);
                [pg_b, ~] = cdl_b{k}();
                hb = squeeze(sum(pg_b, 2)); hb = hb(:);
                H_eff(:, k) = sqrt((1 / PL_lin) * Nt) * (hb / norm(hb));
            end
            
            G_eff = zeros(Nt, num_eve);
            for e = 1:num_eve
                cdl_e{e} = setup_matlab_cdl(cdl_e{e}, Nt, fc, theta_eves(e));
                cdl_e{e}.Seed = randi([0 2^31-1]);
                [pg_e, ~] = cdl_e{e}();
                he = squeeze(sum(pg_e, 2)); he = he(:);
                G_eff(:, e) = sqrt((1 / PL_lin) * Nt) * (he / norm(he));
            end
            
            % Prekoder ZF na dane
            W_raw = H_eff * pinv(H_eff' * H_eff);
            W = W_raw / norm(W_raw, 'fro');
            
            % --- BAZA PRZESTRZENI ZEROWEJ DLA AN ---
            % Z znajduje wektory ortogonalne do kanałów Bobów
            Z = null(H_eff'); 
            
            % --- PĘTLA TRYBÓW (1: Bez AN, 2: Z AN) ---
            for an_mode = 1:2
                if an_mode == 1
                    phi = 1.0;          % 100% mocy na dane
                    Q_AN = zeros(Nt, Nt); % Brak szumu
                else
                    phi = phi_AN;       % 70% mocy na dane
                    P_AN = (1 - phi) * P_tx; % 30% mocy na szum
                    % Kowariancja szumu rozłożona równo w null-space
                    Q_AN = (P_AN / (Nt - K)) * (Z * Z'); 
                end
                
                R_b = zeros(K, 1);
                R_e = zeros(K, 1);
                
                % Macierz kowariancji zakłóceń u Ew (Szum własny + AN stacji bazowej)
                R_AN_eve = G_eff' * Q_AN * G_eff + eye(num_eve);
                R_AN_eve_inv = inv(R_AN_eve);
                
                for k = 1:K
                    % Pojemność Boba (Sygnał przeskalowany przez phi)
                    S_b = phi * P_tx * abs(H_eff(:, k)' * W(:, k))^2;
                    I_b = 0;
                    for j = 1:K
                        if j ~= k
                            I_b = I_b + phi * P_tx * abs(H_eff(:, k)' * W(:, j))^2;
                        end
                    end
                    R_b(k) = log2(1 + S_b / (I_b + 1)); % Bob nie odbiera AN
                    
                    % Pojemność Ewy pod atakiem MRC (Uwzględnia kolorowy szum AN)
                    h_ek = sqrt(phi * P_tx) * G_eff' * W(:, k);
                    R_e(k) = log2(1 + real(h_ek' * R_AN_eve_inv * h_ek));
                end
                
                % Zapis wyników dla danego trybu
                R_s = max(0, R_b - R_e);
                SR_acc(an_mode) = SR_acc(an_mode) + sum(R_s);
                if sum(R_s.^2) == 0
                    F_acc(an_mode) = F_acc(an_mode) + 0;
                else
                    F_acc(an_mode) = F_acc(an_mode) + (sum(R_s))^2 / (K * sum(R_s.^2));
                end
            end
        end
        
        results_SR(b, :, l_idx)   = SR_acc / numIter;
        results_Fair(b, :, l_idx) = F_acc  / numIter;
    end
end

% --- WIZUALIZACJA ---
fig = figure('Color', 'w', 'Position', [100 100 1100 450]);

subplot(1, 2, 1);
plot(L_values, squeeze(results_SR(1,1,:)), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, squeeze(results_SR(1,2,:)), '--b^', 'LineWidth', 2, 'MarkerFaceColor', 'b');
plot(L_values, squeeze(results_SR(2,1,:)), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
plot(L_values, squeeze(results_SR(2,2,:)), '--r^', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Number of colluding eavesdroppers (L)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Impact of Artificial Noise (AN)');
legend('6 GHz (No AN)', '6 GHz (With 30% AN)', '28 GHz (No AN)', '28 GHz (With 30% AN)', 'Location', 'NorthEast');

subplot(1, 2, 2);
plot(L_values, squeeze(results_Fair(1,1,:)), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, squeeze(results_Fair(1,2,:)), '--b^', 'LineWidth', 2, 'MarkerFaceColor', 'b');
plot(L_values, squeeze(results_Fair(2,1,:)), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
plot(L_values, squeeze(results_Fair(2,2,:)), '--r^', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Number of colluding eavesdroppers (L)');
ylabel("Jain's fairness index");
title('Fairness across Bobs');
legend('6 GHz (No AN)', '6 GHz (With AN)', '28 GHz (No AN)', '28 GHz (With AN)', 'Location', 'SouthWest');

sgtitle(sprintf('Targeted Attack Countermeasures: Baseline vs Artificial Noise (K = %d)', K));
save_figure(fig, 'fig_artificial_noise_comparison');

% =========================================================================
% FUNKCJA POMOCNICZA
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