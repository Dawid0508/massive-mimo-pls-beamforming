% =========================================================================
% SCENARIO: 6 GHz vs 28 GHz under targeted/colluding Eves (Spatial Proximity)
% -------------------------------------------------------------------------
% Eves cluster tightly around Bobs (within +/- 6 deg) to exploit beam width.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- PARAMETRY FIZYCZNE I SYSTEMOWE ---
dist        = 50;                     % Dystans do użytkowników i podsłuchiwaczy (m)
L_values    = 1:2:15;                 % Liczba współpracujących Ew (L)
K           = 4;                      % Liczba legalnych użytkowników (Bobów)
numIter     = 10;                     
SNR_tx_dB   = 100;                    % Transmit SNR stacji bazowej (dB)
P_tx        = 10^(SNR_tx_dB / 10);

bands = struct( ...
    'name', {'6 GHz (Massive MIMO)', '28 GHz (Ultra-Massive MIMO)'}, ...
    'fc',   {p.fc_sub6,              p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6,              p.Nt_mmwave}, ...
    'cdl',  {p.cdl_sub6,             p.cdl_mmwave});

results_SR    = zeros(2, length(L_values));
results_Fair  = zeros(2, length(L_values));

max_L = max(L_values);

for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;  cdl_tag = bands(b).cdl;
    
    % Obliczanie bezwzględnego tłumienia przestrzennego (FSPL)
    [PL_lin, PL_dB] = compute_fspl(dist, fc);
    
    fprintf('\n--- Scenariusz: Zacieśniony Atak Kątowy Ew @ %s ---\n', bands(b).name);
    fprintf('  Liczba anten stacji (Nt): %d\n', Nt);
    fprintf('  Dystans: %g m, Tłumienie trasy: %.2f dB\n', dist, PL_dB);
    
    % Wstępna alokacja obiektów kanałowych
    cdl_b = cell(K, 1);
    for k = 1:K, cdl_b{k} = nrCDLChannel; end
    cdl_e = cell(max_L, 1);
    for e = 1:max_L, cdl_e{e} = nrCDLChannel; end

    for l_idx = 1:length(L_values)
        num_eve = L_values(l_idx);
        SR_acc = 0; F_acc = 0;
        
        for it = 1:numIter
            % Bobowie rozproszeni losowo w przestrzeni sektora
            theta_bobs = -60 + 120*rand(1, K);
            
            % --- KLUCZOWA MODYFIKACJA: ATAK GEOMETRYCZNY ---
            % Ewy nie stoją już losowo. Każda Ewa wybiera jednego z Bobów 
            % i podkrada się ekstremalnie blisko pod jego kąt (w zakresie +/- 6 stopni)
            target_bob_idx = randi(K, 1, num_eve); 
            angular_error  = 0 + 12 * rand(1, num_eve); % Odchyłka od wprost w Boba
            theta_eves     = theta_bobs(target_bob_idx) + angular_error;
            
            % Budowanie kanału efektywnego dla Bobów
            H_eff = zeros(Nt, K);
            for k = 1:K
                cdl_b{k} = setup_matlab_cdl(cdl_b{k}, Nt, fc, theta_bobs(k));
                cdl_b{k}.Seed = randi([0 2^31-1]);
                [pg_b, ~] = cdl_b{k}();
                hb = squeeze(sum(pg_b, 2)); hb = hb(:);
                H_eff(:, k) = sqrt((1 / PL_lin) * Nt) * (hb / norm(hb));
            end
            
            % Budowanie kanału efektywnego dla Ew (śledzących Bobów)
            G_eff = zeros(Nt, num_eve);
            for e = 1:num_eve
                cdl_e{e} = setup_matlab_cdl(cdl_e{e}, Nt, fc, theta_eves(e));
                cdl_e{e}.Seed = randi([0 2^31-1]);
                [pg_e, ~] = cdl_e{e}();
                he = squeeze(sum(pg_e, 2)); he = he(:);
                G_eff(:, e) = sqrt((1 / PL_lin) * Nt) * (he / norm(he));
            end
            
            % Wyliczanie prekodera Zero-Forcing (ZF) z normalizacją macierzową
            W = H_eff * pinv(H_eff' * H_eff);
            W = W / norm(W, 'fro');
            
            R_b = zeros(K, 1);
            R_e = zeros(K, 1);
            
            for k = 1:K
                % Legalny użytkownik (Bob)
                S_b = P_tx * abs(H_eff(:, k)' * W(:, k))^2;
                I_b = 0;
                for j = 1:K
                    if j ~= k
                        I_b = I_b + P_tx * abs(H_eff(:, k)' * W(:, j))^2;
                    end
                end
                R_b(k) = log2(1 + S_b / (I_b + 1));
                
                % Współpracujące Ewy (Worst-Case MRC Attack)
                S_e = P_tx * sum(abs(G_eff' * W(:, k)).^2);
                R_e(k) = log2(1 + S_e / 1); 
            end
            
            % Obliczanie Secrecy Rate i indeksu Fairness Jaina
            R_s = max(0, R_b - R_e);
            SR_acc = SR_acc + sum(R_s);
            
            if sum(R_s.^2) == 0
                jains_val = 0;
            else
                jains_val = (sum(R_s))^2 / (K * sum(R_s.^2));
            end
            F_acc = F_acc + jains_val;
        end
        
        results_SR(b, l_idx)   = SR_acc / numIter;
        results_Fair(b, l_idx) = F_acc  / numIter;
    end
end

% --- WIZUALIZACJA METRYK ---
fig = figure('Color', 'w', 'Position', [100 100 1100 420]);

subplot(1, 2, 1);
plot(L_values, results_SR(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, results_SR(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Number of colluding eavesdroppers (L)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Secrecy Rate under Close-Range Attack');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

subplot(1, 2, 2);
plot(L_values, results_Fair(1,:), '--bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, results_Fair(2,:), '--rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
ylim([0 1.05]);
xlabel('Number of colluding eavesdroppers (L)');
ylabel("Jain's fairness index");
title('Fairness across Bobs');
legend(bands(1).name, bands(2).name, 'Location', 'SouthWest');

sgtitle(sprintf('6 GHz vs 28 GHz: Targeted Angular Proximity Attack (\\Delta\\theta \\le 6^\\circ, K = %d)', K));
save_figure(fig, 'fig_colluding_eavesdroppers_proximity');

% =========================================================================
% FUNKCJA POMOCNICZA: Dynamiczna re-konfiguracja istniejącego kanału 3GPP
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