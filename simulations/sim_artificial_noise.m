% =========================================================================
% SCENARIO: Artificial Noise (AN) injection in the null space of H'
% -------------------------------------------------------------------------
% Zaktualizowano do: nrCDLChannel, fizycznego FSPL, Transmit SNR
% oraz Ataku Kątowego (Ewy blisko Bobów).
% Wykres A: Odporność systemu na rosnącą liczbę Ew (przy stałym phi).
% Wykres B: Optymalizacja podziału mocy (phi) pomiędzy Dane a Szum.
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- PARAMETRY FIZYCZNE I SYSTEMOWE --------------------------------------
dist        = 50;                     % Dystans do użytkowników i Ew [m]
L_values    = 1:2:15;                 % Liczba współpracujących Ew (Sweep A)
phi_values  = 0.1:0.1:1.0;            % Frakcja mocy na dane (Sweep B)
L_fixed     = 7;                      % Stała liczba Ew dla Sweepu B
phi_fixed   = 0.7;                    % Stała frakcja phi dla Sweepu A
K           = 4;                      % Liczba legalnych użytkowników (Bobów)
numIter     = 20;                     % Liczba iteracji Monte Carlo

SNR_tx_dB   = 100;                    % Transmit SNR stacji bazowej [dB]
P_tx        = 10^(SNR_tx_dB / 10);
fc          = p.fc_sub6;              % Analizujemy pasmo 6 GHz
Nt          = p.Nt_sub6;              % 32 anteny (Massive MIMO)

[~, PL_dB] = compute_fspl(dist, fc);
fprintf('\n--- Sztuczny Szum (AN) @ 6 GHz ---\n');
fprintf('  Dystans: %g m (FSPL: %.2f dB)\n', dist, PL_dB);
fprintf('  Transmit SNR: %d dB\n', SNR_tx_dB);

% Wstępna alokacja obiektów 3GPP (aby uniknąć narzutu w pętli)
max_L = max([max(L_values), L_fixed]);
cdl_b = cell(K, 1);
for k = 1:K, cdl_b{k} = nrCDLChannel; end
cdl_e = cell(max_L, 1);
for e = 1:max_L, cdl_e{e} = nrCDLChannel; end

% =========================================================================
% --- SWEEP A: Wpływ liczby Ew (L) przy stałym podziale mocy (phi = 0.7)
% =========================================================================
fprintf('Rozpoczynam Sweep A (Liczba Ew)...\n');
SR_no_AN_A   = zeros(1, length(L_values));
SR_with_AN_A = zeros(1, length(L_values));

for l_idx = 1:length(L_values)
    num_eve = L_values(l_idx);
    acc_no = 0; acc_an = 0;
    
    for it = 1:numIter
        [SR_no, SR_an] = run_AN_trial(cdl_b, cdl_e, dist, fc, Nt, K, ...
                                      num_eve, P_tx, phi_fixed);
        acc_no = acc_no + SR_no;
        acc_an = acc_an + SR_an;
    end
    SR_no_AN_A(l_idx)   = acc_no / numIter;
    SR_with_AN_A(l_idx) = acc_an / numIter;
end

% =========================================================================
% --- SWEEP B: Wpływ podziału mocy (phi) przy stałej liczbie Ew (L = 7)
% =========================================================================
fprintf('Rozpoczynam Sweep B (Podział mocy phi)...\n');
SR_with_AN_B = zeros(1, length(phi_values));

for ph_idx = 1:length(phi_values)
    acc = 0;
    for it = 1:numIter
        [~, SR_an] = run_AN_trial(cdl_b, cdl_e, dist, fc, Nt, K, ...
                                  L_fixed, P_tx, phi_values(ph_idx));
        acc = acc + SR_an;
    end
    SR_with_AN_B(ph_idx) = acc / numIter;
end

% =========================================================================
% --- WIZUALIZACJA --------------------------------------------------------
% =========================================================================
fig = figure('Color', 'w', 'Position', [100 100 1100 450]);

subplot(1, 2, 1);
plot(L_values, SR_no_AN_A,   '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(L_values, SR_with_AN_A, '-g^', 'LineWidth', 2, 'MarkerFaceColor', 'g');
grid on; box on;
xlabel('Number of colluding eavesdroppers (L)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Robustness vs L  (\\phi = %.1f)', phi_fixed));
legend('Plain ZF', 'ZF + Artificial Noise', 'Location', 'SouthWest');

subplot(1, 2, 2);
plot(phi_values, SR_with_AN_B, '-mh', 'LineWidth', 2, 'MarkerFaceColor', 'm');
grid on; box on;
xlabel('Data-power fraction \phi (1.0 = No AN)');
ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Power-allocation trade-off  (L = %d Eves)', L_fixed));
xline(phi_fixed, 'k:', sprintf('\\phi = %.1f', phi_fixed), 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');

sgtitle(sprintf('Artificial Noise @ 6 GHz: Targeted Proximity Attack (Nt = %d, SNR_{tx} = %d dB)', Nt, SNR_tx_dB));
save_figure(fig, 'fig_artificial_noise');
plot_an_topology(dist, K, 8)

% =========================================================================
% FUNKCJA WYKONAWCZA: Generuje kanały i wylicza pojemność dla jednej iteracji
% =========================================================================
function [SR_no_AN, SR_with_AN] = run_AN_trial(cdl_b, cdl_e, dist, fc, Nt, K, num_eve, P_tx, phi)
    % 1. Fizyka tłumienia (FSPL)
    [PL_lin, ~] = compute_fspl(dist, fc);

    % 2. Geometria (Targeted Proximity Attack)
    % Ewy skradają się blisko Bobów (+/- 6 stopni)
    theta_bobs = -60 + 120*rand(1, K);
    target_bob_idx = randi(K, 1, num_eve);
    angular_error  = -6 + 12 * rand(1, num_eve);
    theta_eves     = theta_bobs(target_bob_idx) + angular_error;

    % 3. Generowanie kanałów efektywnych 3GPP (Bobowie)
    H_eff = zeros(Nt, K);
    for k = 1:K
        cdl_b{k} = setup_matlab_cdl(cdl_b{k}, Nt, fc, theta_bobs(k));
        cdl_b{k}.Seed = randi([0 2^31-1]);
        [pg_b, ~] = cdl_b{k}();
        hb = squeeze(sum(pg_b, 2)); hb = hb(:);
        H_eff(:, k) = sqrt((1 / PL_lin) * Nt) * (hb / norm(hb));
    end

    % 4. Generowanie kanałów efektywnych 3GPP (Ewy)
    G_eff = zeros(Nt, num_eve);
    for e = 1:num_eve
        cdl_e{e} = setup_matlab_cdl(cdl_e{e}, Nt, fc, theta_eves(e));
        cdl_e{e}.Seed = randi([0 2^31-1]);
        [pg_e, ~] = cdl_e{e}();
        he = squeeze(sum(pg_e, 2)); he = he(:);
        G_eff(:, e) = sqrt((1 / PL_lin) * Nt) * (he / norm(he));
    end

    % 5. Prekodowanie ZF i Baza Przestrzeni Zerowej (Null-space)
    W_raw = H_eff * pinv(H_eff' * H_eff);
    W = W_raw / norm(W_raw, 'fro');
    Z = null(H_eff'); % Wektory ortogonalne do Bobów (tutaj ładujemy szum)

    % 6. Macierze Kowariancji Zakłóceń dla Ew (Optymalny odbiornik MRC)
    % Tryb 1: Brak AN (tylko naturalny szum termiczny znormalizowany do 1)
    R_no_eve_inv = eye(num_eve); 
    
    % Tryb 2: ZF + AN
    P_AN = (1 - phi) * P_tx; % Moc przeznaczona na sztuczny szum
    Q_AN = (P_AN / (Nt - K)) * (Z * Z'); % Macierz kowariancji nadawanego szumu
    R_an_eve = G_eff' * Q_AN * G_eff + eye(num_eve); % Szum odebrany + termiczny
    R_an_eve_inv = inv(R_an_eve);

    R_b_no = zeros(K, 1); R_e_no = zeros(K, 1);
    R_b_an = zeros(K, 1); R_e_an = zeros(K, 1);

    % 7. Obliczenia przepływności (Shannon Capacity)
    for k = 1:K
        % --- Odbiór u Boba ---
        % Bob bez AN (phi = 1.0)
        S_b_no = P_tx * abs(H_eff(:, k)' * W(:, k))^2;
        I_b_no = 0;
        for j = 1:K, if j ~= k, I_b_no = I_b_no + P_tx * abs(H_eff(:, k)' * W(:, j))^2; end; end
        R_b_no(k) = log2(1 + S_b_no / (I_b_no + 1));

        % Bob z AN (Moc użyteczna obcięta do phi * P_tx)
        S_b_an = phi * P_tx * abs(H_eff(:, k)' * W(:, k))^2;
        I_b_an = 0;
        for j = 1:K, if j ~= k, I_b_an = I_b_an + phi * P_tx * abs(H_eff(:, k)' * W(:, j))^2; end; end
        R_b_an(k) = log2(1 + S_b_an / (I_b_an + 1)); % Bob nie odbiera AN, bo Z jest ortogonalne

        % --- Odbiór u Ewy (Kooperatywne MRC) ---
        % Ewa bez AN
        h_ek_no = sqrt(P_tx) * G_eff' * W(:, k);
        R_e_no(k) = log2(1 + real(h_ek_no' * R_no_eve_inv * h_ek_no));

        % Ewa zaatakowana przez AN (Używa optymalnego filtra R_inv)
        h_ek_an = sqrt(phi * P_tx) * G_eff' * W(:, k);
        R_e_an(k) = log2(1 + real(h_ek_an' * R_an_eve_inv * h_ek_an));
    end

    SR_no_AN   = sum(max(0, R_b_no - R_e_no));
    SR_with_AN = sum(max(0, R_b_an - R_e_an));
end

% =========================================================================
% FUNKCJA POMOCNICZA: Szybka rekonfiguracja kątów
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

% =========================================================================
% FUNKCJA POMOCNICZA: Generowanie topologii scenariusza (Targeted Proximity)
% =========================================================================
function plot_an_topology(dist, K, num_eve)
    fig_top = figure('Color', 'w', 'Position', [150 150 700 700]);
    hold on; grid on; box on;
    
    % --- Odtworzenie logiki rozmieszczenia ze skryptu ---
    % Bobs są rozrzuceni w sektorze [-60, 60]
    theta_bobs = -60 + 120 * rand(1, K);
    
    % Ewy celują w losowych Bobów i stoją bardzo blisko nich (+/- 6 stopni)
    target_bob_idx = randi(K, 1, num_eve);
    angular_error  = -6 + 12 * rand(1, num_eve);
    theta_eves     = theta_bobs(target_bob_idx) + angular_error;
    
    % Konwersja na kartezjańskie (BS w 0,0)
    x_bs = 0; y_bs = 0;
    max_d = dist + 15;
    
    % Rysowanie BS
    p_bs = plot(x_bs, y_bs, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Base Station (BS)');
    text(x_bs, y_bs - 3, 'BS (0,0)', 'HorizontalAlignment', 'center', 'Color', 'k');
    
    % Rysowanie Ew (Targeted Proximity)
    x_eves = dist * sind(theta_eves);
    y_eves = dist * cosd(theta_eves);
    
    p_e = [];
    for e = 1:num_eve
        p_e = plot(x_eves(e), y_eves(e), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        text(x_eves(e) + 1.5, y_eves(e) + 1.5, sprintf('E_{%d}', e), 'Color', 'r', 'FontSize', 9);
    end
    if ~isempty(p_e)
        set(p_e, 'DisplayName', sprintf('Eves (L=%d)', num_eve));
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
    title(sprintf('Scenario 4: Artificial Noise (K=%d, L=%d, d=%gm)', K, num_eve, dist));
    legend([p_bs, p_b, p_e], 'Location', 'NorthWest');
    
    % Zapis do pliku
    try
        save_figure(fig_top, '../topology/topology_artificial_noise');
    catch
        warning('Funkcja save_figure nie jest dostępna.');
    end
end