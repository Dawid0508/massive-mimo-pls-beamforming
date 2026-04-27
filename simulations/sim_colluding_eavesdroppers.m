%% 6G PLS: Colluding Eavesdroppers vs Band Frequency (6GHz vs 28GHz)
clear; clc; close all;

% --- Parametry Konfiguracyjne (zgodnie z plikiem kolegi) ---
dist = 30;                  % Dystans w metrach [cite: 2903]
L_values = 1:2:15;          % Liczba współpracujących podsłuchiwaczy (L) [cite: 2881, 2883]
K = 4;                      % Liczba legalnych użytkowników (Bobs)
numIter = 50;               % Liczba powtórzeń Monte Carlo
c = physconst('LightSpeed'); % Prędkość światła [cite: 2904]

% Częstotliwości i anteny (Ultra-Massive MIMO)
fc_6 = 6e9;   Nt_6 = 32;    % 6 GHz [cite: 2906]
fc_mm = 28e9; Nt_mm = 512;  % 28 GHz [cite: 2950]

% --- Obliczanie Tłumienia (FSPL) [cite: 2905] ---
PL_6_lin = 10^((20*log10(dist) + 20*log10(fc_6) - 147.55)/10);
PL_mm_lin = 10^((20*log10(dist) + 20*log10(fc_mm) - 147.55)/10);

% Inicjalizacja tablic na wyniki
results_SR = zeros(2, length(L_values)); % [1: 6GHz, 2: 28GHz]
results_Fairness = zeros(2, length(L_values));

% --- Pętla Symulacji ---
for f_idx = 1:2
    if f_idx == 1
        fc = fc_6; Nt = Nt_6; L_path = PL_6_lin;
    else
        fc = fc_mm; Nt = Nt_mm; L_path = PL_mm_lin;
    end

    % Definicja szyku antenowego 
    ula = phased.ULA('NumElements', Nt, 'ElementSpacing', c/fc/2);
    sv = phased.SteeringVector('SensorArray', ula);

    for l_idx = 1:length(L_values)
        num_eve = L_values(l_idx);
        temp_SR = zeros(numIter, 1);
        temp_Fair = zeros(numIter, 1);

        for iter = 1:numIter
            % 1. Kanały Bobów (różne kąty)
            theta_bobs = -60 + 120*rand(1, K);
            H = step(sv, fc, theta_bobs); % <--- USUNIĘTO TRANSPOZYCJĘ (')

            % 2. Kanały Podsłuchiwaczy (rozsiani losowo)
            theta_eves = -90 + 180*rand(1, num_eve);
            G = step(sv, fc, theta_eves); % <--- USUNIĘTO TRANSPOZYCJĘ (')

            % 3. Precoding Zero-Forcing (ZF)
            W = H * pinv(H' * H);
            W = W / norm(W, 'fro'); % Normalizacja macierzowa

            % 4. Obliczanie SNR u Bobów (uwzględniając tłumienie)
            SNR_transmit = 1000; % Przykładowa moc
            P_eff = SNR_transmit / L_path; 
            
            R_bobs = zeros(K, 1);
            for k = 1:K
                signal = P_eff * abs(H(:,k)' * W(:,k))^2;
                R_bobs(k) = log2(1 + signal);
            end
            
            % 5. Atak: Współpracujący Podsłuchiwacze (MRC) - INDYWIDUALNIE
            R_eve = zeros(K, 1);
            for k = 1:K
                snr_eve_k = 0;
                for e = 1:num_eve
                    % Podsłuchiwacze (Eve) sumują siłę sygnału tylko z wiązki k-tego Boba
                    leakage_k = P_eff * abs(G(:,e)' * W(:,k))^2;
                    snr_eve_k = snr_eve_k + leakage_k;
                end
                R_eve(k) = log2(1 + snr_eve_k);
            end
            
            % 6. Secrecy Rate i Sprawiedliwość Jaina
            SR_vector = max(0, R_bobs - R_eve);
            
            temp_SR(iter) = sum(SR_vector);
            
            % Zabezpieczenie przed dzieleniem przez zero, gdy SR_vector to same zera
            if sum(SR_vector) == 0
                temp_Fair(iter) = 0;
            else
                temp_Fair(iter) = (sum(SR_vector)^2) / (K * sum(SR_vector.^2));
            end
        end
        results_SR(f_idx, l_idx) = mean(temp_SR);
        results_Fairness(f_idx, l_idx) = mean(temp_Fair, 'omitnan');
    end
end

% --- Profesjonalna Wizualizacja ---
figure( 'Position', [100 100 1000 400]);
subplot(1,2,1);
plot(L_values, results_SR(1,:), '-bo', L_values, results_SR(2,:), '-rs', 'LineWidth', 2);
grid on; xlabel('Liczba podsłuchiwaczy (L)'); ylabel('Secrecy Sum-Rate (bps/Hz)');
legend('6 GHz (Massive MIMO)', '28 GHz (Ultra-Massive)', 'Location', 'SouthWest');
title('Wydajność Bezpieczeństwa');

subplot(1,2,2);
plot(L_values, results_Fairness(1,:), '--bo', L_values, results_Fairness(2,:), '--rs', 'LineWidth', 2);
grid on; xlabel('Liczba podsłuchiwaczy (L)'); ylabel('Jain''s Fairness Index');
title('Sprawiedliwość Ochrony');

% --- Zapisywanie wyników ---
% Zapisuje wykres pokazujący różnicę między 6GHz a 28GHz przy ataku grupy Eve
saveas(gcf, '../results/wykres_6GHz_vs_28GHz.png');