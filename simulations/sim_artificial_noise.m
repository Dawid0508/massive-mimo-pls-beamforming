% =========================================================================
% SYMULACJA 3: Ochrona pasma 6 GHz za pomocą Sztucznego Szumu (Artificial Noise)
% =========================================================================
clear; clc; close all;

% --- Parametry Konfiguracyjne ---
dist = 30;                  % Dystans w metrach
L_values = 1:2:15;          % Liczba współpracujących podsłuchiwaczy (L)
K = 4;                      % Liczba legalnych użytkowników (Bobs)
numIter = 50;               % Pętle Monte Carlo
c = physconst('LightSpeed'); 

% Parametry tylko dla pasma 6 GHz (bo tu mieliśmy problem)
fc = 6e9;   Nt = 32;    
PL_lin = 10^((20*log10(dist) + 20*log10(fc) - 147.55)/10);

% Alokacja Mocy
SNR_received_dB = 30; 
P_eff = 10^(SNR_received_dB/10); 
phi = 0.7; % Ułamek mocy na DANE (70%), reszta (30%) idzie na SZUM (AN)

% Inicjalizacja tablic
results_SR_no_AN = zeros(1, length(L_values));
results_SR_with_AN = zeros(1, length(L_values));

% Definicja szyku antenowego
ula = phased.ULA('NumElements', Nt, 'ElementSpacing', c/fc/2);
sv = phased.SteeringVector('SensorArray', ula);

for l_idx = 1:length(L_values)
    num_eve = L_values(l_idx);
    temp_SR_no_AN = zeros(numIter, 1);
    temp_SR_with_AN = zeros(numIter, 1);
    
    for iter = 1:numIter
        % 1. Kanały
        theta_bobs = -60 + 120*rand(1, K);
        H = step(sv, fc, theta_bobs); % Bobowie
        
        theta_eves = -90 + 180*rand(1, num_eve);
        G = step(sv, fc, theta_eves); % Podsłuchiwacze
        
        % 2. Precoding ZF dla Bobów
        W = H * pinv(H' * H);
        W = W / norm(W, 'fro'); 
        
        % 3. Generowanie Przestrzeni Zerowej (Null Space) dla Sztucznego Szumu
        % Funkcja null() znajduje wektory, które po pomnożeniu przez H dają 0
        Z = null(H'); 
        num_null_dims = size(Z, 2); % Liczba wymiarów dla szumu
        
        % 4. Obliczanie SNR dla Bobów
        R_bobs_no_AN = zeros(K, 1);
        R_bobs_with_AN = zeros(K, 1);
        
        for k = 1:K
            % Bez AN (100% mocy na dane)
            sig_no_AN = P_eff * abs(H(:,k)' * W(:,k))^2;
            R_bobs_no_AN(k) = log2(1 + sig_no_AN);
            
            % Z AN (70% mocy na dane, 30% na szum). Bob NIE SŁYSZY szumu!
            sig_with_AN = (phi * P_eff) * abs(H(:,k)' * W(:,k))^2;
            R_bobs_with_AN(k) = log2(1 + sig_with_AN);
        end
        
        % 5. Atak Podsłuchiwaczy (MRC)
        R_eve_no_AN = zeros(K, 1);
        R_eve_with_AN = zeros(K, 1);
        
        for k = 1:K
            snr_eve_no_AN = 0;
            snr_eve_with_AN = 0;
            
            for e = 1:num_eve
                % Wyciek samych danych (dla obu przypadków)
                leakage_data_no_AN = P_eff * abs(G(:,e)' * W(:,k))^2;
                leakage_data_with_AN = (phi * P_eff) * abs(G(:,e)' * W(:,k))^2;
                
                % EVE SŁYSZY SZUM! Uderza on w jej mianownik (zakłóca odbiór)
                % Szum rozkłada się równo na wszystkie wymiary przestrzeni zerowej
                interference_AN = ((1 - phi) * P_eff / num_null_dims) * norm(G(:,e)' * Z)^2;
                
                % MRC sumuje SINR z każdego podsłuchiwacza
                snr_eve_no_AN = snr_eve_no_AN + leakage_data_no_AN;
                snr_eve_with_AN = snr_eve_with_AN + (leakage_data_with_AN / (1 + interference_AN));
            end
            
            R_eve_no_AN(k) = log2(1 + snr_eve_no_AN);
            R_eve_with_AN(k) = log2(1 + snr_eve_with_AN);
        end
        
        % 6. Secrecy Rate
        temp_SR_no_AN(iter) = sum(max(0, R_bobs_no_AN - R_eve_no_AN));
        temp_SR_with_AN(iter) = sum(max(0, R_bobs_with_AN - R_eve_with_AN));
    end
    results_SR_no_AN(l_idx) = mean(temp_SR_no_AN);
    results_SR_with_AN(l_idx) = mean(temp_SR_with_AN);
end

% --- Wykres ---
%--figure('Color', 'w');
plot(L_values, results_SR_no_AN, '-bo', 'LineWidth', 2); hold on;
plot(L_values, results_SR_with_AN, '-g^', 'LineWidth', 2);
grid on;
xlabel('Liczba współpracujących podsłuchiwaczy (L)');
ylabel('Secrecy Sum-Rate (bps/Hz)');
title('Ochrona 6 GHz: Brak AN vs Sztuczny Szum (AN)');
legend('6 GHz (Zwykłe ZF)', '6 GHz (ZF + Sztuczny Szum)', 'Location', 'SouthWest');

% --- Zapisywanie wyników ---
saveas(gcf, '../results/wykres_sztuczny_szum.png');