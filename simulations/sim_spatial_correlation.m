% =========================================================================
% SYMULACJA 1: Wpływ korelacji przestrzennej na Secrecy Rate i Fairness
% =========================================================================
clear; clc; close all;

% --- Parametry Systemu ---
Nt = 64;                % Liczba anten na stacji bazowej (BS)
K = 50;                 % Liczba legalnych użytkowników (Bobs)
rho_values = 0:0.1:0.9; % Wektor wartości korelacji do przebadania [cite: 13]
num_monte_carlo = 100;  % Liczba pętli uśredniających (Monte Carlo)

% Przygotowanie tablic na wyniki (wypełnione zerami dla optymalizacji)
results_sum_rate_zf = zeros(length(rho_values), 1);
results_jains_idx_zf = zeros(length(rho_values), 1);

% --- Główna Pętla Symulacji ---
for r_idx = 1:length(rho_values)
    rho = rho_values(r_idx);

    % Zmienne tymczasowe dla metody Monte Carlo
    temp_sum_rate = 0;
    temp_jains = 0;

    for mc = 1:num_monte_carlo
        % 1. Definicja macierzy korelacji R [cite: 11]
        R = zeros(Nt, Nt);
        for i = 1:Nt
            for j = 1:Nt
                R(i,j) = rho^(abs(i-j));
            end
        end

        % 2. Generowanie nieskorelowanych kanałów (idealny Rayleigh) [cite: 7]
        H_uncorr = (randn(Nt, K) + 1i*randn(Nt, K)) / sqrt(2);

        % 3. Nałożenie korelacji przestrzennej [cite: 12]
        H_corr = sqrtm(R) * H_uncorr; 

        % ===============================================================
        % ===============================================================
        % PRAWDZIWA MATEMATYKA PLS (Precoding Zero-Forcing)
        % ===============================================================
        
        % 1. Obliczenie surowej macierzy ZF (Pseudoodwrotność Moore'a-Penrose'a)
        % Używamy skorelowanego kanału H_corr!
        W_zf_unnorm = H_corr * inv(H_corr' * H_corr);
        
        % ===============================================================
        % 2. Normalizacja Wektorowa (Vector Normalization)
        % Zapewnia, że każda antena/wiązka ma równy limit mocy nadawczej
        % ===============================================================
        W_zf = zeros(Nt, K);
        for k = 1:K
            % Dzielimy każdą kolumnę przez jej własną normę
            W_zf(:,k) = W_zf_unnorm(:,k) / norm(W_zf_unnorm(:,k));
        end
        
        % 3. Obliczenie mocy sygnału odebranego dla każdego Boba
        % Ponieważ to idealne ZF, interferencje wewnątrz komórki znikają (są bliskie 0)
        R_sec = zeros(K, 1);
        noise_var = 1; % Wariancja szumu (AWGN)
        
        for k = 1:K
            % Sygnał użyteczny dla k-tego użytkownika
            signal_power = abs(H_corr(:,k)' * W_zf(:,k))^2;
            
            % Przepustowość Boba (z wzoru Shannona)
            R_bob = log2(1 + (signal_power / noise_var));
            
            % UWAGA: Tutaj na razie zakładamy, że Eve (podsłuchiwacz) ma kanał = 0.
            % W kolejnym kroku dodamy model Eve, żeby policzyć prawdziwy Secrecy Rate!
            R_eve = 0; 
            
            % Wzór na Secrecy Rate: max(0, R_bob - R_eve)
            R_sec(k) = max(0, R_bob - R_eve);
        end
        % ===============================================================


        % 4. Obliczenia wynikowe
        sum_rate = sum(R_sec);

        % Wskaźnik Sprawiedliwości Jaina (Jain's Fairness Index) [cite: 27]
        % Wzór: (Suma R)^2 / (K * Suma R^2) [cite: 28]
        jains_index = (sum(R_sec)^2) / (K * sum(R_sec.^2)); 

        % Dodawanie do średniej Monte Carlo
        temp_sum_rate = temp_sum_rate + sum_rate;
        temp_jains = temp_jains + jains_index;
    end

    % Uśrednianie wyników dla danego rho i zapis do głównej tablicy
    results_sum_rate_zf(r_idx) = temp_sum_rate / num_monte_carlo;
    results_jains_idx_zf(r_idx) = temp_jains / num_monte_carlo;
end

% --- Rysowanie Wykresów ---
figure;
yyaxis left;
plot(rho_values, results_sum_rate_zf, '-bo', 'LineWidth', 2);
ylabel('Sumaryczny Secrecy Rate (bits/s/Hz)');
xlabel('Współczynnik Korelacji Przestrzennej (\rho)');

yyaxis right;
plot(rho_values, results_jains_idx_zf, '-r^', 'LineWidth', 2);
ylabel('Wskaźnik Jaina (Jain''s Fairness Index)');
ylim([0 1.1]); % Wskaźnik Jaina zawsze jest od 0 do 1 [cite: 28]

title('Wpływ korelacji przestrzennej na system Massive MIMO PLS');
grid on;
legend('Secrecy Sum-Rate (ZF)', 'Jain''s Fairness (ZF)', 'Location', 'southwest');

% --- Zapisywanie wyników ---
% Zapisuje wykres udowadniający zapaść algorytmu ZF przy wysokiej korelacji
saveas(gcf, '../results/wykres_korelacja_przestrzenna.png');