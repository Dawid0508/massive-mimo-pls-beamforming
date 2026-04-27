% =========================================================================
% SYMULACJA 4: Aktywny atak Eve - Pilot Contamination (Porwanie Wiązki)
% =========================================================================
clear; clc; close all;

% --- Parametry Systemu ---
fc = 6e9;                   % Pasmo 6 GHz
Nt = 64;                    % 64 anteny (Massive MIMO)
c = physconst('LightSpeed');
numIter = 50;
beta_values = 0:0.1:1;      % Względna moc fałszywego pilota Eve (0% do 100%)

% Pozycje w przestrzeni
theta_bob = -30;            % Bob stoi na -30 stopniach
theta_eve = 20;             % Eve stoi na +20 stopniach

% Modele anten
ula = phased.ULA('NumElements', Nt, 'ElementSpacing', c/fc/2);
sv = phased.SteeringVector('SensorArray', ula);

% Inicjalizacja tablic
results_SR = zeros(length(beta_values), 1);
beam_pattern_clean = zeros(361, 1);
beam_pattern_hacked = zeros(361, 1);
angles = -180:180;

%% Pętla Symulacji
for b_idx = 1:length(beta_values)
    beta = beta_values(b_idx); % Siła ataku Eve w danej pętli
    temp_SR = 0;

    for iter = 1:numIter
        % 1. Prawdziwe Kanały
        h_bob = step(sv, fc, theta_bob)';
        h_eve = step(sv, fc, theta_eve)';

        % Dodanie drobnego rozproszenia (Rayleigh)
        h_bob_c = 0.8*h_bob + 0.2*(randn(1,Nt) + 1j*randn(1,Nt))/sqrt(2);
        h_eve_c = 0.8*h_eve + 0.2*(randn(1,Nt) + 1j*randn(1,Nt))/sqrt(2);

        % ===============================================================
        % 2. FAZA TRENINGOWA (ESTYMACJA KANAŁU) - Tu dzieje się magia!
        % Stacja bazowa szacuje kanał jako sumę sygnału Boba i fałszywego sygnału Eve
        noise_est = 0.05 * (randn(1,Nt) + 1j*randn(1,Nt))/sqrt(2);
        h_est = h_bob_c + sqrt(beta) * h_eve_c + noise_est;
        % ===============================================================

        % 3. Precoding (MRT oparty na SKAŻONYM kanale)
        w = h_est' / norm(h_est);

        % 4. Obliczanie przepustowości
        SNR_tx = 1000; % Przykładowa moc
        R_bob = log2(1 + SNR_tx * abs(h_bob_c * w)^2);
        R_eve = log2(1 + SNR_tx * abs(h_eve_c * w)^2);

        temp_SR = temp_SR + max(0, R_bob - R_eve);

        % Zapisanie kształtu wiązki do wykresu (dla czystego i zhakowanego)
        if iter == 1 
            a_sweep = step(sv, fc, angles)';
            if beta == 0
                beam_pattern_clean = 10*log10(abs(a_sweep * w).^2);
            elseif beta == 1
                beam_pattern_hacked = 10*log10(abs(a_sweep * w).^2);
            end
        end
    end
    results_SR(b_idx) = temp_SR / numIter;
end

%% --- Rysowanie Wykresów ---
figure('Position', [100 100 1000 400]);
%--figure('Color', 'w', 'Position', [100 100 1000 400]);


% Wykres 1: Spadek Secrecy Rate
subplot(1,2,1);
plot(beta_values * 100, results_SR, '-ro', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; xlabel('Moc fałszywego pilota Eve (%)'); ylabel('Secrecy Rate (bps/Hz)');
title('Zapaść systemu przez skażenie pilotów');

% Wykres 2: Porównanie kształtu wiązki
subplot(1,2,2);
plot(angles, beam_pattern_clean - max(beam_pattern_clean), 'b-', 'LineWidth', 1.5); hold on;
plot(angles, beam_pattern_hacked - max(beam_pattern_hacked), 'r--', 'LineWidth', 2);
xline(theta_bob, 'g:', 'Bob (-30^{\circ})', 'LineWidth', 2);
xline(theta_eve, 'k:', 'Eve (+20^{\circ})', 'LineWidth', 2);
grid on; box on;
xlim([-90 90]); ylim([-30 5]);
xlabel('Kąt (Stopnie)'); ylabel('Znormalizowane Wzmocnienie (dB)');
title('Fizyczne "Porwanie Wiązki"');
legend('Wiązka Czysta (Bez Ataku)', 'Wiązka Zhakowana (Beta = 100%)', 'Location', 'South');

% --- Zapisywanie wyników ---
saveas(gcf, '../results/wykres_pilot_contamination.png');