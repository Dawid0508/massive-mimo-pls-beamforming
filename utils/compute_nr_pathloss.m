function [PL_lin, PL_dB] = compute_nr_pathloss(dist_m, fc_Hz)
% COMPUTE_NR_PATHLOSS Oblicza tłumienie trasy zgodnie z 3GPP TR 38.901
% wykorzystując oficjalne funkcje 5G Toolbox.

    % 1. Konfiguracja obiektu nrPathLossConfig
    cfgPL = nrPathLossConfig;
    cfgPL.Scenario = 'UMi'; % Urban Micro - Street Canyon

    % 2. Definicja współrzędnych 3D [x; y; z]
    h_bs = 10.0; % Wysokość stacji bazowej (m)
    h_ut = 1.5;  % Wysokość użytkownika (m)
    
    num_points = length(dist_m);
    pos_bs = repmat([0; 0; h_bs], 1, num_points); 
    pos_ue = [dist_m; zeros(1, num_points); repmat(h_ut, 1, num_points)];
    
    % 3. Warunek LOS / NLOS
    % Funkcja przyjmuje tablicę logiczną: false dla NLOS, true dla LOS.
    is_los = false(1, num_points); 
    
    % 4. Obliczenie tłumienia (w dB)
    % POPRAWNA KOLEJNOŚĆ: config, częstotliwość, los, poz_BS, poz_UE
    PL_dB = nrPathLoss(cfgPL, fc_Hz, is_los, pos_bs, pos_ue);
    
    % 5. Przeliczenie na skalę liniową
    PL_lin = 10.^(PL_dB / 10);
end