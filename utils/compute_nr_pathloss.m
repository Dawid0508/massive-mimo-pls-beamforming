function [PL_lin, PL_dB] = compute_nr_pathloss(dist_m, fc_Hz)
    cfgPL = nrPathLossConfig;
    cfgPL.Scenario = 'UMi'; 
    h_bs = 10.0; 
    h_ut = 1.5;  
    
    num_points = length(dist_m);
    PL_dB = zeros(1, num_points);
    PL_lin = zeros(1, num_points);
    
    pos_bs = [0; 0; h_bs]; % Stacja bazowa twardo w (0,0,10)
    
    % Liczymy path loss niezależnie dla każdego użytkownika (kuloodporne)
    for i = 1:num_points
        pos_ue = [dist_m(i); 0; h_ut];
        PL_dB(i) = nrPathLoss(cfgPL, fc_Hz, false, pos_bs, pos_ue);
        PL_lin(i) = 10^(PL_dB(i) / 10);
    end
end