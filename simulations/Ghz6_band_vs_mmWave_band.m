%% 6G PLS: Massive MIMO vs Ultra-Massive MIMO with Path Loss
clear; clc; close all;

% --- Configuration ---
Nt_6 = 32;               % 6 GHz Massive MIMO
Nt_mm = 512;             % 28 GHz Ultra-Massive MIMO
dist = 1;               % Distance in meters (Alice to Bob/Eve)
SNR_dB = 40:2:80;
SNR_lin = 10.^(SNR_dB/10);
numIter = 200;
c = physconst('LightSpeed');

% Frequency Bands
fc_6GHz = 6e9;          
fc_mmWave = 28e9;       

% --- Path Loss Calculation (FSPL) ---
% FSPL (dB) = 20log10(d) + 20log10(f) - 147.55
PL_6_dB = 20*log10(dist) + 20*log10(fc_6GHz) - 147.55;
PL_mm_dB = 20*log10(dist) + 20*log10(fc_mmWave) - 147.55;

% Convert to Linear Loss (Loss factor > 1)
L_6 = 10^(PL_6_dB/10);
L_mm = 10^(PL_mm_dB/10);

fprintf('Path Loss at 6 GHz: %.2f dB\n', PL_6_dB);
fprintf('Path Loss at 28 GHz: %.2f dB\n', PL_mm_dB);
fprintf('Difference (Penalty for mmWave): %.2f dB\n', PL_mm_dB - PL_6_dB);

% Define Separate Antenna Arrays & Steering Vectors
ula_6 = phased.ULA('NumElements', Nt_6, 'ElementSpacing', c/fc_6GHz/2);
ula_mm = phased.ULA('NumElements', Nt_mm, 'ElementSpacing', c/fc_mmWave/2);
sv_6 = phased.SteeringVector('SensorArray', ula_6);
sv_mm = phased.SteeringVector('SensorArray', ula_mm);

% Results Storage
SecrecyCap = zeros(2, 2, length(SNR_dB)); 

%% PART 1: Monte Carlo Simulation
for f_idx = 1:2
    if f_idx == 1
        fc = fc_6GHz; sv_obj = sv_6; Nt = Nt_6; L = L_6;
    else
        fc = fc_mmWave; sv_obj = sv_mm; Nt = Nt_mm; L = L_mm;
    end
    
    for iter = 1:numIter
        theta_b = -30; 
        theta_e = -40; 
        
        hb = step(sv_obj, fc, theta_b)'; 
        he = step(sv_obj, fc, theta_e)';
        
        if f_idx == 1 % 6GHz Scattering
            hb = 0.8*hb + 0.2*(randn(1,Nt) + 1j*randn(1,Nt))/sqrt(2);
            he = 0.8*he + 0.2*(randn(1,Nt) + 1j*randn(1,Nt))/sqrt(2);
        end

        % Precoding (MRT & ZF)
        w_mrt = hb' / norm(hb);
        P_null = eye(Nt) - (he' * pinv(he * he') * he);
        w_zf = P_null * hb';
        if norm(w_zf) > 1e-6, w_zf = w_zf / norm(w_zf); else, w_zf = zeros(Nt, 1); end

        for s = 1:length(SNR_lin)
            % Scale Transmit Power by the Path Loss
            % P_received = P_transmit / Loss
            P_eff = SNR_lin(s) / L; 
            
            % MRT Capacity
            Rb_mrt = log2(1 + P_eff * abs(hb * w_mrt)^2);
            Re_mrt = log2(1 + P_eff * abs(he * w_mrt)^2);
            SecrecyCap(f_idx, 1, s) = SecrecyCap(f_idx, 1, s) + max(0, Rb_mrt - Re_mrt);
            
            % ZF Capacity
            Rb_zf = log2(1 + P_eff * abs(hb * w_zf)^2);
            Re_zf = log2(1 + P_eff * abs(he * w_zf)^2);
            SecrecyCap(f_idx, 2, s) = SecrecyCap(f_idx, 2, s) + max(0, Rb_zf - Re_zf);
        end
    end
end
SecrecyCap = SecrecyCap / numIter;

%% PART 2: Visualization
snap_theta_b = -30; snap_theta_e = 20;
angles = -90:0.05:90; 
figure('Color', 'w', 'Position', [100 100 1100 800]);

for f_idx = 1:2
    if f_idx == 1
        fc = fc_6GHz; sv_obj = sv_6; Nt = Nt_6; bandStr = '6 GHz (PL included)';
    else
        fc = fc_mmWave; sv_obj = sv_mm; Nt = Nt_mm; bandStr = '28 GHz mmWave (PL included)';
    end
    
    % Re-calculate snapshot patterns (Patterns are independent of Path Loss)
    hb_s = step(sv_obj, fc, snap_theta_b)';
    he_s = step(sv_obj, fc, snap_theta_e)';
    w_mrt_s = hb_s' / norm(hb_s);
    P_null_s = eye(Nt) - (he_s' * pinv(he_s * he_s') * he_s);
    w_zf_s = P_null_s * hb_s'; w_zf_s = w_zf_s / norm(w_zf_s);
    
    a_sweep = step(sv_obj, fc, angles);
    pat_mrt = 10*log10(abs(w_mrt_s' * a_sweep).^2);
    pat_zf  = 10*log10(abs(w_zf_s' * a_sweep).^2);
    
    % Plot Secrecy Capacity
    subplot(2, 2, f_idx);
    plot(SNR_dB, squeeze(SecrecyCap(f_idx,1,:)), 'b-o', 'LineWidth', 1.2); hold on;
    plot(SNR_dB, squeeze(SecrecyCap(f_idx,2,:)), 'r--s', 'LineWidth', 1.2);
    grid on; title(['Secrecy: ', bandStr]); xlabel('SNR (dB)'); ylabel('bps/Hz');
    legend('MRT', 'ZF', 'Location', 'NorthWest');
    
    % Plot Patterns
    subplot(2, 2, f_idx + 2);
    plot(angles, pat_mrt - max(pat_mrt), 'b', 'LineWidth', 1.5); hold on;
    plot(angles, pat_zf - max(pat_zf), 'r--', 'LineWidth', 1.2);
    xline(snap_theta_b, 'g:', 'Bob'); xline(snap_theta_e, 'k:', 'Eve');
    grid on; title(['Normalized Beam Pattern: ', bandStr]); ylim([-40 5]);
end
sgtitle(['6G PLS Analysis at Distance = ' num2str(dist) 'm']);