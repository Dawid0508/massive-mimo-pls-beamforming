% =========================================================================
% SCENARIO: Fixed beam, moving Bob (6 GHz vs 28 GHz, 3GPP CDL)
% Symulacja zmienna w czasie z efektem Dopplera i twardą fizyką (Path Loss)
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- FIZYKA DYSTANSU I MOCY ---
dist_b     = 50;                        % Bob na 50 m
dist_e     = 40;                        % Ewa na 40 m
SNR_tx_dB  = 110;                       % ZMIANA: Transmit SNR (Moc stacji bazowej)
P_tx       = 10^(SNR_tx_dB / 10);       % Moc liniowa stacji

% --- PARAMETRY RUCHU I CZASU (DOPPLER) ---
v_bob_kmh   = 50;      % Bob jedzie 50 km/h
v_eve_kmh   = 0;       % Ewa stoi w miejscu
sample_rate = 10000;   % Próbkowanie 10 kHz (krok co 0.1 ms)
time_span   = 0.01;    % Badamy zachowanie systemu przez 10 milisekund
num_samples = time_span * sample_rate; 
time_ms     = (0:num_samples-1) * (1/sample_rate) * 1000; % Oś X (czas w ms)

bands = struct( ...
    'name', {'6 GHz', '28 GHz'}, ...
    'fc',   {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',   {p.Nt_sub6, p.Nt_mmwave}, ...
    'cdl',  {p.cdl_sub6, p.cdl_mmwave});

theta_beam = -30; % Kąt Boba w t=0 (cel wiązki)
theta_eve  = 20;  % Kąt Ewy

SR_mrt = zeros(2, num_samples);
SR_zf  = zeros(2, num_samples);

for f_idx = 1:2
    fc = bands(f_idx).fc; Nt = bands(f_idx).Nt; cdl_tag = bands(f_idx).cdl;
    
    % Wyliczamy fizyczne tłumienie dla aktualnego pasma
    [PL_b, ~] = compute_fspl(dist_b, fc);
    [PL_e, ~] = compute_fspl(dist_e, fc);
    
    % Inicjalizacja kanałów
    cdl_b = setup_matlab_cdl(Nt, fc, theta_beam, cdl_tag, v_bob_kmh, sample_rate, num_samples);
    cdl_e = setup_matlab_cdl(Nt, fc, theta_eve, cdl_tag, v_eve_kmh, sample_rate, num_samples);
    
    cdl_b.Seed = 101;
    cdl_e.Seed = 202; % Ewa ma swój własny niezależny seed!
    
    [pg_b, ~] = cdl_b();
    [pg_e, ~] = cdl_e();
    
    % --- KROK 1: TRENING I ZAMROŻENIE WIĄZKI (t = 0 ms) ---
    h_b_0 = squeeze(sum(pg_b(1, :, :, :), 2)); h_b_0 = h_b_0(:);
    h_e_0 = squeeze(sum(pg_e(1, :, :, :), 2)); h_e_0 = h_e_0(:);
    
    % Reset normy i wstrzyknięcie zysku szyku oraz tłumienia (Fizyka!)
    h_b_0_eff = sqrt((1/PL_b) * Nt) * (h_b_0 / norm(h_b_0));
    h_e_0_eff = sqrt((1/PL_e) * Nt) * (h_e_0 / norm(h_e_0));
    
    w_mrt_fixed = h_b_0_eff / norm(h_b_0_eff);
    
    P_null = eye(Nt) - (h_e_0_eff * (h_e_0_eff' / (h_e_0_eff' * h_e_0_eff)));
    w_zf_fixed = P_null * h_b_0_eff;
    if norm(w_zf_fixed) > 1e-9
        w_zf_fixed = w_zf_fixed / norm(w_zf_fixed);
    else
        w_zf_fixed = w_mrt_fixed;
    end
    
    % --- KROK 2: EWOLUCJA W CZASIE (Bieg z czasem) ---
    for t = 1:num_samples
        t_b = min(t, size(pg_b, 1));
        h_b_t = squeeze(sum(pg_b(t_b, :, :, :), 2)); h_b_t = h_b_t(:);
        
        t_e = min(t, size(pg_e, 1));
        h_e_t = squeeze(sum(pg_e(t_e, :, :, :), 2)); h_e_t = h_e_t(:);
        
        % Kanał w ułamku sekundy 't' Z FIZYKĄ (reset + zysk + tłumienie)
        h_b_t_eff = sqrt((1/PL_b) * Nt) * (h_b_t / norm(h_b_t));
        h_e_t_eff = sqrt((1/PL_e) * Nt) * (h_e_t / norm(h_e_t));
        
        % Obliczenia pojemności z ustaloną mocą stacji (P_tx)
        R_b_mrt = log2(1 + P_tx * abs(h_b_t_eff' * w_mrt_fixed)^2);
        R_e_mrt = log2(1 + P_tx * abs(h_e_t_eff' * w_mrt_fixed)^2);
        SR_mrt(f_idx, t) = max(0, R_b_mrt - R_e_mrt);
        
        R_b_zf = log2(1 + P_tx * abs(h_b_t_eff' * w_zf_fixed)^2);
        R_e_zf = log2(1 + P_tx * abs(h_e_t_eff' * w_zf_fixed)^2);
        SR_zf(f_idx, t) = max(0, R_b_zf - R_e_zf);
    end
end

% --- TWORZENIE WYKRESÓW ---
angles = -90:0.05:90;
fig = figure('Color', 'w', 'Position', [100 100 1200 520]);

subplot(1, 2, 1);
plot(time_ms, SR_mrt(1,:), '-b', 'LineWidth', 2); hold on;
plot(time_ms, SR_mrt(2,:), '-r', 'LineWidth', 2);
plot(time_ms, SR_zf(1,:),  '--b', 'LineWidth', 2);
plot(time_ms, SR_zf(2,:),  '--r', 'LineWidth', 2);
grid on; box on;
xlabel('Czas od pomiaru kanału (ms)');
ylabel('Secrecy Rate (bits/s/Hz)');
title(sprintf('Degradacja przez Dopplera (Bob: %d km/h)', v_bob_kmh));
legend('6 GHz MRT', '28 GHz MRT', '6 GHz ZF', '28 GHz ZF', 'Location', 'SouthWest');
ylim([0 max(SR_zf(:))*1.1]);

fc2 = bands(2).fc; Nt2 = bands(2).Nt;
[~, sv2] = setup_ula(Nt2, fc2);
cdl_28 = setup_matlab_cdl(Nt2, fc2, theta_beam, bands(2).cdl, 0, sample_rate, 1);
cdl_28.Seed = 101;
[pg_28, ~] = cdl_28();
h_beam2 = squeeze(sum(pg_28, 2)); h_beam2 = h_beam2(:);
w_fix2 = h_beam2 / norm(h_beam2);

a_sweep = step(sv2, fc2, angles);
pat = 10*log10(abs(w_fix2' * a_sweep).^2);
pat = pat - max(pat);

subplot(1, 2, 2);
plot(angles, pat, '-r', 'LineWidth', 1.5); hold on;
xline(theta_beam, 'g--', 'Bob (t=0)', 'LineWidth', 1.5);
xline(theta_eve, 'k:', 'Eve', 'LineWidth', 1.5);
grid on; box on;
xlabel('Angle (deg)'); ylabel('Gain (dB)');
title('28 GHz fixed beam snapshot (t = 0 ms)');
ylim([-40 5]);

sgtitle(sprintf('Channel Aging: Beam @ %d^\\circ, Eve @ %d^\\circ (Transmit SNR = %d dB)', ...
    theta_beam, theta_eve, SNR_tx_dB));
save_figure(fig, 'fig_moving_bob_doppler');

% =========================================================================
function cdl = setup_matlab_cdl(Nt, fc, theta, band_tag, v_kmh, sample_rate, num_samples)
    cdl = nrCDLChannel;
    cdl.DelayProfile = 'CDL-A'; 
    if fc < 10e9  
        cdl.DelaySpread = 30e-9;  
    else          
        cdl.DelaySpread = 10e-9;  
    end
    cdl.CarrierFrequency = fc;
    c = 299792458;                 
    v_ms = v_kmh / 3.6;            
    cdl.MaximumDopplerShift = (v_ms * fc) / c; 
    cdl.SampleRate = sample_rate;
    cdl.NumTimeSamples = num_samples;
    cdl.TransmitAntennaArray.Size = [1 Nt 1 1 1]; 
    cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1]; 
    cdl.TransmitArrayOrientation = [-theta; 0; 0];
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
    cdl.ChannelFiltering = false; 
end