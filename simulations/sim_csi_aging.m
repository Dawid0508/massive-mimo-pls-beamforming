% =========================================================================
% SCENARIO: CSI aging (Doppler) and L-tap Wiener prediction as defence
% (Zaktualizowano do oficjalnego modelu nrCDLChannel z 5G Toolbox)
% Poprawiono logikę Asymetrii Dystansu, Tłumienia i Prawdziwego Dopplera
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Configuration -------------------------------------------------------
fc          = p.fc_sub6;              % use 6 GHz for tractable Nt = 32
Nt          = p.Nt_sub6;
K           = 4;
cdl_tag     = p.cdl_sub6;

v_kmh_vec   = 1:10:120;               % UE velocity sweep [km/h]
v_fixed     = 60;                     % km/h, used for the L sweep
L_vec       = 1:6;                    % Wiener predictor order
L_fixed     = 4;
delta_t     = 1e-3;                   % time slot duration [s]
tau         = 1;                      % prediction horizon (slots ahead)

% --- FIZYKA DYSTANSU I MOCY (Skopiowane z Baseline) ---
dist_b      = 50;                     % Bobowie na 50 m
dist_e      = 40;                     % Ewa ukrywa się na 40 m
SNR_tx_dB   = 110;                    % Moc stacji bazowej [dB]
P_tx        = 10^(SNR_tx_dB / 10);
noise_var   = 1;                      % Znormalizowany szum

% Rozmieszczenie przestrzenne użytkowników (dla CDL)
thetas_b    = linspace(-40, 40, K);   % Bobowie rozsiani sektorowo
theta_e     = 20;                     % Ewa pod określonym kątem

fprintf('\n--- CSI Aging with nrCDLChannel @ %.0f GHz ---\n', fc/1e9);
fprintf('  Transmit SNR: %d dB\n', SNR_tx_dB);
fprintf('  Bob (K=%d, d = %g m), Eve (d = %g m)\n', K, dist_b, dist_e);

numIter     = 20;

% Result storage
SR_no       = zeros(size(v_kmh_vec));
SR_wiener   = zeros(size(v_kmh_vec));
SR_perfect  = zeros(size(v_kmh_vec));
rho_curve   = zeros(size(v_kmh_vec)); % single-lag Jakes correlation

SR_vs_L_no      = zeros(size(L_vec));
SR_vs_L_wiener  = zeros(size(L_vec));
SR_vs_L_perfect = zeros(size(L_vec));

% --- Sweep A: velocity at L = L_fixed -----------------------------------
fprintf('Running Velocity Sweep...\n');
for v_idx = 1:length(v_kmh_vec)
    v = v_kmh_vec(v_idx);
    rho_curve(v_idx) = jakes_correlation(v, fc, delta_t);
    [SR_no(v_idx), SR_wiener(v_idx), SR_perfect(v_idx)] = ...
        run_aging_trial_cdl(Nt, K, P_tx, noise_var, numIter, ...
                            v, fc, delta_t, tau, L_fixed, ...
                            dist_b, dist_e, thetas_b, theta_e, cdl_tag);
end

% --- Sweep B: predictor order L at v = v_fixed --------------------------
fprintf('Running Predictor Order (L) Sweep...\n');
for l_idx = 1:length(L_vec)
    [SR_vs_L_no(l_idx), SR_vs_L_wiener(l_idx), SR_vs_L_perfect(l_idx)] = ...
        run_aging_trial_cdl(Nt, K, P_tx, noise_var, numIter, ...
                            v_fixed, fc, delta_t, tau, L_vec(l_idx), ...
                            dist_b, dist_e, thetas_b, theta_e, cdl_tag);
end

% --- Visualisation -------------------------------------------------------
c = pls_colors();
fig = figure('Color', c.bg, 'Position', [100 100 1200 760]);

% Top-left: Jakes correlation vs velocity
subplot(2, 2, 1);
plot(v_kmh_vec, rho_curve, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b');
grid on; box on;
xlabel('UE velocity (km/h)'); ylabel('J_0(2\pi f_d \Delta t)');
title(sprintf('Jakes correlation @ %.0f GHz, \\Delta t = %.0f ms', fc/1e9, delta_t*1e3));

% Top-right: Secrecy vs velocity for the three strategies
subplot(2, 2, 2);
plot(v_kmh_vec, SR_no,      '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r'); hold on;
plot(v_kmh_vec, SR_wiener,  '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
plot(v_kmh_vec, SR_perfect, '-^', 'Color', c.perfect, 'LineWidth', 2, ...
    'MarkerFaceColor', c.perfect, 'DisplayName', 'Perfect CSI');
grid on; box on;
xlabel('UE velocity (km/h)'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title('Three strategies under CSI aging (CDL Channel)');
legend('No prediction', sprintf('Wiener (L=%d)', L_fixed), 'Perfect CSI', ...
       'Location', 'East');

% Bottom-left: Secrecy vs predictor order L at v_fixed
subplot(2, 2, 3);
plot(L_vec, SR_vs_L_no,      '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r'); hold on;
plot(L_vec, SR_vs_L_wiener,  '-go', 'LineWidth', 2, 'MarkerFaceColor', 'g');
plot(L_vec, SR_vs_L_perfect, '-^', 'Color', c.perfect, 'LineWidth', 2, ...
    'MarkerFaceColor', c.perfect, 'DisplayName', 'Perfect CSI');
grid on; box on;
xlabel('Wiener predictor order L'); ylabel('Secrecy Sum-Rate (bits/s/Hz)');
title(sprintf('Diminishing returns of predictor order  (v = %d km/h)', v_fixed));
legend('No prediction', 'Wiener', 'Perfect CSI', 'Location', 'West');

% Bottom-right: relative gain of Wiener over No-prediction
subplot(2, 2, 4);
gain = (SR_wiener - SR_no) ./ max(SR_perfect - SR_no, eps);
plot(v_kmh_vec, gain*100, '-mo', 'LineWidth', 2, 'MarkerFaceColor', 'm');
grid on; box on; ylim([0 110]);
xlabel('UE velocity (km/h)'); ylabel('Recovered fraction (%)');
title('Wiener as a "shield": fraction of perfect-CSI gap closed');

sgtitle(sprintf('CSI aging at %.0f GHz: Wiener vs Doppler (3GPP nrCDL)', fc/1e9));

save_figure(fig, 'fig_csi_aging_cdl');


% =========================================================================
%                         Local helper functions
% =========================================================================
function [SR_no, SR_w, SR_p] = run_aging_trial_cdl(Nt, K, P_tx, noise_var, ...
                                                   numIter, v, fc, dt, tau, L, ...
                                                   dist_b, dist_e, thetas_b, theta_e, cdl_tag)
% Runs Monte-Carlo trials using the actual nrCDLChannel physics engine.

    SR_no = 0; SR_w = 0; SR_p = 0;

    % Pre-compute the Wiener weights once per call (Using theoretical Jakes)
    % This reflects reality: the base station calculates weights via
    % theory/math, but applies them to a real-world physical channel.
    rr = zeros(L, 1);          
    R  = zeros(L, L);          
    for i = 1:L
        for j = 1:L
            R(i, j) = jakes_correlation(v, fc, (i-j) * dt);
        end
        rr(i) = jakes_correlation(v, fc, (tau + (i-1)) * dt);
    end
    R = R + 1e-9 * eye(L);     
    w_pred = R \ rr;           

    % Obliczenia Tłumienia Przestrzennego (Path Loss)
    [PL_lin_b, ~] = compute_fspl(dist_b, fc);
    [PL_lin_e, ~] = compute_fspl(dist_e, fc);

    % Fizyka Dopplera dla nrCDLChannel
    T = L + tau;
    c = physconst('LightSpeed');
    fd = (v * 1e3 / 3600) * (fc / c); 

    % Inicjalizacja kanałów
    cdl_b = cell(K, 1);
    for k = 1:K
        cdl_b{k} = setup_matlab_cdl_aging(Nt, fc, thetas_b(k), cdl_tag, fd, dt, T);
    end
    cdl_e = setup_matlab_cdl_aging(Nt, fc, theta_e, cdl_tag, fd, dt, T);

    for it = 1:numIter
        % H_time(:, :, t) is the K-user effective channel at slot t
        H_time = zeros(Nt, K, T);
        
        % Generujemy kanał dla Bobów w dziedzinie czasu
        for k = 1:K
            release(cdl_b{k});
            cdl_b{k}.Seed = randi([0 2^31-1]);
            [pg_b, ~] = cdl_b{k}(); 
            % pg_b: [T, NumPaths, 1, Nt]. Redukcja ścieżek z zachowaniem czasu.
            h_raw = squeeze(sum(pg_b, 2)).'; % [Nt, T]
            
            for t_idx = 1:T
                ht = h_raw(:, t_idx);
                % Aplikacja tłumienia fizycznego do kanału efektywnego
                H_time(:, k, t_idx) = sqrt((1/PL_lin_b) * Nt) * (ht / norm(ht));
            end
        end
        
        % Generujemy kanał dla Ewy w dziedzinie czasu
        release(cdl_e);
        cdl_e.Seed = randi([0 2^31-1]);
        [pg_e, ~] = cdl_e();
        he_raw = squeeze(sum(pg_e, 2)).';
        h_eve_time = zeros(Nt, T);
        for t_idx = 1:T
            ht = he_raw(:, t_idx);
            h_eve_time(:, t_idx) = sqrt((1/PL_lin_e) * Nt) * (ht / norm(ht));
        end

        % Past observations at t = 1..L,  prediction horizon at t = L+tau
        H_past = H_time(:, :, 1:L);     
        H_now  = squeeze(H_time(:, :, L+tau)); % The *true* channel right now
        h_eve  = h_eve_time(:, L+tau);

        % --- Strategy 1: no prediction (use the most recent sample) ----
        H_used_no = squeeze(H_past(:, :, L));   

        % --- Strategy 2: L-tap Wiener prediction ---------------------
        H_used_w = zeros(Nt, K);
        for k = 1:K
            % Past samples for user k stacked oldest -> newest
            samples = squeeze(H_past(:, k, end:-1:1)); % Nt x L
            H_used_w(:, k) = samples * w_pred;
        end

        % --- Strategy 3: perfect CSI ---------------------------------
        H_used_p = H_now;

        SR_no = SR_no + secrecy_with_estimate(H_used_no, H_now, h_eve, P_tx, noise_var);
        SR_w  = SR_w  + secrecy_with_estimate(H_used_w,  H_now, h_eve, P_tx, noise_var);
        SR_p  = SR_p  + secrecy_with_estimate(H_used_p,  H_now, h_eve, P_tx, noise_var);
    end

    SR_no = SR_no / numIter;
    SR_w  = SR_w  / numIter;
    SR_p  = SR_p  / numIter;
end


function SR = secrecy_with_estimate(H_est, H_true, h_eve, P_tx, noise_var)
% Builds a ZF precoder from H_est, evaluates it on the true channel
% and uses P_tx for Transmit Power logic.

    K = size(H_est, 2);
    W_raw = H_est * pinv(H_est' * H_est + 1e-9*eye(K));
    % Normalizacja prekodera 
    W = W_raw / norm(W_raw, 'fro'); 

    R_b = zeros(K, 1); R_e = zeros(K, 1);
    for k = 1:K
        sig    = P_tx * abs(H_true(:,k)' * W(:,k))^2;
        intf   = P_tx * (sum(abs(H_true(:,k)' * W).^2) - abs(H_true(:,k)' * W(:,k))^2);
        R_b(k) = log2(1 + sig / (intf + noise_var));

        sig_e  = P_tx * abs(h_eve' * W(:,k))^2;
        intf_e = P_tx * (sum(abs(h_eve' * W).^2) - abs(h_eve' * W(:,k))^2);
        R_e(k) = log2(1 + sig_e / (intf_e + noise_var));
    end
    SR = sum(max(0, R_b - R_e));
end


% =========================================================================
% FUNKCJA POMOCNICZA: Konfiguracja kanału z obsługą upływu czasu (Doppler)
% =========================================================================
function cdl = setup_matlab_cdl_aging(Nt, fc, theta, band_tag, fd, dt, T)
    cdl = nrCDLChannel;
    
    cdl.DelayProfile = 'CDL-A';
    
    if fc < 10e9  
        cdl.DelaySpread = 30e-9;  
    else          
        cdl.DelaySpread = 10e-9;  
    end
    
    cdl.CarrierFrequency = fc;
    
    % KLUCZOWY PARAMETR DLA EFEKTU AGING!
    cdl.MaximumDopplerShift = fd;         
    
    cdl.TransmitAntennaArray.Size = [1 Nt 1 1 1]; 
    cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1]; 
    
    cdl.TransmitArrayOrientation = [-theta; 0; 0];
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
    
    % Ustawienia czasu (dt to nasz delta_t, T to ilosc wymaganych probek)
    cdl.SampleRate = 1/dt;
    cdl.NumTimeSamples = T;
    cdl.ChannelFiltering = false; 
end