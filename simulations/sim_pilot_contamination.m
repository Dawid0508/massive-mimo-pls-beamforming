% =========================================================================
% SCENARIO: Active pilot contamination ("Beam Hijacking")
% -------------------------------------------------------------------------
% Zaktualizowano: nrCDLChannel, fizyczny FSPL, Transmit SNR.
% TOPOLOGIA: Asymetryczny Atak z Bliska. Bob jest daleko (80m), Ewa jest
% blisko (30m). Dzięki mniejszemu tłumieniu trasy, Ewa drastycznie szybko
% przejmuje wiązkę, nawet nadając z ułamkiem mocy Boba (małe beta).
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- PARAMETRY SYSTEMU I TOPOLOGIA ---------------------------------------
beta_values = 0:0.05:1;               % Moc nadawania Ewy (% mocy Boba)
theta_bob   = -20;
dist_b      = 80;                     % Bob jest daleko!
theta_eve   = 40;
dist_e      = 60;                     % Ewa jest bardzo blisko stacji!

numIter     = 20;
SNR_tx_dB   = 100;                    % Transmit SNR [dB]
P_tx        = 10^(SNR_tx_dB / 10);
noise_var   = 1;

bands = struct( ...
    'name',  {'6 GHz (Nt=32)', '28 GHz (Nt=256)'}, ...
    'fc',    {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',    {p.Nt_sub6, 256}, ...    % Ustawiono 256 anten dla 28 GHz
    'cdl',   {p.cdl_sub6, p.cdl_mmwave});

angles = -90:0.5:90;

% Macierze wynikowe
SR_curves   = zeros(2, length(beta_values));
beam_clean  = zeros(2, length(angles));
beam_hijack = zeros(2, length(angles));

% Inicjalizacja modeli kanału
cdl_b = setup_matlab_cdl(nrCDLChannel, 1, bands(1).fc, theta_bob); % Dummy init
cdl_e = setup_matlab_cdl(nrCDLChannel, 1, bands(1).fc, theta_eve); % Dummy init

for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;  cdl_tag = bands(b).cdl;
    
    [PL_lin_b, PL_dB_b] = compute_fspl(dist_b, fc);
    [PL_lin_e, PL_dB_e] = compute_fspl(dist_e, fc);
    
    fprintf('\n--- Pilot Contamination @ %s ---\n', bands(b).name);
    fprintf('  Bob: %g m (FSPL: %.2f dB)\n', dist_b, PL_dB_b);
    fprintf('  Eve: %g m (FSPL: %.2f dB)\n', dist_e, PL_dB_e);
    
    [~, sv] = setup_ula(Nt, fc);
    a_sweep = step(sv, fc, angles);
    
    % Re-konfiguracja dla danego pasma i liczby anten
    cdl_b = setup_matlab_cdl(cdl_b, Nt, fc, theta_bob);
    cdl_e = setup_matlab_cdl(cdl_e, Nt, fc, theta_eve);
    
    for be_idx = 1:length(beta_values)
        beta = beta_values(be_idx);
        SR_acc = 0;
        bp_acc = zeros(length(angles), 1);
        
        for it = 1:numIter
            release(cdl_b);
            cdl_b.Seed = randi([0 2^31-1]);
            [pg_b, ~] = cdl_b();
            hb = squeeze(sum(pg_b, 2)); hb = hb(:);
            % Fizyczny kanał Boba
            h_eff_b = sqrt((1 / PL_lin_b) * Nt) * (hb / norm(hb));
            
            release(cdl_e);
            cdl_e.Seed = randi([0 2^31-1]);
            [pg_e, ~] = cdl_e();
            he = squeeze(sum(pg_e, 2)); he = he(:);
            % Fizyczny kanał Ewy
            h_eff_e = sqrt((1 / PL_lin_e) * Nt) * (he / norm(he));
            
            % Stacja bazowa szacuje kanał w fazie uplink.
            % Odbiera nałożone na siebie piloty (Ewa wysyła z ułamkiem beta mocy Boba)
            n_est = sqrt(1 / PL_lin_b) * 0.1 * (randn(Nt,1) + 1j*randn(Nt,1))/sqrt(2);
            h_est = h_eff_b + sqrt(beta) * h_eff_e + n_est;
            
            % Stacja używa zatrutej estymaty do stworzenia prekodera MRT
            w = h_est / norm(h_est);
            
            % Obliczanie pojemności dla stacji transmitującej P_tx do Boba
            R_b = log2(1 + P_tx * abs(h_eff_b' * w)^2 / noise_var);
            R_e = log2(1 + P_tx * abs(h_eff_e' * w)^2 / noise_var);
            
            SR_acc = SR_acc + max(0, R_b - R_e);
            bp_acc = bp_acc + abs(a_sweep' * w).^2;
        end
        
        SR_curves(b, be_idx) = SR_acc / numIter;
        bp_avg = bp_acc / numIter;
        
        % Zapisywanie uśrednionego kształtu wiązki dla skrajnych przypadków
        if abs(beta) < 1e-6
            beam_clean(b, :) = 10*log10(bp_avg).';
        elseif abs(beta - 1) < 1e-6
            beam_hijack(b, :) = 10*log10(bp_avg).';
        end
    end
end

% --- WIZUALIZACJA --------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

% Górny Lewy: Spadek pojemności
subplot(2, 2, 1);
plot(beta_values*100, SR_curves(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(beta_values*100, SR_curves(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel("Eve's pilot transmit power \beta (% of Bob)");
ylabel('Secrecy Rate (bits/s/Hz)');
title('Secrecy Collapse under Pilot Contamination');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

% Górny Prawy: Relatywny spadek (Pokazuje jak szybko pęka zabezpieczenie)
subplot(2, 2, 2);
ratio = SR_curves ./ (SR_curves(:,1) + eps);
plot(beta_values*100, ratio(1,:), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(beta_values*100, ratio(2,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel("Eve's pilot transmit power \beta (% of Bob)");
ylabel('Normalised Secrecy Rate');
title('Relative Collapse Speed');
legend(bands(1).name, bands(2).name, 'Location', 'NorthEast');

% Dolne wykresy: Kradzież wiązki (Beam Hijacking) dla 6 GHz i 28 GHz
for b = 1:2
    subplot(2, 2, 2 + b);
    plot(angles, beam_clean(b,:)  - max(beam_clean(b,:)),  'b-',  'LineWidth', 1.7); hold on;
    plot(angles, beam_hijack(b,:) - max(beam_hijack(b,:)), 'r--', 'LineWidth', 2.0);
    mark_bob(theta_bob, sprintf('Bob (%d^{\\circ}, %dm)', theta_bob, dist_b), 'right');
    mark_eve(theta_eve, sprintf('Eve (+%d^{\\circ}, %dm)', theta_eve, dist_e), 'left');
    pls_axis_prefs(gca, 'refLabelV', 'top');
    grid on; box on;
    xlim([-60 60]); ylim([-30 5]);
    xlabel('Angle (deg)'); ylabel('Normalised gain (dB)');
    title(sprintf('Ergodic Beam Pattern @ %s', bands(b).name));
    legend('Clean (\beta=0)', 'Hijacked (\beta=1)', 'Location', 'South');
end

sgtitle(sprintf('Active Pilot Contamination (Target: Bob %dm, Attacker: Eve %dm)', dist_b, dist_e));
save_figure(fig, 'fig_pilot_contamination');

% =========================================================================
% FUNKCJA POMOCNICZA
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