% =========================================================================
% SCENARIO: Hardware phase noise impact on Massive MIMO PLS
% -------------------------------------------------------------------------
% (Zaktualizowano: FSPL, nrCDLChannel, stały budżet P_tx).
% Celowe wymuszenie bliskości kątowej, aby uwypuklić efekt "Beam Smearing".
% Real RF chains add per-antenna phase jitter, modeled as:
%       W_err = W .* exp(1j * sigma_phi * randn(size(W)))
% =========================================================================
pls_startup();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
p = default_params();
rng(p.rng_seed);

% --- Konfiguracja -------------------------------------------------------
sigma_deg_vec = 0:1:15;                     % Odchylenie std. błędu fazy [deg]
K             = 4;                          % Liczba użytkowników (Bobów)
numIter       = 20;
dist          = 50;                         % Dystans do użytkowników [m]
SNR_tx_dB     = 100;                        % Moc stacji bazowej (SNR_tx) [dB]
P_tx          = 10^(SNR_tx_dB / 10);

bands = struct( ...
    'name',  {'6 GHz (Massive MIMO)', '28 GHz (Ultra-Massive MIMO)'}, ...
    'fc',    {p.fc_sub6, p.fc_mmwave}, ...
    'Nt',    {p.Nt_sub6, p.Nt_mmwave}, ...
    'cdl',   {p.cdl_sub6, p.cdl_mmwave});

SR_results = zeros(2, 2, length(sigma_deg_vec));   % {6/28} x {Matrix/Vector} x sigma

% Konfiguracja migawki (snapshot) charakterystyki promieniowania (tylko 28 GHz)
sigma_snap_deg = [0, 5, 10];
angles    = -90:0.25:90;
bp_snap   = zeros(length(sigma_snap_deg), length(angles));

% Sztywne kąty referencyjne dla wykresu promieniowania (i dla Ewy w symulacji)
theta_b_ref = -10;
theta_e_ref =  -5; % Ewa czai się 5 stopni obok Boba referencyjnego!

% --- Wstępna alokacja obiektów 3GPP ---
cdl_b = cell(K, 1);
for k = 1:K, cdl_b{k} = nrCDLChannel; end
cdl_e = nrCDLChannel;

% --- Pętla Symulacji -----------------------------------------------------
for b = 1:2
    fc = bands(b).fc;  Nt = bands(b).Nt;  cdl_tag = bands(b).cdl;
    [PL_lin, PL_dB] = compute_fspl(dist, fc);
    
    fprintf('\n--- Symulacja Szumu Fazowego @ %s ---\n', bands(b).name);
    fprintf('  Dystans: %g m (FSPL: %.2f dB)\n', dist, PL_dB);
    
    [~, sv] = setup_ula(Nt, fc);
    
    for s_idx = 1:length(sigma_deg_vec)
        sigma_phi = deg2rad(sigma_deg_vec(s_idx));
        SR_acc = zeros(2, 1);
        
        for it = 1:numIter
            % Bobowie rozstawieni w sektorze, ale zawsze "zostawiamy"
            % referencyjnego Boba i celującą w niego Ewę
            theta_bobs = -60 + 120*rand(1, K);
            theta_bobs(1) = theta_b_ref; 
            
            H_eff = zeros(Nt, K);
            for k = 1:K
                cdl_b{k} = setup_matlab_cdl(cdl_b{k}, Nt, fc, theta_bobs(k));
                cdl_b{k}.Seed = randi([0 2^31-1]);
                [pg_b, ~] = cdl_b{k}();
                hb = squeeze(sum(pg_b, 2)); hb = hb(:);
                H_eff(:, k) = sqrt((1 / PL_lin) * Nt) * (hb / norm(hb));
            end
            
            % Jedna Ewa podkradająca się pod referencyjnego Boba
            cdl_e = setup_matlab_cdl(cdl_e, Nt, fc, theta_e_ref);
            cdl_e.Seed = randi([0 2^31-1]);
            [pg_e, ~] = cdl_e();
            he = squeeze(sum(pg_e, 2)); he = he(:);
            he_eff = sqrt((1 / PL_lin) * Nt) * (he / norm(he));
            
            % Obliczanie nominalnych prekoderów (ZF)
            W_raw = H_eff * pinv(H_eff' * H_eff);
            W_mat = W_raw / norm(W_raw, 'fro');
            
            W_vec = zeros(Nt, K);
            for k = 1:K
                if norm(W_raw(:, k)) > 1e-9
                    W_vec(:, k) = W_raw(:, k) / norm(W_raw(:, k)) * sqrt(1/K);
                end
            end
            
            % Aplikacja zniekształceń fazowych układów RF!
            err_mat = exp(1j * sigma_phi * randn(Nt, K));
            W_mat_e = W_mat .* err_mat;
            W_vec_e = W_vec .* err_mat;
            
            for n_idx = 1:2
                if n_idx == 1, W = W_mat_e; else, W = W_vec_e; end
                R_b = zeros(K, 1); R_e = zeros(K, 1);
                
                for k = 1:K
                    % Pojemność Boba 
                    sig_b  = P_tx * abs(H_eff(:,k)' * W(:,k))^2;
                    intf_b = 0;
                    for j = 1:K, if j~=k, intf_b = intf_b + P_tx * abs(H_eff(:,k)' * W(:,j))^2; end; end
                    R_b(k) = log2(1 + sig_b / (intf_b + 1)); % Szum znormalizowany do 1
                    
                    % Ewa podsłuchująca ten sam strumień
                    sig_e  = P_tx * abs(he_eff' * W(:,k))^2;
                    intf_e = 0;
                    for j = 1:K, if j~=k, intf_e = intf_e + P_tx * abs(he_eff' * W(:,j))^2; end; end
                    R_e(k) = log2(1 + sig_e / (intf_e + 1));
                end
                SR_acc(n_idx) = SR_acc(n_idx) + sum(max(0, R_b - R_e));
            end
        end
        SR_results(b, :, s_idx) = SR_acc / numIter;
    end
    
    % Generowanie snapshotów wiązki (tylko dla Ultra-Massive MIMO 28 GHz)
    if b == 2
        % Tu liczymy z idealnego (czystego geometrycznie) środowiska, 
        % aby wykres był czytelny.
        cdl_b_snap = setup_matlab_cdl(nrCDLChannel, Nt, fc, theta_b_ref);
        cdl_b_snap.Seed = 101;
        [pg_b_snap, ~] = cdl_b_snap();
        hb_snap = squeeze(sum(pg_b_snap, 2)); hb_snap = hb_snap(:);
        w0 = hb_snap / norm(hb_snap); 
        
        a_sweep = step(sv, fc, angles);
        for ss = 1:length(sigma_snap_deg)
            sigma_phi = deg2rad(sigma_snap_deg(ss));
            bp_acc = zeros(length(angles), 1);
            for it = 1:200
                w = w0 .* exp(1j * sigma_phi * randn(Nt, 1));
                bp_acc = bp_acc + abs(a_sweep' * w).^2;
            end
            bp_snap(ss, :) = 10*log10(bp_acc.' / 200);
        end
    end
end

% --- WIZUALIZACJA --------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1200 760]);

% Górny rząd: Porównanie SR względem błędu fazy dla 6 GHz i 28 GHz
for b = 1:2
    subplot(2, 2, b);
    plot(sigma_deg_vec, squeeze(SR_results(b,1,:)), '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
    plot(sigma_deg_vec, squeeze(SR_results(b,2,:)), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    grid on; box on;
    xlabel('Phase-error std \sigma_\phi (deg)');
    ylabel('Secrecy Sum-Rate (bits/s/Hz)');
    title(['Robustness to Phase Noise: ', bands(b).name]);
    legend('Matrix Normalization', 'Vector Normalization', 'Location', 'SouthWest');
end

% Dolny Lewy: Relatywna degradacja (Matrix norm)
subplot(2, 2, 3);
deg6  = squeeze(SR_results(1, 1, :)) / max(SR_results(1, 1, 1), eps);
deg28 = squeeze(SR_results(2, 1, :)) / max(SR_results(2, 1, 1), eps);
plot(sigma_deg_vec, deg6,  '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
plot(sigma_deg_vec, deg28, '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; box on; ylim([0 1.05]);
xlabel('Phase-error std \sigma_\phi (deg)');
ylabel('Normalised Secrecy Capacity');
title('Relative Degradation (Matrix Norm.)');
legend('6 GHz (Nt=32)', '28 GHz (Nt=512)', 'Location', 'SouthWest');

% Dolny Prawy: Zjawisko Beam Smearing na mmWave
subplot(2, 2, 4);
colors = lines(length(sigma_snap_deg));
for ss = 1:length(sigma_snap_deg)
    plot(angles, bp_snap(ss,:) - max(bp_snap(ss,:)), 'LineWidth', 2, 'Color', colors(ss,:)); hold on;
end
mark_bob(theta_b_ref, sprintf('Bob (%d^{\\circ})', theta_b_ref), 'right');
mark_eve(theta_e_ref, sprintf('Eve (+%d^{\\circ})', theta_e_ref), 'left');
pls_axis_prefs(gca, 'refLabelV', 'top');
grid on; box on;
xlim([-30 10]); ylim([-40 5]); % Przybliżenie na główny listek!
xlabel('Angle (deg)'); ylabel('Normalised gain (dB)');
title('28 GHz: Beam Smearing Effect');
legend(arrayfun(@(s) sprintf('\\sigma_\\phi = %d^{\\circ}', s), sigma_snap_deg, 'UniformOutput', false), ...
       'Location', 'SouthEast');

sgtitle(sprintf('Hardware Phase Noise Impact (K = %d, Target Proximity Attack)', K));
save_figure(fig, 'fig_phase_noise');

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