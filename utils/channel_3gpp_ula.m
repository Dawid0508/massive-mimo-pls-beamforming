function h = channel_3gpp_ula(sv, fc, theta_los_deg, band_tag)
% CHANNEL_3GPP_ULA  TR 38.901-inspired CDL cluster model on a ULA.
%
%   h = channel_3gpp_ula(sv, fc, theta_los_deg, band_tag)
%
%   band_tag : 'sub6'  -> CDL-A-like (rich NLOS, K = 7 dB)
%              'mmwave' -> CDL-D-like (LOS street canyon, K = 18 dB)
%
%   Returns an Nt x 1 column vector. Random phase per cluster is drawn
%   on each call (Monte-Carlo fading).

    a_los = step(sv, fc, theta_los_deg);
    a_los = a_los(:) / norm(a_los);

    if nargin < 4 || isempty(band_tag)
        if fc >= 20e9
            band_tag = 'mmwave';
        else
            band_tag = 'sub6';
        end
    end

    switch lower(band_tag)
        case {'sub6', '6ghz', 'sub-6'}
            prof  = cdl_profile_sub6();
            K_lin = 10^(7/10);    % CDL-A typical NLOS urban
        case {'mmwave', '28ghz', 'mmw'}
            prof  = cdl_profile_mmwave();
            K_lin = 10^(18/10);   % CDL-D strong LOS
        otherwise
            error('channel_3gpp_ula:band', 'band_tag must be ''sub6'' or ''mmwave''.');
    end

    h_nlos = zeros(numel(a_los), 1);
    for c = 1:numel(prof.pwr_dB)
        theta_c = theta_los_deg + prof.aoa_deg(c);
        a_c = step(sv, fc, theta_c);
        a_c = a_c(:) / norm(a_c);
        h_nlos = h_nlos + 10^(prof.pwr_dB(c)/20) * exp(1j*2*pi*rand) * a_c;
    end
    if norm(h_nlos) > 1e-12
        h_nlos = h_nlos / norm(h_nlos);
    end

    h = sqrt(K_lin/(K_lin+1)) * a_los + sqrt(1/(K_lin+1)) * h_nlos;
    h = h(:);
end

% -------------------------------------------------------------------------
function prof = cdl_profile_sub6()
% Reduced CDL-A (38.901 Table 7.7.1): wider angular spread at sub-6 GHz.
    prof.aoa_deg = [-22.8, -5.2, 8.6, 18.4, -31.5, 27.1, -12.4, 35.2, ...
                     -41.0, 11.3, -8.9, 22.7, -18.6, 42.5, -28.3, 5.8, ...
                     -15.1, 31.9, -36.4, 14.2, -25.7, 19.6, -9.4, 38.1];
    prof.pwr_dB  = [-13.4, -18.8, -21.0, -22.8, -17.9, -20.1, -21.9, ...
                     -22.9, -18.6, -20.0, -21.8, -19.2, -21.7, -22.6, ...
                     -17.8, -19.9, -21.5, -22.4, -18.2, -20.3, -21.6, ...
                     -22.3, -19.0, -22.1];
end

function prof = cdl_profile_mmwave()
% Reduced CDL-D (38.901): dominant LOS + few strong clusters at mmWave.
    prof.aoa_deg = [0.0, 9.3, -9.3, 18.7, -18.7, 4.6, -4.6];
    prof.pwr_dB  = [0.0, -13.5, -13.5, -18.8, -18.8, -21.5, -21.5];
end
