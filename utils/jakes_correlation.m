function rho = jakes_correlation(velocity_kmh, fc_Hz, delta_t_s)
% JAKES_CORRELATION  Temporal correlation coefficient for a flat-fading
% Rayleigh channel under the classical Jakes (uniform-scattering) model.
%
%   rho = jakes_correlation(v_kmh, fc, dt)
%
%   rho = J_0(2 * pi * f_d * dt)        with  f_d = v / lambda
%
%   v_kmh : user velocity [km/h]
%   fc_Hz : carrier frequency [Hz]
%   dt    : time lag [s]
%
%   Output rho in [-1, 1]; rho = 1 at dt = 0, oscillates and decays
%   as the lag grows. besselj(0, .) is the standard MATLAB form.

    c = physconst('LightSpeed');
    v = velocity_kmh / 3.6;             % km/h -> m/s
    f_d = v * fc_Hz / c;                % maximum Doppler shift [Hz]
    rho = besselj(0, 2*pi*f_d*delta_t_s);
end
