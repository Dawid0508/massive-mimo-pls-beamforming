function [ula, sv, lambda] = setup_ula(Nt, fc)
% SETUP_ULA  Build a half-wavelength Uniform Linear Array and its
% steering-vector object for a given carrier.
%
%   [ula, sv, lambda] = setup_ula(Nt, fc)
%
%   Nt     - number of antenna elements
%   fc     - carrier frequency [Hz]
%
%   ula    - phased.ULA object (half-wavelength spacing)
%   sv     - phased.SteeringVector bound to ula
%   lambda - wavelength [m]

    c      = physconst('LightSpeed');
    lambda = c / fc;
    ula    = phased.ULA('NumElements', Nt, 'ElementSpacing', lambda/2);
    sv     = phased.SteeringVector('SensorArray', ula);
end
