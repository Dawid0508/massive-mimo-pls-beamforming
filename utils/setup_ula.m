function [ula, sv, lambda] = setup_ula(Nt, fc)
% SETUP_ULA  Build a half-wavelength Uniform Linear Array and its
% steering-vector object for a given carrier.
%
%   [ula, sv, lambda] = setup_ula(Nt, fc)
%
%   Requires: Phased Array System Toolbox (see assert_requirements).

    assert_requirements('Phased_Array_System_Toolbox');

    c      = physconst('LightSpeed');
    lambda = c / fc;
    ula    = phased.ULA('NumElements', Nt, 'ElementSpacing', lambda/2);
    sv     = phased.SteeringVector('SensorArray', ula);
end
