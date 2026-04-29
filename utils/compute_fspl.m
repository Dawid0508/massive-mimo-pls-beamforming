function [PL_lin, PL_dB] = compute_fspl(dist_m, fc_Hz)
% COMPUTE_FSPL  Free-Space Path Loss (Friis), linear and dB.
%
%   FSPL_dB = 20*log10(d) + 20*log10(f) - 147.55
%
%   The constant -147.55 collapses the 4*pi/c term in the Friis equation
%   when d is in metres and f is in Hz.

    PL_dB  = 20*log10(dist_m) + 20*log10(fc_Hz) - 147.55;
    PL_lin = 10.^(PL_dB/10);
end
