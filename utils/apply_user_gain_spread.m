function H = apply_user_gain_spread(H, spread_dB)
% APPLY_USER_GAIN_SPREAD  Per-user large-scale fading on channel columns.
%
%   H = apply_user_gain_spread(H)
%   H = apply_user_gain_spread(H, spread_dB)
%
%   Without gain differences, i.i.d. Rayleigh ZF columns often have
%   nearly equal norms, so matrix and vector normalization collapse to
%   the same precoder. A fixed spread across users breaks that symmetry.

    if nargin < 2 || isempty(spread_dB)
        spread_dB = 12;
    end

    K = size(H, 2);
    beta_dB = linspace(-spread_dB/2, spread_dB/2, K);
    beta = 10.^(beta_dB(:)' / 20);
    H = H .* beta;
end
