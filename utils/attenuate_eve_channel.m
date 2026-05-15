function h_e = attenuate_eve_channel(h_e, atten_dB)
% ATTENUATE_EVE_CHANNEL  Scale Eve's channel for abstract Rayleigh scenarios.
%
%   A single-antenna Eve with unit-variance Rayleigh fading otherwise
%   often decodes every stream better than Bob, driving all secrecy rates
%   to zero. atten_dB weakens Eve so PLS metrics remain informative.

    if nargin < 2 || isempty(atten_dB)
        p = default_params();
        atten_dB = p.eve_attn_dB;
    end
    h_e = h_e * 10.^(-atten_dB / 20);
end
