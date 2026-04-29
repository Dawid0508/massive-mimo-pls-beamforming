function R_s = secrecy_rate(R_bob, R_eve)
% SECRECY_RATE  Element-wise non-negative Secrecy Rate.
%
%   R_s(k) = max(0, R_bob(k) - R_eve(k))
%
%   Accepts scalars or vectors. The clamp at 0 reflects the standard PLS
%   definition: if the eavesdropper's channel is better than the legitimate
%   one, no positive-rate secrecy code exists (Wyner 1975).

    R_s = max(0, R_bob - R_eve);
end
