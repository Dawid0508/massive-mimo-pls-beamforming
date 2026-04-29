function J = jains_fairness(rates)
% JAINS_FAIRNESS  Jain's fairness index for a vector of per-user rates.
%
%   J = (sum r)^2 / (K * sum r^2)
%
%   Returns a value in [1/K, 1]. J = 1 means perfectly equal allocation.
%   When the input is identically zero we return 0 (degenerate; standard
%   convention to avoid 0/0).

    rates = rates(:);
    K = numel(rates);
    s1 = sum(rates);
    s2 = sum(rates.^2);
    if s2 < eps
        J = 0;
    else
        J = (s1^2) / (K * s2);
    end
end
