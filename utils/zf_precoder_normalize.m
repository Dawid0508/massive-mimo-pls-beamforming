function W = zf_precoder_normalize(W_raw, P_rx, mode)
% ZF_PRECODER_NORMALIZE  Scale a raw ZF precoder under two power constraints.
%
%   W = zf_precoder_normalize(W_raw, P_rx, mode)
%
%   'matrix' : ||W||_F^2 = P_rx  (joint budget — favors weak streams)
%   'vector' : ||W(:,k)||^2 = P_rx/K for each k  (equal per-user power)

    switch lower(mode)
        case 'matrix'
            nf = norm(W_raw, 'fro');
            if nf < 1e-12
                W = zeros(size(W_raw));
            else
                W = W_raw / nf * sqrt(P_rx);
            end

        case 'vector'
            [Nt, K] = size(W_raw);
            W = zeros(Nt, K);
            pk = P_rx / K;
            for k = 1:K
                col = W_raw(:, k);
                nk = norm(col);
                if nk > 1e-12
                    W(:, k) = col / nk * sqrt(pk);
                end
            end

        otherwise
            error('zf_precoder_normalize:mode', ...
                'mode must be ''matrix'' or ''vector''.');
    end
end
