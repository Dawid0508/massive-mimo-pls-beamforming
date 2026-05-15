function [SR_sum, J] = compute_zf_secrecy_metrics(H, h_eve, W, noise_var)
% COMPUTE_ZF_SECRECY_METRICS  Sum-rate and Jain index for ZF precoding.
%
%   SR_sum : sum of per-user secrecy rates  max(0, R_b(k) - R_e(k))
%   J      : Jain fairness over *Bob* rates R_b(k). Normalization studies
%            compare fairness among legitimate users; secrecy rates are
%            often all zero if Jain is taken on R_s with a strong Eve.

    K = size(H, 2);
    R_b = zeros(K, 1);
    R_e = zeros(K, 1);
    for k = 1:K
        sig   = abs(H(:,k)' * W(:,k))^2;
        intf  = sum(abs(H(:,k)' * W).^2) - sig;
        R_b(k) = log2(1 + sig / (intf + noise_var));

        sig_e  = abs(h_eve' * W(:,k))^2;
        intf_e = sum(abs(h_eve' * W).^2) - sig_e;
        R_e(k) = log2(1 + sig_e / (intf_e + noise_var));
    end

    R_s = secrecy_rate(R_b, R_e);
    SR_sum = sum(R_s);
    J = jains_fairness(R_b);
end
