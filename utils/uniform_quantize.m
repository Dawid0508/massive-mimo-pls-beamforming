function y = uniform_quantize(x, b)
% UNIFORM_QUANTIZE  Mid-rise uniform quantizer for complex baseband signal.
%
%   y = uniform_quantize(x, b)
%
%   x : complex matrix; real and imaginary parts are quantised
%       independently (typical I/Q DAC architecture).
%   b : number of bits per real branch. b = Inf returns x unchanged
%       (ideal DAC reference).
%
%   The clipping range is set to 3 * std(x) which is the textbook
%   choice for Gaussian-distributed signals (less than 0.3% clip
%   probability per branch).

    if isinf(b)
        y = x;
        return
    end

    L = 2^b;
    re = real(x); im = imag(x);

    % per-branch dynamic range from the empirical std (3-sigma loading)
    sigma_x = std([re(:); im(:)]);
    A     = 3 * sigma_x;
    delta = 2 * A / L;

    qr = delta * (floor(re/delta) + 0.5);
    qi = delta * (floor(im/delta) + 0.5);
    qr = max(min(qr, A - delta/2), -A + delta/2);
    qi = max(min(qi, A - delta/2), -A + delta/2);

    y = qr + 1j*qi;
end
