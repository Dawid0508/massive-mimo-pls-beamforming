function c = pls_colors()
% PLS_COLORS  Consistent palette for dark-background PLS figures.
%
%   c = pls_colors()
%
%   Use for Bob/Eve markers, band curves, and reference lines so nothing
%   is drawn in black on a black figure background.

    c.bg      = [0.07 0.07 0.09];
    c.fg      = [0.93 0.93 0.96];
    c.grid    = [0.38 0.38 0.44];

    c.bob     = [0.25 0.88 0.48];
    c.eve     = [1.00 0.48 0.58];   % salmon — distinct from Bob green
    c.beam    = [0.98 0.78 0.22];
    c.ref     = [0.62 0.82 1.00];   % reference / theory lines
    c.perfect = [0.92 0.88 0.45];   % "ideal" baselines

    c.sub6    = [0.35 0.65 1.00];
    c.mmwave  = [1.00 0.42 0.35];
    c.bs      = [0.75 0.78 0.85];
end
