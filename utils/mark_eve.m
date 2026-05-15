function h = mark_eve(x, label, labelSideOfLine)
% MARK_EVE  Vertical marker on beam-pattern plots (dark-theme safe).
%
%   mark_eve(x, label, side) — side 'left' | 'right' positions the label text.
    c = pls_colors();
    if nargin < 2 || isempty(label)
        label = 'Eve';
    end
    h = xline(x, ':', label, 'Color', c.eve, 'LineWidth', 1.5);
    if nargin >= 3 && ~isempty(labelSideOfLine)
        try
            setappdata(h, 'plsConstLabelSide', lower(char(labelSideOfLine)));
        catch %#ok<*CTCH>
        end
    end
end
