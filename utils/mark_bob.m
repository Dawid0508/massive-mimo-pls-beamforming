function h = mark_bob(x, label, labelSideOfLine)
% MARK_BOB  Vertical marker on beam-pattern plots (dark-theme safe).
%
%   mark_bob(x)
%   mark_bob(x, label)
%   mark_bob(x, label, side)    side 'left' | 'right': label sits on that side
%                               text is horizontally aligned accordingly.
    c = pls_colors();
    if nargin < 2 || isempty(label)
        label = 'Bob';
    end
    h = xline(x, ':', label, 'Color', c.bob, 'LineWidth', 1.5);
    if nargin >= 3 && ~isempty(labelSideOfLine)
        try
            setappdata(h, 'plsConstLabelSide', lower(char(labelSideOfLine)));
        catch %#ok<*CTCH>
        end
    end
end
