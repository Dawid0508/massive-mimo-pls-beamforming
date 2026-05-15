function apply_plot_style(fig)
% APPLY_PLOT_STYLE  Dark theme, palette fixes, reference-line label layout.
%
%   apply_plot_style()      % style gcf
%   apply_plot_style(fig)   % style a figure handle
%
%   Called from save_figure / save_topology before export. Recolors black
%   strokes, styles legends, and places xline/yline labels away from legends
%   without resizing axes (labels opposite the legend side).

    if nargin < 1 || isempty(fig)
        fig = gcf;
    end
    if ~isgraphics(fig)
        return;
    end

    c = pls_colors();
    set(fig, 'Color', c.bg, 'InvertHardcopy', 'off');

    axs = findall(fig, 'Type', 'axes');
    for k = 1:numel(axs)
        ax = axs(k);
        set(ax, 'Color', c.bg, 'XColor', c.fg, 'YColor', c.fg, 'ZColor', c.fg, ...
            'GridColor', c.grid, 'MinorGridColor', c.grid);
        style_axis_label(ax.Title,  c.fg);
        style_axis_label(ax.XLabel, c.fg);
        style_axis_label(ax.YLabel, c.fg);
        style_axis_label(ax.ZLabel, c.fg);
        lg = ax.Legend;
        if ~isempty(lg) && isgraphics(lg)
            set(lg, 'Color', c.bg, 'TextColor', c.fg, 'EdgeColor', c.grid);
        end
    end

    txts = findall(fig, 'Type', 'text');
    for k = 1:numel(txts)
        if is_black(get_color_rgb(txts(k).Color))
            set(txts(k), 'Color', c.fg);
        end
    end

    lines = findall(fig, 'Type', 'line');
    for k = 1:numel(lines)
        ln = lines(k);
        col = get_color_rgb(ln.Color);
        if is_black(col)
            set(ln, 'Color', infer_line_color(ln, c));
        end
        mfc = get_color_rgb(ln.MarkerFaceColor);
        if is_black(mfc)
            set(ln, 'MarkerFaceColor', ln.Color);
        end
        mec = get_color_rgb(ln.MarkerEdgeColor);
        if is_black(mec)
            set(ln, 'MarkerEdgeColor', ln.Color);
        end
    end

    consts = findall(fig, 'Type', 'ConstantLine');
    for k = 1:numel(consts)
        cl = consts(k);
        col = get_color_rgb(cl.Color);
        if is_black(col)
            cl.Color = infer_constant_color(cl, c);
        end
        style_constant_line_label(cl);
    end

    for k = 1:numel(axs)
        layout_axis_ref_labels(axs(k));
    end
end

% -------------------------------------------------------------------------
function layout_axis_ref_labels(ax)
    cls = findall(ax, 'Type', 'ConstantLine');
    vPref = ref_label_pref(ax);
    doStagger = ref_stagger_pref(ax, cls);

    for c = 1:numel(cls)
        if has_constant_label(cls(c))
            place_reference_line_label(cls(c), ax, vPref, constant_line_ref_side(cls(c), ax));
        end
    end
    if doStagger
        stagger_axis_ref_labels(ax);
    end
end

function vPref = ref_label_pref(ax)
    vPref = 'auto';
    if isstruct(ax.UserData) && isfield(ax.UserData, 'plsRefLabelV')
        vPref = ax.UserData.plsRefLabelV;
    end
end

function side = constant_line_ref_side(cl, ax)
% Per-line appdata (mark_bob/mark_eve) overrides axes-wide default.
    side = 'auto';
    try
        if isprop(cl, 'Parent') && ~isempty(cl.Parent)
            s = getappdata(cl, 'plsConstLabelSide');
            if ~isempty(s) && ischar(s)
                side = lower(s);
                return;
            end
        end
    catch %#ok<*CTCH>
    end
    if isstruct(ax.UserData) && isfield(ax.UserData, 'plsRefSideOfLine')
        s = ax.UserData.plsRefSideOfLine;
        if ~strcmpi(s, 'auto')
            side = lower(s);
        end
    end
end

function orient = ref_label_orient(ax)
    orient = 'horizontal';
    if isstruct(ax.UserData) && isfield(ax.UserData, 'plsRefLabelOrient')
        orient = ax.UserData.plsRefLabelOrient;
    end
end

function tf = ref_stagger_pref(ax, cls)
    tf = false;
    nLab = 0;
    for c = 1:numel(cls)
        if has_constant_label(cls(c))
            nLab = nLab + 1;
        end
    end
    if nLab >= 2
        tf = true;
    end
    if isstruct(ax.UserData) && isfield(ax.UserData, 'plsStaggerRef')
        tf = logical(ax.UserData.plsStaggerRef);
    end
end

function place_reference_line_label(cl, ax, vPref, sidePref)
% Place xline/yline label; vPref is 'top', 'bottom', or 'auto'.
% sidePref for xline: 'left' (text left of line, HA=right), 'right', or 'auto'.

    if nargin < 3 || isempty(vPref)
        vPref = 'auto';
    end
    if nargin < 4 || isempty(sidePref)
        sidePref = 'auto';
    end
    cl.LabelOrientation = ref_label_orient(ax);

    if isprop(cl, 'Orientation') && strcmpi(cl.Orientation, 'horizontal')
        cl.LabelHorizontalAlignment = 'right';
        if strcmpi(vPref, 'auto')
            cl.LabelVerticalAlignment = label_vertical_vs_legend(ax);
        else
            cl.LabelVerticalAlignment = vPref;
        end
        return;
    end

    if strcmpi(vPref, 'auto')
        vAlign = label_vertical_vs_legend(ax);
    else
        vAlign = vPref;
    end
    cl.LabelVerticalAlignment = vAlign;
    if ~isprop(cl, 'LabelHorizontalAlignment')
        return;
    end
    xv = cl.Value;
    if isnumeric(xv) && isscalar(xv)
        if strcmpi(sidePref, 'left')
            cl.LabelHorizontalAlignment = 'right';
        elseif strcmpi(sidePref, 'right')
            cl.LabelHorizontalAlignment = 'left';
        else
            xl = ax.XLim;
            if xv <= mean(xl)
                cl.LabelHorizontalAlignment = 'right';
            else
                cl.LabelHorizontalAlignment = 'left';
            end
        end
    end
end

function stagger_axis_ref_labels(ax)
% Separate overlapping reference-line labels on the same axes.

    vert = gobjects(0);
    horiz = gobjects(0);
    cls = findall(ax, 'Type', 'ConstantLine');
    for c = 1:numel(cls)
        if ~has_constant_label(cls(c))
            continue;
        end
        if isprop(cls(c), 'Orientation') && strcmpi(cls(c).Orientation, 'horizontal')
            horiz(end+1) = cls(c); %#ok<AGROW>
        else
            vert(end+1) = cls(c); %#ok<AGROW>
        end
    end

    if numel(vert) >= 2
        [~, ord] = sort([vert.Value]);
        vert = vert(ord);
        yl = ax.YLim;
        dy = 0.06 * (yl(2) - yl(1));
        nudgeIdx = 0;
        vPrefVert = ref_label_pref(ax);
        for i = 1:numel(vert)
            if strcmpi(vPrefVert, 'auto')
                vert(i).LabelVerticalAlignment = 'top';
            else
                vert(i).LabelVerticalAlignment = vPrefVert;
            end
            side = constant_line_ref_side(vert(i), ax);
            if strcmpi(side, 'left')
                vert(i).LabelHorizontalAlignment = 'right';
            elseif strcmpi(side, 'right')
                vert(i).LabelHorizontalAlignment = 'left';
            elseif mod(i, 2) == 1
                vert(i).LabelHorizontalAlignment = 'right';
            else
                vert(i).LabelHorizontalAlignment = 'left';
            end
            if strcmpi(ref_label_orient(ax), 'aligned')
                nudge_constant_line_label_y(vert(i), ax, nudgeIdx * dy);
                nudgeIdx = nudgeIdx + 1;
            else
                nudge_constant_line_label_y(vert(i), ax, (i - 1) * dy);
            end
        end
    end

    if numel(horiz) >= 2
        [~, ord] = sort([horiz.Value]);
        horiz = horiz(ord);
        xl = ax.XLim;
        dx = 0.025 * (xl(2) - xl(1));
        for i = 1:numel(horiz)
            horiz(i).LabelHorizontalAlignment = 'right';
            if mod(i, 2) == 1
                horiz(i).LabelVerticalAlignment = 'top';
            else
                horiz(i).LabelVerticalAlignment = 'bottom';
            end
            nudge_constant_line_label_x(horiz(i), ax, (i - 1) * dx);
        end
    elseif numel(horiz) == 1
        horiz(1).LabelHorizontalAlignment = 'right';
        horiz(1).LabelVerticalAlignment = 'bottom';
    end
end

function nudge_constant_line_label_y(cl, ax, dy)
% Nudge ConstantLine label vertically in data units if supported.
    if ~isscalar(dy) || dy == 0
        return;
    end
    if ~isprop(cl, 'Label')
        return;
    end
    lab = cl.Label;
    if numel(lab) < 1
        return;
    end
    lab = lab(1);
    if ~isgraphics(lab)
        return;
    end

    if ~isprop(lab, 'Position') || ~isprop(lab, 'Units')
        return;
    end
    try
        was = lab.Units;
        lab.Units = 'data';
        p = lab.Position;
        if numel(p) >= 2
            lab.Position = [p(1), p(2) + dy, p(3:min(3, end))];
        end
        lab.Units = was;
    catch
    end
end

function nudge_constant_line_label_x(cl, ax, dx)
% Nudge ConstantLine label horizontally in data units if supported.
    if ~isscalar(dx) || dx == 0
        return;
    end
    if ~isprop(cl, 'Label')
        return;
    end
    lab = cl.Label;
    if numel(lab) < 1
        return;
    end
    lab = lab(1);
    if ~isgraphics(lab)
        return;
    end

    if ~isprop(lab, 'Position') || ~isprop(lab, 'Units')
        return;
    end
    try
        was = lab.Units;
        lab.Units = 'data';
        p = lab.Position;
        if numel(p) >= 2
            lab.Position = [p(1) + dx, p(2), p(3:min(3, end))];
        end
        lab.Units = was;
    catch
    end
end

function vAlign = label_vertical_vs_legend(ax)
    vAlign = 'bottom';
    lg = ax.Legend;
    if isempty(lg) || ~all(isgraphics(lg(:)))
        return;
    end
    loc = lower(char(lg.Location));
    if any(strcmp(loc, {'south', 'southwest', 'southeast'}))
        vAlign = 'top';
    end
end

function tf = has_constant_label(cl)
    tf = false;
    if isempty(cl.Label)
        return;
    end
    lab = cl.Label;
    if ischar(lab) || isstring(lab)
        tf = strlength(strtrim(lab)) > 0;
    elseif numel(lab) >= 1 && isgraphics(lab(1)) && isprop(lab(1), 'String')
        s = lab(1).String;
        if iscell(s), s = s{1}; end
        tf = strlength(strtrim(string(s))) > 0;
    end
end

function rgb = get_color_rgb(col)
    if ischar(col) || isstring(col)
        switch char(col)
            case 'k', rgb = [0 0 0];
            case 'w', rgb = [1 1 1];
            case 'r', rgb = [1 0 0];
            case 'g', rgb = [0 1 0];
            case 'b', rgb = [0 0 1];
            case 'c', rgb = [0 1 1];
            case 'm', rgb = [1 0 1];
            case 'y', rgb = [1 1 0];
            otherwise, rgb = [0 0 0];
        end
    else
        rgb = col(1:min(3, numel(col)));
        if numel(rgb) < 3
            rgb = [0 0 0];
        end
    end
end

function tf = is_black(rgb)
    tf = all(rgb < 0.15);
end

function col = infer_line_color(ln, c)
    lbl = '';
    if ~isempty(ln.DisplayName)
        lbl = lower(ln.DisplayName);
    end
    col = pick_from_label(lbl, c);
    if is_black(col)
        col = c.ref;
    end
end

function col = infer_constant_color(cl, c)
    col = pick_from_label(constant_line_label_string(cl), c);
    if is_black(col)
        col = c.ref;
    end
end

function style_axis_label(h, fg)
    if ~isempty(h) && isgraphics(h)
        set(h, 'Color', fg);
    end
end

function style_constant_line_label(cl)
    lab = cl.Label;
    if isempty(lab)
        return;
    end
    if isgraphics(lab)
        set(lab, 'Color', cl.Color, 'FontWeight', 'bold');
    end
end

function lbl = constant_line_label_string(cl)
    lbl = '';
    if isempty(cl.Label)
        return;
    end
    lab = cl.Label;
    if ischar(lab) || isstring(lab)
        lbl = lower(char(lab));
    elseif isgraphics(lab)
        lbl = lower(lab.String);
    end
end

function col = pick_from_label(lbl, c)
    if contains(lbl, 'eve')
        col = c.eve;
    elseif contains(lbl, 'bob')
        col = c.bob;
    elseif contains(lbl, 'beam') || contains(lbl, 'boresight')
        col = c.beam;
    elseif contains(lbl, 'perfect') || contains(lbl, 'ideal')
        col = c.perfect;
    elseif contains(lbl, 'phi') || contains(lbl, 'fixed')
        col = c.beam;
    else
        col = [0 0 0];
    end
end
