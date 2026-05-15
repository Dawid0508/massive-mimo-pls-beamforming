function pls_axis_prefs(ax, varargin)
% PLS_AXIS_PREFS  Per-axes layout hints for apply_plot_style.
%
%   pls_axis_prefs(ax, 'refLabelV', 'top')     % xline/yline labels at top
%   pls_axis_prefs(ax, 'refLabelV', 'bottom')
%   pls_axis_prefs(ax, 'refSideOfLine', 'left')   % xline text to left of line (HA=right)
%   pls_axis_prefs(ax, 'refSideOfLine', 'right')  % xline text to right of line (HA=left)
%   pls_axis_prefs(ax, 'refLabelOrient', 'aligned')  % text parallel to vertical xline
%   pls_axis_prefs(ax, 'staggerRef', true)    % separate overlapping ref labels

    if nargin < 1 || isempty(ax)
        ax = gca;
    end

    p = inputParser;
    addParameter(p, 'refLabelV', '', @(s) isempty(s) || any(strcmpi(s, {'top', 'bottom', 'auto'})));
    addParameter(p, 'refSideOfLine', '', @(s) isempty(s) || any(strcmpi(s, {'left', 'right', 'auto'})));
    addParameter(p, 'refLabelOrient', '', @(s) isempty(s) || any(strcmpi(s, {'horizontal', 'aligned'})));
    addParameter(p, 'staggerRef', [], @(x) isempty(x) || islogical(x));
    parse(p, varargin{:});

    if isstruct(ax.UserData)
        uds = ax.UserData;
    else
        uds = struct();
    end

    if ~isempty(p.Results.refLabelV)
        uds.plsRefLabelV = lower(char(p.Results.refLabelV));
    end
    if ~isempty(p.Results.refSideOfLine)
        uds.plsRefSideOfLine = lower(char(p.Results.refSideOfLine));
    end
    if ~isempty(p.Results.refLabelOrient)
        uds.plsRefLabelOrient = lower(char(p.Results.refLabelOrient));
    end
    if ~isempty(p.Results.staggerRef)
        uds.plsStaggerRef = logical(p.Results.staggerRef);
    end
    ax.UserData = uds;
end
