function save_figure(fig_handle, filename)
% SAVE_FIGURE  Persist a figure with white background to /results.
%
%   save_figure(fig_handle, filename)
%
%   fig_handle - handle returned by figure(); if empty, uses gcf
%   filename   - basename without extension; resolved against /results

    if nargin < 1 || isempty(fig_handle)
        fig_handle = gcf;
    end

    set(fig_handle, 'Color', 'w');
    set(fig_handle, 'InvertHardcopy', 'off');   % preserve white bg on print

    p = default_params();
    if ~exist(p.results_dir, 'dir')
        mkdir(p.results_dir);
    end
    out_path = fullfile(p.results_dir, [filename, '.png']);
    exportgraphics(fig_handle, out_path, 'Resolution', 200, ...
                   'BackgroundColor', 'white');
    fprintf('[save_figure] saved -> %s\n', out_path);
end
