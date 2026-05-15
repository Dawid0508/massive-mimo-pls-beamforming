function save_figure(fig_handle, filename, output_dir)
% SAVE_FIGURE  Persist a figure with dark background to /results.
%
%   save_figure(fig_handle, filename)
%   save_figure(fig_handle, filename, output_dir)   % optional override
%
%   Applies apply_plot_style (black background, light axes, no black strokes)
%   before export.

    if nargin < 1 || isempty(fig_handle)
        fig_handle = gcf;
    end

    apply_plot_style(fig_handle);

    p = default_params();
    if nargin >= 3 && ~isempty(output_dir)
        out_dir = output_dir;
    else
        out_dir = p.results_dir;
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    out_path = fullfile(out_dir, [filename, '.png']);
    c = pls_colors();
    exportgraphics(fig_handle, out_path, 'Resolution', 200, ...
                   'BackgroundColor', c.bg);
    fprintf('[save_figure] saved -> %s\n', out_path);
end
