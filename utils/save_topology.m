function save_topology(fig_handle, filename)
% SAVE_TOPOLOGY  Export topology figure to /topology (dark theme).

    if nargin < 1 || isempty(fig_handle)
        fig_handle = gcf;
    end

    apply_plot_style(fig_handle);

    repo = fileparts(fileparts(mfilename('fullpath')));
    topo_dir = fullfile(repo, 'topology');
    if ~exist(topo_dir, 'dir')
        mkdir(topo_dir);
    end
    out_path = fullfile(topo_dir, [filename, '.png']);
    c = pls_colors();
    exportgraphics(fig_handle, out_path, 'Resolution', 200, ...
                   'BackgroundColor', c.bg);
    fprintf('[save_topology] saved -> %s\n', out_path);
    close(fig_handle);
end
