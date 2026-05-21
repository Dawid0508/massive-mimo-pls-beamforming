function save_figure(fig_handle, filename, output_dir)
% SAVE_FIGURE  Persist a figure with white background to /results.
%
%   save_figure(fig_handle, filename)
%   save_figure(fig_handle, filename, output_dir)   % optional override

    if nargin < 1 || isempty(fig_handle)
        fig_handle = gcf;
    end
    
    % ZAKOMENTOWANE: Ta funkcja niszczyła Twoje białe tło i psuła osie!
    % apply_plot_style(fig_handle); 
    
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
    
    % ZMIANA: Twarde ustawienie eksportu na całkowicie białe tło ('white')
    exportgraphics(fig_handle, out_path, 'Resolution', 300, ...
                   'BackgroundColor', 'white');
                   
    fprintf('[save_figure] saved -> %s\n', out_path);
end