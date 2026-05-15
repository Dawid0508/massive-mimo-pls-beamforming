function fig = plot_topology(cfg)
% PLOT_TOPOLOGY  Top-down link geometry for a PLS scenario (dark theme).

    if numel(cfg) ~= 1
        error('plot_topology:cfg', ...
            'Expected one scalar config struct (got %d elements).', numel(cfg));
    end

    if isfield(cfg, 'bs_angle_deg') && ~isempty(cfg.bs_angle_deg)
        bs_angle_deg = cfg.bs_angle_deg;
    else
        bs_angle_deg = 90;
    end

    c = pls_colors();
    fig = figure('Color', c.bg, 'Position', [80 80 560 520], 'Visible', 'off');
    hold on; axis equal; grid on; box on;
    set(gca, 'Color', c.bg, 'XColor', c.fg, 'YColor', c.fg, 'GridColor', c.grid);

    plot(0, 0, 's', 'MarkerSize', 12, 'MarkerFaceColor', c.bs, ...
         'MarkerEdgeColor', c.fg, 'LineWidth', 1.2);
    text(0.3, -0.8, 'BS (ULA)', 'FontSize', 9, 'FontWeight', 'bold', 'Color', c.fg);

    br = deg2rad(bs_angle_deg);
    quiver(0, 0, 2*cos(br), 2*sin(br), 0, 'Color', c.beam, ...
           'LineWidth', 1.2, 'MaxHeadSize', 0.8);
    text(2.2*cos(br), 2.2*sin(br), 'boresight', 'FontSize', 8, 'Color', c.beam);

    max_r = 5;
    nodes = cfg.nodes;
    for n = 1:numel(nodes)
        nd = nodes(n);
        th = deg2rad(nd.angle_deg);
        x  = nd.dist_m * cos(th);
        y  = nd.dist_m * sin(th);
        if isfield(nd, 'color') && ~isempty(nd.color)
            col = nd.color;
        else
            col = c.sub6;
        end
        plot(x, y, 'o', 'MarkerSize', 10, 'MarkerFaceColor', col, ...
             'MarkerEdgeColor', min(col + 0.25, 1));
        text(x + 0.6, y + 0.4, nd.name, 'FontSize', 9, 'Color', col);
        max_r = max(max_r, nd.dist_m);
    end

    if isfield(cfg, 'beam_angle_deg') && ~isempty(cfg.beam_angle_deg)
        bw = 6;
        if isfield(cfg, 'beam_width_deg') && ~isempty(cfg.beam_width_deg)
            bw = cfg.beam_width_deg;
        end
        ba = deg2rad(cfg.beam_angle_deg);
        r_beam = max_r * 0.85;
        th1 = ba - deg2rad(bw/2);
        th2 = ba + deg2rad(bw/2);
        fill([0, r_beam*cos(th1), r_beam*cos(th2)], ...
             [0, r_beam*sin(th1), r_beam*sin(th2)], ...
             c.beam, 'FaceAlpha', 0.22, 'EdgeColor', c.beam, 'LineStyle', '--');
        text(r_beam*0.55*cos(ba), r_beam*0.55*sin(ba), 'fixed beam', ...
             'FontSize', 8, 'Color', c.beam);
    end

    xlim([-max_r-2, max_r+2]);
    ylim([-max_r-2, max_r+2]);
    xlabel('x [m]', 'Color', c.fg);
    ylabel('y [m]', 'Color', c.fg);
    title(cfg.title, 'Interpreter', 'none', 'Color', c.fg);
    legend({'BS'}, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
           'TextColor', c.fg, 'Color', c.bg, 'EdgeColor', c.grid);

    save_topology(fig, cfg.filename);
end
