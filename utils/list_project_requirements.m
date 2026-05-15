function list_project_requirements()
% LIST_PROJECT_REQUIREMENTS  Print MATLAB / toolbox dependencies for this repo.
%
%   >> list_project_requirements

    fprintf('=== PLS project dependencies ===\n\n');

    fprintf('MATLAB\n');
    fprintf('  Minimum release : R2016b (9.1) — strings in run_all, isgraphics in plots\n');
    fprintf('  Your release    : %s\n\n', version);

    fprintf('Required add-on\n');
    print_toolbox('Phased Array System Toolbox', 'phased', ...
        {'physconst', 'phased.ULA', 'phased.SteeringVector', 'step(SteeringVector,...)'} );

    fprintf('Base MATLAB only (no extra toolbox)\n');
    fprintf('  Linear algebra : pinv, toeplitz, sqrtm, \\ \n');
    fprintf('  Special funcs  : besselj (Jakes model in jakes_correlation.m)\n');
    fprintf('  Plotting       : figure, plot, saveas (export via save_figure.m)\n\n');

    fprintf('Not used\n');
    fprintf('  Communications Toolbox, WLAN Toolbox, 5G Toolbox, Antenna Toolbox\n\n');

    fprintf('Entry points that enforce requirements\n');
    fprintf('  run_all.m, pls_startup.m, setup_ula.m, default_params.m\n\n');

    fprintf('Scenarios using Phased Array (ULA + 3GPP CDL steering)\n');
    tags = { ...
        'Ghz6_band_vs_mmWave_band', ...
        'sim_moving_bob', ...
        'sim_pilot_contamination', ...
        'sim_colluding_eavesdroppers', ...
        'sim_artificial_noise', ...
        'sim_phase_noise', ...
        'sim_location_error'};
    for k = 1:numel(tags)
        fprintf('  - %s\n', tags{k});
    end
    fprintf('\nScenarios using physconst via default_params only (Rayleigh / abstract)\n');
    tags2 = { ...
        'sim_spatial_correlation', ...
        'sim_fairness_normalization', ...
        'sim_channel_hardening', ...
        'sim_low_res_dac', ...
        'sim_pilot_jamming', ...
        'sim_csi_aging'};
    for k = 1:numel(tags2)
        fprintf('  - %s\n', tags2{k});
    end
    fprintf('\n(generate_topologies.m needs default_params → Phased Array)\n');
end

function print_toolbox(displayName, shortName, symbols)
    lic = license('test', 'Phased_Array_System_Toolbox');
    v = ver(shortName);
    installed = lic && ~isempty(v);
    if installed
        fprintf('  [OK]   %s (%s)\n', displayName, v(1).Version);
    else
        fprintf('  [MISS] %s\n', displayName);
    end
    for k = 1:numel(symbols)
        fprintf('         %s\n', symbols{k});
    end
    fprintf('\n');
end
