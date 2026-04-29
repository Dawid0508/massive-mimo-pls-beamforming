function run_all()
% =========================================================================
% RUN_ALL  Execute every PLS scenario and refresh the figures in /results.
% -------------------------------------------------------------------------
% Usage: from the repository root, just type
%        >> run_all
% Each scenario is run in its own try/catch so a single failure does not
% abort the whole batch.
% =========================================================================
    clc; close all;
    
    repo_root = fileparts(mfilename('fullpath'));
    addpath(fullfile(repo_root, 'utils'));
    addpath(fullfile(repo_root, 'simulations'));
    
    scripts = { ...
        'Ghz6_band_vs_mmWave_band', ...
        'sim_spatial_correlation', ...
        'sim_colluding_eavesdroppers', ...
        'sim_artificial_noise', ...
        'sim_pilot_contamination', ...
        'sim_fairness_normalization', ...
        'sim_phase_noise', ...
        'sim_channel_hardening', ...
        'sim_low_res_dac', ...
        'sim_pilot_jamming', ...
        'sim_csi_aging', ...
        'sim_location_error'};
        
    t0 = tic;
    for k = 1:numel(scripts)
        name = scripts{k};
        fprintf('\n=== [%d/%d] %s ===\n', k, numel(scripts), name);
        t = tic;
        try
            script_path = fullfile(repo_root, 'simulations', [name, '.m']);
            evalin('base', sprintf('run(''%s'');', script_path));
            
            fprintf('[OK] %s  (%.1f s)\n', name, toc(t));
        catch ME
            fprintf(2, '[FAIL] %s : %s\n', name, ME.message);
        end
        evalin('base', 'close all;');
    end
    fprintf('\nAll scenarios completed in %.1f s.\n', toc(t0));
end