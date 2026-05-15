function run_all()
% =========================================================================
% RUN_ALL  Execute every PLS scenario and refresh the figures in /results.
% -------------------------------------------------------------------------
% Usage: from the repository root, just type
%        >> run_all
%
% Console output is preserved during the batch (scenarios skip clc).
% A full transcript is also written to results/run_all.log.
% Topology maps: run generate_topologies separately.
% =========================================================================

    repo_root = fileparts(mfilename('fullpath'));
    utils_dir = fullfile(repo_root, 'utils');
    sim_dir   = fullfile(repo_root, 'simulations');
    addpath(utils_dir, sim_dir);
    evalin('base', sprintf('addpath(''%s''); addpath(''%s'');', ...
        escape_for_matlab_str(utils_dir), escape_for_matlab_str(sim_dir)));

    assert_requirements();

    p = default_params();
    if ~exist(p.results_dir, 'dir')
        mkdir(p.results_dir);
    end
    log_path = fullfile(p.results_dir, 'run_all.log');

    assignin('base', 'PLS_BATCH_RUN', true);
    diary(log_path);
    diary on;

    fprintf('=== PLS batch run started %s ===\n', datestr(now));
    fprintf('Log file: %s\n\n', log_path);

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
        'sim_location_error', ...
        'sim_moving_bob'};

    n_ok   = 0;
    n_fail = 0;
    status = strings(numel(scripts), 1);
    elapsed_s = zeros(numel(scripts), 1);

    t0 = tic;
    for scIdx = 1:numel(scripts)
        name = scripts{scIdx};
        fprintf('\n=== [%d/%d] %s ===\n', scIdx, numel(scripts), name);
        t_script = tic;
        try
            script_path = fullfile(repo_root, 'simulations', [name, '.m']);
            % Run in base workspace so scenario variables cannot clobber
            % this function's loop counters or summary arrays.
            evalin('base', sprintf('run(''%s'');', script_path));
            elapsed_s(scIdx) = toc(t_script);
            status(scIdx) = "OK";
            n_ok = n_ok + 1;
            fprintf('[OK] %s  (%.1f s)\n', name, elapsed_s(scIdx));
        catch ME
            elapsed_s(scIdx) = toc(t_script);
            status(scIdx) = "FAIL";
            n_fail = n_fail + 1;
            fprintf(2, '[FAIL] %s  (%.1f s)\n', name, elapsed_s(scIdx));
            fprintf(2, '       %s\n', ME.message);
            if ~isempty(ME.stack)
                st = ME.stack(1);
                fprintf(2, '       %s (line %d)\n', st.file, st.line);
            end
        end
        evalin('base', 'close all;');
        drawnow;
    end

    total_s = toc(t0);
    fprintf('\n=== Summary (%d OK, %d FAIL, %.1f s total) ===\n', n_ok, n_fail, total_s);
    for scIdx = 1:numel(scripts)
        fprintf('  %-32s  %4s  %6.1f s\n', scripts{scIdx}, char(status(scIdx)), ...
            elapsed_s(scIdx));
    end
    fprintf('\nFull transcript: %s\n', log_path);
    if n_fail > 0
        fprintf(2, '\nSome scenarios failed — see log for details.\n');
    end

    diary off;
    assignin('base', 'PLS_BATCH_RUN', false);
end

function s = escape_for_matlab_str(p)
    s = strrep(p, '\', '\\');
    s = strrep(s, '''', '''''');
end
