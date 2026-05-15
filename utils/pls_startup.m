function pls_startup()
% PLS_STARTUP  Per-scenario init: clear workspace, optional clc, close figures.
%
%   When run_all sets PLS_BATCH_RUN in the base workspace, clc is skipped so
%   console output from the batch run remains visible.

    utils_dir = fileparts(mfilename('fullpath'));
    addpath(utils_dir);
    assert_requirements();

    batch = pls_batch_mode();
    evalin('base', 'clear');
    if batch
        assignin('base', 'PLS_BATCH_RUN', true);
    end
    if ~batch
        clc;
    end
    close all;
end
