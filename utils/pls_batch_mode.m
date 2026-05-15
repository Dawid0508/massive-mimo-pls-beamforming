function tf = pls_batch_mode()
% PLS_BATCH_MODE  True while run_all is executing scenarios.
    tf = false;
    if evalin('base', 'exist(''PLS_BATCH_RUN'', ''var'')')
        tf = logical(evalin('base', 'PLS_BATCH_RUN'));
    end
end
