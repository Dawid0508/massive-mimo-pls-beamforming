function print_scenario_snr(varargin)
% PRINT_SCENARIO_SNR  Consistent TX/RX SNR report for all link actors.
%
%   print_scenario_snr('title', 'Pilot contamination @ 6 GHz', ...
%       'SNR_rx_dB', 20, 'dist_m', 30, 'fc_Hz', 6e9, ...
%       'actors', {'Bob', 'Eve'})
%
%   print_scenario_snr('title', 'Fairness ZF', 'SNR_rx_dB', 20, ...
%       'actors', {'Bob (K=8)', 'Eve'})     % no dist/fc -> normalized channel
%
%   SNR_rx_dB is the configured received SNR at Bob (after path loss when
%   dist_m and fc_Hz are set). Other actors use the same link budget unless
%   an actor struct overrides .SNR_rx_dB, .dist_m, or .fc_Hz.
%
%   Optional name-value pairs:
%     'title'       - scenario / band header
%     'SNR_rx_dB'   - scalar reference received SNR at Bob [dB]
%     'dist_m'      - link distance [m] (omit for abstract / normalized)
%     'fc_Hz'       - carrier [Hz] (omit with dist_m for abstract)
%     'actors'      - cellstr of labels, or struct array with fields:
%                     .name (required), optional .SNR_rx_dB, .dist_m, .fc_Hz
%     'notes'       - extra line printed after actors (char/string)

    opts = parse_inputs(varargin{:});

    fprintf('\n--- %s ---\n', opts.title);
    if ~isempty(opts.notes)
        fprintf('  %s\n', opts.notes);
    end

    for a = 1:numel(opts.actors)
        act = opts.actors(a);
        snr_rx = act.SNR_rx_dB;
        if is_fspl_link(act)
            [PL_lin, PL_dB] = compute_fspl(act.dist_m, act.fc_Hz);
            snr_tx = snr_rx + PL_dB;
            fprintf(['  %-22s  SNR_tx = %6.2f dB, SNR_rx = %6.2f dB', ...
                '  (d = %.0f m, fc = %.2f GHz, PL = %.2f dB)\n'], ...
                act.name, snr_tx, snr_rx, act.dist_m, act.fc_Hz/1e9, PL_dB);
        else
            fprintf('  %-22s  SNR_tx = %6.2f dB, SNR_rx = %6.2f dB', ...
                act.name, snr_rx, snr_rx);
            fprintf('  (normalized Rayleigh, no FSPL)\n');
        end
    end
end

% -------------------------------------------------------------------------
function opts = parse_inputs(varargin)
    opts.title     = 'Scenario';
    opts.SNR_rx_dB = 20;
    opts.dist_m    = [];
    opts.fc_Hz     = [];
    opts.actors    = struct('name', {}, 'SNR_rx_dB', {}, 'dist_m', {}, 'fc_Hz', {});
    opts.notes     = '';

    if mod(nargin, 2) ~= 0
        error('print_scenario_snr:args', 'Use name-value pairs only.');
    end

    for k = 1:2:nargin
        key = lower(varargin{k});
        val = varargin{k+1};
        switch key
            case 'title'
                opts.title = char(val);
            case {'snr_rx_db', 'snr_db'}
                opts.SNR_rx_dB = val;
            case 'dist_m'
                opts.dist_m = val;
            case 'fc_hz'
                opts.fc_Hz = val;
            case 'actors'
                opts.actors = normalize_actors(val, opts.SNR_rx_dB, opts.dist_m, opts.fc_Hz);
            case 'notes'
                opts.notes = char(val);
            otherwise
                error('print_scenario_snr:arg', 'Unknown option "%s".', varargin{k});
        end
    end

    if isempty(opts.actors)
        opts.actors = normalize_actors({'Bob', 'Eve'}, opts.SNR_rx_dB, opts.dist_m, opts.fc_Hz);
    end
end

function actors = normalize_actors(val, SNR_rx_dB, dist_m, fc_Hz)
    if isstruct(val)
        actors = val(:);
        for i = 1:numel(actors)
            if ~isfield(actors(i), 'name') || isempty(actors(i).name)
                error('print_scenario_snr:actor', 'Each actor needs a .name field.');
            end
            if ~isfield(actors(i), 'SNR_rx_dB') || isempty(actors(i).SNR_rx_dB)
                actors(i).SNR_rx_dB = SNR_rx_dB;
            end
            if ~isfield(actors(i), 'dist_m')
                actors(i).dist_m = dist_m;
            end
            if ~isfield(actors(i), 'fc_Hz')
                actors(i).fc_Hz = fc_Hz;
            end
        end
        return;
    end

    if iscell(val)
        if numel(val) == 1 && iscell(val{1})
            names = val{1}(:);
        else
            names = val(:);
        end
    elseif ischar(val) || isstring(val)
        names = cellstr(val);
    else
        error('print_scenario_snr:actors', ...
            'actors must be a cell array of names, a char vector, or a struct array.');
    end

    names = cellfun(@char, names, 'UniformOutput', false);
    actors = repmat(struct('name', '', 'SNR_rx_dB', SNR_rx_dB, ...
        'dist_m', dist_m, 'fc_Hz', fc_Hz), numel(names), 1);
    for i = 1:numel(names)
        actors(i).name = names{i};
    end
end

function tf = is_fspl_link(act)
    tf = ~isempty(act.dist_m) && ~isempty(act.fc_Hz) ...
        && isfinite(act.dist_m) && isfinite(act.fc_Hz) && act.dist_m > 0 && act.fc_Hz > 0;
end
