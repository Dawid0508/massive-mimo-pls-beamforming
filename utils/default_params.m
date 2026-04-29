function p = default_params()
% DEFAULT_PARAMS  Returns the canonical parameter set used across all
% simulations in this project. Keeping a single source of truth avoids
% silent inconsistencies (e.g. different SNR conventions) between scripts.
%
%   p = default_params() returns a struct with the fields below. Override
%   any field in the calling script after retrieving the struct.

    p.c          = physconst('LightSpeed');

    % --- Carriers and arrays (matches Ghz6_band_vs_mmWave_band.m) ---
    p.fc_sub6    = 6e9;          % Sub-6 GHz Massive MIMO carrier [Hz]
    p.fc_mmwave  = 28e9;         % mmWave Ultra-Massive MIMO carrier [Hz]
    p.Nt_sub6    = 32;           % Antennas at sub-6 GHz BS
    p.Nt_mmwave  = 512;          % Antennas at mmWave BS

    % --- Power / noise convention --------------------------------------
    %   Noise variance is fixed at 1 (linear). All "SNR" values are
    %   transmit-side SNR in dB; effective received SNR equals
    %   10^(SNR_dB/10) / FSPL_lin.
    p.noise_var  = 1;
    p.SNR_dB     = 30;           % default operating point [dB]

    % --- Monte-Carlo --------------------------------------------------
    p.numIter    = 200;          % default Monte-Carlo iterations
    p.rng_seed   = 2026;         % fixed seed for reproducibility

    % --- Output --------------------------------------------------------
    p.results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
end
