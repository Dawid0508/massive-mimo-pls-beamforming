function p = default_params()
% DEFAULT_PARAMS  Returns the canonical parameter set used across all
% simulations in this project. Keeping a single source of truth avoids
% silent inconsistencies (e.g. different SNR conventions) between scripts.
%
%   p = default_params() returns a struct with the fields below. Override
%   any field in the calling script after retrieving the struct.
%
%   Requires Phased Array System Toolbox (physconst); see assert_requirements.

    assert_requirements('Phased_Array_System_Toolbox');

    p.c          = physconst('LightSpeed');

    % --- Carriers and arrays (matches Ghz6_band_vs_mmWave_band.m) ---
    p.fc_sub6    = 6e9;          % Sub-6 GHz Massive MIMO carrier [Hz]
    p.fc_mmwave  = 28e9;         % mmWave Ultra-Massive MIMO carrier [Hz]
    p.Nt_sub6    = 32;           % Antennas at sub-6 GHz BS
    p.Nt_mmwave  = 512;          % Antennas at mmWave BS

    % --- Power / noise convention --------------------------------------
    %   Noise variance is fixed at 1 (linear). SNR_dB is always the
    %   *received* SNR at Bob after path loss. Transmit power:
    %       P_tx = rx_snr_power('tx_for_rx', SNR_dB, PL_lin)
    %   Abstract scenarios (no FSPL) use unit-variance Rayleigh channels;
    %   there received SNR equals the configured linear power in the rate.
    p.noise_var  = 1;
    p.SNR_rx_dB  = 30;           % default received SNR at Bob [dB]
    p.SNR_dB     = p.SNR_rx_dB;  % alias for backward compatibility
    p.link_dist_m = 30;          % default Alice–user distance [m]
    p.eve_attn_dB = 15;          % Eve channel attenuation (abstract Rayleigh)

    % --- 3GPP TR 38.901 CDL tags (used by channel_3gpp_ula) ------------
    p.cdl_sub6    = 'sub6';      % CDL-A-like at 6 GHz
    p.cdl_mmwave  = 'mmwave';    % CDL-D-like at 28 GHz

    % --- Monte-Carlo --------------------------------------------------
    p.numIter    = 200;          % default Monte-Carlo iterations
    p.rng_seed   = 2026;         % fixed seed for reproducibility

    % --- Output --------------------------------------------------------
    p.results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
end
