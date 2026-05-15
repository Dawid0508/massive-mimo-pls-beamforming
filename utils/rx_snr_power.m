function varargout = rx_snr_power(varargin)
% RX_SNR_POWER  Received-SNR power helpers (project-wide convention).
%
%   P_rx = rx_snr_power('linear', SNR_rx_dB)
%   P_tx = rx_snr_power('tx_for_rx', SNR_rx_dB, PL_lin)
%   Legacy 'print' mode -> use print_scenario_snr instead.
%
%   SNR_rx_dB always denotes the desired SNR at the legitimate receiver
%   *after* path loss. Transmit power is P_tx = P_rx * PL_lin.

    if nargin < 1
        error('rx_snr_power:mode', 'Specify mode: ''linear'', ''tx_for_rx'', or ''print''.');
    end

    mode = lower(varargin{1});

    switch mode
        case 'linear'
            SNR_rx_dB = varargin{2};
            varargout{1} = 10.^(SNR_rx_dB / 10);

        case 'tx_for_rx'
            SNR_rx_dB = varargin{2};
            PL_lin    = varargin{3};
            P_rx      = 10.^(SNR_rx_dB / 10);
            varargout{1} = P_rx * PL_lin;

        case 'print'
            SNR_rx_dB = varargin{2};
            dist_m    = varargin{3};
            fc_Hz     = varargin{4};
            if nargin >= 5 && ~isempty(varargin{5})
                title = varargin{5};
            else
                title = 'Link budget';
            end
            print_scenario_snr('title', title, 'SNR_rx_dB', SNR_rx_dB, ...
                'dist_m', dist_m, 'fc_Hz', fc_Hz, 'actors', {'Bob', 'Eve'});

        otherwise
            error('rx_snr_power:mode', 'Unknown mode "%s".', mode);
    end
end
