MC = 10000;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)

snr_db = 10;                                                               % SNR in dB
channel_model = {'ur-los','rayleigh'};                                     % Channel models

for chn_idx = [1 2]
    channel_type = channel_model{chn_idx};
    for M = [64 256]                                                       % Number of antennas at the base station
        for K = [18 36 72]                                                 % Number of users at the cell
            % for snr_db = [-20 -15 -10 -5 0 5 10]                           % SNR in dB
                run automatic_selection.m 
            % end
        end
    end
end