MC = 10000;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
L = 13;                                                                    % Number of selected users

channel_type = 'rayleigh';

for M = [64]                                                           % Number of antennas at the base station
    for K = [36]                                                     % Number of users at the cell
        for snr_db = [10]                               % SNR in dB
            run user_selection.m 
        end
    end
end
