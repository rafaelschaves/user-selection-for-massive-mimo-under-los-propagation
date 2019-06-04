MC = 10000;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
L = 31;                                                                    % Number of selected users

channel_type = 'rayleigh';

for M = [256]                                                           % Number of antennas at the base station
    for K = [72]                                                     % Number of users at the cell
        for snr_db = [5 10]                               % SNR in dB
            run user_selection.m 
        end
    end
end
