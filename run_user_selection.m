MC = 10000;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
L = 13;                                                                    % Number of selected users

for M = [64 256]                                                           % Number of antennas at the base station
    for K = [18 36 72]                                                     % Number of users at the cell
        for snr_db = [-20 -15 -10 -5 0 5 10]                               % SNR in dB
            channel_type = 'ur-los';
    
            run user_selection.m 
        end
    end
end