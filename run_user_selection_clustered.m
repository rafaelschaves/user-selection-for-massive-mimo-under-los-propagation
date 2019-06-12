MC = 10000;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
L = 67;                                                                    % Number of selected users
snr_db = 10;                                                               % SNR in dB

for M = [64 256]                                                           % Number of antennas at the base station
    for K = [72]                                                           % Number of users at the cell
        for theta_mid = [0 pi/4 pi/2]
            for theta_step = [pi/180 pi/18 pi/9 pi/6 pi/2] 
                run user_selection_clustered.m
            end
        end
    end
end