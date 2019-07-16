MC = 10000;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
                                                                 
channel_type = 'ur-los';

for M = [64 128 256]                                                       % Number of antennas at the base station
    for r_k = 0.25:0.25:1.25
        K = M*r_k;                                                         % Number of users at the cell
        for r_l = 0.25:0.25:0.75
            L = K*r_l;                                                     % Number of selected users
            run user_selection.m 
        end
    end
end
