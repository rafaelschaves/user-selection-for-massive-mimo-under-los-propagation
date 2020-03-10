MC = 1000;                                                                 % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)

channel_type = 'ur-los';

M = 50;                                                                    % Number of antennas at the base station

for K = [10 30 50 70]                                                      % Number of users at the cell
    if K == 70
        for L = 1:M-1                                                      % Number of selected users
            run user_selection.m
        end
    else
        for L = 1:M-1                                                      % Number of selected users
            run user_selection.m
        end
    end
end