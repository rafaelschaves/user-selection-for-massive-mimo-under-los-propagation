MC = 20;                                                                   % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)

channel_type = 'ur-los';

M = 50;                                                                    % Number of antennas at the base station

for K = [10 30 50 70]                                                      % Number of users at the cell
    K
    if K > M
        for L = 1:M-1                                                      % Number of selected users
            L
            run user_selection.m
        end
    else
        for L = 1:K-1                                                      % Number of selected users
            L
            run user_selection.m
        end
    end
end