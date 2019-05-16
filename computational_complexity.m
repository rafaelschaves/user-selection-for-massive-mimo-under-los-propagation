clear;
close all;
clc;

addpath('./functions/')

root = './results/comp_time_';

MC    = 10000;                                                              % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
N_ALG = 4;                                                                 % Number of algorithms user to perform user scheduling
N_CHN = 3;                                                                 % Number of channel models tested

M = 64;                                                                   % Number of antennas at the base station
K = 18;                                                                    % Number of users at the cell
L = K - 1;                                                                 % Maximum number of selected users

commcell.nAntennas       = M;                                              % Number of Antennas
commcell.nUsers          = K;                                              % Number of Users
commcell.radius          = 200;                                            % Cell's raidus (circumradius) in meters
commcell.bsHeight        = 32;                                             % Height of base station in meters
commcell.userHeight      = [1 2];                                          % Height of user terminals in meters ([min max])
commcell.nPaths          = 30;                                             % Number of Multipaths
commcell.frequency       = 1.9e9;                                          % Carrier frequency in Hz
commcell.meanShadowFad   = 0;                                              % Shadow fading mean in dB
commcell.stdDevShadowFad = 8;                                              % Shadow fading standard deviation in dB
commcell.city            = 'large';                                        % Type of city

% Initialization

func = cell(N_ALG,1);

comp_time = zeros(MC,L,N_ALG,N_CHN);

channel_type = {'ur-los','sparse','rayleigh'};

for chn_idx = 1:N_CHN
    for l = 1:L
        l
        
        for mc = 1:MC
            [G,~] = massiveMIMOChannel(commcell,channel_type{chn_idx});
            
            func{1} = @() userScheduling(G,'random selection','fixed',l);
            func{2} = @() userScheduling(G,'semi-orthogonal selection','fixed',l);
            func{3} = @() userScheduling(G,'correlation-based selection','fixed',l);
            func{4} = @() userScheduling(G,'ici-based selection','fixed',l);
            
            comp_time(mc,l,1,chn_idx) = timeit(func{1},2);                 % Random Selection
            comp_time(mc,l,2,chn_idx) = timeit(func{2},2);                 % Semi-orthogonal Selection
            comp_time(mc,l,3,chn_idx) = timeit(func{3},2);                 % Correlation-based Selection
            comp_time(mc,l,4,chn_idx) = timeit(func{4},2);                 % ICI-based Selection
        end
    end
end

save([root 'M_' num2str(M) '_K_' num2str(K) '_MC_' num2str(MC) '.mat'],'comp_time');