addpath('./functions/');

% Cheking deirectory

dir_save  = './results/power_allocation/downlink/';
root_save = [dir_save 'results_upper_bound_'];

if ~exist(dir_save,'dir')
    mkdir(dir_save);
end

% Checking variables

if ~exist('MC','var')
    MC = 1000;                                                             % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
end

if ~exist('M','var')
    M = 64;                                                                % Number of antennas at the base station
end

if ~exist('K','var')
    K = 8;                                                                 % Number of users at the cell
end

if ~exist('radius','var')
    radius = 500;
end

if ~exist('bs_power','var')
    bs_power = 1;
end

if ~exist('channel_type','var')
    channel_type = 'ur-los';
end

commcell.nAntennas       = M;                                              % Number of Antennas
commcell.nUsers          = K;                                              % Number of Users
commcell.radius          = radius;                                         % Cell's raidus (circumradius) in meters
commcell.bsHeight        = 32;                                             % Height of base station in meters
commcell.userHeight      = [1 2];                                          % Height of user terminals in meters ([min max])
commcell.frequency       = 1.9e9;                                          % Carrier frequency in Hz
commcell.meanShadowFad   = 0;                                              % Shadow fading mean in dB
commcell.stdDevShadowFad = 8;                                              % Shadow fading standard deviation in dB
commcell.city            = 'large';                                        % Type of city

linkprop.bsPower         = bs_power;                                       % in Watts
linkprop.userPower       = 0.2;                                            % in Watts
linkprop.AntennaGainBS   = 0;                                              % in dBi
linkprop.AntennaGainUser = 0;                                              % in dBi
linkprop.noiseFigureBS   = 9;                                              % in dB
linkprop.noiseFigureUser = 9 ;                                             % in dB
linkprop.bandwidth       = 20e6;                                           % in Hz

[~,snr_db] = linkBudgetCalculation(linkprop);                              % SNR in dB
snr        = 10.^(snr_db/10);

gamma      = zeros(MC,1);
gamma_u    = zeros(MC,1);
max_lambda = zeros(MC,1);

for mc = 1:MC
    mc
    
    [G,beta]   = massiveMIMOChannel(commcell,channel_type);
    
    [gamma(mc),~,gamma_u(mc),~,R] = maxMinFairness(G,beta,snr,'algorithm 1');
    
    max_lambda(mc) = max(sum(R,2));
end

result = sum(gamma > 1./max_lambda);

save([root_save strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_cell_radius_' num2str(radius) '_m_BS_power_' ...
      num2str(bs_power) '_W_MC_' num2str(MC) '.mat'],'result');