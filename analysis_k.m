addpath('./functions/')

% Cheking deirectory

dir_save  = './results/users/';
root_save = [dir_save 'throughput_outdoors_pedestrian_mf_'];

if ~exist(dir_save,'dir')
    mkdir(dir_save);
end

% Checking variables

if ~exist('MC','var')
    MC = 10000;                                                            % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
end

if ~exist('M','var')
    M = 64;                                                                % Number of antennas at the base station
end

if ~exist('K','var')
    K = 1;                                                                 % Number of users at the cell
end

if ~exist('snr_db','var')
    snr_db = 10;                                                           % SNR in dB
end

if ~exist('channel_type','var')
    channel_type = 'ur-los';
end

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

settings.coherenceTime           = 15000;                                  % Coherence time in samples
settings.PilotTime               = K;                                      % Pilot time in samples
settings.uplinkDownlinkTimeRatio = 0.5;                                    % Ratio between the uplink and downlink payload time
settings.bandwidth               = 20e6;                                   % Sytem bandwidth in Hz
settings.cellArea                = 1;                                      % Cell area in km^2
settings.snr                     = 10.^((snr_db)/10);                      % SNR

% Initialization

thrput_u = zeros(K,MC);
thrput_d = zeros(K,MC);

pow_upl = ones(K,1);
pow_dow = ones(K,1)/K;
              
for mc = 1:MC
    mc
    
    [G,~] = massiveMIMOChannel(commcell,channel_type);
        
    [Q,W] = decoderMatrix(G,'mf');
    
    settings.linkType = 'uplink';
    
    thrput_u(:,mc) = throughput(G,Q,pow_upl,settings);
   
    settings.linkType = 'downlink';
    
    thrput_d(:,mc) = throughput(G,W,pow_dow,settings);
end

save([root_save strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_SNR_' num2str(snr_db) '_dB_MC_' num2str(MC) '.mat'], ...
      'thrput_d','thrput_u');