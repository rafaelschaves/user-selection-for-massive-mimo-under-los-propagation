addpath('./functions/')

% Cheking deirectory

dir_save_dow = './results/scheduling/downlink/';
dir_save_upl = './results/scheduling/uplink/';

root_save_dow = [dir_save_dow 'rate_mf_'];
root_save_upl = [dir_save_upl 'rate_mf_'];

if ~exist(dir_save_dow,'dir')
    mkdir(dir_save_dow);
end

if ~exist(dir_save_upl,'dir')
    mkdir(dir_save_upl);
end

% Checking variables

if ~exist('MC','var')
    MC = 10000;                                                            % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
end

if ~exist('M','var')
    M = 64;                                                                % Number of antennas at the base station
end

if ~exist('K','var')
    K = 72;                                                                % Number of users at the cell
end

if ~exist('L','var')
    L = 13;                                                                % Number of selected users
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

N_ALG = 4;

settings.coherenceTime           = 15000;                                  % Coherence time in samples
settings.PilotTime               = K;                                      % Pilot time in samples
settings.uplinkDownlinkTimeRatio = 0.5;                                    % Ratio between the uplink and downlink payload time
settings.bandwidth               = 20e6;                                   % Sytem bandwidth in Hz
settings.cellArea                = 1;                                      % Cell area in km^2
settings.snr                     = 10.^((snr_db)/10);                      % SNR

% Initialization

thrput_u     = zeros(K,MC);
thrput_d     = zeros(K,MC);
thrput_u_sel = zeros(L,MC,N_ALG);
thrput_d_sel = zeros(L,MC,N_ALG);

psi     = zeros(K,MC);
psi_sel = zeros(L,MC,N_ALG);

algorithm_type = {'random selection', ...
                  'semi-orthogonal selection', ...
                  'correlation-based selection', ...
                  'ici-based selection'};

pow_upl = ones(K,1);
pow_dow = ones(K,1)/K;

pow_upl_sel = ones(L,1);
pow_dow_sel = ones(L,1)/L;
              
for mc = 1:MC
    mc
    
    [G,~] = massiveMIMOChannel(commcell,channel_type);
    
    psi(:,mc)   = ici(G);
    
    % No Selection
    
    [Q,W] = decoderMatrix(G,'mf');
    
    settings.linkType = 'uplink';
    
    thrput_u(:,mc) = throughput(G,Q,pow_upl,settings);
   
    settings.linkType = 'downlink';
    
    thrput_d(:,mc) = throughput(G,W,pow_dow,settings);
          
    for alg_idx = 1:N_ALG
        [~,H_sel] = userScheduling(G,algorithm_type{alg_idx},'fixed',L,[]);
        
        psi_sel(:,mc,alg_idx) = ici(H_sel);

        [Q,W] = decoderMatrix(H_sel,'mf');
                                                      
        settings.linkType = 'uplink';                                      

        thrput_u_sel(:,mc,alg_idx) = throughput(H_sel,Q,pow_upl_sel,settings);
        
        settings.linkType = 'downlink';                                      
               
        thrput_d_sel(:,mc,alg_idx) = throughput(H_sel,W,pow_dow_sel,settings);    
    end
end

save([root_save_dow 'M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) ...
      '_SNR_' num2str(snr_db) '_dB_MC_' num2str(MC) '.mat'], ...
      'thrput_d','psi','thrput_d_sel','psi_sel');

save([root_save_upl 'M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) ...
      '_SNR_' num2str(snr_db) '_dB_MC_' num2str(MC) '.mat'], ...
      'thrput_u','psi','thrput_u_sel','psi_sel');