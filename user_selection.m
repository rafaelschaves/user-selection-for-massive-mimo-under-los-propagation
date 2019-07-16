addpath('./functions/')

% Cheking deirectory

dir_save_dl = './results/scheduling/downlink/';
dir_save_ul = './results/scheduling/uplink/';

% root_save_dl = [dir_save_dow 'throughput_outdoors_pedestrian_mf_'];
% root_save_ul = [dir_save_upl 'throughput_outdoors_pedestrian_mf_'];

root_save_dl = [dir_save_dl 'spectral_efficiency_mf_'];
root_save_ul = [dir_save_ul 'spectral_efficiency_mf_'];

if ~exist(dir_save_dl,'dir')
    mkdir(dir_save_dl);
end

if ~exist(dir_save_ul,'dir')
    mkdir(dir_save_ul);
end

% Checking variables

% if ~exist('MC','var')
%     MC = 10000;                                                            % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
% end
% 
% if ~exist('M','var')
%     M = 64;                                                                % Number of antennas at the base station
% end
% 
% if ~exist('K','var')
%     K = 72;                                                                % Number of users at the cell
% end
% 
% if ~exist('L','var')
%     L = 13;                                                                % Number of selected users
% end
% 
% if ~exist('snr_db','var')
%     snr_db = 10;                                                           % SNR in dB
% end
% 
% if ~exist('channel_type','var')
%     channel_type = 'ur-los';
% end

N_ALG = 4;

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

settings.coherenceTime           = 1;                                      % Coherence time in samples
settings.PilotTime               = 0;                                      % Pilot time in samples
settings.uplinkDownlinkTimeRatio = 0.5;                                    % Ratio between the uplink and downlink payload time
settings.bandwidth               = 20e6;                                   % Sytem bandwidth in Hz
settings.cellArea                = 1;                                      % Cell area in km^2

linkprop.bsPower         = 10;                                             % in Watts
linkprop.userPower       = 0.2;                                            % in Watts
linkprop.AntennaGainBS   = 0;                                              % in dBi
linkprop.AntennaGainUser = 0;                                              % in dBi
linkprop.noiseFigureBS   = 9;                                              % in dB
linkprop.noiseFigureUser = 9 ;                                             % in dB
linkprop.bandwidth       = 20e6;                                           % in Hz

[snr_u_db,snr_d_db] = linkBudgetCalculation(linkprop);                     % SNR in dB
                
beta_db = -148 - 37.6*log10(commcell.radius/1000);

snr_u_eff = round(snr_u_db + beta_db);
snr_d_eff = round(snr_d_db + beta_db);

snr_u = 10.^(snr_u_eff/10);
snr_d = 10.^(snr_d_eff/10);

% Initialization

se_u     = zeros(K,MC);
se_d     = zeros(K,MC);
se_u_sel = zeros(L,MC,N_ALG);
se_d_sel = zeros(L,MC,N_ALG);

psi     = zeros(K,MC);
psi_sel = zeros(L,MC,N_ALG);

algorithm_type = {'random selection', ...
                  'semi-orthogonal selection', ...
                  'correlation-based selection', ...
                  'ici-based selection'};

pow_ul = ones(K,1);
pow_dl = ones(K,1)/K;

pow_ul_sel = ones(L,1);
pow_dl_sel = ones(L,1)/L;
              
for mc = 1:MC
    mc
    
    [G,~] = massiveMIMOChannel(commcell,channel_type);
    
    psi(:,mc)   = ici(G);
    
    % No Selection
    
    [Q,W] = decoderMatrix(G,'mf');
    
    [~,se_u(:,mc)] = throughput(G,Q,pow_ul,'uplink',snr_u,settings); 
    [~,se_d(:,mc)] = throughput(G,W,pow_dl,'downlink',snr_d,settings);
          
    for alg_idx = 1:N_ALG
        [~,H_sel] = userSelector(G,algorithm_type{alg_idx},'fixed',L,[]);
        
        psi_sel(:,mc,alg_idx) = ici(H_sel);

        [Q,W] = decoderMatrix(H_sel,'mf');
                                                      
        [~,se_u_sel(:,mc,alg_idx)] = throughput(H_sel,Q,pow_ul_sel,'uplink',snr_u,settings);               
        [~,se_d_sel(:,mc,alg_idx)] = throughput(H_sel,W,pow_dl_sel,'downlink',snr_d,settings);    
    end
end

save([root_save_dl strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_d_eff) '_dB_MC_' ...
      num2str(MC) '.mat'],'se_d','psi','se_d_sel','psi_sel');

save([root_save_ul strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_u_eff) '_dB_MC_' ...
      num2str(MC) '.mat'],'se_u','psi','se_u_sel','psi_sel');