clear;
close all;
clc;

addpath('./functions/')

% Cheking deirectory

% dir_save_dl = './results/scheduling/downlink/';
% dir_save_ul = './results/scheduling/uplink/';

% root_save_dl = [dir_save_dow 'throughput_outdoors_pedestrian_mf_'];
% root_save_ul = [dir_save_upl 'throughput_outdoors_pedestrian_mf_'];

% root_save_dl = [dir_save_dl 'spectral_efficiency_mf_'];
% root_save_ul = [dir_save_ul 'spectral_efficiency_mf_'];
% 
% if ~exist(dir_save_dl,'dir')
%     mkdir(dir_save_dl);
% end
% 
% if ~exist(dir_save_ul,'dir')
%     mkdir(dir_save_ul);
% end

% Checking variables

MC = 100;                                                                  % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
M = 64;                                                                    % Number of antennas at the base station
K = 80;                                                                    % Number of users at the cell
channel_type = 'ur-los';

commcell.nAntennas       = M;                                              % Number of Antennas
commcell.nUsers          = K;                                              % Number of Users
commcell.radius          = 2000;                                            % Cell's raidus (circumradius) in meters
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

[~,snr_db] = linkBudgetCalculation(linkprop);                              % SNR in dB
                
snr = 10.^(snr_db/10);

% Initialization

se_ep     = zeros(K,MC);
gamma_ep  = zeros(K,MC);
eta_max   = zeros(K,MC);
gamma_max = zeros(MC,1);

for mc = 1:MC
    mc
    
    [G,beta] = massiveMIMOChannel(commcell,channel_type);
         
    [~,W] = decoderMatrix(G,'mf');
    
    [~,se_ep(:,mc),gamma_ep(:,mc)] = throughput(G,W,1/K,'downlink',snr*beta,settings);
    
    [gamma_max(mc),eta_max(:,mc)] = maxMinFairness(G,beta,snr);
end

aux   = repmat(10*log10(gamma_max),1,K)';
aux_2 = repmat(0.5*log(1 + gamma_max),1,K)';

figure;

ecdf(10*log10(min(gamma_ep)));
hold on;
ecdf(10*log10(gamma_max));

% save([root_save_dl strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
%       num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_d_eff) '_dB_MC_' ...
%       num2str(MC) '.mat'],'se_d','psi','se_d_sel','psi_sel');