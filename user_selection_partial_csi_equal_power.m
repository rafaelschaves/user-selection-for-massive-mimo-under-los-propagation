clear;
close all;
clc;

addpath('./functions/')

% Cheking deirectory

dir_save_dl  = './results/scheduling/downlink/partial_csi';
root_save_dl = [dir_save_dl 'spectral_efficiency_mf_'];

if ~exist(dir_save_dl,'dir')
    mkdir(dir_save_dl);
end

% Checking variables

if ~exist('MC_OUTTER','var')
    MC_OUTTER = 1000;                                                      % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
end

if ~exist('MC_INNER','var')
    MC_INNER = 100;                                                        % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
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

if ~exist('channel_type','var')
    channel_type = 'ur-los';
end

commcell.nAntennas       = M;                                              % Number of Antennas
commcell.nUsers          = K;                                              % Number of Users
commcell.radius          = 100;                                            % Cell's raidus (circumradius) in meters
commcell.bsHeight        = 32;                                             % Height of base station in meters
commcell.userHeight      = [1 2];                                          % Height of user terminals in meters ([min max])
commcell.frequency       = 1.9e9;                                          % Carrier frequency in Hz
commcell.meanShadowFad   = 0;                                              % Shadow fading mean in dB
commcell.stdDevShadowFad = 8;                                              % Shadow fading standard deviation in dB
commcell.city            = 'large';                                        % Type of city

settings.coherenceTime           = 1;                                      % Coherence time in samples
settings.PilotTime               = 0;                                      % Pilot time in samples
settings.uplinkDownlinkTimeRatio = 0.5;                                    % Ratio between the uplink and downlink payload time
settings.bandwidth               = 20e6;                                   % Sytem bandwidth in Hz
settings.cellArea                = 1;                                      % Cell area in km^2

linkprop.bsPower         = 1;                                              % in Watts
linkprop.userPower       = 0.2;                                            % in Watts
linkprop.AntennaGainBS   = 0;                                              % in dBi
linkprop.AntennaGainUser = 0;                                              % in dBi
linkprop.noiseFigureBS   = 9;                                              % in dB
linkprop.noiseFigureUser = 9 ;                                             % in dB
linkprop.bandwidth       = 20e6;                                           % in Hz

[~,snr_dl_db] = linkBudgetCalculation(linkprop);                           % SNR in dB
                
beta_db = -148 - 37.6*log10(commcell.radius/1000);
beta    = 10^(beta_db/10);

snr_dl_eff = round(snr_dl_db + beta_db);
snr_dl     = 10.^(snr_dl_eff/10);

xi = 0:0.1:1;

N_XI  = length(xi);
N_ALG = 4;

% Initialization

psi     = zeros(K,MC_OUTTER);
psi_sel = zeros(L,MC_INNER,N_XI,MC_OUTTER,N_ALG);

se_dl     = zeros(K,MC_OUTTER);
se_dl_sel = zeros(L,MC_INNER,N_XI,MC_OUTTER,N_ALG);

algorithm_type = {'random selection', ...
                  'semi-orthogonal selection', ...
                  'correlation-based selection', ...
                  'ici-based selection'};

for mc_outter = 1:MC_OUTTER
    mc_outter
    
    [G,~] = massiveMIMOChannel(commcell,channel_type);
    
    [~,se_dl(:,mc_outter)] = throughput(G,precoderMatrix(G,'mf'),1/K,'downlink',snr_dl,settings);
    
    for n_xi = 1:N_XI
        for mc_inner = 1:MC_INNER
            E = (randn(M,K) + 1i*randn(M,K))/sqrt(2);
            
            G_hat = xi(n_xi)*G + sqrt(1 - xi(n_xi)^2)*E;
            
            for alg_idx = 1:N_ALG
                [~,G_sel] = userSelector(G_hat,algorithm_type{alg_idx},'fixed',L,[]);
                
                [~,se_dl_sel(:,mc_inner,n_xi,mc_outter,alg_idx)] = throughput(G_sel,precoderMatrix(G_sel,'mf'),1/L,'downlink',snr_dl,settings);
            end
        end
    end
end

save([root_save_dl strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_dl_eff) '_dB_radius' ...
      num2str(commcell.radius) 'BS_power' num2str(linkprop.bsPower) '_MC_' ...
      num2str(MC_OUTTER) '.mat'],'se_dl','psi','se_dl_sel','psi_sel');