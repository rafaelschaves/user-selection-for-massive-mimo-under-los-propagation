clear;
close all;
clc;

addpath('./functions/')

% Cheking deirectory

dir_save  = './results/scheduling/downlink/partial_csi/';
root_save = [dir_save 'spectral_efficiency_mf_'];

if ~exist(dir_save,'dir')
    mkdir(dir_save);
end

% Checking variables

if ~exist('MC_1','var')
    MC_1 = 10;                                                           % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
end

if ~exist('MC_2','var')
    MC_2 = 100;                                                            % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
end

if ~exist('M','var')
    M = 64;                                                                % Number of antennas at the base station
end

if ~exist('K','var')
    K = 72;                                                                % Number of users at the cell
end

if ~exist('L','var')
    L = 8;                                                                 % Number of selected users
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

[~,snr_db] = linkBudgetCalculation(linkprop);                              % SNR in dB
                
beta_db = -148 - 37.6*log10(commcell.radius/1000);
% beta    = 10^(beta_db/10);
beta = 1;

snr_eff = round(snr_db + beta_db);
snr     = 10.^(snr_eff/10);

xi = 0:0.1:1;

N_XI  = length(xi);
N_ALG = 4;

% Initialization

psi     = zeros(K,MC_2,N_XI,MC_1);
psi_sel = zeros(L,MC_2,N_XI,MC_1,N_ALG);

se_ep  = zeros(K,MC_2,N_XI,MC_1);
se_max = zeros(K,MC_2,N_XI,MC_1);

se_ep_sel  = zeros(L,MC_2,N_XI,MC_1,N_ALG);
se_max_sel = zeros(L,MC_2,N_XI,MC_1,N_ALG);

gamma_max_0     = zeros(MC_2,N_XI,MC_1);
gamma_max_sel_0 = zeros(MC_2,N_XI,MC_1,N_ALG);


algorithm_type = {'random selection', ...
                  'semi-orthogonal selection', ...
                  'correlation-based selection', ...
                  'ici-based selection'};

eta_max     = zeros(K,MC_2,N_XI,MC_1);
eta_max_sel = zeros(L,MC_2,N_XI,MC_1,N_ALG);

for mc_1 = 1:MC_1
    mc_1
    
    [G,~] = massiveMIMOChannel(commcell,channel_type);
        
    for n_xi = 1:N_XI
        for mc_2 = 1:MC_2
            E = (randn(M,K) + 1i*randn(M,K))/sqrt(2);
            
            G_hat = xi(n_xi)*G + sqrt(1 - xi(n_xi)^2)*E;
            
            psi(:,mc_2,n_xi,mc_1) = ici(G_hat);
            
            [gamma_max_0(mc_2,n_xi,mc_1),eta_max(:,mc_2,n_xi,mc_1)] = maxMinFairness(G_hat,beta,snr);
    
            W = precoderMatrix(G_hat,'mf');
                
            [~,se_ep(:,mc_2,n_xi,mc_1)]  = throughput(G,W,1/K                      ,'downlink',snr,settings);
            [~,se_max(:,mc_2,n_xi,mc_1)] = throughput(G,W,eta_max(:,mc_2,n_xi,mc_1),'downlink',snr,settings);
            
            for alg_idx = 1:N_ALG
                [user_sel,G_sel_hat] = userSelector(G_hat,algorithm_type{alg_idx},'fixed',L,[]);
                
                G_sel = G(:,user_sel);
                
                psi_sel(:,mc_2,n_xi,mc_1,alg_idx) = ici(G_sel_hat);
                
                [gamma_max_sel_0(mc_2,n_xi,mc_1,alg_idx),eta_max_sel(:,mc_2,n_xi,mc_1,alg_idx)] = maxMinFairness(G_sel_hat,beta,snr);
                
                W = precoderMatrix(G_sel_hat,'mf');
                
                [~,se_ep_sel(:,mc_2,n_xi,mc_1,alg_idx)]  = throughput(G_sel,W,1/L                                  ,'downlink',snr,settings);
                [~,se_max_sel(:,mc_2,n_xi,mc_1,alg_idx)] = throughput(G_sel,W,eta_max_sel(:,mc_2,n_xi,mc_1,alg_idx),'downlink',snr,settings);
            end
        end
    end
end

save([root_save strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_L_' num2str(L) '_radius_' num2str(commcell.radius) ...
      '_m_BS_power_' num2str(linkprop.bsPower) '_W_MC_' num2str(MC_1) ...
      '.mat'],'se_ep','se_max','gamma_max_0','eta_max','psi','se_ep_sel', ...
      'se_max_sel','gamma_max_sel_0','eta_max_sel','psi_sel');
