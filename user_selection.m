clear;
close all;
clc;

addpath('./functions/')

% Cheking deirectory

dir_save_dow = './results/scheduling/downlink/';
dir_save_upl = './results/scheduling/uplink';

root_save_dow = [dir_save_dow 'rate_mf_'];
root_save_upl = [dir_save_pul 'rate_mf_'];

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
    K = 18;                                                                % Number of users at the cell
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

snr = 10.^((snr_db)/10);                                                   % SNR

% Initialization

rate_u     = zeros(K,MC);
rate_d     = zeros(K,MC);
rate_u_sel = zeros(L,MC,N_ALG);
rate_d_sel = zeros(L,MC,N_ALG);

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
    
    rate_u(:,mc) = rateCalculation(G,Q,pow_upl,snr,'uplink');
    rate_d(:,mc) = rateCalculation(G,W,pow_dow,snr,'downlink');
        
    for alg_idx = 1:N_ALG
        [~,H_sel] = userScheduling(G,algorithm_type{alg_idx},'fixed',L,[]);
        
        psi_sel(:,mc,alg_idx) = ici(H_sel);

        [Q,W] = decoderMatrix(H_sel,'mf');
                                                      
        rate_u_sel(:,mc,alg_idx) = rateCalculation(H_sel,Q,pow_upl_sel,snr,'uplink');
        rate_d_sel(:,mc,alg_idx) = rateCalculation(H_sel,W,pow_dow_sel,snr,'downlink');    
    end
end

save([root_save_dow strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_db) '_dB_MC_' ...
      num2str(MC) '.mat'],'rate_d','psi','rate_d_sel','psi_sel');

save([root_save_upl strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_db) '_dB_MC_' ...
      num2str(MC) '.mat'],'rate_u','psi','rate_u_sel','psi_sel');