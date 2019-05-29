clear;
close all;
clc;

addpath('./functions/')

root_downlink = './results/auto_scheduling/downlink/rate_downlink_mf_';
root_uplink   = './results/auto_scheduling/uplink/rate_uplink_mf_';

MC    = 10000;                                                             % Size of the Monte Carlo ensemble (Varies the channel realizarions)
N_ALG = 2;

M = 64;                                                                    % Number of antennas at the base station
K = 18;                                                                    % Number of users at the cell

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

% linkprop.bsPower         = 10;                                             % in Watts
% linkprop.userPower       = 0.2;                                            % in Watts
% linkprop.AntennaGainBS   = 0;                                              % in dBi
% linkprop.AntennaGainUser = 0;                                              % in dBi
% linkprop.noiseFigureBS   = 9;                                              % in dB
% linkprop.noiseFigureUser = 9 ;                                             % in dB
% linkprop.bandwidth       = 20e6;                                           % in Hz

% [snr_u_db,snr_d_db] = linkBudgetCalculation(linkprop);                     % SNR in dB

% beta_db = -135;

% snr_u_eff = round(snr_u_db + beta_db);
% snr_d_eff = round(snr_d_db + beta_db);

snr_u_eff = 20;
snr_d_eff = 20;

snr_u = 10.^((snr_u_eff)/10);                                              % Uplink SNR
snr_d = 10.^((snr_d_eff)/10);                                              % Downlink SNR

% Normalized threshold

tau_norm_min = 0;
tau_norm_max = 1;

tau_norm_step = 0.01;

tau_norm = tau_norm_min:tau_norm_step:tau_norm_max;

N_TAU = length(tau_norm);

tau_icibs_max = 0.25;

tau(:,1) = tau_norm';                                                      % Threshold - CBS
tau(:,2) = tau_icibs_max*tau_norm';                                        % Threshold - ICIBS

% Initialization

L = zeros(MC,N_TAU,N_ALG);

user_sel = cell(MC,N_TAU,N_ALG);
rate_u   = cell(MC,N_TAU,N_ALG);
rate_d   = cell(MC,N_TAU,N_ALG);
psi      = cell(MC,N_TAU,N_ALG);

channel_type   = 'ur-los';
algorithm_type = {'correlation-based selection','ici-based selection'};

% Correlation-based Selection

% ICI-based Selection

for mc = 1:MC
    mc
    
    [G,~] = massiveMIMOChannel(commcell,channel_type);
    
    for alg_idx = 1:N_ALG
        for tau_idx = 1:N_TAU
            [user_sel{mc,tau_idx,alg_idx},H_sel] = userScheduling(G,algorithm_type{alg_idx},'automatic',[],tau(tau_idx,alg_idx));
            
            L(mc,tau_idx,alg_idx) = size(H_sel,2);
            
            [Q,W] = decoderMatrix(H_sel,'mf');
                        
            pow_upl = ones(L(mc,tau_idx,alg_idx),1);
            pow_dow = ones(L(mc,tau_idx,alg_idx),1)/L(mc,tau_idx,alg_idx);
            
            rate_u{mc,tau_idx,alg_idx} = rateCalculation(H_sel,Q,pow_upl,snr_u,'uplink');
            rate_d{mc,tau_idx,alg_idx} = rateCalculation(H_sel,W,pow_dow,snr_d,'downlink');
            
            psi{mc,tau_idx,alg_idx} = ici(H_sel);
        end
    end
end

save([root_downlink strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' num2str(K) '_SNR_' ...
    num2str(snr_u_eff) '_dB_MC_' num2str(MC) '.mat'],'user_sel','rate_d','psi','L');

save([root_uplink strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' num2str(K) '_SNR_' ...
    num2str(snr_u_eff) '_dB_MC_' num2str(MC) '.mat'],'user_sel','rate_u','psi','L');