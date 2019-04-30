clear;
close all;
clc;

addpath('./functions/')

MC = 10000;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)

M = 500;                                                                   % Number of antennas at the base station
K = 30;                                                                     % Number of users at the cell
L = 25;

commcell.nAntennas       = M;                                              % Number of Antennas
commcell.nUsers          = K;                                              % Number of Users
commcell.radius          = 500;                                            % Cell's raidus (circumradius) in meters
commcell.bsHeight        = 30;                                             % Height of base station in meters
commcell.userHeight      = [1 2];                                          % Height of user terminals in meters ([min max])
commcell.nPaths          = 30;                                              % Number of Multipaths
commcell.frequency       = 2e9;                                            % Carrier frequency in Hz
commcell.meanShadowFad   = 0;                                              % Shadow fading mean in dB
commcell.stdDevShadowFad = 8;                                              % Shadow fading standard deviation in dB
commcell.city            = 'large';                                        % Type of city

linkprop.bsPower         = 10;                                             % in Watts
linkprop.userPower       = 1;                                         % in Watts
linkprop.AntennaGainBS   = 0;                                              % in dBi
linkprop.AntennaGainUser = 0;                                              % in dBi
linkprop.noiseFigureBS   = 9;                                              % in dB
linkprop.noiseFigureUser = 9 ;                                             % in dB
linkprop.bandwidth       = 20e6;                                           % in Hz

[snr_u_db,snr_d_db] = linkBudgetCalculation(linkprop);                     % SNR in dB
                
beta_db = -129.1419;

snr_u_eff = round(snr_u_db + beta_db);
snr_d_eff = round(snr_d_db + beta_db);

%snr_u_eff = -20;
%snr_d_eff = -3;

snr_u = 10.^((snr_u_eff)/10);                                              % Uplink SNR
snr_d = 10.^((snr_d_eff)/10);                                              % Downlink SNR

% Initialization

H = zeros(M,K,MC);                                                         % Channel matrix

user_set_rs    = zeros(L,MC);
user_set_sos   = zeros(L,MC);
user_set_icibs = zeros(L,MC);

gamma_u = zeros(K,MC);
gamma_d = zeros(K,MC);

gamma_rs_u = zeros(L,MC);
gamma_rs_d = zeros(L,MC);

gamma_sos_u = zeros(L,MC);
gamma_sos_d = zeros(L,MC);

gamma_icibs_u = zeros(L,MC);
gamma_icibs_d = zeros(L,MC);

rate_u  = zeros(K,MC);
rate_d  = zeros(K,MC);

rate_rs_u  = zeros(L,MC);
rate_rs_d  = zeros(L,MC);

rate_sos_u  = zeros(L,MC);
rate_sos_d  = zeros(L,MC);

rate_icibs_u  = zeros(L,MC);
rate_icibs_d  = zeros(L,MC);

psi       = zeros(K,MC);
psi_rs    = zeros(L,MC);
psi_sos   = zeros(L,MC);
psi_icibs = zeros(L,MC);

for out_mc = 1:MC
    out_mc
    
    [H(:,:,out_mc),beta] = massiveMIMOChannel(commcell,'ur-los');
    
    H(:,:,out_mc) = H(:,:,out_mc)*sqrt(diag(1./beta));
    
    % No Selection
    
    [rate_u(:,out_mc),gamma_u(:,out_mc)] = rateCalculation(H(:,:,out_mc),snr_u,'uplink');
    [rate_d(:,out_mc),gamma_d(:,out_mc)] = rateCalculation(H(:,:,out_mc),snr_d,'downlink');
    
    psi(:,out_mc)   = ici(H(:,:,out_mc));
    
    % Random Selection
    
    [user_set_rs(:,out_mc),H_rs] = userScheduling(H(:,:,out_mc),L, ...
                                                  'random selection');
    
    [rate_rs_u(:,out_mc),gamma_rs_u(:,out_mc)] = rateCalculation(H_rs,snr_u,'uplink');
    [rate_rs_d(:,out_mc),gamma_rs_d(:,out_mc)] = rateCalculation(H_rs,snr_d,'downlink');
    
    psi_rs(:,out_mc) = ici(H_rs);

    % Semi-orthogonal Selection
    
    [user_set_sos(:,out_mc),H_sos] = userScheduling(H(:,:,out_mc),L, ...
                                                    'semi-orthogonal selection');
    
    [rate_sos_u(:,out_mc),gamma_sos_u(:,out_mc)] = rateCalculation(H_sos,snr_u,'uplink');
    [rate_sos_d(:,out_mc),gamma_sos_d(:,out_mc)] = rateCalculation(H_sos,snr_d,'downlink');

    psi_sos(:,out_mc)   = ici(H_sos);
    
    % ICI-based Selection
    
    [user_set_icibs(:,out_mc),H_icibs] = userScheduling(H(:,:,out_mc),L, ...
                                                      'ici-based selection');
    
    [rate_icibs_u(:,out_mc),gamma_icibs_u(:,out_mc)] = rateCalculation(H_icibs,snr_u,'uplink');
    [rate_icibs_d(:,out_mc),gamma_icibs_d(:,out_mc)] = rateCalculation(H_icibs,snr_d,'downlink');

    psi_icibs(:,out_mc)   = ici(H_icibs);
end

save(['./results/rate_mf_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_u_eff) '_dB_MC_' num2str(MC) '.mat'], ...
      'H', ...
      'gamma_u','rate_u','gamma_d','rate_d','psi', ...
      'gamma_rs_u','rate_rs_u','gamma_rs_d','rate_rs_d','psi_rs','user_set_rs', ...
      'gamma_sos_u','rate_sos_u','gamma_sos_d','rate_sos_d','psi_sos','user_set_sos', ...
      'gamma_icibs_u','rate_icibs_u','gamma_icibs_d','rate_icibs_d','psi_icibs','user_set_icibs');