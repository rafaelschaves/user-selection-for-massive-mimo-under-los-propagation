clear;
close all;
clc;

addpath('./functions/')

MC = 10000;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)

M = 256;                                                                   % Number of antennas at the base station
K = 5;                                                                    % Number of users at the cell
L = 3;

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

linkprop.bsPower         = 10;                                             % in Watts
linkprop.userPower       = 0.2;                                            % in Watts
linkprop.AntennaGainBS   = 0;                                              % in dBi
linkprop.AntennaGainUser = 0;                                              % in dBi
linkprop.noiseFigureBS   = 9;                                              % in dB
linkprop.noiseFigureUser = 9 ;                                             % in dB
linkprop.bandwidth       = 20e6;                                           % in Hz

[snr_u_db,snr_d_db] = linkBudgetCalculation(linkprop);                     % SNR in dB
                
beta_db = -135;

snr_u_eff = round(snr_u_db + beta_db);
snr_d_eff = round(snr_d_db + beta_db);

snr_u = 10.^((snr_u_eff)/10);                                              % Uplink SNR
snr_d = 10.^((snr_d_eff)/10);                                              % Downlink SNR

% Initialization

H = zeros(M,K,MC);                                                         % Channel matrix

user_set_rs    = zeros(L,MC);
user_set_sos   = zeros(L,MC);
user_set_cbs   = zeros(L,MC);
user_set_icibs = zeros(L,MC);

gamma_u = zeros(K,MC);
gamma_d = zeros(K,MC);

gamma_rs_u = zeros(L,MC);
gamma_rs_d = zeros(L,MC);

gamma_sos_u = zeros(L,MC);
gamma_sos_d = zeros(L,MC);

gamma_cbs_u = zeros(L,MC);
gamma_cbs_d = zeros(L,MC);

gamma_icibs_u = zeros(L,MC);
gamma_icibs_d = zeros(L,MC);

rate_u  = zeros(K,MC);
rate_d  = zeros(K,MC);

rate_rs_u  = zeros(L,MC);
rate_rs_d  = zeros(L,MC);

rate_sos_u  = zeros(L,MC);
rate_sos_d  = zeros(L,MC);

rate_cbs_u  = zeros(L,MC);
rate_cbs_d  = zeros(L,MC);

rate_icibs_u  = zeros(L,MC);
rate_icibs_d  = zeros(L,MC);

psi       = zeros(K,MC);
psi_rs    = zeros(L,MC);
psi_sos   = zeros(L,MC);
psi_cbs   = zeros(L,MC);
psi_icibs = zeros(L,MC);

for out_mc = 1:MC
    out_mc
    
    [H(:,:,out_mc),beta] = massiveMIMOChannel(commcell,'ur-los');
    
    H(:,:,out_mc) = H(:,:,out_mc)*sqrt(diag(1./beta));
        
    % No Selection
    
    h_norm_ns     = vecnorm(H(:,:,out_mc));
    h_norm_ns_mtx = repmat(h_norm_ns,M,1);
    
    H_norm_ns = H(:,:,out_mc)./h_norm_ns_mtx;    
    
    Q_mf_ns = H_norm_ns;
    W_mf_ns = conj(H_norm_ns);

    pow_upl_ns = ones(K,1);
    pow_dow_ns = ones(K,1)/K;
    
    [rate_u(:,out_mc),gamma_u(:,out_mc)] = rateCalculation(H(:,:,out_mc),Q_mf_ns,pow_upl_ns,snr_u,'uplink');
    [rate_d(:,out_mc),gamma_d(:,out_mc)] = rateCalculation(H(:,:,out_mc),W_mf_ns,pow_dow_ns,snr_d,'downlink');
    
    psi(:,out_mc)   = ici(H(:,:,out_mc));
    
    % Random Selection
    
    [user_set_rs(:,out_mc),H_rs] = userScheduling(H(:,:,out_mc),L,'random selection');
    
    h_norm_rs     = vecnorm(H_rs);
    h_norm_rs_mtx = repmat(h_norm_rs,M,1);
    
    H_norm_rs = H_rs./h_norm_rs_mtx;    
    
    Q_mf_rs = H_norm_rs;
    W_mf_rs = conj(H_norm_rs);

    pow_upl_rs = ones(L,1);
    pow_dow_rs = ones(L,1)/L;
                                              
    [rate_rs_u(:,out_mc),gamma_rs_u(:,out_mc)] = rateCalculation(H_rs,Q_mf_rs,pow_upl_rs,snr_u,'uplink');
    [rate_rs_d(:,out_mc),gamma_rs_d(:,out_mc)] = rateCalculation(H_rs,W_mf_rs,pow_dow_rs,snr_d,'downlink');
    
    psi_rs(:,out_mc) = ici(H_rs);

    % Semi-orthogonal Selection
    
    [user_set_sos(:,out_mc),H_sos] = userScheduling(H(:,:,out_mc),L,'semi-orthogonal selection');
    
    h_norm_sos     = vecnorm(H_sos);
    h_norm_sos_mtx = repmat(h_norm_sos,M,1);
    
    H_norm_sos = H_sos./h_norm_sos_mtx;    
    
    Q_mf_sos = H_norm_sos;
    W_mf_sos = conj(H_norm_sos);

    pow_upl_sos = ones(L,1);
    pow_dow_sos = ones(L,1)/L;
                                                
    [rate_sos_u(:,out_mc),gamma_sos_u(:,out_mc)] = rateCalculation(H_sos,Q_mf_sos,pow_upl_sos,snr_u,'uplink');
    [rate_sos_d(:,out_mc),gamma_sos_d(:,out_mc)] = rateCalculation(H_sos,W_mf_sos,pow_dow_sos,snr_d,'downlink');

    psi_sos(:,out_mc)   = ici(H_sos);
    
    % Correlation-based Selection
    
    [user_set_cbs(:,out_mc),H_cbs] = userScheduling(H(:,:,out_mc),L,'correlation-based selection');
    
    h_norm_cbs     = vecnorm(H_cbs);
    h_norm_cbs_mtx = repmat(h_norm_cbs,M,1);
    
    H_norm_cbs = H_cbs./h_norm_cbs_mtx;    
    
    Q_mf_cbs = H_norm_cbs;
    W_mf_cbs = conj(H_norm_cbs);

    pow_upl_cbs = ones(L,1);
    pow_dow_cbs = ones(L,1)/L;
                                                
    [rate_cbs_u(:,out_mc),gamma_cbs_u(:,out_mc)] = rateCalculation(H_cbs,Q_mf_cbs,pow_upl_cbs,snr_u,'uplink');
    [rate_cbs_d(:,out_mc),gamma_cbs_d(:,out_mc)] = rateCalculation(H_cbs,W_mf_cbs,pow_dow_cbs,snr_d,'downlink');

    psi_cbs(:,out_mc)   = ici(H_cbs);
    
    % ICI-based Selection
    
    [user_set_icibs(:,out_mc),H_icibs] = userScheduling(H(:,:,out_mc),L,'ici-based selection');
    
    h_norm_icibs     = vecnorm(H_icibs);
    h_norm_icibs_mtx = repmat(h_norm_icibs,M,1);
    
    H_norm_icibs = H_icibs./h_norm_icibs_mtx;    
    
    Q_mf_icibs = H_norm_icibs;
    W_mf_icibs = conj(H_norm_icibs);

    pow_upl_icibs = ones(L,1);
    pow_dow_icibs = ones(L,1)/L;
    
    [rate_icibs_u(:,out_mc),gamma_icibs_u(:,out_mc)] = rateCalculation(H_icibs,Q_mf_icibs,pow_upl_icibs,snr_u,'uplink');
    [rate_icibs_d(:,out_mc),gamma_icibs_d(:,out_mc)] = rateCalculation(H_icibs,W_mf_icibs,pow_dow_icibs,snr_d,'downlink');

    psi_icibs(:,out_mc)   = ici(H_icibs);
end

save(['./results/rate_mf_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_u_eff) '_dB_MC_' num2str(MC) '.mat'], ...
      'H', ...
      'gamma_u','rate_u','gamma_d','rate_d','psi', ...
      'gamma_rs_u','rate_rs_u','gamma_rs_d','rate_rs_d','psi_rs','user_set_rs', ...
      'gamma_sos_u','rate_sos_u','gamma_sos_d','rate_sos_d','psi_sos','user_set_sos', ...
      'gamma_icibs_u','rate_icibs_u','gamma_icibs_d','rate_icibs_d','psi_icibs','user_set_icibs');