clear;
close all;
clc;

addpath('./functions/')

root_downlink = './results/scheduling/downlink/rate_downlink_mf_';
root_uplink   = './results/scheduling/uplink/rate_uplink_mf_';

MC    = 10000;                                                             % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
N_ALG = 4;

M = 256;                                                                    % Number of antennas at the base station
K = 18;                                                                    % Number of users at the cell
L = 13;

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

% snr_u_eff = round(snr_u_db + beta_db);
% snr_d_eff = round(snr_d_db + beta_db);

snr_u_eff = 10;
snr_d_eff = 10;

snr_u = 10.^((snr_u_eff)/10);                                              % Uplink SNR
snr_d = 10.^((snr_d_eff)/10);                                              % Downlink SNR

% Initialization

user_set = zeros(L,MC,N_ALG);

gamma_u = zeros(K,MC);
gamma_d = zeros(K,MC);

rate_u = zeros(K,MC);
rate_d = zeros(K,MC);

psi = zeros(K,MC);

gamma_u_alg = zeros(L,MC,N_ALG);
gamma_d_alg = zeros(L,MC,N_ALG);

rate_u_alg = zeros(L,MC,N_ALG);
rate_d_alg = zeros(L,MC,N_ALG);

psi_alg = zeros(L,MC,N_ALG);

channel_type = 'rayleigh';

for out_mc = 1:MC
    out_mc
    
    [G,beta] = massiveMIMOChannel(commcell,channel_type);
    
    % No Selection
    
    h_norm_ns     = vecnorm(G);
    h_norm_ns_mtx = repmat(h_norm_ns,M,1);
    
    H_norm_ns = G./h_norm_ns_mtx;    
    
    Q_mf_ns = H_norm_ns;
    W_mf_ns = conj(H_norm_ns);

    pow_upl_ns = ones(K,1);
    pow_dow_ns = ones(K,1)/K;
    
    [rate_u(:,out_mc),gamma_u(:,out_mc)] = rateCalculation(G,Q_mf_ns,pow_upl_ns,snr_u,'uplink');
    [rate_d(:,out_mc),gamma_d(:,out_mc)] = rateCalculation(G,W_mf_ns,pow_dow_ns,snr_d,'downlink');
    
    psi(:,out_mc)   = ici(G);
    
    % Random Selection
    
    [user_set(:,out_mc,1),H_rs] = userScheduling(G,L,'random selection');
    
    h_norm_rs     = vecnorm(H_rs);
    h_norm_rs_mtx = repmat(h_norm_rs,M,1);
    
    H_norm_rs = H_rs./h_norm_rs_mtx;    
    
    Q_mf_rs = H_norm_rs;
    W_mf_rs = conj(H_norm_rs);

    pow_upl_rs = ones(L,1);
    pow_dow_rs = ones(L,1)/L;
                                              
    [rate_u_alg(:,out_mc,1),gamma_u_alg(:,out_mc,1)] = rateCalculation(H_rs,Q_mf_rs,pow_upl_rs,snr_u,'uplink');
    [rate_d_alg(:,out_mc,1),gamma_d_alg(:,out_mc,1)] = rateCalculation(H_rs,W_mf_rs,pow_dow_rs,snr_d,'downlink');
    
    psi_alg(:,out_mc,1) = ici(H_rs);

    % Semi-orthogonal Selection
    
    [user_set(:,out_mc,2),H_sos] = userScheduling(G,L,'semi-orthogonal selection');
    
    h_norm_sos     = vecnorm(H_sos);
    h_norm_sos_mtx = repmat(h_norm_sos,M,1);
    
    H_norm_sos = H_sos./h_norm_sos_mtx;    
    
    Q_mf_sos = H_norm_sos;
    W_mf_sos = conj(H_norm_sos);

    pow_upl_sos = ones(L,1);
    pow_dow_sos = ones(L,1)/L;
                                                
    [rate_u_alg(:,out_mc,2),gamma_u_alg(:,out_mc,2)] = rateCalculation(H_sos,Q_mf_sos,pow_upl_sos,snr_u,'uplink');
    [rate_d_alg(:,out_mc,2),gamma_d_alg(:,out_mc,2)] = rateCalculation(H_sos,W_mf_sos,pow_dow_sos,snr_d,'downlink');

    psi_alg(:,out_mc,2) = ici(H_sos);
    
    % Correlation-based Selection
    
    [user_set(:,out_mc,3),H_cbs] = userScheduling(G,L,'correlation-based selection');
    
    h_norm_cbs     = vecnorm(H_cbs);
    h_norm_cbs_mtx = repmat(h_norm_cbs,M,1);
    
    H_norm_cbs = H_cbs./h_norm_cbs_mtx;    
    
    Q_mf_cbs = H_norm_cbs;
    W_mf_cbs = conj(H_norm_cbs);

    pow_upl_cbs = ones(L,1);
    pow_dow_cbs = ones(L,1)/L;
                                                
    [rate_u_alg(:,out_mc,3),gamma_u_alg(:,out_mc,3)] = rateCalculation(H_cbs,Q_mf_cbs,pow_upl_cbs,snr_u,'uplink');
    [rate_d_alg(:,out_mc,3),gamma_d_alg(:,out_mc,3)] = rateCalculation(H_cbs,W_mf_cbs,pow_dow_cbs,snr_d,'downlink');

    psi_alg(:,out_mc,3) = ici(H_cbs);
    
    % ICI-based Selection
    
    [user_set(:,out_mc,4),H_icibs] = userScheduling(G,L,'ici-based selection');
    
    h_norm_icibs     = vecnorm(H_icibs);
    h_norm_icibs_mtx = repmat(h_norm_icibs,M,1);
    
    H_norm_icibs = H_icibs./h_norm_icibs_mtx;    
    
    Q_mf_icibs = H_norm_icibs;
    W_mf_icibs = conj(H_norm_icibs);

    pow_upl_icibs = ones(L,1);
    pow_dow_icibs = ones(L,1)/L;
    
    [rate_u_alg(:,out_mc,4),gamma_u_alg(:,out_mc,4)] = rateCalculation(H_icibs,Q_mf_icibs,pow_upl_icibs,snr_u,'uplink');
    [rate_d_alg(:,out_mc,4),gamma_d_alg(:,out_mc,4)] = rateCalculation(H_icibs,W_mf_icibs,pow_dow_icibs,snr_d,'downlink');

    psi_alg(:,out_mc,4) = ici(H_icibs);
end

save([root_downlink channel_type '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_u_eff) '_dB_MC_' num2str(MC) '.mat'], ...
      'rate_d','psi','rate_d_alg','user_set','psi_alg');

save([root_uplink channel_type '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_u_eff) '_dB_MC_' num2str(MC) '.mat'], ...
      'rate_u','psi','rate_u_alg','user_set','psi_alg');