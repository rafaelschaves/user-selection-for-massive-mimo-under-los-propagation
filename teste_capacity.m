clear;
close all;
clc;

addpath('./functions/')

MC = 10000;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)

M = 500;                                                                   % Number of antennas at the base station
K = 5;                                                                     % Number of users at the cell

snr_db = 0;                                                                % SNR in dB
snr    = 10.^(snr_db/10);                                                  % SNR

commcell.nAntennas       = M;                                              % Number of Antennas
commcell.nUsers          = K;                                              % Number of Users
commcell.radius          = 500;                                            % Cell's raidus (circumradius) in meters
commcell.bsHeight        = 30;                                             % Height of base station in meters
commcell.userHeight      = [1 2];                                          % Height of user terminals in meters ([min max])
commcell.nPaths          = 1;                                              % Number of Multipaths
commcell.frequency       = 2e9;                                            % Carrier frequency in Hz
commcell.meanShadowFad   = 0;                                              % Shadow fading mean in dB
commcell.stdDevShadowFad = 8;                                              % Shadow fading standard deviation in dB
commcell.city            = 'large';                                        % Type of city

% Initialization

H = zeros(M,K,MC);                                                         % Channel matrix

idx_aux = zeros(K-1,MC);

user_idx_rs    = zeros(MC,1);
user_idx_sos   = zeros(MC,1);
user_idx_icibs = zeros(MC,1);

gamma_u = zeros(K,MC);
gamma_d = zeros(K,MC);

gamma_rs_u = zeros(K-1,MC);
gamma_rs_d = zeros(K-1,MC);

gamma_sos_u = zeros(K-1,MC);
gamma_sos_d = zeros(K-1,MC);

gamma_icibs_u = zeros(K-1,MC);
gamma_icibs_d = zeros(K-1,MC);

rate_u  = zeros(K,MC);
rate_d  = zeros(K,MC);

rate_rs_u  = zeros(K-1,MC);
rate_rs_d  = zeros(K-1,MC);

rate_sos_u  = zeros(K-1,MC);
rate_sos_d  = zeros(K-1,MC);

rate_icibs_u  = zeros(K-1,MC);
rate_icibs_d  = zeros(K-1,MC);

psi       = zeros(K,MC);
psi_rs    = zeros(K-1,MC);
psi_sos   = zeros(K-1,MC);
psi_icibs = zeros(K-1,MC);

for out_mc = 1:MC
    out_mc
    
    [H(:,:,out_mc),beta] = massiveMIMOChannel(commcell,'rayleigh');
    
    H(:,:,out_mc) = H(:,:,out_mc)*sqrt(diag(1./beta));
    
    % No Selection
    
    [rate_u(:,out_mc),gamma_u(:,out_mc)] = rateCalculation(H(:,:,out_mc),snr,'uplink');
    [rate_d(:,out_mc),gamma_d(:,out_mc)] = rateCalculation(H(:,:,out_mc),snr,'downlink');
    
    psi(:,out_mc)   = ici(H(:,:,out_mc));
    
    % Random Selection
    
    user_idx_rs(out_mc) = randi([1 K]);
    
    H_rs = H(:,:,out_mc);
    H_rs(:,user_idx_rs(out_mc)) = [];
    
    [rate_rs_u(:,out_mc),gamma_rs_u(:,out_mc)] = rateCalculation(H_rs,snr,'uplink');
    [rate_rs_d(:,out_mc),gamma_rs_d(:,out_mc)] = rateCalculation(H_rs,snr,'downlink');
    
    psi_rs(:,out_mc)   = ici(H_rs);
    
    % Semi-orthogonal Selection
    
    g_i = [];
    
    I_M = eye(M);
    
     H_aux = H(:,:,out_mc);
     
    for i = 1:K-1
        if(i == 1)
            G = H(:,:,out_mc);
            g_norm = vecnorm(G,2);
            
            [~,idx_aux(i,out_mc)] = max(g_norm);
            
            g_i = [g_i G(:,idx_aux(i,out_mc))/norm(G(:,idx_aux(i,out_mc)),2)];
            
            H_aux(:,idx_aux(i,out_mc)) = zeros(M,1);
        else
            P = I_M - sum(g_i*g_i.',2);
            G = P*H_aux;
            
            g_norm = vecnorm(G,2);
            
            [~,idx_aux(i,out_mc)] = max(g_norm);
            
            g_i = [g_i G(:,idx_aux(i,out_mc))/norm(G(:,idx_aux(i,out_mc)),2)];
            
            H_aux(:,idx_aux(i,out_mc)) = zeros(M,1);
        end
    end
    
    aux = 1:K;
    aux(idx_aux(:,out_mc)) = [];
    
    user_idx_sos(out_mc) = aux;
    
    H_sos = H(:,:,out_mc);
    H_sos(:,user_idx_sos(out_mc)) = [];
    
    [rate_sos_u(:,out_mc),gamma_sos_u(:,out_mc)] = rateCalculation(H_sos,snr,'uplink');
    [rate_sos_d(:,out_mc),gamma_sos_d(:,out_mc)] = rateCalculation(H_sos,snr,'downlink');

    psi_sos(:,out_mc)   = ici(H_sos);
    
    % ICI-based Selection
    
    [~,user_idx_icibs(out_mc)] = max(psi(:,out_mc));
    
    H_icibs = H(:,:,out_mc);
    H_icibs(:,user_idx_icibs(out_mc)) = [];
    
    [rate_icibs_u(:,out_mc),gamma_icibs_u(:,out_mc)] = rateCalculation(H_icibs,snr,'uplink');
    [rate_icibs_d(:,out_mc),gamma_icibs_d(:,out_mc)] = rateCalculation(H_icibs,snr,'downlink');

    psi_icibs(:,out_mc)   = ici(H_icibs);
end

save(['./results/rate_mf_M_' num2str(M) '_K_' num2str(K) '_MC_' num2str(MC) '.mat'], ...
      'H', ...
      'gamma_u','rate_u','gamma_d','rate_d','psi', ...
      'gamma_rs_u','rate_rs_u','gamma_rs_d','rate_rs_d','psi_rs','user_idx_rs', ...
      'gamma_sos_u','rate_sos_u','gamma_sos_d','rate_sos_d','psi_sos','user_idx_sos', ...
      'gamma_icibs_u','rate_icibs_u','gamma_icibs_d','rate_icibs_d','psi_icibs','user_idx_icibs');