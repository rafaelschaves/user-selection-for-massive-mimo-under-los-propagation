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

user_idx_ici_based = zeros(MC,1);
user_idx_random = zeros(MC,1);

gamma = zeros(K,MC);
rate  = zeros(K,MC);
psi   = zeros(K,MC);

gamma_ici_based = zeros(K-1,MC);
rate_ici_based  = zeros(K-1,MC);
psi_ici_based   = zeros(K-1,MC);

gamma_random = zeros(K-1,MC);
rate_random  = zeros(K-1,MC);
psi_random   = zeros(K-1,MC);

for out_mc = 1:MC
    out_mc
    
    [H(:,:,out_mc),beta] = massiveMIMOChannel(commcell,'rayleigh');
    
    H(:,:,out_mc) = H(:,:,out_mc)*sqrt(diag(1./beta));
        
    gamma(:,out_mc) = sinr(H(:,:,out_mc),snr); 
    rate(:,out_mc)  = log2(1 + gamma(:,out_mc));
    psi(:,out_mc)   = ici(H(:,:,out_mc));
    
    % Removing user randomly
    
    user_idx_random(out_mc) = randi([1 K]);
    
    H_random = H(:,:,out_mc);
    H_random(:,user_idx_random(out_mc)) = [];
    
    gamma_random(:,out_mc) = sinr(H_random,snr); 
    rate_random(:,out_mc)  = log2(1 + gamma_random(:,out_mc));
    psi_random(:,out_mc)   = ici(H_random);
    
    % Removing user based in ICI
    
    [~,user_idx_ici_based(out_mc)] = max(psi(:,out_mc));
    
    H_ici_based = H(:,:,out_mc);
    H_ici_based(:,user_idx_ici_based(out_mc)) = [];
    
    gamma_ici_based(:,out_mc) = sinr(H_ici_based,snr); 
    rate_ici_based(:,out_mc)  = log2(1 + gamma_ici_based(:,out_mc));
    psi_ici_based(:,out_mc)   = ici(H_ici_based);
end

save(['./results/rate_mf_M_' num2str(M) '_K_' num2str(K) '_MC_' num2str(MC) '.mat'], ...
      'gamma','rate','psi', ...
      'gamma_random','rate_random','psi_random', ...
      'gamma_ici_based','rate_ici_based','psi_ici_based','user_idx_ici_based');