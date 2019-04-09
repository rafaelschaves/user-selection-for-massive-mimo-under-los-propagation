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

H = zeros(M,K,OUTER_MC);                                                   % Channel matrix

user_idx = zeros(OUTER_MC,1);

gamma = zeros(K,OUTER_MC);
rate  = zeros(K,OUTER_MC);
psi   = zeros(K,OUTER_MC);

gamma_prime = zeros(K-1,OUTER_MC);
rate_prime  = zeros(K-1,OUTER_MC);
psi_prime   = zeros(K-1,OUTER_MC);

for out_mc = 1:OUTER_MC
    out_mc
    
    [H(:,:,out_mc),beta] = massiveMIMOChannel(commcell,'rayleigh');
    
    H(:,:,out_mc) = H(:,:,out_mc)*sqrt(diag(1./beta));
        
    gamma(:,out_mc) = sinr(H(:,:,out_mc),snr); 
    rate(:,out_mc)  = log2(1 + gamma(:,out_mc));
    psi(:,out_mc)   = ici(H(:,:,out_mc));
    
    % Removing user based in ICI
    
    [~,user_idx(out_mc)] = max(psi(:,out_mc));
    
    H_aux = H(:,:,out_mc);
    H_aux(:,user_idx(out_mc)) = [];
    
    gamma_prime(:,out_mc) = sinr(H_aux,snr); 
    rate_prime(:,out_mc)  = log2(1 + gamma_prime(:,out_mc));
    psi_prime(:,out_mc)   = ici(H_aux);
end

save(['./results/rate_mf_M_' num2str(M) '_K_' num2str(K) '_MC_' num2str(MC) '.mat'], ...
      'gamma','rate','psi','gamma_prime','rate_prime','psi_prime','user_idx');