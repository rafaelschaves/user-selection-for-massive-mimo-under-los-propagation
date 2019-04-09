clear;
close all;
clc;

addpath('./functions/')

OUTER_MC = 10000;                                                          % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
INNER_MC = 200;                                                            % Size of the inner Monte Carlo ensemble (Varies the messages for one channel realization)

B = 4;                                                                     % Number of bits in each symbol
N = 750;                                                                   % Number of blocks in the transmission

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

sinr = zeros(K,OUTER_MC);
rate = zeros(K,OUTER_MC);

for out_mc = 1:OUTER_MC
    out_mc
    
    [H(:,:,out_mc),beta] = massiveMIMOChannel(commcell,'rayleigh');
    
    H(:,:,out_mc) = H(:,:,out_mc)*sqrt(diag(1./beta));
    
    for k = 1:K
        H_aux      = H(:,:,out_mc);
        H_aux(:,k) = [];
        
        sinr(k,out_mc) = real(norm(H(:,k,out_mc),2)^2/(1/snr + sum(H(:,k,out_mc)'*H_aux/norm(H(:,k,out_mc),2)^2)));
        rate(k,out_mc) = log(1+sinr(k,out_mc));
    end
end

% save(['ber_' decpar.decoder '_M_'  num2str(M) '_K_' num2str(K) '_N_' num2str(N) '_MC_' num2str(MONTE_CARLO) '.mat'],'ber','H');