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

snr_db = -20:2:10;                                                         % SNR in dB
snr    = 10.^(snr_db/10);                                                  % SNR
n_snr  = length(snr_db);

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

b_hat = zeros(B*N,K);                                                      % Estimated message in bits

H = zeros(M,K,OUTER_MC);                                                   % Channel matrix

ber = zeros(n_snr,K,INNER_MC,OUTER_MC);                                    % Bit-error rate

for out_mc = 1:OUTER_MC
    out_mc
    
    [H(:,:,out_mc),beta] = massiveMIMOChannel(commcell,'rayleigh');
    
    H(:,:,out_mc) = H(:,:,out_mc)*sqrt(diag(1./beta));
    
    for inn_mc = 1:INNER_MC        
        [s,Ps,b] = userTX(K,N,B);                                          % Transmitted signals by each user
        
        % Decoder parameters
            
        decpar = struct('decoder','mf','power',Ps);
        
        for snr_idx = 1:n_snr
            y = channel(s,Ps,H(:,:,out_mc),snr(snr_idx),'uplink');         % Channel
            
            b_hat = baseStationRX(y,H(:,:,out_mc),decpar,B);               % Base station receiver
            
            [~,ber(snr_idx,:,inn_mc,out_mc)] = biterr(b_hat, ...           % Estimated bits at the base station
                                                      b, ...               % Transmitted bits by users
                                                      [], ...
                                                      'column-wise');      % Bit-error rate calculation
        end
    end
end

save(['ber_mf_M_'  num2str(M) '_K_' num2str(K) '_N_' num2str(N) '_MC_' num2str(OUTER_MC) '.mat'],'ber','H');