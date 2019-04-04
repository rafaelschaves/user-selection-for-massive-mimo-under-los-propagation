clear;
close all;
clc;

addpath('./functions/')

MONTE_CARLO = 1;                                                       % Size of the Monte Carlo ensemble

B = 4;                                                                     % Number of bits in each symbol
N = 750;                                                                   % Number of blocks in the transmission

M = 500;                                                                   % Number of antennas at the base station
K = 5;                                                                     % Number of users at the cell

snr_db = -20:1:10;                                                         % SNR in dB
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

H = zeros(M,K,MONTE_CARLO);                                                % Channel matrix

ber = zeros(n_snr,K,MONTE_CARLO);                                          % Bit-error rate

for mc_idx = 1:MONTE_CARLO
    
    mc_idx
    
    [s,Ps,b] = userTX(K,N,B);                                              % Signal generation for each user
    
    % Decoder parameters
    
    decpar.decoder = 'mf';
    decpar.power    = Ps;
    
    for snr_idx = 1:n_snr
        
        [y,H(:,:,mc_idx)] = channel(s,Ps,snr(snr_idx),commcell,'rayleigh');
                                                           
        s_hat = decoder(y,H(:,:,mc_idx),decpar);                           % Decoding received signal
        
        % Signal decodification for each user
        
        for k = 1:K
            b_hat(:,k) = qamdemod(s_hat(k,:).',2^B,'OutputType','bit');
        end
        
        [~,ber(snr_idx,:,mc_idx)] = biterr(b_hat,b,[],'column-wise');
        
    end
    
end

% save(['ber_' decpar.precoder '_M_'  num2str(M) '_K_' num2str(K) '_N_' num2str(N) '_MC_' num2str(MONTE_CARLO) '.mat'],'ber','H','D');