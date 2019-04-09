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

gamma = zeros(K,OUTER_MC);
gamma_prime = zeros(K,OUTER_MC);
rate  = zeros(K,OUTER_MC);
psi   = zeros(K,OUTER_MC);

for out_mc = 1:OUTER_MC
    out_mc
    
    [H(:,:,out_mc),beta] = massiveMIMOChannel(commcell,'rayleigh');
    
    H(:,:,out_mc) = H(:,:,out_mc)*sqrt(diag(1./beta));
    
    psi(:,out_mc) = ici(H(:,:,out_mc));
    
    gamma_prime(:,out_mc) = sinr_prime(H(:,:,out_mc),snr); 
    
    for k = 1:K
        gamma(k,out_mc) = sinr(H(:,:,out_mc),k,snr);
        rate(k,out_mc) = log(1+gamma(k,out_mc));
    end
end

% Ploting Figures

linewidth = 2;
fontname  = 'Times New Roman';
fontsize  = 20;

BIN_WIDTH_CDF  = 0.005;

% root_erg_cap  = './Figures/Capacity/erg_cap';
% root_out_prob = './Figures/Capacity/out_prob';

[values, edges] = histcounts(sum(rate),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');

figure;

set(gcf,'position',[0 0 800 600]);

plot(edges,[values 1],'linewidth',linewidth);

xlabel('Capacity (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Outage Probability','fontname',fontname,'fontsize',fontsize);

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 1]);

% saveas(gcf,[root_out_prob '_M_' num2str(M(m))],'fig');
% saveas(gcf,[root_out_prob '_M_' num2str(M(m))],'png');
% saveas(gcf,[root_out_prob '_M_' num2str(M(m))],'epsc2');

% save(['ber_' decpar.decoder '_M_'  num2str(M) '_K_' num2str(K) '_N_' num2str(N) '_MC_' num2str(MONTE_CARLO) '.mat'],'ber','H');

function [gamma] = sinr(H,k,snr)

h_k = H(:,k);

H(:,k) = [];

gamma = norm(h_k,2)^2/(1/snr + sum(abs(h_k'*H).^2/norm(h_k,2)^2));

end

function [gamma] = sinr_prime(H,snr)

M = size(H,1);

h_norm = vecnorm(H,2)';
H_norm = repmat(h_norm,M,1);

H_n = H./H_norm;

D = H_n'*H;

int = sum(abs(D),2) - h_norm';

gamma = (h_nrom.^2/(1/snr + sum(abs(int).^2))';

end

function [psi] = ici(H)

M = size(H,1);

h_norm = vecnorm(H,2);
H_norm = repmat(h_norm,M,1);

H_n = H./H_norm;

D_n = H_n'*H_n;

psi = sum(abs(D_n),2) - 1;
end