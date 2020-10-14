clear;
close all;
clc;

% Macros

MC   = 1000;                                                              % Size of the monte-carlo ensemble
N_MC = 5;

M = 50;                                                                   % Number of antennas at base station
K = [50 75];                                                        % Number of users at the cell 

% M = 50  & K = [10 25 50 75]
% M = 100 & K = [10 25 50 75 100 150] 
% M = 200 & K = [10 25 50 75 100 150 200 250]

snr = -5;

err = [0 pi/(6*M) pi/(3*M) pi/(2*M)];

N_K = length(K);
N_ALG = 3;                                                                 % Number of algorithms for perform user scheduling
N_PRE = 3;
N_ERR = length(err);

err = 0:ERR_STE:MAX_ERR;

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Downlink/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Selection/Downlink/';

zero_pad_1 = '%03d';
zero_pad_2 = '%02d';

chn_type = 'ur_los';

avg_sum_se = zeros(N_PRE,N_ERR,N_K);

for k = 1:N_K
    se_all_mc = zeros(K(k),N_PRE,N_ERR,MC*N_MC);
    
    for n_mc = 1:N_MC
        load([root_load 'spectral_efficiency_all_L_' chn_type '_partial_csi_2_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K(k)) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC) '_' sprintf(zero_pad_2,n_mc) '.mat']);
        
        idx_ini = (n_mc - 1)*MC + 1;
        idx_end = n_mc*MC;
        
        se_all_mc(:,:,:,idx_ini:idx_end) = se;
        
        clear se se_s_all_L;
    end
    
    sum_se = reshape(sum(se_all_mc,1),N_PRE,N_ERR,MC*N_MC);
    
    avg_sum_se(:,:,k) = mean(sum_se,3);
end

% N_BIN = 15;

% N_BIN = 100;
% 
% cdf_sum_se = cell(N_ALG+1,N_PRE,N_ERR);
% edg_sum_se = cell(N_ALG+1,N_PRE,N_ERR);
% 
% for n_err = 1:N_ERR
%     for n_pre = 1:N_PRE
%         [cdf_sum_se{1,n_pre},edg_sum_se{1,n_pre}] = histcounts(sum_se(n_pre,n_err,:),N_BIN,'normalization','cdf');
%         
%         for n_alg = 1:N_ALG
%             [cdf_sum_se{n_alg+1,n_pre,n_err},edg_sum_se{n_alg+1,n_pre,n_err}] = histcounts(sum_se_s_star(:,n_pre,n_alg,n_err),N_BIN,'normalization','cdf');
%         end
%     end
% end

% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 30;

marker = {'o','s','^'};

linestyle = {'-','--',':'};

savefig = 1;

% NS - No selection
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_K      = {'$K = 10$','$K = 25$','$K = 50$','$K = 75$'};
legend_prec   = {'MRT','ZF','MMSE'};
legend_algo   = {'NS','SOS','CBS','ICIBS'};
legend_algo_2 = {'SOS','CBS','ICIBS'};

location_1 = 'northwest';
location_2 = 'northeast';
location_3 = 'southwest';
location_4 = 'southeast';

colours = [0.0000 0.4470 0.7410;
           0.8500 0.3250 0.0980;
           0.9290 0.6940 0.1250;
           0.4940 0.1840 0.5560;
           0.4660 0.6740 0.1880;
           0.3010 0.7450 0.9330;
           0.6350 0.0780 0.1840];

figure;
       
set(gcf,'position',[0 0 800 600]);

plot(err,avg_sum_se(1,:,1),linestyle{1},'color',colours(1,:),'linewidth',linewidth);
hold on;
plot(err,avg_sum_se(1,:,2),linestyle{1},'color',colours(2,:),'linewidth',linewidth);
plot(err,avg_sum_se(1,:,3),linestyle{1},'color',colours(3,:),'linewidth',linewidth);
plot(err,avg_sum_se(1,:,4),linestyle{1},'color',colours(4,:),'linewidth',linewidth);
plot(err,avg_sum_se(2,:,1),linestyle{2},'color',colours(1,:),'linewidth',linewidth);
plot(err,avg_sum_se(2,:,2),linestyle{2},'color',colours(2,:),'linewidth',linewidth);
plot(err,avg_sum_se(2,:,3),linestyle{2},'color',colours(3,:),'linewidth',linewidth);
plot(err,avg_sum_se(3,:,1),linestyle{3},'color',colours(1,:),'linewidth',linewidth);
plot(err,avg_sum_se(3,:,2),linestyle{3},'color',colours(2,:),'linewidth',linewidth);
plot(err,avg_sum_se(3,:,3),linestyle{3},'color',colours(3,:),'linewidth',linewidth);

xticks([0 pi/144 pi/72]);
xticklabels({'0','1.25','2.5'})

xlabel('Maximum estimation error','fontname',fontname,'fontsize',fontsize);
ylabel('Sum-spectral efficiency','fontname',fontname,'fontsize',fontsize);

legend(legend_K,'fontname',fontname,'fontsize',fontsize,'location',location_2,'interpreter','latex');
legend box off;

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([0 pi/72]);
ylim([0 20]);

if (savefig == 1)
    saveas(gcf,[root_save 'se_partial_csi_' chn_type '_M_' sprintf(zero_pad_1,M) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
    saveas(gcf,[root_save 'se_partial_csi_' chn_type '_M_' sprintf(zero_pad_1,M) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
    saveas(gcf,[root_save 'se_partial_csi_' chn_type '_M_' sprintf(zero_pad_1,M) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
end