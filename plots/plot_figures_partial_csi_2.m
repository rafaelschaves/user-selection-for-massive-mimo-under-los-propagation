clear;
close all;
clc;

% Macros

MC   = 1000;                                                              % Size of the monte-carlo ensemble
N_MC = 5;

M = 50;                                                                   % Number of antennas at base station
K = 10;                                                                   % Number of users at the cell 

% M = 50  & K = [10 25 50 75]
% M = 100 & K = [10 25 50 75 100 150] 
% M = 200 & K = [10 25 50 75 100 150 200 250]

if K > M
    L_max = M;
else
    L_max = K-1;
end

snr = -5;

MAX_ERR = pi/72;
ERR_STE = pi/720;

N_ALG = 3;                                                                 % Number of algorithms for perform user scheduling
N_PRE = 3;
N_ERR = 1 + MAX_ERR/ERR_STE;

err = 0:ERR_STE:MAX_ERR;

% Roots

<<<<<<< HEAD
%root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Downlink/';
root_load = 'D:\PhD\user-selection\Partial CSI\';
=======
root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Downlink/';
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Selection/Downlink/';

zero_pad_1 = '%03d';
zero_pad_2 = '%02d';

chn_type = 'ur_los';

% Loading data

se_all_mc     = zeros(K,N_PRE,N_ERR,MC*N_MC);
se_s_L_all_mc = zeros(L_max,L_max,N_PRE,N_ALG,N_ERR,MC*N_MC);
sum_se_s      = zeros(L_max,N_PRE,N_ALG,N_ERR,MC*N_MC);

for n_mc = 1:N_MC
    load([root_load 'spectral_efficiency_all_L_' chn_type '_partial_csi_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC) '_' sprintf(zero_pad_2,n_mc) '.mat']);
    
    idx_ini = (n_mc - 1)*MC + 1;
    idx_end = n_mc*MC;
        
    se_all_mc(:,:,:,idx_ini:idx_end)         = se;
    se_s_L_all_mc(:,:,:,:,:,idx_ini:idx_end) = se_s_all_L;
    
    clear se se_s_all_L;
end

sum_se = reshape(sum(se_all_mc,1),N_PRE,N_ERR,MC*N_MC);

for L = 1:L_max
    sum_se_s(L,:,:,:,:) = sum(se_s_L_all_mc(:,L,:,:,:,:),1);
end

avg_sum_se   = mean(sum_se,3);
avg_sum_se_s = mean(sum_se_s,5);

% std_sum_se = std(sum_se,[],2);
% std_sum_se_s = std(sum_se_s,[],4);

[max_sum_se_s,L_star] = max(avg_sum_se_s,[],1);

max_sum_se_s = reshape(max_sum_se_s,N_PRE,N_ALG,N_ERR);
L_star       = reshape(L_star,N_PRE,N_ALG,N_ERR);

% L_star = 80*ones(3,3);

for n_err = 1:N_ERR
    for n_alg = 1:N_ALG
        for n_pre = 1:N_PRE
            sum_se_s_star(:,n_pre,n_alg,n_err) = sum_se_s(L_star(n_alg,n_pre),n_pre,n_alg,n_err,:);
        end
    end
end

% N_BIN = 15;

N_BIN = 100;

cdf_sum_se = cell(N_ALG+1,N_PRE,N_ERR);
edg_sum_se = cell(N_ALG+1,N_PRE,N_ERR);

for n_err = 1:N_ERR
    for n_pre = 1:N_PRE
        [cdf_sum_se{1,n_pre},edg_sum_se{1,n_pre}] = histcounts(sum_se(n_pre,n_err,:),N_BIN,'normalization','cdf');
        
        for n_alg = 1:N_ALG
            [cdf_sum_se{n_alg+1,n_pre,n_err},edg_sum_se{n_alg+1,n_pre,n_err}] = histcounts(sum_se_s_star(:,n_pre,n_alg,n_err),N_BIN,'normalization','cdf');
        end
    end
end

% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 30;

marker = {'o','s','^'};

linestyle = {'-','--',':'};

savefig = 0;

% NS - No selection
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

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

if K <= M
    plot(err,avg_sum_se(1,:),linestyle{1},'color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(err,avg_sum_se(2,:),linestyle{1},'color',colours(2,:),'linewidth',linewidth);
    plot(err,avg_sum_se(3,:),linestyle{1},'color',colours(3,:),'linewidth',linewidth);
else
    plot(err,avg_sum_se(1,:),linestyle{1},'color',colours(1,:),'linewidth',linewidth);
end

xticks([0 pi/180 pi/72]);
xticklabels({'0','5','10'})

xlabel('Maximum estimation error','fontname',fontname,'fontsize',fontsize);
ylabel('Sum-spectral efficiency','fontname',fontname,'fontsize',fontsize);

if K == 10
    legend(legend_prec,'fontname',fontname,'fontsize',fontsize,'location',location_2);
    legend box off;
end

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([0 pi/72]);

% figure;
%        
% set(gcf,'position',[0 0 800 600]);
% 
% if K <= M
%     plot(1:L_max,avg_sum_se(1)*ones(L_max,1),'-k','linewidth',linewidth);
%     hold on;
%     plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se(2)*ones(L_max,1),'--k','linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se(3)*ones(L_max,1),':k','linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
% else
%     plot(1:L_max,avg_sum_se(1)*ones(L_max,1),'-k','linewidth',linewidth);
%     hold on;
%     plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
% end
% 
% xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
% ylabel('Sum-spectral efficiency','fontname',fontname,'fontsize',fontsize);
% 
% if K == 10
%     legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_4);
%     legend box off;
% end
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% xlim([1 L_max]);
% 
% switch chn_type
%     case 'rayleigh'
%     case 'ur_los'
%         switch M
%             case 50
%                 switch K
%                     case 10
%                         ylim([4 13]);
%                         
%                         dim = [0.7 0.825 0.2 0.07];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.375 .275 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([7 9]);
%                         ylim([12 12.5]);
%                     case 25
%                         ylim([4 16]);
%                         
%                         dim = [0.56 0.85 0.2 0.07];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.45 .3 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([15 18]);
%                         ylim([15.5 15.9]);
%                     case 50
%                         ylim([0 18]);
%                         
%                         dim = [0.45 0.85 0.175 0.07];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.3 .375 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([20 30]);
%                         ylim([17 17.8]);
%                     case 75
%                         ylim([0 19]);
%                         
%                         dim = [0.5 0.85 0.2 0.07];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.35 .35 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([25 36]);
%                         ylim([17.5 18.55]);
%                 end
%             case 100
%                 switch K
%                     case 10
%                         ylim([5 20]);
%                         
%                         dim = [0.825 0.825 0.075 0.075];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.375 .275 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se(1)*ones(L_max,1),'-k','linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se(2)*ones(L_max,1),'--k','linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se(3)*ones(L_max,1),':k','linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([8 9]);
%                         ylim([18 19.5]);
%                     case 25
%                         ylim([5 27]);
%                         
%                         dim = [0.7 0.875 0.15 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.45 .3 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([17 22]);
%                         ylim([26 27]);
%                     case 50
%                         ylim([5 32]);
%                         
%                         dim = [0.6 0.875 0.15 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.45 .275 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([31 37]);
%                         ylim([31 32]);
%                     case 75
%                         ylim([5 35]);
%                         
%                         dim = [0.525 0.875 0.15 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.4 .265 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([35 50]);
%                         ylim([33 34.3]);
%                     case 100
%                         ylim([5 36]);
%                         
%                         dim = [0.475 0.875 0.15 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.325 .3 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([45 60]);
%                         ylim([34.2 35.6]);
%                     case 150
%                         ylim([0 37.1]);
%                         
%                         dim = [0.5 0.875 0.225 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.4 .3 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([55 76]);
%                         ylim([35 37.1]);
%                 end
%             case 200
%                 switch K
%                     case 10
%                     case 25
%                         ylim([5 43]);
%                         
%                         dim = [0.75 0.875 0.15 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.475 .3 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([20 23]);
%                         ylim([41.5 42.6]);
%                     case 50
%                         ylim([5 55]);
%                         
%                         dim = [0.7 0.875 0.1 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.45 .275 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([36 42]);
%                         ylim([52.5 54]);
%                     case 75
%                         ylim([5 60]);
%                         
%                         dim = [0.655 0.875 0.125 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.45 .275 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([52 58]);
%                         ylim([59.5 60.1]);
%                     case 100
%                         ylim([5 64]);
%                         
%                         dim = [0.625 0.875 0.125 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.45 .275 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([64 72]);
%                         ylim([62 64]);
%                     case 150
%                         ylim([5 69]);
%                         
%                         dim = [0.55 0.85 0.125 0.075];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.4 .3 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([80 95]);
%                         ylim([66 69]);
%                     case 200
%                     case 250
%                 end
%         end
%     otherwise
% end
% 
% if (savefig == 1)
%     saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
%     saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
%     saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
% end
% 
% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
% hold on;
% plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
% plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
% plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
% plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
% plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
% plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
% plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
% plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
% 
% xlabel('Sum-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_algo_2,'fontname',fontname,'fontsize',fontsize,'location',location_1);
% legend box off;
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% switch chn_type
%     case 'rayleigh'
%     case 'ur_los'
%         switch M
%             case 50
%                 switch K
%                     case 10
%                         xlim([9 13])
%                         ylim([0 1]);
%                         
%                         dim = [0.2 0.2 0.4 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([9.5 11]);
%                         ylim([0.01 0.05]);
%                     case 25
%                         xlim([13 16.5]);
%                         ylim([0 1]);
%                         
%                         dim = [0.225 0.2 0.35 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([13.5 15]);
%                         ylim([0.03 0.05]);
%                     case 50
%                         xlim([15 18.2]);
%                         ylim([0 1]);
%                         
%                         dim = [0.375 0.2 0.3 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.275 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([16 17.3]);
%                         ylim([0.03 0.05]);
%                     case 75
%                         xlim([15.5 19]);
%                         ylim([0 1]);
%                         
%                         dim = [0.425 0.2 0.275 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([17.2 18]);
%                         ylim([0.03 0.05]);
%                 end
%             case 100
%                 switch K
%                     case 10
%                         xlim([14 20])
%                         ylim([0 1]);
%                         
%                         dim = [0.3 0.2 0.2 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([15 16.4]);
%                         ylim([0.03 0.05]);
%                     case 25
%                         xlim([22 28]);
%                         ylim([0 1]);
%                         
%                         dim = [0.225 0.2 0.35 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([23 24.6]);
%                         ylim([0.03 0.05]);
%                     case 50
%                         xlim([28 33]);
%                         ylim([0 1]);
%                         
%                         dim = [0.275 0.2 0.275 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([29 30.4]);
%                         ylim([0.03 0.05]);
%                     case 75
%                         xlim([30 35]);
%                         ylim([0 1]);
%                         
%                         dim = [0.425 0.2 0.2 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([31.8 33.2]);
%                         ylim([0.03 0.05]);
%                     case 100
%                         xlim([31.5 36.3]);
%                         ylim([0 1]);
%                         
%                         dim = [0.425 0.2 0.225 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([33.3 34.7]);
%                         ylim([0.03 0.05]);
%                     case 150
%                         xlim([32 37.7]);
%                         ylim([0 1]);
%                         
%                         dim = [0.475 0.2 0.275 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([34.6 36.4]);
%                         ylim([0.03 0.05]);
%                 end
%             case 200
%                 switch K
%                     case 10
%                         xlim([22 27.5]);
%                         ylim([0 1]);
%                         
%                         dim = [0.275 0.2 0.35 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([23 25.5]);
%                         ylim([0.015 0.05]);
%                     case 25
%                         xlim([35 44]);
%                         ylim([0 1]);
%                         
%                         dim = [0.225 0.2 0.2 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .375 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([36 38]);
%                         ylim([0.015 0.05]);
%                     case 50
%                         xlim([47 55.5]);
%                         ylim([0 1]);
%                         
%                         dim = [0.225 0.175 0.225 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([48 51.5]);
%                         ylim([0.02 0.05]);
%                     case 75
%                         xlim([55 61.8]);
%                         ylim([0 1]);
%                         
%                         dim = [0.15 0.175 0.225 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([55 57.5]);
%                         ylim([0.02 0.05]);
%                     case 100
%                         xlim([59 65.5]);
%                         ylim([0 1]);
%                         
%                         dim = [0.25 0.175 0.225 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([60 61.8]);
%                         ylim([0.02 0.05]);
%                     case 150
%                         xlim([62 70]);
%                         ylim([0 1]);
%                         
%                         dim = [0.3 0.175 0.3 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .375 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([64 67]);
%                         ylim([0.02 0.05]);
%                     case 200
%                     case 250
%                 end
%         end
%     otherwise
% end
% 
% if (savefig == 1)
%     saveas(gcf,[root_save 'cdf_sum_se_star_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
%     saveas(gcf,[root_save 'cdf_sum_se_star_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
%     saveas(gcf,[root_save 'cdf_sum_se_star_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
% end