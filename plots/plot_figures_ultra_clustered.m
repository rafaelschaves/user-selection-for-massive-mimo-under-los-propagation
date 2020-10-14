clear;
% close all;
clc;

% Macros

MC   = 500;                                                              % Size of the monte-carlo ensemble
N_MC = 1;

M = 50;                                                                  % Number of antennas at base station
K = 75;                                                                   % Number of users at the cell 

% M = 50  & K = [10 25 50 75]
% M = 100 & K = [10 25 50 75 100 150] 
% M = 200 & K = [10 25 50 75 100 150 200 250]

if K > M
    L_max = M;
else
    L_max = K-1;
end

snr = -5;

theta_step = 10:5:50;

N_ALG = 3;                                                                 % Number of algorithms for perform user scheduling
N_PRE = 3;
N_STP = length(theta_step);

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Downlink/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Selection/Downlink/';

zero_pad_1 = '%03d';
zero_pad_2 = '%02d';

chn_type = 'ur_los';

% Loading data

se_all_mc     = zeros(K,N_PRE,N_STP,MC*N_MC);
se_s_L_all_mc = zeros(L_max,L_max,N_PRE,N_ALG,N_STP,MC*N_MC);

se_users_s    = zeros(L_max*(L_max+1)*MC*N_MC/2,N_PRE,N_ALG,N_STP);
sum_se_s      = zeros(L_max,N_PRE,N_ALG,N_STP,MC*N_MC);

for stp_idx = 1:N_STP
    for n_mc = 1:N_MC
        load([root_load 'spectral_efficiency_all_L_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_00_theta_step_' ...
              num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC) '_' sprintf(zero_pad_2,n_mc) '.mat']);
        
        idx_ini = (n_mc - 1)*MC + 1;
        idx_end = n_mc*MC;
        
        se_all_mc(:,:,stp_idx,idx_ini:idx_end)         = se;
        se_s_L_all_mc(:,:,:,:,stp_idx,idx_ini:idx_end) = se_s_all_L; 
        
        clear se se_s_all_L;
    end
end

for stp_idx = 1:N_STP
    for alg_idx = 1:N_ALG
        for pre_idx = 1:N_PRE
            se_users_s(:,pre_idx,alg_idx,stp_idx) = nonzeros(se_s_L_all_mc(:,:,pre_idx,alg_idx,stp_idx,:));
        end
    end
end

sum_se = reshape(sum(se_all_mc,1),N_PRE,N_STP,MC*N_MC);

for L = 1:L_max
    sum_se_s(L,:,:,:,:) = sum(se_s_L_all_mc(:,L,:,:,:,:),1);
end

avg_sum_se   = mean(sum_se,3);
avg_sum_se_s = mean(sum_se_s,5);

L_aux = repmat((1:L_max)',1,N_PRE,N_ALG,N_STP);

per_user_se    = avg_sum_se/K;
per_user_se_s = reshape(mean(avg_sum_se_s./L_max,1),N_PRE,N_ALG,N_STP);
per_user_se_s = permute(per_user_se_s,[3 1 2]);

N_BIN = 100;
BIN_EDG = linspace(0,1,N_BIN);

cdf_per_user_se = cell(N_PRE,N_ALG,N_STP);
edg_per_user_se = cell(N_PRE,N_ALG,N_STP);

% cdf_sum_se = cell(N_ALG+1,N_PRE);
% edg_sum_se = cell(N_ALG+1,N_PRE);
% 
% for n_pre = 1:N_PRE
%     [cdf_sum_se{1,n_pre},edg_sum_se{1,n_pre}] = histcounts(sum_se(n_pre,:),N_BIN,'normalization','cdf');
% 
%     for n_alg = 1:N_ALG
%         [cdf_sum_se{n_alg+1,n_pre},edg_sum_se{n_alg+1,n_pre}] = histcounts(sum_se_s_star(:,n_pre,n_alg),N_BIN,'normalization','cdf');
%     end
% end

for stp_idx = 1:N_STP
    for alg_idx = 1:N_ALG
        for pre_idx = 1:N_PRE
            [cdf_per_user_se{pre_idx,alg_idx,stp_idx},edg_per_user_se{pre_idx,alg_idx,stp_idx}] = histcounts(se_users_s(:,pre_idx,alg_idx,stp_idx),N_BIN,'binedges',BIN_EDG,'normalization','cdf');
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

% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_alg_plus_pre     = {'SOS','CBS','ICIBS','MRT','ZF','MMSE'};
legend_section_plus_alg = {'$\Delta\theta = 1^{\circ}$','$\Delta\theta = 2.5^{\circ}$','$\Delta\theta = 5^{\circ}$','SOS','CBS','ICIBS'};

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
           0.6350 0.0780 0.1840;
           0.0000 0.0000 0.0000];


figure;

set(gcf,'position',[0 0 800 600]);

% Plots for user selection algorihtm legends
plot(theta_step/10,per_user_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
hold on;
plot(theta_step/10,per_user_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
plot(theta_step/10,per_user_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
% Plots for precoding algorithm legends
plot(theta_step/10,per_user_se_s(:,1,1),'-' ,'color',colours(8,:),'linewidth',linewidth);
plot(theta_step/10,per_user_se_s(:,2,1),'--','color',colours(8,:),'linewidth',linewidth);
plot(theta_step/10,per_user_se_s(:,3,1),':' ,'color',colours(8,:),'linewidth',linewidth);
% Plots for results
plot(theta_step/10,per_user_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
plot(theta_step/10,per_user_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
plot(theta_step/10,per_user_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
plot(theta_step/10,per_user_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
plot(theta_step/10,per_user_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
plot(theta_step/10,per_user_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
plot(theta_step/10,per_user_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
plot(theta_step/10,per_user_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
plot(theta_step/10,per_user_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);

xlabel('$\Delta\theta$ (in degrees)','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
ylabel('Average per-user SE','fontname',fontname,'fontsize',fontsize);

legend(legend_alg_plus_pre,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',2);
legend box off;
    
set(gca,'fontname',fontname,'fontsize',fontsize);

% xticks([7 10]);
% xticklabels({'0','1.25','2.5'})                       

 for stp_idx = 1:N_STP
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    if K <= M
        % Plots for user selection algorithm legends
        plot(1:K,[avg_sum_se_s(:,1,1,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:K,[avg_sum_se_s(:,1,2,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,1,3,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(3,:),'linewidth',linewidth);
        % Plots for precoding algorithm legends
        plot(1:K,[avg_sum_se_s(:,1,1,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(8,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,2,1,stp_idx); avg_sum_se(2,stp_idx)],'--','color',colours(8,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,3,1,stp_idx); avg_sum_se(3,stp_idx)],':' ,'color',colours(8,:),'linewidth',linewidth);
        % Plots for results
        plot(1:K,[avg_sum_se_s(:,1,1,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,1,2,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,1,3,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,2,1,stp_idx); avg_sum_se(2,stp_idx)],'--','color',colours(1,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,2,2,stp_idx); avg_sum_se(2,stp_idx)],'--','color',colours(2,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,2,3,stp_idx); avg_sum_se(2,stp_idx)],'--','color',colours(3,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,3,1,stp_idx); avg_sum_se(3,stp_idx)],':' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,3,2,stp_idx); avg_sum_se(3,stp_idx)],':' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,3,3,stp_idx); avg_sum_se(3,stp_idx)],':' ,'color',colours(3,:),'linewidth',linewidth);
    else
        % Plots for user selection algorithms legend
        plot(1:L_max,avg_sum_se_s(:,1,1,stp_idx),'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:L_max,avg_sum_se_s(:,1,2,stp_idx),'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,1,3,stp_idx),'-' ,'color',colours(3,:),'linewidth',linewidth);
        % Plots for precoding algorithms legend
        plot(1:L_max,avg_sum_se_s(:,1,1,stp_idx),'-' ,'color',colours(8,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,2,1,stp_idx),'--','color',colours(8,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,3,1,stp_idx),':' ,'color',colours(8,:),'linewidth',linewidth); 
        % Plots for result
        plot(1:L_max,avg_sum_se_s(:,1,1,stp_idx),'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,1,2,stp_idx),'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,1,3,stp_idx),'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,2,1,stp_idx),'--','color',colours(1,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,2,2,stp_idx),'--','color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,2,3,stp_idx),'--','color',colours(3,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,3,1,stp_idx),':' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,3,2,stp_idx),':' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,3,3,stp_idx),':' ,'color',colours(3,:),'linewidth',linewidth);
    end
    
    xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
    ylabel('Sum-spectral efficiency','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_alg_plus_pre,'fontname',fontname,'fontsize',fontsize,'location',location_2,'numcolumns',2);
    legend box off;
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    if K <= M
        xlim([1 K]);
    else
        xlim([1 L_max]);
    end
end

for pre_idx = 1:N_PRE
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    if K <= M
        % Plots for cell section legends
        plot(1:K,[avg_sum_se_s(:,pre_idx,1,1); avg_sum_se(pre_idx,1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:K,[avg_sum_se_s(:,pre_idx,1,4); avg_sum_se(pre_idx,4)],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,pre_idx,1,9); avg_sum_se(pre_idx,9)],'-' ,'color',colours(3,:),'linewidth',linewidth);
        % Plots for user selection algorithm legends
        plot(1:K,[avg_sum_se_s(:,pre_idx,1,1); avg_sum_se(pre_idx,1)],'-' ,'color',colours(8,:),'linewidth',linewidth);     
        plot(1:K,[avg_sum_se_s(:,pre_idx,2,1); avg_sum_se(pre_idx,1)],'--','color',colours(8,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,pre_idx,3,1); avg_sum_se(pre_idx,1)],':' ,'color',colours(8,:),'linewidth',linewidth);
        % Plots for results
        plot(1:K,[avg_sum_se_s(:,pre_idx,1,1); avg_sum_se(pre_idx,1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,pre_idx,2,1); avg_sum_se(pre_idx,1)],'--','color',colours(1,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,pre_idx,1,1); avg_sum_se(pre_idx,1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,pre_idx,1,4); avg_sum_se(pre_idx,4)],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,pre_idx,2,4); avg_sum_se(pre_idx,4)],'--','color',colours(2,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,pre_idx,3,4); avg_sum_se(pre_idx,4)],':' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,pre_idx,1,9); avg_sum_se(pre_idx,9)],'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,pre_idx,2,9); avg_sum_se(pre_idx,9)],'--','color',colours(3,:),'linewidth',linewidth);
        plot(1:K,[avg_sum_se_s(:,pre_idx,3,9); avg_sum_se(pre_idx,9)],':' ,'color',colours(3,:),'linewidth',linewidth);
    else
        % Plots for cell section legends
        plot(1:L_max,avg_sum_se_s(:,pre_idx,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:L_max,avg_sum_se_s(:,pre_idx,1,4),'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,pre_idx,1,9),'-' ,'color',colours(3,:),'linewidth',linewidth);
        % Plots for user selection algorithm legends
        plot(1:L_max,avg_sum_se_s(:,pre_idx,1,1),'-' ,'color',colours(8,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,pre_idx,2,1),'--','color',colours(8,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,pre_idx,3,1),':' ,'color',colours(8,:),'linewidth',linewidth);
        % Plots for resuls
        plot(1:L_max,avg_sum_se_s(:,pre_idx,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,pre_idx,2,1),'--','color',colours(1,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,pre_idx,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,pre_idx,1,4),'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,pre_idx,2,4),'--','color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,pre_idx,3,4),':' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,pre_idx,1,9),'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,pre_idx,2,9),'--','color',colours(3,:),'linewidth',linewidth);
        plot(1:L_max,avg_sum_se_s(:,pre_idx,3,9),':' ,'color',colours(3,:),'linewidth',linewidth);
    end
    
    xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
    ylabel('Sum-spectral efficiency','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_section_plus_alg,'fontname',fontname,'fontsize',fontsize,'location',location_2,'numcolumns',2,'interpreter','latex');
    legend box off;
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    ylim([0 10]);
    
    if K <= M
        xlim([1 K]);
    else
        xlim([1 L_max]);
    end
end

for stp_idx = 1:N_STP
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot(edg_per_user_se{1,1,stp_idx},[cdf_per_user_se{1,1,stp_idx} 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(edg_per_user_se{1,2,stp_idx},[cdf_per_user_se{1,2,stp_idx} 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
    plot(edg_per_user_se{1,3,stp_idx},[cdf_per_user_se{1,3,stp_idx} 1],'-' ,'color',colours(3,:),'linewidth',linewidth);
    plot(edg_per_user_se{2,1,stp_idx},[cdf_per_user_se{2,1,stp_idx} 1],'--','color',colours(1,:),'linewidth',linewidth);
    plot(edg_per_user_se{2,2,stp_idx},[cdf_per_user_se{2,2,stp_idx} 1],'--','color',colours(2,:),'linewidth',linewidth);
    plot(edg_per_user_se{2,3,stp_idx},[cdf_per_user_se{2,3,stp_idx} 1],'--','color',colours(3,:),'linewidth',linewidth);
    plot(edg_per_user_se{3,1,stp_idx},[cdf_per_user_se{3,1,stp_idx} 1],':' ,'color',colours(1,:),'linewidth',linewidth);
    plot(edg_per_user_se{3,2,stp_idx},[cdf_per_user_se{3,2,stp_idx} 1],':' ,'color',colours(2,:),'linewidth',linewidth);
    plot(edg_per_user_se{3,3,stp_idx},[cdf_per_user_se{3,3,stp_idx} 1],':' ,'color',colours(3,:),'linewidth',linewidth);
    
    xlabel('per-user SE (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
    
%     legend(legend_alg_plus_pre,'fontname',fontname,'fontsize',fontsize,'location',location_1);
%     legend box off;
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    xlim([0 1]);
    ylim([0 1]);
end

% switch chn_type
%     case 'rayleigh'
%     case 'ur_los'
%         switch M
%             case 50
%                 ylim([0 20]);
%                 
%                 switch K
%                     case 10                        
%                         dim = [0.7 0.6 0.2 0.07];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.24 .585 .3 .3]);
%                         box on;
%                             
%                         plot(1:K,[avg_sum_se_s(:,1,1); avg_sum_se(1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:K,[avg_sum_se_s(:,1,2); avg_sum_se(1)],'-' ,'color',colours(2,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,1,3); avg_sum_se(1)],'-' ,'color',colours(3,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,2,1); avg_sum_se(2)],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,2,2); avg_sum_se(2)],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,2,3); avg_sum_se(2)],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,3,1); avg_sum_se(3)],':' ,'color',colours(1,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,3,2); avg_sum_se(3)],':' ,'color',colours(2,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,3,3); avg_sum_se(3)],':' ,'color',colours(3,:),'linewidth',linewidth);
%                         
%                         % xticks([7 10]);
%                         % xticklabels({'0','1.25','2.5'})
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([7 10]);
%                         ylim([12 12.5]);
%                     case 25
%                         dim = [0.56 0.725 0.2 0.07];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.35 .275 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([15 20]);
%                         ylim([15.5 16]);
%                     case 50
%                         dim = [0.425 0.8 0.175 0.07];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.275 .275 .3 .3]);
%                         box on;
%                         
%                         plot(1:L_max,avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:L_max,avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
%                         plot(1:L_max,avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([20 30]);
%                         ylim([17 18]);
%                     case 75
%                         dim = [0.5 0.825 0.2 0.07];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.3 .275 .3 .3]);
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
%                 ylim([0 40]);
%                 
%                 switch K
%                     case 10                        
%                         dim = [0.75 0.49 0.15 0.075];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .575 .3 .3]);
%                         box on;
%                         
%                         plot(1:K,[avg_sum_se_s(:,1,1); avg_sum_se(1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(1:K,[avg_sum_se_s(:,1,2); avg_sum_se(1)],'-' ,'color',colours(2,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,1,3); avg_sum_se(1)],'-' ,'color',colours(3,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,2,1); avg_sum_se(2)],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,2,2); avg_sum_se(2)],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,2,3); avg_sum_se(2)],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,3,1); avg_sum_se(3)],':' ,'color',colours(1,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,3,2); avg_sum_se(3)],':' ,'color',colours(2,:),'linewidth',linewidth);
%                         plot(1:K,[avg_sum_se_s(:,3,3); avg_sum_se(3)],':' ,'color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([7 10]);
%                         ylim([18 19.5]);
%                     case 25                 
%                         dim = [0.7 0.65 0.15 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.55 .275 .3 .3]);
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
%                         xlim([18 22]);
%                         ylim([26 27]);
%                     case 50   
%                         dim = [0.6 0.75 0.15 0.05];
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
%                         dim = [0.525 0.775 0.15 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.275 .275 .3 .3]);
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
%                         dim = [0.45 0.8 0.15 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.325 .275 .3 .3]);
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
%                         dim = [0.5 0.825 0.225 0.05];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.3 .3 .3 .3]);
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

if (savefig == 1)
    saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
    saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
    saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
end

% figure;

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