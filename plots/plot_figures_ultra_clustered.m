clear;
close all;
clc;

% Macros

MC   = 500;                                                                                                                                           % Size of the monte-carlo ensemble
N_MC = 4;

M = 50;                                                                                                                                               % Number of antennas at base station
K = 75;                                                                                                                                               % Number of users at the cell

L = [ceil(K/5) ceil(K/2)];

% M = 50  & K = [10 25 50 75]
% M = 100 & K = [10 25 50 75 100 150] 
% M = 200 & K = [10 25 50 75 100 150 200 250]

if K > M
    L_max = M;
else
    L_max = K-1;
end

snr = -5;

theta_mid   = 45;
theta_step  = 10:5:50;

bandwidth   = 20e6;
dl_ul_ratio = 0.5;

N_ALG = 3;                                                                                                                                            % Number of algorithms for perform user scheduling
N_PRE = 3;
N_STP = length(theta_step);
N_L   = length(L);

% Roots

% root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Downlink/';
root_load = 'D:\PhD\user-selection\Ultra Clustered\';
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
        load([root_load 'spectral_efficiency_all_L_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' ...
              sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC) '_' ...
              sprintf(zero_pad_2,n_mc) '.mat']);
        
        idx_ini = (n_mc - 1)*MC + 1;
        idx_end = n_mc*MC;
        
        se_all_mc(:,:,stp_idx,idx_ini:idx_end)         = bandwidth*dl_ul_ratio*se;
        se_s_L_all_mc(:,:,:,:,stp_idx,idx_ini:idx_end) = bandwidth*dl_ul_ratio*se_s_all_L; 
        
        clear se se_s_all_L;
    end
end

sum_se = reshape(sum(se_all_mc,1),N_PRE,N_STP,MC*N_MC);

for l = 1:L_max
    sum_se_s(l,:,:,:,:) = sum(se_s_L_all_mc(:,l,:,:,:,:),1);
end

avg_sum_se   = mean(sum_se,3);
avg_sum_se_s = mean(sum_se_s,5);

L_aux = repmat((1:L_max)',1,N_PRE,N_ALG,N_STP);

per_user_se   = avg_sum_se/K;
per_user_se_s = reshape(mean(avg_sum_se_s./L_aux,1),N_PRE,N_ALG,N_STP);
per_user_se_s = permute(per_user_se_s,[3 1 2]);

N_BIN = 100;

cdf_sum_se_s = cell(N_PRE,N_ALG,N_STP,2);
edg_sum_se_s = cell(N_PRE,N_ALG,N_STP,2);

for l_idx = 1:N_L
    for stp_idx = 1:N_STP
        for alg_idx = 1:N_ALG
            for pre_idx = 1:N_PRE
                [cdf_sum_se_s{pre_idx,alg_idx,stp_idx,l_idx},edg_sum_se_s{pre_idx,alg_idx,stp_idx,l_idx}] = histcounts(sum_se_s(L(l_idx),pre_idx,alg_idx,stp_idx,:),N_BIN,'normalization','cdf');
            end
        end
    end
end

% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 30;

marker    = {'o','s','^'};
linestyle = {'-','--',':'};

savefig = 1;
plotse  = 0;

if M == 50
    OM = 1e-6;
elseif M == 100
    OM = 1e-6;
end

% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_alg_plus_pre     = {'SOS','CBS','ICIBS','MRT','ZF','MMSE'};
legend_alg_plus_pre_2   = {'SOS','CBS','ICIBS','MRT','MMSE'};
legend_section_plus_alg = {'$\Delta\theta = 1^{\circ}$','$\Delta\theta = 2.5^{\circ}$','$\Delta\theta = 5^{\circ}$','SOS','CBS','ICIBS'};

um = {'(kbps)','(Mbps)','(Gbps)'};

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
plot(theta_step/10,OM*per_user_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
hold on;
plot(theta_step/10,OM*per_user_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
plot(theta_step/10,OM*per_user_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
% Plots for precoding algorithm legends
plot(theta_step/10,OM*per_user_se_s(:,1,1),'-' ,'color',colours(8,:),'linewidth',linewidth);
plot(theta_step/10,OM*per_user_se_s(:,2,1),'--','color',colours(8,:),'linewidth',linewidth);
plot(theta_step/10,OM*per_user_se_s(:,3,1),':' ,'color',colours(8,:),'linewidth',linewidth);
% Plots for results
plot(theta_step/10,OM*per_user_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
plot(theta_step/10,OM*per_user_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
plot(theta_step/10,OM*per_user_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
plot(theta_step/10,OM*per_user_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
plot(theta_step/10,OM*per_user_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
plot(theta_step/10,OM*per_user_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
plot(theta_step/10,OM*per_user_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
plot(theta_step/10,OM*per_user_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
plot(theta_step/10,OM*per_user_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);

xlabel('$\Delta\theta$ (in degrees)','fontname',fontname,'fontsize',fontsize,'interpreter','latex');

if plotse == 1
    ylabel('Per-user SE','fontname',fontname,'fontsize',fontsize);
else
    ylabel(['Per-user throughput ' um{abs(log10(OM))/3}],'fontname',fontname,'fontsize',fontsize);
    % ylabel('Per-user throughput','fontname',fontname,'fontsize',fontsize);
end

legend(legend_alg_plus_pre,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',2);
legend box off;
    
% xticks([1 2 3 4 5]);
% xticklabels({'1^{\circ}','2^{\circ}','3^{\circ}','4^{\circ}','5^{\circ}'})                       

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 5]);

if savefig == 1
    if plotse == 1
        saveas(gcf,[root_save 'per_user_se_dtheta_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
        saveas(gcf,[root_save 'per_user_se_dtheta_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
        saveas(gcf,[root_save 'per_user_se_dtheta_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
    else
        saveas(gcf,[root_save 'per_user_throughput_dtheta_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
        saveas(gcf,[root_save 'per_user_throughput_dtheta_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
        saveas(gcf,[root_save 'per_user_throughput_dtheta_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');     
    end
end

for stp_idx = 2
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    if K <= M
        % Plots for user selection algorithm legends
        plot(1:K,OM*[avg_sum_se_s(:,1,1,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:K,OM*[avg_sum_se_s(:,1,2,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,1,3,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(3,:),'linewidth',linewidth);
        % Plots for precoding algorithm legends
        plot(1:K,OM*[avg_sum_se_s(:,1,1,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(8,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,2,1,stp_idx); avg_sum_se(2,stp_idx)],'--','color',colours(8,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,3,1,stp_idx); avg_sum_se(3,stp_idx)],':' ,'color',colours(8,:),'linewidth',linewidth);
        % Plots for results
        plot(1:K,OM*[avg_sum_se_s(:,1,1,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,1,2,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,1,3,stp_idx); avg_sum_se(1,stp_idx)],'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,2,1,stp_idx); avg_sum_se(2,stp_idx)],'--','color',colours(1,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,2,2,stp_idx); avg_sum_se(2,stp_idx)],'--','color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,2,3,stp_idx); avg_sum_se(2,stp_idx)],'--','color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,3,1,stp_idx); avg_sum_se(3,stp_idx)],':' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,3,2,stp_idx); avg_sum_se(3,stp_idx)],':' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,3,3,stp_idx); avg_sum_se(3,stp_idx)],':' ,'color',colours(3,:),'linewidth',linewidth);
    else
        % Plots for user selection algorithms legend
        plot(1:L_max,OM*avg_sum_se_s(:,1,1,stp_idx),'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:L_max,OM*avg_sum_se_s(:,1,2,stp_idx),'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_se_s(:,1,3,stp_idx),'-' ,'color',colours(3,:),'linewidth',linewidth);
        % Plots for precoding algorithms legend
        plot(1:L_max,OM*avg_sum_se_s(:,1,1,stp_idx),'-' ,'color',colours(8,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_se_s(:,2,1,stp_idx),'--','color',colours(8,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_se_s(:,3,1,stp_idx),':' ,'color',colours(8,:),'linewidth',linewidth);
        % Plots for result
        plot(1:L_max,OM*avg_sum_se_s(:,1,1,stp_idx),'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_se_s(:,1,2,stp_idx),'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_se_s(:,1,3,stp_idx),'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_se_s(:,2,1,stp_idx),'--','color',colours(1,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_se_s(:,2,2,stp_idx),'--','color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_se_s(:,2,3,stp_idx),'--','color',colours(3,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_se_s(:,3,1,stp_idx),':' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_se_s(:,3,2,stp_idx),':' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_se_s(:,3,3,stp_idx),':' ,'color',colours(3,:),'linewidth',linewidth);
    end
    
    xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
    
    if plotse == 1
        ylabel('Sum-spectral efficiency','fontname',fontname,'fontsize',fontsize);
    else
        ylabel(['Throughput ' um{abs(log10(OM))/3}],'fontname',fontname,'fontsize',fontsize);
    end
    
    legend(legend_alg_plus_pre,'fontname',fontname,'fontsize',fontsize,'location',location_2,'numcolumns',2);
    legend box off;
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    if K <= M
        xlim([1 K]);
    else
        xlim([1 L_max]);
    end
    
    ylim([0 80]);
    
    if savefig == 1
        if plotse == 1
            saveas(gcf,[root_save 'sum_se_all_L_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
            saveas(gcf,[root_save 'sum_se_all_L_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
            saveas(gcf,[root_save 'sum_se_all_L_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
        else
            saveas(gcf,[root_save 'throughput_all_L_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
            saveas(gcf,[root_save 'throughput_all_L_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
            saveas(gcf,[root_save 'throughput_all_L_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
        end
    end
    
    for l_idx = 1:N_L
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        % Plots for user selection algorihtm legends
        plot(OM*edg_sum_se_s{1,1,stp_idx,l_idx},[cdf_sum_se_s{1,1,stp_idx,l_idx} 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(OM*edg_sum_se_s{1,2,stp_idx,l_idx},[cdf_sum_se_s{1,2,stp_idx,l_idx} 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(OM*edg_sum_se_s{1,3,stp_idx,l_idx},[cdf_sum_se_s{1,3,stp_idx,l_idx} 1],'-' ,'color',colours(3,:),'linewidth',linewidth);
        % Plots for precoding algorithm legends
        plot(OM*edg_sum_se_s{1,1,stp_idx,l_idx},[cdf_sum_se_s{1,1,stp_idx,l_idx} 1],'-' ,'color',colours(8,:),'linewidth',linewidth);
        % plot(OM*edg_sum_se_s{2,1,stp_idx,l_idx},[cdf_sum_se_s{2,1,stp_idx,l_idx} 1],'--','color',colours(8,:),'linewidth',linewidth);
        plot(OM*edg_sum_se_s{3,1,stp_idx,l_idx},[cdf_sum_se_s{3,1,stp_idx,l_idx} 1],':' ,'color',colours(8,:),'linewidth',linewidth);
        % Plots for results
        plot(OM*edg_sum_se_s{1,1,stp_idx,l_idx},[cdf_sum_se_s{1,1,stp_idx,l_idx} 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(OM*edg_sum_se_s{1,2,stp_idx,l_idx},[cdf_sum_se_s{1,2,stp_idx,l_idx} 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(OM*edg_sum_se_s{1,3,stp_idx,l_idx},[cdf_sum_se_s{1,3,stp_idx,l_idx} 1],'-' ,'color',colours(3,:),'linewidth',linewidth);
        % plot(OM*edg_sum_se_s{2,1,stp_idx,l_idx},[cdf_sum_se_s{2,1,stp_idx,l_idx} 1],'--','color',colours(1,:),'linewidth',linewidth);
        % plot(OM*edg_sum_se_s{2,2,stp_idx,l_idx},[cdf_sum_se_s{2,2,stp_idx,l_idx} 1],'--','color',colours(2,:),'linewidth',linewidth);
        % plot(OM*edg_sum_se_s{2,3,stp_idx,l_idx},[cdf_sum_se_s{2,3,stp_idx,l_idx} 1],'--','color',colours(3,:),'linewidth',linewidth);
        plot(OM*edg_sum_se_s{3,1,stp_idx,l_idx},[cdf_sum_se_s{3,1,stp_idx,l_idx} 1],':' ,'color',colours(1,:),'linewidth',linewidth);
        plot(OM*edg_sum_se_s{3,2,stp_idx,l_idx},[cdf_sum_se_s{3,2,stp_idx,l_idx} 1],':' ,'color',colours(2,:),'linewidth',linewidth);
        plot(OM*edg_sum_se_s{3,3,stp_idx,l_idx},[cdf_sum_se_s{3,3,stp_idx,l_idx} 1],':' ,'color',colours(3,:),'linewidth',linewidth);
        
        if plotse == 1
            xlabel('Sum-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        else
            xlabel(['Throughput ' um{abs(log10(OM))/3}],'fontname',fontname,'fontsize',fontsize);
        end
        
        ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
        
        %if theta_mid == 0 && stp_idx == 2 && l_idx == 1
            legend(legend_alg_plus_pre_2,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',1);
            legend box off;
        %end
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        if stp_idx == 2
            if l_idx == 1
                xlim([0 30]);
            elseif l_idx == 2
                xlim([0 25]);
            end
        end
        
        % xlim([0 45]);
        ylim([0 1]);
        
        if savefig == 1
            if plotse == 1
                saveas(gcf,[root_save 'cdf_sum_se_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L(l_idx)) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
                saveas(gcf,[root_save 'cdf_sum_se_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L(l_idx)) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
                saveas(gcf,[root_save 'cdf_sum_se_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L(l_idx)) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
            else
                saveas(gcf,[root_save 'cdf_throughput_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L(l_idx)) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
                saveas(gcf,[root_save 'cdf_throughput_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L(l_idx)) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
                saveas(gcf,[root_save 'cdf_throughput_clustered_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L(l_idx)) '_theta_mid_' sprintf(zero_pad_2,theta_mid) '_theta_step_' num2str(theta_step(stp_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
            end
        end
    end
end