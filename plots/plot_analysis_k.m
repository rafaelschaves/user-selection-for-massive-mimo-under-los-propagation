clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = [64 128 256];                                                          % Number of antennas at base station
K = [1 2 4 8 16 32 64 128 256 512 1024];                                   % Number of mobile users

snr = 10;

M_SIZ = length(M);                                                         % Size of the antenna set
K_SIZ = length(K);                                                         % Size of the user set
N_CHN = 2;                                                                 % Number of channel models simulated

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Users/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Users/';

chn_type = {'ur_los','rayleigh'};

% Loading data

throughput_u = cell(K_SIZ,M_SIZ,N_CHN);                                    % Uplink throughput
throughput_d = cell(K_SIZ,M_SIZ,N_CHN);

for chn_idx = 1:N_CHN
    for m = 1:M_SIZ
        for k = 1:K_SIZ
            load([root_load 'throughput_outdoors_pedestrian_mf_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K(k)) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC) '.mat']);
            
            throughput_u{k,m,chn_idx} = thrput_u;
            throughput_d{k,m,chn_idx} = thrput_d;
        end
    end
end

% Post processing - Calculating the CDF

avg_sum_throughput_u = zeros(K_SIZ,M_SIZ,N_CHN);
avg_sum_throughput_d = zeros(K_SIZ,M_SIZ,N_CHN);

avg_throughput_per_user_u = zeros(K_SIZ,M_SIZ,N_CHN);
avg_throughput_per_user_d = zeros(K_SIZ,M_SIZ,N_CHN);

for chn_idx = 1:N_CHN
    for m = 1:M_SIZ
        for k = 1:K_SIZ
            avg_sum_throughput_u(k,m,chn_idx) = mean(sum(throughput_u{k,m,chn_idx},1));
            avg_sum_throughput_d(k,m,chn_idx) = mean(sum(throughput_d{k,m,chn_idx},1));
            
            avg_throughput_per_user_u(k,m,chn_idx) = mean(throughput_u{k,m,chn_idx}(:));
            avg_throughput_per_user_d(k,m,chn_idx) = mean(throughput_d{k,m,chn_idx}(:));
        end
    end
end

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

savefig = 0;

legend_M = {'M = 64','M = 128','M = 256'};

location = 'northwest';

colours = get(gca,'colororder');
close;

figure;

set(gcf,'position',[0 0 800 600]);

plot(K,1e-9*avg_sum_throughput_d(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot(K,1e-9*avg_sum_throughput_d(:,2,1),'-','color',colours(2,:),'linewidth',linewidth);
plot(K,1e-9*avg_sum_throughput_d(:,3,1),'-','color',colours(3,:),'linewidth',linewidth);
plot(K,1e-9*avg_sum_throughput_d(:,1,2),'--','color',colours(1,:),'linewidth',linewidth);
plot(K,1e-9*avg_sum_throughput_d(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
plot(K,1e-9*avg_sum_throughput_d(:,3,2),'--','color',colours(3,:),'linewidth',linewidth);

xlabel('Number of users','fontname',fontname,'fontsize',fontsize);
ylabel('Average sum-throughput (Gbps)','fontname',fontname,'fontsize',fontsize);

legend(legend_M,'fontname',fontname,'fontsize',fontsize,'location',location);
legend box off;

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([1 1024]);

if (savefig == 1)
    saveas(gcf,[root_save 'out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K(k)) '_SNR_' num2str(snr) '_dB'],'fig');
    saveas(gcf,[root_save 'out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K(k)) '_SNR_' num2str(snr) '_dB'],'png');
    saveas(gcf,[root_save 'out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K(k)) '_SNR_' num2str(snr) '_dB'],'epsc2');
end

figure;

set(gcf,'position',[0 0 800 600]);

plot(K,1e-6*avg_throughput_per_user_d(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot(K,1e-6*avg_throughput_per_user_d(:,2,1),'-','color',colours(2,:),'linewidth',linewidth);
plot(K,1e-6*avg_throughput_per_user_d(:,3,1),'-','color',colours(3,:),'linewidth',linewidth);
plot(K,1e-6*avg_throughput_per_user_d(:,1,2),'--','color',colours(1,:),'linewidth',linewidth);
plot(K,1e-6*avg_throughput_per_user_d(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
plot(K,1e-6*avg_throughput_per_user_d(:,3,2),'--','color',colours(3,:),'linewidth',linewidth);

xlabel('Number of users','fontname',fontname,'fontsize',fontsize);
ylabel('Average throughput per user (Mbps)','fontname',fontname,'fontsize',fontsize);

legend(legend_M,'fontname',fontname,'fontsize',fontsize,'location',location);
legend box off;

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([1 1024]);

if (savefig == 1)
    saveas(gcf,[root_save 'out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K(k)) '_SNR_' num2str(snr) '_dB'],'fig');
    saveas(gcf,[root_save 'out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K(k)) '_SNR_' num2str(snr) '_dB'],'png');
    saveas(gcf,[root_save 'out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K(k)) '_SNR_' num2str(snr) '_dB'],'epsc2');
end