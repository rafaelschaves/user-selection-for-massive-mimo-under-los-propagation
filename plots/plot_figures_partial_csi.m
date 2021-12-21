clear;
close all;
clc;

% Macros

MC = 5;                                                              % Size of the monte-carlo ensemble
MC_ERR = 10;

M = 50;                                                                   % Number of antennas at base station
K = 75;                                                                   % Number of users at the cell 

if K > M
    L_max = M;
else
    L_max = K-1;
end

snr = -5;

bandwidth   = 20e6;
dl_ul_ratio = 0.5;

var_err = [0 1e-3 1e-2 1e-1 1];

N_ALG = 3;                                                                 % Number of algorithms for perform user scheduling
N_PRE = 3;
N_ERR = length(var_err);

% Roots

root_load = 'G:\My Drive\UFRJ\PhD\Codes\user-scheduling-massive-mimo\Results\Selection\Downlink\';
root_save = 'G:\My Drive\UFRJ\PhD\Codes\user-scheduling-massive-mimo\Figures\Selection\Downlink\';

zero_pad_1 = '%03d';
zero_pad_2 = '%02d';

% Loading data

sum_se_s = zeros(L_max,N_PRE,N_ALG,N_ERR,MC_ERR,MC);

load([root_load 'se_all_L_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC) '.mat']);
    
sum_se = reshape(sum(se,1),N_PRE,N_ERR,MC_ERR,MC);

for l = 1:L_max
    sum_se_s(l,:,:,:,:,:) = sum(se_s_all_L(1:l,l,:,:,:,:,:),1);
end

avg_sum_thrgpt_eps   = bandwidth*dl_ul_ratio*mean(sum_se,3);
avg_sum_thrgpt_s_eps = bandwidth*dl_ul_ratio*mean(sum_se_s,5);

avg_sum_thrgpt_eps   = reshape(avg_sum_thrgpt_eps,N_PRE,N_ERR,MC);
avg_sum_thrgpt_s_eps = reshape(avg_sum_thrgpt_s_eps,L_max,N_PRE,N_ALG,N_ERR,MC);

[max_sum_thrgpt_s,L_star] = max(avg_sum_thrgpt_s_eps,[],1);

max_sum_thrgpt_s = reshape(max_sum_thrgpt_s,N_PRE,N_ALG,N_ERR,MC);
L_star           = reshape(L_star,N_PRE,N_ALG,N_ERR,MC);

for mc = 1:MC
    for n_err = 1:N_ERR
        for n_alg = 1:N_ALG
            for n_pre = 1:N_PRE
                sum_thrgpt_s_star(:,n_pre,n_alg,n_err) = avg_sum_thrgpt_s_eps(L_star(n_alg,n_pre),n_pre,n_alg,n_err,:);
            end
        end
    end
end

avg_sum_thrgpt   = mean(avg_sum_thrgpt_eps,3);
avg_sum_thrgpt_s = mean(avg_sum_thrgpt_s_eps,5);

avg_max_sum_thrgpt_s = mean(max_sum_thrgpt_s,4);
avg_L_star = mean(L_star,4);

% N_BIN = 100;
% 
% cdf_thrgpt_se = cell(N_ALG+1,N_PRE,N_ERR);
% edg_thrgpt_se = cell(N_ALG+1,N_PRE,N_ERR);
% 
% for n_err = 1:N_ERR
%     for n_pre = 1:N_PRE
%         [cdf_thrgpt_se{1,n_pre},edg_thrgpt_se{1,n_pre}] = histcounts(sum_se(n_pre,n_err,:),N_BIN,'normalization','cdf');
%         
%         for n_alg = 1:N_ALG
%             [cdf_thrgpt_se{n_alg+1,n_pre,n_err},edg_thrgpt_se{n_alg+1,n_pre,n_err}] = histcounts(sum_thrgpt_s_star(:,n_pre,n_alg,n_err),N_BIN,'normalization','cdf');
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

if M == 50
    OM = 1e-6;
elseif M == 100
    OM = 1e-6;
end

% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_pre = {'MRT','ZF','MMSE'};
legend_alg = {'SOS','CBS','ICIBS'};
legend_err = {'$\sigma_{\varepsilon}^{2} = 0$', ...
              '$\sigma_{\varepsilon}^{2} = 10^{-3}$', ...
              '$\sigma_{\varepsilon}^{2} = 10^{-2}$', ...
              '$\sigma_{\varepsilon}^{2} = 10^{-1}$', ...
              '$\sigma_{\varepsilon}^{2} = 1$'};

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
  
for alg = 1:N_ALG
    for pre = 1:N_PRE
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        plot(1:L_max,OM*avg_sum_thrgpt_s(:,pre,alg,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:L_max,OM*avg_sum_thrgpt_s(:,pre,alg,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_thrgpt_s(:,pre,alg,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_thrgpt_s(:,pre,alg,4),'-' ,'color',colours(4,:),'linewidth',linewidth);
        plot(1:L_max,OM*avg_sum_thrgpt_s(:,pre,alg,5),'-' ,'color',colours(5,:),'linewidth',linewidth);
        
        xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
        ylabel(['Throughput ' um{abs(log10(OM))/3}],'fontname',fontname,'fontsize',fontsize);
        
        if pre == 1
            legend(legend_err,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_4,'numcolumns',2);
            legend box off;
        end

        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        if K <= M
            xlim([1 K]);
        else
            xlim([1 L_max]);
        end
        
        if savefig == 1
            saveas(gcf,[root_save 'throughput_pcsi_all_L_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_' legend_pre{pre} '_' legend_alg{alg} '_MC_' num2str(MC)],'fig');
            saveas(gcf,[root_save 'throughput_pcsi_all_L_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_' legend_pre{pre} '_' legend_alg{alg} '_MC_' num2str(MC)],'png');
            saveas(gcf,[root_save 'throughput_pcsi_all_L_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_' legend_pre{pre} '_' legend_alg{alg} '_MC_' num2str(MC)],'epsc2');
        end
    end
end

for pre = 1:N_PRE
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    semilogx(var_err,reshape(OM*avg_max_sum_thrgpt_s(pre,1,:),1,[]),'-' ,'color',colours(1,:),'linewidth',linewidth);
    hold on;
    semilogx(var_err,reshape(OM*avg_max_sum_thrgpt_s(pre,2,:),1,[]),'-' ,'color',colours(2,:),'linewidth',linewidth);
    semilogx(var_err,reshape(OM*avg_max_sum_thrgpt_s(pre,3,:),1,[]),'-' ,'color',colours(3,:),'linewidth',linewidth);
    
    xlabel('$\sigma_{\varepsilon}^{2}$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
    ylabel(['Throughput ' um{abs(log10(OM))/3}],'fontname',fontname,'fontsize',fontsize);
    
    legend(legend_alg,'fontname',fontname,'fontsize',fontsize,'location',location_3,'numcolumns',1);
    legend box off;
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    if savefig == 1
        saveas(gcf,[root_save 'throughput_star_pcsi_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_' legend_pre{pre} '_MC_' num2str(MC)],'fig');
        saveas(gcf,[root_save 'throughput_star_pcsi_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_' legend_pre{pre} '_MC_' num2str(MC)],'png');
        saveas(gcf,[root_save 'throughput_star_pcsi_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_' legend_pre{pre} '_MC_' num2str(MC)],'epsc2');
    end
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    semilogx(var_err,reshape(avg_L_star(pre,1,:),1,[]),'-' ,'color',colours(1,:),'linewidth',linewidth);
    hold on;
    semilogx(var_err,reshape(avg_L_star(pre,2,:),1,[]),'-' ,'color',colours(2,:),'linewidth',linewidth);
    semilogx(var_err,reshape(avg_L_star(pre,3,:),1,[]),'-' ,'color',colours(3,:),'linewidth',linewidth);
    
    xlabel('$\sigma_{\varepsilon}^{2}$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
    ylabel('$L^{\star}$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
    
    if pre == 1
        legend(legend_alg,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',1);
    else
        legend(legend_alg,'fontname',fontname,'fontsize',fontsize,'location',location_3,'numcolumns',1);
    end
    legend box off;
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    if savefig == 1
        saveas(gcf,[root_save 'L_star_pcsi_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_' legend_pre{pre} '_MC_' num2str(MC)],'fig');
        saveas(gcf,[root_save 'L_star_pcsi_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_' legend_pre{pre} '_MC_' num2str(MC)],'png');
        saveas(gcf,[root_save 'L_star_pcsi_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_' legend_pre{pre} '_MC_' num2str(MC)],'epsc2');
    end
end