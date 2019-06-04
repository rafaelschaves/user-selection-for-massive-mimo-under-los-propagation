clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = [64];                                                              % Number of antennas at base station
K = [72];                                                            % Number of mobile users
L = 13;                                                                    % Number of selected users

%snr = (-20:5:10)';                                                         % SNR in dB
snr = 10;

M_SIZ = length(M);                                                         % Size of the antennas set
N_ALG = 4;                                                                 % Number of algorithms for perform user scheduling
N_SNR = length(snr);                                                       % Size of the SNR set 
N_CHN = 2;                                                                 % Number of channel models simulated

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Scheduling/Uplink/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Scheduling/Uplink/';

chn_type = {'ur_los','rayleigh'};

% Loading data

rate     = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                  % Rate using all K users
rate_sel = zeros(L,MC,N_ALG,M_SIZ,N_SNR,N_CHN);                            % Rate using only L users

for chn_idx = 1:N_CHN
    for m = 1:M_SIZ
        for snr_idx = 1:N_SNR
            load([root_load num2str(M(m)) '/rate_mf_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) '_dB_MC_' num2str(MC) '.mat']);
            
            rate(:,:,m,snr_idx,chn_idx)       = rate_u;
            rate_sel(:,:,:,m,snr_idx,chn_idx) = rate_u_sel;
        end
    end
end

% Post processing - Calculating the CDF

bin_width = 0.0005;

prob = cell(N_ALG+1,M_SIZ,N_SNR,N_CHN);
edge = cell(N_ALG+1,M_SIZ,N_SNR,N_CHN); 

for chn_idx = 1:N_CHN
    for snr_idx = 1:N_SNR
        for m = 1:M_SIZ
            [prob{1,m,snr_idx,chn_idx},edge{1,m,snr_idx,chn_idx}] = histcounts(mean(rate(:,:,m,snr_idx,chn_idx)),'binwidth',bin_width,'normalization','cdf');
            [prob{2,m,snr_idx,chn_idx},edge{2,m,snr_idx,chn_idx}] = histcounts(mean(rate_sel(:,:,1,m,snr_idx,chn_idx)),'binwidth',bin_width,'normalization','cdf');
            [prob{3,m,snr_idx,chn_idx},edge{3,m,snr_idx,chn_idx}] = histcounts(mean(rate_sel(:,:,2,m,snr_idx,chn_idx)),'binwidth',bin_width,'normalization','cdf');
            [prob{4,m,snr_idx,chn_idx},edge{4,m,snr_idx,chn_idx}] = histcounts(mean(rate_sel(:,:,3,m,snr_idx,chn_idx)),'binwidth',bin_width,'normalization','cdf');
            [prob{5,m,snr_idx,chn_idx},edge{5,m,snr_idx,chn_idx}] = histcounts(mean(rate_sel(:,:,4,m,snr_idx,chn_idx)),'binwidth',bin_width,'normalization','cdf');
            
            prob{1,m,snr_idx,chn_idx} = [prob{1,m,snr_idx,chn_idx} 1];
            prob{2,m,snr_idx,chn_idx} = [prob{2,m,snr_idx,chn_idx} 1];
            prob{3,m,snr_idx,chn_idx} = [prob{3,m,snr_idx,chn_idx} 1];
            prob{4,m,snr_idx,chn_idx} = [prob{4,m,snr_idx,chn_idx} 1];
            prob{5,m,snr_idx,chn_idx} = [prob{5,m,snr_idx,chn_idx} 1];
            
            edge{1,m,snr_idx,chn_idx} = edge{1,m,snr_idx,chn_idx} + bin_width/2;
            edge{2,m,snr_idx,chn_idx} = edge{2,m,snr_idx,chn_idx} + bin_width/2;
            edge{3,m,snr_idx,chn_idx} = edge{3,m,snr_idx,chn_idx} + bin_width/2;
            edge{4,m,snr_idx,chn_idx} = edge{4,m,snr_idx,chn_idx} + bin_width/2;
            edge{5,m,snr_idx,chn_idx} = edge{5,m,snr_idx,chn_idx} + bin_width/2;
        end
    end
end

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

savefig = 1;

% NS - No selection
% RS - Random selection
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_algo = {'NS','RS','SOS','CBS','ICIBS'};

location = 'northwest';

colours = get(gca,'colororder');
close;

for chn_idx = 1:N_CHN
    for snr_idx = 1:N_SNR        
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        plot(edge{1,1,snr_idx,chn_idx},prob{1,1,snr_idx,chn_idx},'-','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(edge{2,1,snr_idx,chn_idx},prob{2,1,snr_idx,chn_idx},'-','color',colours(2,:),'linewidth',linewidth);
        plot(edge{3,1,snr_idx,chn_idx},prob{3,1,snr_idx,chn_idx},'-','color',colours(3,:),'linewidth',linewidth);
        plot(edge{4,1,snr_idx,chn_idx},prob{4,1,snr_idx,chn_idx},'-','color',colours(4,:),'linewidth',linewidth);
        plot(edge{5,1,snr_idx,chn_idx},prob{5,1,snr_idx,chn_idx},'-','color',colours(5,:),'linewidth',linewidth);
%         plot(edge{1,2,snr_idx,chn_idx},prob{1,2,snr_idx,chn_idx},'--','color',colours(1,:),'linewidth',linewidth);
%         plot(edge{2,2,snr_idx,chn_idx},prob{2,2,snr_idx,chn_idx},'--','color',colours(2,:),'linewidth',linewidth);
%         plot(edge{3,2,snr_idx,chn_idx},prob{3,2,snr_idx,chn_idx},'--','color',colours(3,:),'linewidth',linewidth);
%         plot(edge{4,2,snr_idx,chn_idx},prob{4,2,snr_idx,chn_idx},'--','color',colours(4,:),'linewidth',linewidth);
%         plot(edge{5,2,snr_idx,chn_idx},prob{5,2,snr_idx,chn_idx},'--','color',colours(5,:),'linewidth',linewidth);
        
        if(snr_idx == 7 && chn_idx == 1)
            dim(1,:) = [0.33 0.5 0.225 0.1];
            dim(2,:) = [0.57 0.5 0.225 0.1];
            
            annotation('ellipse',dim(1,:),'linewidth',linewidth);
            annotation('ellipse',dim(2,:),'linewidth',linewidth);
            
            cord(1,:) = [1.5 9.7];
            cord(2,:) = [0.425 0.425];
            
            text(cord(1,1),cord(2,1),'$M = 64$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
            text(cord(1,2),cord(2,2),'$M = 256$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
        end
        
        if(snr_idx == 7 && chn_idx == 3)
            % xlim([0 6]);
            % dim(1,:) = [0.375 0.5 0.2 0.1];
            % dim(2,:) = [0.600 0.5 0.2 0.1];
            
            % xlim([1 6]);
            dim(1,:) = [0.280 0.5 0.2 0.1];
            dim(2,:) = [0.570 0.5 0.2 0.1];
            
            annotation('ellipse',dim(1,:),'linewidth',linewidth);
            annotation('ellipse',dim(2,:),'linewidth',linewidth);
            
            cord(1,:) = [3.0 4.9];
            cord(2,:) = [0.4 0.4];
            
            text(cord(1,1),cord(2,1),'$M = 64$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
            text(cord(1,2),cord(2,2),'$M = 256$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
        end

        xlabel('Normalized sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
        
        legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
        legend box off;
        
        if((snr_idx == 7 && chn_idx == 1) || (snr_idx == 7 && chn_idx == 3))
            legend box off;
        end
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        if(chn_idx == 3)
            xlim([1 6]);
        end
        
        ylim([0 1]);
        
        if (savefig == 1)
            saveas(gcf,[root_save '64/out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) '_dB'],'fig');
            saveas(gcf,[root_save '64/out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) '_dB'],'png');
            saveas(gcf,[root_save '64/out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) '_dB'],'epsc2');
        end
    end
end