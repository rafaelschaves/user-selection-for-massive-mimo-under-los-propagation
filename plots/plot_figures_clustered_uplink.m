%clear;
%close all;
%clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = [256];                                                                 % Number of antennas at base station
K = [36];                                                                  % Number of mobile users
L = 13;                                                                    % Number of selected users

%snr = (-20:5:10)';                                                        % SNR in dB
snr = 10;

M_SIZ        = length(M);                                                         % Size of the antennas set
N_ALG        = 4;                                                                 % Number of algorithms for perform user scheduling
N_THETA_MID  = 3;
N_THETA_STEP = 5;

theta_mid  = [0 45 90];
theta_step = [1 10 20 30 90];

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Scheduling/Clustered/Uplink/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Scheduling/Clustered/Uplink/';

% Loading data

rate     = zeros(K,MC,N_THETA_STEP,N_THETA_MID,M_SIZ);                     % Rate using all K users
rate_sel = zeros(L,MC,N_ALG,N_THETA_STEP,N_THETA_MID,M_SIZ);               % Rate using only L users

for m = 1:M_SIZ
    for mid_idx = 1:N_THETA_MID
        for step_idx = 1:N_THETA_STEP
            load([root_load '/rate_mf_ur_los_M_' num2str(M(m)) '_K_' ...
                  num2str(K) '_L_' num2str(L) '_theta_mid_' ...
                  num2str(theta_mid(mid_idx)) '_theta_step_' ...
                  num2str(theta_step(step_idx)) '_SNR_10_dB_MC_' ...
                  num2str(MC) '.mat']);
    
            rate(:,:,step_idx,mid_idx,m)       = rate_u;
            rate_sel(:,:,:,step_idx,mid_idx,m) = rate_u_sel;
        end
    end
end
% Post processing - Calculating the CDF

bin_width = 0.0005;

prob = cell(N_ALG+1,N_THETA_STEP,N_THETA_MID,M_SIZ);
edge = cell(N_ALG+1,N_THETA_STEP,N_THETA_MID,M_SIZ); 

for m = 1:M_SIZ
    for mid_idx = 1:N_THETA_MID
        for step_idx = 1:N_THETA_STEP
            [prob{1,step_idx,mid_idx,m},edge{1,step_idx,mid_idx,m}] = ...
                histcounts(mean(rate(:,:,step_idx,mid_idx,m)),'binwidth',bin_width,'normalization','cdf');
            
            prob{1,step_idx,mid_idx,m} = [prob{1,step_idx,mid_idx,m} 1];
            edge{1,step_idx,mid_idx,m} = edge{1,step_idx,mid_idx,m} + bin_width/2;
            for alg_idx = 1:N_ALG
                [prob{1+alg_idx,step_idx,mid_idx,m},edge{1+alg_idx,step_idx,mid_idx,m}] = ...
                    histcounts(mean(rate_sel(:,:,alg_idx,step_idx,mid_idx,m)),'binwidth',bin_width,'normalization','cdf');
                
                prob{1+alg_idx,step_idx,mid_idx,m} = [prob{1+alg_idx,step_idx,mid_idx,m} 1];
                edge{1+alg_idx,step_idx,mid_idx,m} = edge{1+alg_idx,step_idx,mid_idx,m} + bin_width/2;
            end
        end
    end
end

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

savefig = 0;

% NS - No selection
% RS - Random selection
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_algo = {'NS','RS','SOS','CBS','ICIBS'};

location = 'northwest';

colours = get(gca,'colororder');
close;

for mid_idx = 1:N_THETA_MID
    for step_idx = 1:N_THETA_STEP
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        plot(edge{1,step_idx,mid_idx,1},prob{1,step_idx,mid_idx,1},'-','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(edge{2,step_idx,mid_idx,1},prob{2,step_idx,mid_idx,1},'-','color',colours(2,:),'linewidth',linewidth);
        plot(edge{3,step_idx,mid_idx,1},prob{3,step_idx,mid_idx,1},'-','color',colours(3,:),'linewidth',linewidth);
        plot(edge{4,step_idx,mid_idx,1},prob{4,step_idx,mid_idx,1},'-','color',colours(4,:),'linewidth',linewidth);
        plot(edge{5,step_idx,mid_idx,1},prob{5,step_idx,mid_idx,1},'-','color',colours(5,:),'linewidth',linewidth);
        
        xlabel('Normalized sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
        
        legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
        legend box off;
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        ylim([0 1]);
        
        if (savefig == 1)
            saveas(gcf,[root_save '64/out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) '_dB'],'fig');
            saveas(gcf,[root_save '64/out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) '_dB'],'png');
            saveas(gcf,[root_save '64/out_prob_norm_sum_rate_' chn_type{chn_idx} '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) '_dB'],'epsc2');
        end
    end
end