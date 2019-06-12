clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = 64;                                                                    % Number of antennas at base station
K = 18;                                                                    % Number of mobile users
L = 13;                                                                    % Number of selected users

snr = 10;

N_ALG        = 4;                                                          % Number of algorithms for perform user scheduling
N_THETA_MID  = 3;
N_THETA_STEP = 5;

theta_mid  = [0 45 90];
theta_step = [1 10 20 30 90];

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Clustered/Uplink/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Clustered/Uplink/';

if ~exist(root_save,'dir')
    mkdir(root_save);
end

% Loading data

thrput     = zeros(K,MC,N_THETA_STEP,N_THETA_MID);                         % Rate using all K users
thrput_sel = zeros(L,MC,N_ALG,N_THETA_STEP,N_THETA_MID);                   % Rate using only L users

for mid_idx = 1:N_THETA_MID
    for step_idx = 1:N_THETA_STEP
        load([root_load 'throughput_outdoors_pedestrian_mf_ur_los_M_' num2str(M) '_K_' ...
            num2str(K) '_L_' num2str(L) '_theta_mid_' ...
            num2str(theta_mid(mid_idx)) '_theta_step_' ...
            num2str(theta_step(step_idx)) '_SNR_10_dB_MC_' ...
            num2str(MC) '.mat']);
        
        thrput(:,:,step_idx,mid_idx)       = thrput_u;
        thrput_sel(:,:,:,step_idx,mid_idx) = thrput_u_sel;
    end
end

% Post processing - Calculating the CDF

bin_width = 0.0005;

prob = cell(N_ALG+1,N_THETA_STEP,N_THETA_MID);
edge = cell(N_ALG+1,N_THETA_STEP,N_THETA_MID); 

for mid_idx = 1:N_THETA_MID
    for step_idx = 1:N_THETA_STEP
        [prob{1,step_idx,mid_idx},edge{1,step_idx,mid_idx}] = ...
            histcounts(mean(thrput(:,:,step_idx,mid_idx)),'binwidth',bin_width,'normalization','cdf');
        
        prob{1,step_idx,mid_idx} = [prob{1,step_idx,mid_idx} 1];
        edge{1,step_idx,mid_idx} = edge{1,step_idx,mid_idx} + bin_width/2;
        for alg_idx = 1:N_ALG
            [prob{1+alg_idx,step_idx,mid_idx},edge{1+alg_idx,step_idx,mid_idx}] = ...
                histcounts(mean(thrput_sel(:,:,alg_idx,step_idx,mid_idx)),'binwidth',bin_width,'normalization','cdf');
            
            prob{1+alg_idx,step_idx,mid_idx} = [prob{1+alg_idx,step_idx,mid_idx} 1];
            edge{1+alg_idx,step_idx,mid_idx} = edge{1+alg_idx,step_idx,mid_idx} + bin_width/2;
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

for mid_idx = 1:N_THETA_MID
    for step_idx = 1:N_THETA_STEP
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        plot(1e-6*edge{1,step_idx,mid_idx},prob{1,step_idx,mid_idx},'-','color',colours(1,:),'linewidth',linewidth);
        hold on;
        
        for alg_idx = 1:N_ALG
            plot(1e-6*edge{1 + alg_idx,step_idx,mid_idx},prob{1 + alg_idx,step_idx,mid_idx},'-','color',colours(1 + alg_idx,:),'linewidth',linewidth);
        end
                
        xlabel('Normalized sum-throughput (Mbps/cell)','fontname',fontname,'fontsize',fontsize);
        ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
        
        legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
        legend box off;
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([0 inf]);
        ylim([0 1]);
        
        if (savefig == 1)
            saveas(gcf,[root_save 'throughput_out_ped_mf_ur_los_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_theta_mid_' num2str(theta_mid(mid_idx)) '_theta_step_' num2str(theta_step(step_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC)],'fig');
            saveas(gcf,[root_save 'throughput_out_ped_mf_ur_los_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_theta_mid_' num2str(theta_mid(mid_idx)) '_theta_step_' num2str(theta_step(step_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC)],'png');
            saveas(gcf,[root_save 'throughput_out_ped_mf_ur_los_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_theta_mid_' num2str(theta_mid(mid_idx)) '_theta_step_' num2str(theta_step(step_idx)) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC)],'epsc2');
        end
    end
end