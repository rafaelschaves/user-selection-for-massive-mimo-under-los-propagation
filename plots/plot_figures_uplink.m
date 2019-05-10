clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = [64 256];                                                              % Number of antennas at base station
K = 18;                                                                    % Number of mobile users
L = 13;                                                                    % Number of selected users

snr = (-20:5:10)';                                                         % SNR in dB

M_SIZ = length(M);                                                         % Size of the antennas set
N_ALG = 4;                                                                 % Number of algorithms for perform user scheduling
N_SNR = length(snr);                                                       % Size of the SNR set 
N_CHN = 3;                                                                 % Number of channel models simulated

% Loading data

rate = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                      % Rate using all K users
psi  = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                      % Interchannel interference for all users

rate_sel = zeros(L,MC,N_ALG,M_SIZ,N_SNR,N_CHN);                            % Rate using only L users
% psi_sel  = zeros(L,MC,N_ALG,M_SIZ,N_SNR,N_CHN);                            % Interchannel interference for L users

for m_idx = 1:length(M)
    for snr_idx = 1:N_SNR
        load(['../results/uplink/rate_uplink_mf_ur-los_M_' num2str(M(m_idx)) ...
            '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) ...
            '_dB_MC_' num2str(MC) '.mat']);
        
        rate(:,:,m_idx,snr_idx,1) = rate_u;
        psi(:,:,m_idx,snr_idx,1)  = psi;
        
        rate_sel(:,:,:,m_idx,snr_idx,1) = rate_u_alg;
        % psi_sel(:,:,:,m_idx,snr_idx,1)  = psi_alg;
        
        load(['../results/uplink/rate_uplink_mf_sparse_M_' num2str(M(m_idx)) ...
            '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) ...
            '_dB_MC_' num2str(MC) '.mat']);
        
        rate(:,:,m_idx,snr_idx,2) = rate_u;
        psi(:,:,m_idx,snr_idx,2)  = psi;
        
        rate_sel(:,:,:,m_idx,snr_idx,2) = rate_u_alg;
        % psi_sel(:,:,:,m_idx,snr_idx,2)  = psi_alg;
        
        load(['../results/uplink/rate_uplink_mf_rayleigh_M_' num2str(M(m_idx)) ...
            '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) ...
            '_dB_MC_' num2str(MC) '.mat']);
        
        rate(:,:,m_idx,snr_idx,3) = rate_u;
        psi(:,:,m_idx,snr_idx,3)  = psi;
        
        rate_sel(:,:,:,m_idx,snr_idx,3) = rate_u_alg;
        % psi_sel(:,:,:,m_idx,snr_idx,3)  = psi_alg;
    end
end

% Post processing - Calculating the CDF

bin_width = 0.0005;

prob = cell(N_ALG+1,M_SIZ,N_SNR,N_CHN);
edge = cell(N_ALG+1,M_SIZ,N_SNR,N_CHN); 

for chn_idx = 1:N_CHN
    for snr_idx = 1:N_SNR
        for m_idx = 1:M_SIZ
            [prob{1,m_idx,snr_idx,chn_idx},edge{1,m_idx,snr_idx,chn_idx}] = histcounts(mean(rate(:,:,m_idx,snr_idx,chn_idx)),'binwidth',bin_width,'normalization','cdf');
            [prob{2,m_idx,snr_idx,chn_idx},edge{2,m_idx,snr_idx,chn_idx}] = histcounts(mean(rate_sel(:,:,1,m_idx,snr_idx,chn_idx)),'binwidth',bin_width,'normalization','cdf');
            [prob{3,m_idx,snr_idx,chn_idx},edge{3,m_idx,snr_idx,chn_idx}] = histcounts(mean(rate_sel(:,:,2,m_idx,snr_idx,chn_idx)),'binwidth',bin_width,'normalization','cdf');
            [prob{4,m_idx,snr_idx,chn_idx},edge{4,m_idx,snr_idx,chn_idx}] = histcounts(mean(rate_sel(:,:,3,m_idx,snr_idx,chn_idx)),'binwidth',bin_width,'normalization','cdf');
            [prob{5,m_idx,snr_idx,chn_idx},edge{5,m_idx,snr_idx,chn_idx}] = histcounts(mean(rate_sel(:,:,4,m_idx,snr_idx,chn_idx)),'binwidth',bin_width,'normalization','cdf');
            
            prob{1,m_idx,snr_idx,chn_idx} = [prob{1,m_idx,snr_idx,chn_idx} 1];
            prob{2,m_idx,snr_idx,chn_idx} = [prob{2,m_idx,snr_idx,chn_idx} 1];
            prob{3,m_idx,snr_idx,chn_idx} = [prob{3,m_idx,snr_idx,chn_idx} 1];
            prob{4,m_idx,snr_idx,chn_idx} = [prob{4,m_idx,snr_idx,chn_idx} 1];
            prob{5,m_idx,snr_idx,chn_idx} = [prob{5,m_idx,snr_idx,chn_idx} 1];
            
            edge{1,m_idx,snr_idx,chn_idx} = edge{1,m_idx,snr_idx,chn_idx} + bin_width/2;
            edge{2,m_idx,snr_idx,chn_idx} = edge{2,m_idx,snr_idx,chn_idx} + bin_width/2;
            edge{3,m_idx,snr_idx,chn_idx} = edge{3,m_idx,snr_idx,chn_idx} + bin_width/2;
            edge{4,m_idx,snr_idx,chn_idx} = edge{4,m_idx,snr_idx,chn_idx} + bin_width/2;
            edge{5,m_idx,snr_idx,chn_idx} = edge{5,m_idx,snr_idx,chn_idx} + bin_width/2;
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
channel_mod = {'ur_los','sparse','rayleigh'};

location = 'northwest';

root_out_prob = '../figures/rate/out_prob_';

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
        plot(edge{1,2,snr_idx,chn_idx},prob{1,2,snr_idx,chn_idx},'--','color',colours(1,:),'linewidth',linewidth);
        plot(edge{2,2,snr_idx,chn_idx},prob{2,2,snr_idx,chn_idx},'--','color',colours(2,:),'linewidth',linewidth);
        plot(edge{3,2,snr_idx,chn_idx},prob{3,2,snr_idx,chn_idx},'--','color',colours(3,:),'linewidth',linewidth);
        plot(edge{4,2,snr_idx,chn_idx},prob{4,2,snr_idx,chn_idx},'--','color',colours(4,:),'linewidth',linewidth);
        plot(edge{5,2,snr_idx,chn_idx},prob{5,2,snr_idx,chn_idx},'--','color',colours(5,:),'linewidth',linewidth);

        xlabel('Uplink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
        
        legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        ylim([0 1]);
        
        if (savefig == 1)
            saveas(gcf,[root_out_prob 'uplink_' channel_mod{chn_idx} '_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'fig');
            saveas(gcf,[root_out_prob 'uplink_' channel_mod{chn_idx} '_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'png');
            saveas(gcf,[root_out_prob 'uplink_' channel_mod{chn_idx} '_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'epsc2');
        end
    end
end