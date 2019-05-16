clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = 64;                                                              % Number of antennas at base station
K = 18;                                                                    % Number of mobile users

snr = 20;                                                         % SNR in dB

tau = 0.2;

M_SIZ = length(M);                                                         % Size of the antennas set
N_ALG = 2;                                                                 % Number of algorithms for perform user scheduling
N_SNR = length(snr);                                                       % Size of the SNR set 
N_CHN = 1;                                                                 % Number of channel models simulated

% Loading data

% rate = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                      % Rate using all K users
% psi  = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                      % Interchannel interference for all users

% rate_sel = zeros(L,MC,N_ALG,M_SIZ,N_SNR,N_CHN);                            % Rate using only L users
% psi_sel  = zeros(L,MC,N_ALG,M_SIZ,N_SNR,N_CHN);                            % Interchannel interference for L users

for m_idx = 1:length(M)
    for snr_idx = 1:N_SNR
        load(['../results/uplink/rate_uplink_auto_scheduling_mf_ur-los_M_' num2str(M(m_idx)) ...
            '_K_' num2str(K) '_tau_' num2str(tau) '_SNR_' num2str(snr(snr_idx)) ...
            '_dB_MC_' num2str(MC) '.mat']);
                
        % rate_sel(:,:,:,m_idx,snr_idx,1) = rate_u_alg;
        % psi_sel(:,:,:,m_idx,snr_idx,1)  = psi_alg;
        
%         load(['../results/uplink/rate_uplink_mf_sparse_M_' num2str(M(m_idx)) ...
%             '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) ...
%             '_dB_MC_' num2str(MC) '.mat']);
%         
%         rate(:,:,m_idx,snr_idx,2) = rate_u;
%         psi(:,:,m_idx,snr_idx,2)  = psi;
%         
%         rate_sel(:,:,:,m_idx,snr_idx,2) = rate_u_alg;
%         % psi_sel(:,:,:,m_idx,snr_idx,2)  = psi_alg;
%         
%         load(['../results/uplink/rate_uplink_mf_rayleigh_M_' num2str(M(m_idx)) ...
%             '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) ...
%             '_dB_MC_' num2str(MC) '.mat']);
%         
%         rate(:,:,m_idx,snr_idx,3) = rate_u;
%         psi(:,:,m_idx,snr_idx,3)  = psi;
%         
%         rate_sel(:,:,:,m_idx,snr_idx,3) = rate_u_alg;
%         % psi_sel(:,:,:,m_idx,snr_idx,3)  = psi_alg;
    end
end

% Post processing - Calculating the CDF

up_rate_cbs = [];
up_rate_icibs = [];

up_rate_per_terminal = zeros(MC,N_ALG);

for mc = 1:MC
    up_rate_cbs = [up_rate_cbs; rate_u_alg{mc,1}];
    up_rate_icibs = [up_rate_icibs; rate_u_alg{mc,2}];

    up_rate_per_terminal(mc,1) = mean(rate_u_alg{mc,1});
    up_rate_per_terminal(mc,2) = mean(rate_u_alg{mc,2});
end

bin_width = 0.0005;

prob = cell(N_ALG);
edge = cell(N_ALG); 

[prob{1},edge{1}] = histcounts(up_rate_per_terminal(:,1),'binwidth',bin_width,'normalization','cdf');
[prob{2},edge{2}] = histcounts(up_rate_per_terminal(:,2),'binwidth',bin_width,'normalization','cdf');

% [prob{1},edge{1}] = histcounts(up_rate_cbs,'binwidth',bin_width,'normalization','cdf');
% [prob{2},edge{2}] = histcounts(up_rate_icibs,'binwidth',bin_width,'normalization','cdf');

prob{1} = [prob{1} 1];
prob{2} = [prob{2} 1];

edge{1} = edge{1} + bin_width/2;
edge{2} = edge{2} + bin_width/2;

bin_width_L = 1;

[prob_L{1},edge_L{1}] = histcounts(K - L(:,1),'binwidth',bin_width,'normalization','cdf');
[prob_L{2},edge_L{2}] = histcounts(K - L(:,2),'binwidth',bin_width,'normalization','cdf');

prob_L{1} = [prob_L{1} 1];
prob_L{2} = [prob_L{2} 1];

edge_L{1} = edge_L{1};
edge_L{2} = edge_L{2} + bin_width_L/2;


% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

savefig = 0;

% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_algo = {'CBS','ICIBS'};
channel_mod = {'ur_los','sparse','rayleigh'};

location = 'northwest';

root_out_prob = '../figures/rate/out_prob_';

colours = get(gca,'colororder');
close;

for chn_idx = 1:N_CHN
    for snr_idx = 1:N_SNR        
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        plot(edge{1},prob{1},'-','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(edge{2},prob{2},'-','color',colours(2,:),'linewidth',linewidth);
       
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

for chn_idx = 1:N_CHN
    for snr_idx = 1:N_SNR        
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        plot(edge_L{1},prob_L{1},'-','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(edge_L{2},prob_L{2},'-','color',colours(2,:),'linewidth',linewidth);
       
        xlabel('Number of dropped users','fontname',fontname,'fontsize',fontsize);
        ylabel('CDF','fontname',fontname,'fontsize',fontsize);
        
        legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([1 10]);
        ylim([0 1]);
        
        if (savefig == 1)
            saveas(gcf,[root_out_prob 'uplink_' channel_mod{chn_idx} '_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'fig');
            saveas(gcf,[root_out_prob 'uplink_' channel_mod{chn_idx} '_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'png');
            saveas(gcf,[root_out_prob 'uplink_' channel_mod{chn_idx} '_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'epsc2');
        end
    end
end