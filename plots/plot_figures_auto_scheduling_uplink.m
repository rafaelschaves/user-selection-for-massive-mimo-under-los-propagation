clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = 64;                                                                    % Number of antennas at base station
K = 18;                                                                    % Number of mobile users

snr = 20;                                                                  % SNR in dB

tau_min = 0.01;
tau_max = 0.5;

tau_step = 0.005;

tau = tau_min:tau_step:tau_max;

M_SIZ = length(M);                                                         % Size of the antennas set
N_ALG = 2;                                                                 % Number of algorithms for perform user scheduling
N_SNR = length(snr);                                                       % Size of the SNR set 
N_CHN = 1;                                                                 % Number of channel models simulated
N_TAU = length(tau);

% Root

root = '../results/auto_scheduling/uplink/rate_uplink_auto_scheduling_mf_';

% Loading data

% rate = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                      % Rate using all K users
% psi  = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                      % Interchannel interference for all users

% rate_sel = zeros(L,MC,N_ALG,M_SIZ,N_SNR,N_CHN);                            % Rate using only L users
% psi_sel  = zeros(L,MC,N_ALG,M_SIZ,N_SNR,N_CHN);                            % Interchannel interference for L users

for m_idx = 1:length(M)
    for snr_idx = 1:N_SNR
        load([root 'ur-los_M_' num2str(M(m_idx)) '_K_' num2str(K) '_tau_' ...
              num2str(N_TAU) '_SNR_' num2str(snr(snr_idx)) '_dB_MC_' ...
              num2str(MC) '.mat']);
        
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

% up_rate_cbs = [];
% up_rate_icibs = [];

up_norm_sum_rate = zeros(MC,N_TAU,N_ALG);

prob = cell(N_TAU,N_ALG);
edge = cell(N_TAU,N_ALG); 

prob_L = cell(N_TAU,N_ALG);
edge_L = cell(N_TAU,N_ALG); 

bin_width = 0.0005;
    
for tau_idx = 1:N_TAU
    for mc = 1:MC
        % up_rate_cbs = [up_rate_cbs; rate_u_alg{mc,1}];
        % up_rate_icibs = [up_rate_icibs; rate_u_alg{mc,2}];
        
        up_norm_sum_rate(mc,tau_idx,1) = mean(rate_u_alg{mc,tau_idx,1});
        up_norm_sum_rate(mc,tau_idx,2) = mean(rate_u_alg{mc,tau_idx,2});
    end
    
    [prob{tau_idx,1},edge{tau_idx,1}] = histcounts(up_norm_sum_rate(:,tau_idx,1),'binwidth',bin_width,'normalization','cdf');
    [prob{tau_idx,2},edge{tau_idx,2}] = histcounts(up_norm_sum_rate(:,tau_idx,2),'binwidth',bin_width,'normalization','cdf');
    
    % [prob{1},edge{1}] = histcounts(up_rate_cbs,'binwidth',bin_width,'normalization','cdf');
    % [prob{2},edge{2}] = histcounts(up_rate_icibs,'binwidth',bin_width,'normalization','cdf');
    
    prob{tau_idx,1} = [prob{tau_idx,1} 1];
    prob{tau_idx,2} = [prob{tau_idx,2} 1];
    
    edge{tau_idx,1} = edge{tau_idx,1} + bin_width/2;
    edge{tau_idx,2} = edge{tau_idx,2} + bin_width/2;
    
    [prob_L{tau_idx,1},edge_L{tau_idx,1}] = histcounts(K - L(:,tau_idx,1),'binwidth',bin_width,'normalization','cdf');
    [prob_L{tau_idx,2},edge_L{tau_idx,2}] = histcounts(K - L(:,tau_idx,2),'binwidth',bin_width,'normalization','cdf');
    
    prob_L{tau_idx,1} = [prob_L{tau_idx,1} 1];
    prob_L{tau_idx,2} = [prob_L{tau_idx,2} 1];    
end

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

savefig = 1;

% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_algo = {'CBS','ICIBS'};

tau_idx_aux = [1 10 15 20 50 70 99];

legend_tau  = {'$\tau = 0.010$', ...
               '$\tau = 0.055$', ...
               '$\tau = 0.080$', ...
               '$\tau = 0.105$', ...
               '$\tau = 0.255$', ...
               '$\tau = 0.355$', ...
               '$\tau = 0.500$'};
      
channel_mod = {'ur_los','sparse','rayleigh'};
algorithm = {'cbs','icibs'};

location_1 = 'northwest';
location_2 = 'southeast';

root_rate = '../figures/auto_scheduling/uplink/normalized_sum_rate/auto_scheduling_norm_sum_rate_uplink_';
root_drop = '../figures/auto_scheduling/uplink/dropped_users/cdf_dropped_users_uplink_';

colours = get(gca,'colororder');
close;

L_max = 17;

N_TAU = length(tau_idx_aux);

for chn_idx = 1:N_CHN
    for alg_idx = 1:N_ALG
        for snr_idx = 1:N_SNR
            figure;
                
            set(gcf,'position',[0 0 800 600]);
   
            for tau_idx = 1:N_TAU    
                % plot(edge{tau_idx,alg_idx},prob{tau_idx,alg_idx},'-','linewidth',linewidth);
                plot(edge{tau_idx_aux(tau_idx),alg_idx},prob{tau_idx_aux(tau_idx),alg_idx},'-','linewidth',linewidth);
                
                if(tau_idx == 1)
                    hold on;
                end
            end
            
            xlabel('Uplink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
            ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
                
            legend(legend_tau,'fontname',fontname,'fontsize',fontsize,'location',location_1,'interpreter','latex');
            % legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_1);
            % ,'color',colours(1,:)
                
            set(gca,'fontname',fontname,'fontsize',fontsize);
                
            ylim([0 1]);
        end
        
        if (savefig == 1)
            saveas(gcf,[root_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'fig');
            saveas(gcf,[root_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'png');
            saveas(gcf,[root_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'epsc2');
        end
    end
end

for chn_idx = 1:N_CHN
    for alg_idx = 1:N_ALG
        for snr_idx = 1:N_SNR
            figure;
            
            set(gcf,'position',[0 0 800 600]);
            
            for tau_idx = 1:N_TAU
                % plot(edge_L{tau_idx,alg_idx},prob_L{tau_idx,alg_idx},'-','linewidth',linewidth);
                plot(edge_L{tau_idx_aux(tau_idx),alg_idx},prob_L{tau_idx_aux(tau_idx),alg_idx},'-','linewidth',linewidth);
                
                if(tau_idx == 1)
                    hold on;
                end             
            end
            
            xlabel('Number of dropped users','fontname',fontname,'fontsize',fontsize);
            ylabel('CDF','fontname',fontname,'fontsize',fontsize);
            
            legend(legend_tau,'fontname',fontname,'fontsize',fontsize,'location',location_1,'interpreter','latex');
            % legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_2);
            
            set(gca,'fontname',fontname,'fontsize',fontsize);
            
            xlim([0 L_max]);
            ylim([0 1]);
        end
        
        if (savefig == 1)
            saveas(gcf,[root_drop channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'fig');
            saveas(gcf,[root_drop channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'png');
            saveas(gcf,[root_drop channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'epsc2');
        end
    end
end