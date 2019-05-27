clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = 64;                                                                    % Number of antennas at base station
K = 18;                                                                    % Number of mobile users

snr = 20;                                                                  % SNR in dB

% Threshold - CBS

tau_cbs_min = 0;
tau_cbs_max = 1;

tau_cbs_step = 0.01;

tau_cbs = tau_cbs_min:tau_cbs_step:tau_cbs_max;

N_TAU_CBS = length(tau_cbs);

% Threshold - ICIBS

tau_icibs_min = 0;
tau_icibs_max = 0.25;

tau_icibs_step = 0.0025;

tau_icibs = tau_icibs_min:tau_icibs_step:tau_icibs_max;

N_TAU_ICIBS = length(tau_icibs);

M_SIZ = length(M);                                                         % Size of the antennas set
N_ALG = 2;                                                                 % Number of algorithms for perform user scheduling
N_SNR = length(snr);                                                       % Size of the SNR set 
N_CHN = 1;                                                                 % Number of channel models simulated

% Root

% root_load = '../results/auto_scheduling/uplink/rate_uplink_auto_scheduling_mf_';
root_load = '../results/auto_scheduling/uplink/rate_uplink_mf_';

root_save_rate = '../figures/auto_scheduling/uplink/normalized_sum_rate/auto_scheduling_norm_sum_rate_uplink_';
root_save_drop = '../figures/auto_scheduling/uplink/dropped_users/cdf_dropped_users_uplink_';

% Loading data

% rate = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                      % Rate using all K users
% psi  = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                      % Interchannel interference for all users

% rate_sel = zeros(L,MC,N_ALG,M_SIZ,N_SNR,N_CHN);                            % Rate using only L users
% psi_sel  = zeros(L,MC,N_ALG,M_SIZ,N_SNR,N_CHN);                            % Interchannel interference for L users

for m_idx = 1:length(M)
    for snr_idx = 1:N_SNR
        load([root_load 'ur-los_M_' num2str(M(m_idx)) '_K_' num2str(K) ...
              '_SNR_' num2str(snr(snr_idx)) '_dB_MC_' ...
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

max_rate     = 15;
rate_samples = 15000;

L_max = 18;

tau = tau_cbs;
N_TAU = N_TAU_CBS;

up_norm_sum_rate = zeros(MC,N_TAU,N_ALG);

edges_rate = max_rate*linspace(0,1,rate_samples);
edges_sele = 0:L_max;

prob_rate   = zeros(rate_samples - 1,N_TAU,N_ALG);
prob_rate_f = zeros(rate_samples,N_TAU,N_ALG);

prob_sele = zeros(L_max,N_TAU,N_ALG);

edges_rate_grid = repmat(edges_rate,N_TAU,1);
tau_rate_grid   = repmat(tau',1,rate_samples);

edges_sele_grid = repmat(edges_sele(1:end-1),N_TAU,1);
tau_sele_grid   = repmat(tau',1,L_max);
    
for tau_idx = 2:N_TAU_CBS
    for mc = 1:MC
        % up_norm_sum_rate(mc,tau_idx,1) = mean(rate_u_alg{mc,tau_idx,1});
        % up_norm_sum_rate(mc,tau_idx,2) = mean(rate_u_alg{mc,tau_idx,2});
        
        up_norm_sum_rate(mc,tau_idx,1) = mean(rate_u_cbs{mc,tau_idx});
    end
        
    prob_rate(:,tau_idx,1) = histcounts(up_norm_sum_rate(:,tau_idx,1),edges_rate,'normalization','cdf');
    
    prob_sele(:,tau_idx,1) = histcounts(K - L_cbs(:,tau_idx),edges_sele,'normalization','probability');
        
    prob_rate_f(:,tau_idx,1) = [prob_rate(:,tau_idx,1); 1];
end

for tau_idx = 2:N_TAU_ICIBS
    for mc = 1:MC
        % up_norm_sum_rate(mc,tau_idx,1) = mean(rate_u_alg{mc,tau_idx,1});
        % up_norm_sum_rate(mc,tau_idx,2) = mean(rate_u_alg{mc,tau_idx,2});
        
        up_norm_sum_rate(mc,tau_idx,2) = mean(rate_u_icibs{mc,tau_idx});
    end
        
    prob_rate(:,tau_idx,2) = histcounts(up_norm_sum_rate(:,tau_idx,2),edges_rate,'normalization','cdf');
    
    prob_sele(:,tau_idx,2) = histcounts(K - L_icibs(:,tau_idx),edges_sele,'normalization','probability');
        
    prob_rate_f(:,tau_idx,2) = [prob_rate(:,tau_idx,2); 1];
end

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

savefig = 0;
      
channel_mod = {'ur_los','sparse','rayleigh'};
algorithm = {'cbs','icibs'};

location_1 = 'northwest';
location_2 = 'southeast';

n_contour = 15;

for chn_idx = 1:N_CHN
    for alg_idx = 1:N_ALG
        for snr_idx = 1:N_SNR
            figure;
                
            set(gcf,'position',[0 0 800 600]);
   
            surf(edges_rate_grid,tau_rate_grid,prob_rate_f(:,:,alg_idx)','edgecolor','none');

            xlabel('Uplink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
            ylabel('$\tau$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
            zlabel('Outage probability','fontname',fontname,'fontsize',fontsize);
            
            colorbar;
            
            set(gca,'fontname',fontname,'fontsize',fontsize);
            
            xlim([0 15]);
            ylim([0 1]);
            zlim([0 1]);    
        end
        
        if (savefig == 1)
            saveas(gcf,[root_save_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'fig');
            saveas(gcf,[root_save_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'png');
            saveas(gcf,[root_save_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'epsc2');
        end
    end
end

for chn_idx = 1:N_CHN
    for alg_idx = 1:N_ALG
        for snr_idx = 1:N_SNR
            figure;
                
            set(gcf,'position',[0 0 800 600]);
   
            contour(edges_rate_grid,tau_rate_grid,prob_rate_f(:,:,alg_idx)',n_contour,'linewidth',linewidth);
                        
            xlabel('Uplink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
            ylabel('$\tau$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
            
            colorbar;
            
            set(gca,'fontname',fontname,'fontsize',fontsize);         
        end
        
        if (savefig == 1)
            saveas(gcf,[root_save_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'fig');
            saveas(gcf,[root_save_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'png');
            saveas(gcf,[root_save_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'epsc2');
        end
    end
end

for chn_idx = 1:N_CHN
    for alg_idx = 1:N_ALG
        for snr_idx = 1:N_SNR
            figure;
            
            set(gcf,'position',[0 0 800 600]);
                       
            surf(edges_sele_grid,tau_sele_grid,prob_sele(:,:,alg_idx)','edgecolor','none');
            
            view(2);
            
            xlabel('Number of dropped users','fontname',fontname,'fontsize',fontsize);
            ylabel('$\tau$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
            zlabel('CDF','fontname',fontname,'fontsize',fontsize);
            
            colorbar;
                        
            set(gca,'fontname',fontname,'fontsize',fontsize);
            
            xlim([0 L_max]);
            ylim([0 1]);
            zlim([0 1]);
        end
        
        if (savefig == 1)
            saveas(gcf,[root_drop channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'fig');
            saveas(gcf,[root_drop channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'png');
            saveas(gcf,[root_drop channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_idx)],'epsc2');
        end
    end
end