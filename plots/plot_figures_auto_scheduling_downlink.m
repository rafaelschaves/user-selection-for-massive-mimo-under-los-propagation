clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = 256;                                                                    % Number of antennas at base station
K = 18;                                                                    % Number of mobile users

snr = 20;                                                                  % SNR in dB

% Normalized threshold

tau_norm_min = 0;
tau_norm_max = 1;

tau_norm_step = 0.01;

tau_norm = tau_norm_min:tau_norm_step:tau_norm_max;

tau_icibs_max = 0.25;

N_ALG = 2;                                                                 % Number of algorithms for perform user scheduling
N_SNR = length(snr);                                                       % Size of the SNR set 
N_CHN = 1;                                                                 % Number of channel models simulated
N_TAU = length(tau_norm);

T_c = 2.5e-3;
time = T_c*(0:MC-1);

% Root

root_load = '../results/auto_scheduling/downlink/rate_downlink_mf_';

root_save_rate = '../figures/auto_scheduling/downlink/normalized_sum_rate/cdf_';
root_save_drop = '../figures/auto_scheduling/downlink/dropped_users/pdf_';
root_save_time = '../figures/auto_scheduling/downlink/time/';

% Loading data

% rate_sel = zeros(L,MC,N_ALG,M_SIZ,N_SNR,N_CHN);                            % Rate using only L users

for snr_idx = 1:N_SNR
    load([root_load 'ur_los_M_' num2str(M) '_K_' num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr(snr_idx)) '_dB_MC_' num2str(MC) '.mat']);
    
    % rate_sel(:,:,:,m_idx,snr_idx,1) = rate_u_alg;
    
    % load(['../results/uplink/rate_uplink_mf_sparse_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) '_dB_MC_' num2str(MC) '.mat']);
    %
    % rate_sel(:,:,:,m_idx,snr_idx,2) = rate_u_alg;
    %
    % load(['../results/uplink/rate_uplink_mf_rayleigh_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) '_dB_MC_' num2str(MC) '.mat']);
    %
    % rate_sel(:,:,:,m_idx,snr_idx,3) = rate_u_alg;
end

rate = zeros(MC,K,N_TAU,N_ALG);

for alg_idx = 1:N_ALG
    for mc = 1:MC
        for tau_idx = 1:N_TAU
            rate(mc,user_sel{mc,tau_idx,alg_idx},tau_idx,alg_idx) = rate_d{mc,tau_idx,alg_idx};
        end
    end
end

% Post processing - Calculating the CDF

max_rate     = 15;
rate_samples = 15000;

L_max = 19;

do_norm_sum_rate = zeros(MC,N_TAU,N_ALG);

edges_rate = max_rate*linspace(0,1,rate_samples);
edges_sele = 0:L_max;

prob_rate   = zeros(rate_samples - 1,N_TAU,N_ALG);
prob_rate_f = zeros(rate_samples,N_TAU,N_ALG);

prob_sele = zeros(L_max,N_TAU,N_ALG);

edges_rate_grid = repmat(edges_rate,N_TAU,1);
tau_rate_grid   = repmat(tau_norm',1,rate_samples);

edges_sele_grid = repmat(edges_sele(1:end-1),N_TAU,1);
tau_sele_grid   = repmat(tau_norm',1,L_max);
    
for alg_idx = 1:N_ALG
    for tau_idx = 1:N_TAU
        for mc = 1:MC
            do_norm_sum_rate(mc,tau_idx,alg_idx) = mean(rate_d{mc,tau_idx,alg_idx});
        end
        
        prob_rate(:,tau_idx,alg_idx) = histcounts(do_norm_sum_rate(:,tau_idx,alg_idx),edges_rate,'normalization','cdf');
        prob_sele(:,tau_idx,alg_idx) = histcounts(K - L(:,tau_idx,alg_idx),edges_sele,'normalization','probability');
        
        prob_rate_f(:,tau_idx,alg_idx) = [prob_rate(:,tau_idx,alg_idx); 1];
    end
end

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

savefig = 1;
      
channel_mod = {'ur_los','sparse','rayleigh'};
algorithm   = {'cbs','icibs'};

legend_tau_cbs    = {['$ \tau = ' num2str(tau_norm(end)) '$'], ...
                     ['$ \tau = ' num2str(tau_norm(10)) '$'], ...
                     ['$ \tau = ' num2str(tau_norm(5)) '$'], ...
                     ['$ \tau = ' num2str(tau_norm(1)) '$']};

legend_tau_cbs2   = {['$ \tau = ' num2str(tau_norm(end)) '$'], ...
                     ['$ \tau = ' num2str(tau_norm(10)) '$'], ...
                     ['$ \tau = ' num2str(tau_norm(5)) '$']};
                
legend_tau_icibs  = {['$ \tau = ' num2str(tau_icibs_max*tau_norm(end)) '$'], ...
                     ['$ \tau = ' num2str(tau_icibs_max*tau_norm(10)) '$'], ...
                     ['$ \tau = ' num2str(tau_icibs_max*tau_norm(5)) '$'], ...
                     ['$ \tau = ' num2str(tau_icibs_max*tau_norm(1)) '$']};

legend_tau_icibs2 = {['$ \tau = ' num2str(tau_icibs_max*tau_norm(end)) '$'], ...
                     ['$ \tau = ' num2str(tau_icibs_max*tau_norm(10)) '$'], ...
                     ['$ \tau = ' num2str(tau_icibs_max*tau_norm(5)) '$']};

location_1 = 'northwest';
location_2 = 'southeast';
location_3 = 'southwest';

n_contour = 20;

colours = get(gca,'colororder');
close;                        

% Pictures - Rate x tau x outage probability

for chn_idx = 1:N_CHN
    for alg_idx = 1:N_ALG
        for snr_idx = 1:N_SNR
            figure;
                
            set(gcf,'position',[0 0 800 600]);
            
            if(alg_idx == 1)
                surf(edges_rate_grid,tau_rate_grid,prob_rate_f(:,:,alg_idx)','edgecolor','none');
            else
                surf(edges_rate_grid,tau_icibs_max*tau_rate_grid,prob_rate_f(:,:,alg_idx)','edgecolor','none');
            end
            
            xlabel('Uplink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
            ylabel('$\tau$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
            zlabel('Outage probability','fontname',fontname,'fontsize',fontsize);
            
            colorbar;
            
            set(gca,'fontname',fontname,'fontsize',fontsize);
            
            xlim([0 max_rate]);
            
            if(alg_idx == 1)
                ylim([0 tau_norm_max]);
            else
                ylim([0 tau_icibs_max]);
            end
            
            zlim([0 1]);    
        end
        
        if (savefig == 1)
            saveas(gcf,[root_save_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr(snr_idx)) '_dB'],'fig');
            saveas(gcf,[root_save_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr(snr_idx)) '_dB'],'png');
            saveas(gcf,[root_save_rate channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr(snr_idx)) '_dB'],'epsc2');
        end
    end
end

% Pictures - Contour levels for Rate x tau x outage probability

for chn_idx = 1:N_CHN
    for alg_idx = 1:N_ALG
        for snr_idx = 1:N_SNR
            figure;
                
            set(gcf,'position',[0 0 800 600]);
            
            if(alg_idx == 1)
                contour(edges_rate_grid,tau_rate_grid,prob_rate_f(:,:,alg_idx)',n_contour,'linewidth',linewidth);                   
            else
                contour(edges_rate_grid,tau_icibs_max*tau_rate_grid,prob_rate_f(:,:,alg_idx)',n_contour,'linewidth',linewidth);
            end
            xlabel('Uplink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
            ylabel('$\tau$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
            
            colorbar;
            
            set(gca,'fontname',fontname,'fontsize',fontsize);         
        end
        
        if (savefig == 1)
            saveas(gcf,[root_save_rate 'countour_' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr(snr_idx)) '_dB'],'fig');
            saveas(gcf,[root_save_rate 'countour_' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr(snr_idx)) '_dB'],'png');
            saveas(gcf,[root_save_rate 'countour_' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr(snr_idx)) '_dB'],'epsc2');
        end
    end
end

% Pictures - Number of dropped users x tau

for chn_idx = 1:N_CHN
    for alg_idx = 1:N_ALG
        for snr_idx = 1:N_SNR
            figure;
            
            set(gcf,'position',[0 0 800 600]);
            
            if(alg_idx == 1)
                surf(edges_sele_grid,tau_sele_grid,prob_sele(:,:,alg_idx)','edgecolor','none');
            else
                surf(edges_sele_grid,tau_icibs_max*tau_sele_grid,prob_sele(:,:,alg_idx)','edgecolor','none');        
            end
                       
            view(2);
            
            xlabel('Number of dropped users','fontname',fontname,'fontsize',fontsize);
            ylabel('$\tau$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
            zlabel('CDF','fontname',fontname,'fontsize',fontsize);
            
            colorbar;
                        
            set(gca,'fontname',fontname,'fontsize',fontsize);
            
            xlim([0 L_max-1]);
            
            if(alg_idx == 1)
                ylim([0 tau_norm_max]);
            else
                ylim([0 tau_icibs_max]);
            end

            zlim([0 1]);
        end
        
        if (savefig == 1)
            saveas(gcf,[root_save_drop channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr(snr_idx)) '_dB'],'fig');
            saveas(gcf,[root_save_drop channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr(snr_idx)) '_dB'],'png');
            saveas(gcf,[root_save_drop channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr(snr_idx)) '_dB'],'epsc2');
        end
    end
end

% Picture - Selected users x time

for chn_idx = 1:N_CHN
    for alg_idx = 1:N_ALG
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        plot(time,L(:,end,chn_idx),'o-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
        hold on;
        plot(time,L(:,50,chn_idx),'o-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
        plot(time,L(:,10,chn_idx),'o-','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
        
        xlabel('Time (s)','fontname',fontname,'fontsize',fontsize);
        ylabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
        
        if(alg_idx == 1)
            legend(legend_tau_cbs2,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_3);
            legend box off;
        else
            legend(legend_tau_icibs2,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_3);
            legend box off;
        end
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([0 0.1]);
        ylim([0 L_max-1]);
        
        if (savefig == 1)
            saveas(gcf,[root_save_time 'selected_user/' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_SNR_' num2str(snr(snr_idx)) '_dB'],'fig');
            saveas(gcf,[root_save_time 'selected_user/' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_SNR_' num2str(snr(snr_idx)) '_dB'],'png');
            saveas(gcf,[root_save_time 'selected_user/' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_SNR_' num2str(snr(snr_idx)) '_dB'],'epsc2');
        end
    end
end

% Picture - Rate x time

for chn_idx = 1:N_CHN
    for alg_idx = 1:N_ALG
        for k = 1:K
            figure;
            
            set(gcf,'position',[0 0 800 600]);
            
            plot(time,rate(:,k,end,chn_idx),'o-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
            hold on;
            plot(time,rate(:,k,50,chn_idx),'o-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);   
            plot(time,rate(:,k,10,chn_idx),'o-','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
            plot(time,rate(:,k,1,chn_idx),'o-','color',colours(4,:),'linewidth',linewidth,'markersize',markersize);
         
            xlabel('Time (s)','fontname',fontname,'fontsize',fontsize);
            ylabel('Rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
            
            if(alg_idx == 1)
                legend(legend_tau_cbs,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_1);
                legend box off;
            else
                legend(legend_tau_icibs,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_1);
                legend box off;
            end
            
            set(gca,'fontname',fontname,'fontsize',fontsize);
            
            xlim([0 0.1]);
            ylim([0 max_rate + 7]);
            
            if (savefig == 1)
                saveas(gcf,[root_save_time 'rate/' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_' num2str(k) '_user_SNR_' num2str(snr(snr_idx)) '_dB'],'fig');
                saveas(gcf,[root_save_time 'rate/' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_' num2str(k) '_user_SNR_' num2str(snr(snr_idx)) '_dB'],'png');
                saveas(gcf,[root_save_time 'rate/' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_' num2str(k) '_user_SNR_' num2str(snr(snr_idx)) '_dB'],'epsc2');
            end
        end   
    end
end

% Picture - Power per user x time

for chn_idx = 1:N_CHN
    for alg_idx = 1:N_ALG
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        plot(time,1./L(:,end,chn_idx),'o-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
        hold on;
        plot(time,1./L(:,50,chn_idx),'o-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
        plot(time,1./L(:,10,chn_idx),'o-','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
        
        xlabel('Time (s)','fontname',fontname,'fontsize',fontsize);
        ylabel('Power per active user','fontname',fontname,'fontsize',fontsize);
        
        if(alg_idx == 1)
            legend(legend_tau_cbs2,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_3);
            legend box off;
        else
            legend(legend_tau_icibs2,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_3);
            legend box off;
        end
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([0 0.1]);
        ylim([0 0.1]);
        
        if (savefig == 1)
            saveas(gcf,[root_save_time 'power/' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_SNR_' num2str(snr(snr_idx)) '_dB'],'fig');
            saveas(gcf,[root_save_time 'power/' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_SNR_' num2str(snr(snr_idx)) '_dB'],'png');
            saveas(gcf,[root_save_time 'power/' channel_mod{chn_idx} '_' algorithm{alg_idx} '_M_' num2str(M) '_SNR_' num2str(snr(snr_idx)) '_dB'],'epsc2');
        end
    end
end