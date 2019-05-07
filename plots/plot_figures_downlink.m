clear;
close all;
clc;

MC = 10000;

M = 64;
K = 18;
L = 13;

N_ALG = 4;
N_SNR = 7;

snr = [-20 -15 -10 -5 0 5 10]';

rate_d_ur_los = zeros(K,MC,N_SNR);

rate_d_ur_los_alg = zeros(L,MC,N_ALG,N_SNR);

rate_d_sparse = zeros(K,MC,N_SNR);

rate_d_sparse_alg = zeros(L,MC,N_ALG,N_SNR);

rate_d_rayleigh = zeros(K,MC,N_SNR);

rate_d_rayleigh_alg = zeros(L,MC,N_ALG,N_SNR);

load('../results/rate_mf_ur-los_M_64_K_18_L_13_SNR_-20_dB_MC_10000.mat');

rate_d_ur_los(:,:,1) = rate_d;

rate_d_ur_los_alg(:,:,1,1) = rate_rs_d;
rate_d_ur_los_alg(:,:,2,1) = rate_sos_d;
rate_d_ur_los_alg(:,:,3,1) = rate_cbs_d;
rate_d_ur_los_alg(:,:,4,1) = rate_icibs_d;

load('../results/rate_mf_ur-los_M_64_K_18_L_13_SNR_-15_dB_MC_10000.mat');

rate_d_ur_los(:,:,2) = rate_d;

rate_d_ur_los_alg(:,:,1,2) = rate_rs_d;
rate_d_ur_los_alg(:,:,2,2) = rate_sos_d;
rate_d_ur_los_alg(:,:,3,2) = rate_cbs_d;
rate_d_ur_los_alg(:,:,4,2) = rate_icibs_d;

load('../results/rate_mf_ur-los_M_64_K_18_L_13_SNR_-10_dB_MC_10000.mat');

rate_d_ur_los(:,:,3) = rate_d;

rate_d_ur_los_alg(:,:,1,3) = rate_rs_d;
rate_d_ur_los_alg(:,:,2,3) = rate_sos_d;
rate_d_ur_los_alg(:,:,3,3) = rate_cbs_d;
rate_d_ur_los_alg(:,:,4,3) = rate_icibs_d;

load('../results/rate_mf_ur-los_M_64_K_18_L_13_SNR_-5_dB_MC_10000.mat');

rate_d_ur_los(:,:,4) = rate_d;

rate_d_ur_los_alg(:,:,1,4) = rate_rs_d;
rate_d_ur_los_alg(:,:,2,4) = rate_sos_d;
rate_d_ur_los_alg(:,:,3,4) = rate_cbs_d;
rate_d_ur_los_alg(:,:,4,4) = rate_icibs_d;

load('../results/rate_mf_ur-los_M_64_K_18_L_13_SNR_0_dB_MC_10000.mat');

rate_d_ur_los(:,:,5) = rate_d;

rate_d_ur_los_alg(:,:,1,5) = rate_rs_d;
rate_d_ur_los_alg(:,:,2,5) = rate_sos_d;
rate_d_ur_los_alg(:,:,3,5) = rate_cbs_d;
rate_d_ur_los_alg(:,:,4,5) = rate_icibs_d;

load('../results/rate_mf_ur-los_M_64_K_18_L_13_SNR_5_dB_MC_10000.mat');

rate_d_ur_los(:,:,6) = rate_d;

rate_d_ur_los_alg(:,:,1,6) = rate_rs_d;
rate_d_ur_los_alg(:,:,2,6) = rate_sos_d;
rate_d_ur_los_alg(:,:,3,6) = rate_cbs_d;
rate_d_ur_los_alg(:,:,4,6) = rate_icibs_d;

load('../results/rate_mf_ur-los_M_64_K_18_L_13_SNR_10_dB_MC_10000.mat');

rate_d_ur_los(:,:,7) = rate_d;

rate_d_ur_los_alg(:,:,1,7) = rate_rs_d;
rate_d_ur_los_alg(:,:,2,7) = rate_sos_d;
rate_d_ur_los_alg(:,:,3,7) = rate_cbs_d;
rate_d_ur_los_alg(:,:,4,7) = rate_icibs_d;

load('../results/rate_mf_sparse_M_64_K_18_L_13_SNR_-20_dB_MC_10000.mat');

rate_d_sparse(:,:,1) = rate_d;

rate_d_sparse_alg(:,:,1,1) = rate_rs_d;
rate_d_sparse_alg(:,:,2,1) = rate_sos_d;
rate_d_sparse_alg(:,:,3,1) = rate_cbs_d;
rate_d_sparse_alg(:,:,4,1) = rate_icibs_d;

load('../results/rate_mf_sparse_M_64_K_18_L_13_SNR_-15_dB_MC_10000.mat');

rate_d_sparse(:,:,2) = rate_d;

rate_d_sparse_alg(:,:,1,2) = rate_rs_d;
rate_d_sparse_alg(:,:,2,2) = rate_sos_d;
rate_d_sparse_alg(:,:,3,2) = rate_cbs_d;
rate_d_sparse_alg(:,:,4,2) = rate_icibs_d;

load('../results/rate_mf_sparse_M_64_K_18_L_13_SNR_-10_dB_MC_10000.mat');

rate_d_sparse(:,:,3) = rate_d;

rate_d_sparse_alg(:,:,1,3) = rate_rs_d;
rate_d_sparse_alg(:,:,2,3) = rate_sos_d;
rate_d_sparse_alg(:,:,3,3) = rate_cbs_d;
rate_d_sparse_alg(:,:,4,3) = rate_icibs_d;

load('../results/rate_mf_sparse_M_64_K_18_L_13_SNR_-5_dB_MC_10000.mat');

rate_d_sparse(:,:,4) = rate_d;

rate_d_sparse_alg(:,:,1,4) = rate_rs_d;
rate_d_sparse_alg(:,:,2,4) = rate_sos_d;
rate_d_sparse_alg(:,:,3,4) = rate_cbs_d;
rate_d_sparse_alg(:,:,4,4) = rate_icibs_d;

load('../results/rate_mf_sparse_M_64_K_18_L_13_SNR_0_dB_MC_10000.mat');

rate_d_sparse(:,:,5) = rate_d;

rate_d_sparse_alg(:,:,1,5) = rate_rs_d;
rate_d_sparse_alg(:,:,2,5) = rate_sos_d;
rate_d_sparse_alg(:,:,3,5) = rate_cbs_d;
rate_d_sparse_alg(:,:,4,5) = rate_icibs_d;

load('../results/rate_mf_sparse_M_64_K_18_L_13_SNR_5_dB_MC_10000.mat');

rate_d_sparse(:,:,6) = rate_d;

rate_d_sparse_alg(:,:,1,6) = rate_rs_d;
rate_d_sparse_alg(:,:,2,6) = rate_sos_d;
rate_d_sparse_alg(:,:,3,6) = rate_cbs_d;
rate_d_sparse_alg(:,:,4,6) = rate_icibs_d;

load('../results/rate_mf_sparse_M_64_K_18_L_13_SNR_10_dB_MC_10000.mat');

rate_d_sparse(:,:,7) = rate_d;

rate_d_sparse_alg(:,:,1,7) = rate_rs_d;
rate_d_sparse_alg(:,:,2,7) = rate_sos_d;
rate_d_sparse_alg(:,:,3,7) = rate_cbs_d;
rate_d_sparse_alg(:,:,4,7) = rate_icibs_d;

load('../results/rate_mf_rayleigh_M_64_K_18_L_13_SNR_-20_dB_MC_10000.mat');

rate_d_rayleigh(:,:,1) = rate_d;

rate_d_rayleigh_alg(:,:,1,1) = rate_rs_d;
rate_d_rayleigh_alg(:,:,2,1) = rate_sos_d;
rate_d_rayleigh_alg(:,:,3,1) = rate_cbs_d;
rate_d_rayleigh_alg(:,:,4,1) = rate_icibs_d;

load('../results/rate_mf_rayleigh_M_64_K_18_L_13_SNR_-15_dB_MC_10000.mat');

rate_d_rayleigh(:,:,2) = rate_d;

rate_d_rayleigh_alg(:,:,1,2) = rate_rs_d;
rate_d_rayleigh_alg(:,:,2,2) = rate_sos_d;
rate_d_rayleigh_alg(:,:,3,2) = rate_cbs_d;
rate_d_rayleigh_alg(:,:,4,2) = rate_icibs_d;

load('../results/rate_mf_rayleigh_M_64_K_18_L_13_SNR_-10_dB_MC_10000.mat');

rate_d_rayleigh(:,:,1) = rate_d;

rate_d_rayleigh_alg(:,:,1,3) = rate_rs_d;
rate_d_rayleigh_alg(:,:,2,3) = rate_sos_d;
rate_d_rayleigh_alg(:,:,3,3) = rate_cbs_d;
rate_d_rayleigh_alg(:,:,4,3) = rate_icibs_d;

load('../results/rate_mf_rayleigh_M_64_K_18_L_13_SNR_-5_dB_MC_10000.mat');

rate_d_rayleigh(:,:,4) = rate_d;

rate_d_rayleigh_alg(:,:,1,4) = rate_rs_d;
rate_d_rayleigh_alg(:,:,2,4) = rate_sos_d;
rate_d_rayleigh_alg(:,:,3,4) = rate_cbs_d;
rate_d_rayleigh_alg(:,:,4,4) = rate_icibs_d;

load('../results/rate_mf_rayleigh_M_64_K_18_L_13_SNR_0_dB_MC_10000.mat');

rate_d_rayleigh(:,:,5) = rate_d;

rate_d_rayleigh_alg(:,:,1,5) = rate_rs_d;
rate_d_rayleigh_alg(:,:,2,5) = rate_sos_d;
rate_d_rayleigh_alg(:,:,3,5) = rate_cbs_d;
rate_d_rayleigh_alg(:,:,4,5) = rate_icibs_d;

load('../results/rate_mf_rayleigh_M_64_K_18_L_13_SNR_5_dB_MC_10000.mat');

rate_d_rayleigh(:,:,6) = rate_d;

rate_d_rayleigh_alg(:,:,1,6) = rate_rs_d;
rate_d_rayleigh_alg(:,:,2,6) = rate_sos_d;
rate_d_rayleigh_alg(:,:,3,6) = rate_cbs_d;
rate_d_rayleigh_alg(:,:,4,6) = rate_icibs_d;

load('../results/rate_mf_rayleigh_M_64_K_18_L_13_SNR_10_dB_MC_10000.mat');

rate_d_rayleigh(:,:,7) = rate_d;

rate_d_rayleigh_alg(:,:,1,7) = rate_rs_d;
rate_d_rayleigh_alg(:,:,2,7) = rate_sos_d;
rate_d_rayleigh_alg(:,:,3,7) = rate_cbs_d;
rate_d_rayleigh_alg(:,:,4,7) = rate_icibs_d;

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

BIN_WIDTH_CDF  = 0.005;

BAR_SIZE = 0.8;

% NS - No selection
% RS - Random selection
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_algo = {'NS','RS','SOS','CBS','ICIBS'};
legend_link = {'Uplink','Downlink'};

location = 'northwest';

cat = categorical(legend_algo);
cat = reordercats(cat,legend_algo);

root_rate_fit = '../figures/rate/fit_';
root_erg_rate = '../figures/rate/erg_cap_';
root_out_prob = '../figures/rate/out_prob_';

savefig = 0;

% p_o = 0.05;
% 
% R_sum_d_ns    = 31.46;
% R_sum_d_rs    = 26.16;
% R_sum_d_sos   = 26.19;
% R_sum_d_icibs = 28.16;
% 
% R_d_ns    = 6.285;
% R_d_rs    = 6.535;
% R_d_sos   = 6.540;
% R_d_icibs = 7.035;

colours = get(gca,'colororder');
close;

% p_d = zeros(3,K);
% 
% psi_range = (0:0.01:0.4);
% 
% f_d = zeros(K,length(psi_range));
% 
% for k = 1:1
%     curvefit = fit(psi(k,:)',rate_d(k,:)','poly2');
%     
%     p_d(1,k) = curvefit.p1;
%     p_d(2,k) = curvefit.p2;
%     p_d(3,k) = curvefit.p3;
%         
%     f_d(k,:) = p_d(1,k)*psi_range.^2 + p_d(2,k)*psi_range + p_d(3,k);
%     
%     figure;
%     
%     set(gcf,'position',[0 0 800 600]);
%     
%     plot(psi(k,:),rate_d(k,:),'.','color',colours(1,:),'linewidth',linewidth);
%     hold on;
%     plot(psi_range,f_d(k,:),'-','color',colours(2,:),'linewidth',linewidth);
%     
%     xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
%     ylabel('Downlink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
%     
%     legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
%     
%     set(gca,'fontname',fontname,'fontsize',fontsize);
%     
%     xlim([0 0.4]);
%     
%     saveas(gcf,[root_rate_fit 'downlink_M_' num2str(M) '_' num2str(k) '_user'],'fig');
%     saveas(gcf,[root_rate_fit 'downlink_M_' num2str(M) '_' num2str(k) '_user'],'png');
%     saveas(gcf,[root_rate_fit 'downlink_M_' num2str(M) '_' num2str(k) '_user'],'epsc2');
% end

for snr_idx = 1:N_SNR
    [values_d_ur_los,edges_d_ur_los]             = histcounts(sum(rate_d_ur_los(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_ur_los,edges_rs_d_ur_los]       = histcounts(sum(rate_d_ur_los_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_ur_los,edges_sos_d_ur_los]     = histcounts(sum(rate_d_ur_los_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_ur_los,edges_cbs_d_ur_los]     = histcounts(sum(rate_d_ur_los_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_ur_los,edges_icibs_d_ur_los] = histcounts(sum(rate_d_ur_los_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_d_sparse,edges_d_sparse]             = histcounts(sum(rate_d_sparse(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_sparse,edges_rs_d_sparse]       = histcounts(sum(rate_d_sparse_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_sparse,edges_sos_d_sparse]     = histcounts(sum(rate_d_sparse_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_sparse,edges_cbs_d_sparse]     = histcounts(sum(rate_d_sparse_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_sparse,edges_icibs_d_sparse] = histcounts(sum(rate_d_sparse_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_d_rayleigh,edges_d_rayleigh]             = histcounts(sum(rate_d_rayleigh(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_rayleigh,edges_rs_d_rayleigh]       = histcounts(sum(rate_d_rayleigh_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_rayleigh,edges_sos_d_rayleigh]     = histcounts(sum(rate_d_rayleigh_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_rayleigh,edges_cbs_d_rayleigh]     = histcounts(sum(rate_d_rayleigh_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_rayleigh,edges_icibs_d_rayleigh] = histcounts(sum(rate_d_rayleigh_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot([edges_d_ur_los],[values_d_ur_los 1],'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot([edges_rs_d_ur_los],[values_rs_d_ur_los 1],'-','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_ur_los],[values_sos_d_ur_los 1],'-','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_ur_los],[values_cbs_d_ur_los 1],'-','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_ur_los],[values_icibs_d_ur_los 1],'-','color',colours(5,:),'linewidth',linewidth);
    plot([edges_d_sparse],[values_d_sparse 1],'--','color',colours(1,:),'linewidth',linewidth);
    plot([edges_rs_d_sparse],[values_rs_d_sparse 1],'--','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_sparse],[values_sos_d_sparse 1],'--','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_sparse],[values_cbs_d_sparse 1],'--','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_sparse],[values_icibs_d_sparse 1],'--','color',colours(5,:),'linewidth',linewidth);
    plot([edges_d_rayleigh],[values_d_rayleigh 1],':','color',colours(1,:),'linewidth',linewidth);
    plot([edges_rs_d_rayleigh],[values_rs_d_rayleigh 1],':','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_rayleigh],[values_sos_d_rayleigh 1],':','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_rayleigh],[values_cbs_d_rayleigh 1],':','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_rayleigh],[values_icibs_d_rayleigh 1],':','color',colours(5,:),'linewidth',linewidth);
    
    title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);
    
    xlabel('Downlink sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    % xlim([20 45]);
    ylim([0 1]);
    
    if (savefig == 1)
        saveas(gcf,[root_out_prob 'downlink_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'fig');
        saveas(gcf,[root_out_prob 'downlink_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'png');
        saveas(gcf,[root_out_prob 'downlink_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'epsc2');
    end
end

for snr_idx = 1:N_SNR
    [values_d_ur_los,edges_d_ur_los]             = histcounts(sum(rate_d_ur_los(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_ur_los,edges_rs_d_ur_los]       = histcounts(sum(rate_d_ur_los_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_ur_los,edges_sos_d_ur_los]     = histcounts(sum(rate_d_ur_los_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_ur_los,edges_cbs_d_ur_los]     = histcounts(sum(rate_d_ur_los_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_ur_los,edges_icibs_d_ur_los] = histcounts(sum(rate_d_ur_los_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
        
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot([edges_d_ur_los],[values_d_ur_los 1],'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot([edges_rs_d_ur_los],[values_rs_d_ur_los 1],'-','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_ur_los],[values_sos_d_ur_los 1],'-','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_ur_los],[values_cbs_d_ur_los 1],'-','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_ur_los],[values_icibs_d_ur_los 1],'-','color',colours(5,:),'linewidth',linewidth);
    
    title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);
    
    xlabel('Downlink sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    % xlim([20 45]);
    ylim([0 1]);
    
    if (savefig == 1)
        saveas(gcf,[root_out_prob 'downlink_ur_los_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'fig');
        saveas(gcf,[root_out_prob 'downlink_ur_los_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'png');
        saveas(gcf,[root_out_prob 'downlink_ur_los_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'epsc2');
    end
end

for snr_idx = 1:N_SNR
    [values_d_sparse,edges_d_sparse]             = histcounts(sum(rate_d_sparse(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_sparse,edges_rs_d_sparse]       = histcounts(sum(rate_d_sparse_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_sparse,edges_sos_d_sparse]     = histcounts(sum(rate_d_sparse_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_sparse,edges_cbs_d_sparse]     = histcounts(sum(rate_d_sparse_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_sparse,edges_icibs_d_sparse] = histcounts(sum(rate_d_sparse_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
        
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot([edges_d_sparse],[values_d_sparse 1],'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot([edges_rs_d_sparse],[values_rs_d_sparse 1],'-','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_sparse],[values_sos_d_sparse 1],'-','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_sparse],[values_cbs_d_sparse 1],'-','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_sparse],[values_icibs_d_sparse 1],'-','color',colours(5,:),'linewidth',linewidth);
    
    title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);
    
    xlabel('Downlink sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    % xlim([20 45]);
    ylim([0 1]);
    
    if (savefig == 1)
        saveas(gcf,[root_out_prob 'downlink_sparse_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'fig');
        saveas(gcf,[root_out_prob 'downlink_sparse_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'png');
        saveas(gcf,[root_out_prob 'downlink_sparse_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'epsc2');
    end
end

for snr_idx = 1:N_SNR    
    [values_d_rayleigh,edges_d_rayleigh]             = histcounts(sum(rate_d_rayleigh(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_rayleigh,edges_rs_d_rayleigh]       = histcounts(sum(rate_d_rayleigh_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_rayleigh,edges_sos_d_rayleigh]     = histcounts(sum(rate_d_rayleigh_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_rayleigh,edges_cbs_d_rayleigh]     = histcounts(sum(rate_d_rayleigh_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_rayleigh,edges_icibs_d_rayleigh] = histcounts(sum(rate_d_rayleigh_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot([edges_d_rayleigh],[values_d_rayleigh 1],'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot([edges_rs_d_rayleigh],[values_rs_d_rayleigh 1],'-','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_rayleigh],[values_sos_d_rayleigh 1],'-','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_rayleigh],[values_cbs_d_rayleigh 1],'-','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_rayleigh],[values_icibs_d_rayleigh 1],'-','color',colours(5,:),'linewidth',linewidth);
    
    title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);
    
    xlabel('Downlink sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    % xlim([20 45]);
    ylim([0 1]);
    
    if (savefig == 1)
        saveas(gcf,[root_out_prob 'downlink_rayleigh_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'fig');
        saveas(gcf,[root_out_prob 'downlink_rayleigh_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'png');
        saveas(gcf,[root_out_prob 'downlink_rayleigh_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'epsc2');
    end
end

for snr_idx = 1:N_SNR
    [values_d_ur_los,edges_d_ur_los]             = histcounts(mean(rate_d_ur_los(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_ur_los,edges_rs_d_ur_los]       = histcounts(mean(rate_d_ur_los_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_ur_los,edges_sos_d_ur_los]     = histcounts(mean(rate_d_ur_los_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_ur_los,edges_cbs_d_ur_los]     = histcounts(mean(rate_d_ur_los_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_ur_los,edges_icibs_d_ur_los] = histcounts(mean(rate_d_ur_los_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_d_sparse,edges_d_sparse]             = histcounts(mean(rate_d_sparse(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_sparse,edges_rs_d_sparse]       = histcounts(mean(rate_d_sparse_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_sparse,edges_sos_d_sparse]     = histcounts(mean(rate_d_sparse_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_sparse,edges_cbs_d_sparse]     = histcounts(mean(rate_d_sparse_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_sparse,edges_icibs_d_sparse] = histcounts(mean(rate_d_sparse_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_d_rayleigh,edges_d_rayleigh]             = histcounts(mean(rate_d_rayleigh(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_rayleigh,edges_rs_d_rayleigh]       = histcounts(mean(rate_d_rayleigh_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_rayleigh,edges_sos_d_rayleigh]     = histcounts(mean(rate_d_rayleigh_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_rayleigh,edges_cbs_d_rayleigh]     = histcounts(mean(rate_d_rayleigh_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_rayleigh,edges_icibs_d_rayleigh] = histcounts(mean(rate_d_rayleigh_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot([edges_d_ur_los],[values_d_ur_los 1],'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot([edges_rs_d_ur_los],[values_rs_d_ur_los 1],'-','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_ur_los],[values_sos_d_ur_los 1],'-','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_ur_los],[values_cbs_d_ur_los 1],'-','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_ur_los],[values_icibs_d_ur_los 1],'-','color',colours(5,:),'linewidth',linewidth);
    plot([edges_d_sparse],[values_d_sparse 1],'--','color',colours(1,:),'linewidth',linewidth);
    plot([edges_rs_d_sparse],[values_rs_d_sparse 1],'--','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_sparse],[values_sos_d_sparse 1],'--','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_sparse],[values_cbs_d_sparse 1],'--','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_sparse],[values_icibs_d_sparse 1],'--','color',colours(5,:),'linewidth',linewidth);
    plot([edges_d_rayleigh],[values_d_rayleigh 1],':','color',colours(1,:),'linewidth',linewidth);
    plot([edges_rs_d_rayleigh],[values_rs_d_rayleigh 1],':','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_rayleigh],[values_sos_d_rayleigh 1],':','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_rayleigh],[values_cbs_d_rayleigh 1],':','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_rayleigh],[values_icibs_d_rayleigh 1],':','color',colours(5,:),'linewidth',linewidth);
    
    title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);
    
    xlabel('Downlink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    % xlim([5 9.5]);
    ylim([0 1]);
    
    if (savefig == 1)
        saveas(gcf,[root_out_prob 'downlink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'fig');
        saveas(gcf,[root_out_prob 'downlink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'png');
        saveas(gcf,[root_out_prob 'downlink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'epsc2');
    end
end

for snr_idx = 1:N_SNR
    [values_d_ur_los,edges_d_ur_los]             = histcounts(mean(rate_d_ur_los(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_ur_los,edges_rs_d_ur_los]       = histcounts(mean(rate_d_ur_los_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_ur_los,edges_sos_d_ur_los]     = histcounts(mean(rate_d_ur_los_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_ur_los,edges_cbs_d_ur_los]     = histcounts(mean(rate_d_ur_los_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_ur_los,edges_icibs_d_ur_los] = histcounts(mean(rate_d_ur_los_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot([edges_d_ur_los],[values_d_ur_los 1],'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot([edges_rs_d_ur_los],[values_rs_d_ur_los 1],'-','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_ur_los],[values_sos_d_ur_los 1],'-','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_ur_los],[values_cbs_d_ur_los 1],'-','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_ur_los],[values_icibs_d_ur_los 1],'-','color',colours(5,:),'linewidth',linewidth);
    
    title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);
    
    xlabel('Downlink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    % xlim([5 9.5]);
    ylim([0 1]);
    
    if (savefig == 1)
        saveas(gcf,[root_out_prob 'downlink_ur_los_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'fig');
        saveas(gcf,[root_out_prob 'downlink_ur_los_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'png');
        saveas(gcf,[root_out_prob 'downlink_ur_los_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'epsc2');
    end
end

for snr_idx = 1:N_SNR
    [values_d_sparse,edges_d_sparse]             = histcounts(mean(rate_d_sparse(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_sparse,edges_rs_d_sparse]       = histcounts(mean(rate_d_sparse_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_sparse,edges_sos_d_sparse]     = histcounts(mean(rate_d_sparse_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_sparse,edges_cbs_d_sparse]     = histcounts(mean(rate_d_sparse_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_sparse,edges_icibs_d_sparse] = histcounts(mean(rate_d_sparse_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
        
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot([edges_d_sparse],[values_d_sparse 1],'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot([edges_rs_d_sparse],[values_rs_d_sparse 1],'-','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_sparse],[values_sos_d_sparse 1],'-','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_sparse],[values_cbs_d_sparse 1],'-','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_sparse],[values_icibs_d_sparse 1],'-','color',colours(5,:),'linewidth',linewidth);
    
    title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);
    
    xlabel('Downlink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    % xlim([5 9.5]);
    ylim([0 1]);
    
    if (savefig == 1)
        saveas(gcf,[root_out_prob 'downlink_sparse_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'fig');
        saveas(gcf,[root_out_prob 'downlink_sparse_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'png');
        saveas(gcf,[root_out_prob 'downlink_sparse_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'epsc2');
    end
end

for snr_idx = 1:N_SNR
    [values_d_rayleigh,edges_d_rayleigh]             = histcounts(mean(rate_d_rayleigh(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    [values_rs_d_rayleigh,edges_rs_d_rayleigh]       = histcounts(mean(rate_d_rayleigh_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_sos_d_rayleigh,edges_sos_d_rayleigh]     = histcounts(mean(rate_d_rayleigh_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_cbs_d_rayleigh,edges_cbs_d_rayleigh]     = histcounts(mean(rate_d_rayleigh_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    [values_icibs_d_rayleigh,edges_icibs_d_rayleigh] = histcounts(mean(rate_d_rayleigh_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot([edges_d_rayleigh],[values_d_rayleigh 1],'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot([edges_rs_d_rayleigh],[values_rs_d_rayleigh 1],'-','color',colours(2,:),'linewidth',linewidth);
    plot([edges_sos_d_rayleigh],[values_sos_d_rayleigh 1],'-','color',colours(3,:),'linewidth',linewidth);
    plot([edges_cbs_d_rayleigh],[values_cbs_d_rayleigh 1],'-','color',colours(4,:),'linewidth',linewidth);
    plot([edges_icibs_d_rayleigh],[values_icibs_d_rayleigh 1],'-','color',colours(5,:),'linewidth',linewidth);
    
    title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);
    
    xlabel('Downlink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    % xlim([5 9.5]);
    ylim([0 1]);
    
    if (savefig == 1)
        saveas(gcf,[root_out_prob 'downlink_rayleigh_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'fig');
        saveas(gcf,[root_out_prob 'downlink_rayleigh_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'png');
        saveas(gcf,[root_out_prob 'downlink_rayleigh_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' snr(snr_idx)],'epsc2');
    end
end