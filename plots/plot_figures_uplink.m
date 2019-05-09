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

rate_u_ur_los = zeros(K,MC,N_SNR);
psi_ur_los    = zeros(K,MC,N_SNR);

rate_u_ur_los_alg = zeros(L,MC,N_ALG,N_SNR);
psi_ur_los_alg    = zeros(L,MC,N_ALG,N_SNR);

rate_u_sparse = zeros(K,MC,N_SNR);
psi_sparse    = zeros(K,MC,N_SNR);

rate_u_sparse_alg = zeros(L,MC,N_ALG,N_SNR);
psi_sparse_alg    = zeros(L,MC,N_ALG,N_SNR);

rate_u_rayleigh = zeros(K,MC,N_SNR);
psi_rayleigh    = zeros(K,MC,N_SNR);

rate_u_rayleigh_alg = zeros(L,MC,N_ALG,N_SNR);
psi_rayleigh_alg    = zeros(L,MC,N_ALG,N_SNR);

load('../results/uplink/rate_uplink_mf_ur-los_M_64_K_18_L_13_SNR_-20_dB_MC_10000.mat');

rate_u_ur_los(:,:,1) = rate_u;
psi_ur_los(:,:,1)  = psi;

rate_u_ur_los_alg(:,:,:,1) = rate_u_alg;
psi_ur_los_alg(:,:,:,1)  = psi_alg;

load('../results/uplink/rate_uplink_mf_ur-los_M_64_K_18_L_13_SNR_-15_dB_MC_10000.mat');

rate_u_ur_los(:,:,2) = rate_u;
psi_ur_los(:,:,2)  = psi;

rate_u_ur_los_alg(:,:,:,2) = rate_u_alg;
psi_ur_los_alg(:,:,:,2)  = psi_alg;

load('../results/uplink/rate_uplink_mf_ur-los_M_64_K_18_L_13_SNR_-10_dB_MC_10000.mat');

rate_u_ur_los(:,:,3) = rate_u;
psi_ur_los(:,:,3)  = psi;

rate_u_ur_los_alg(:,:,:,3) = rate_u_alg;
psi_ur_los_alg(:,:,:,3)  = psi_alg;

load('../results/uplink/rate_uplink_mf_ur-los_M_64_K_18_L_13_SNR_-5_dB_MC_10000.mat');

rate_u_ur_los(:,:,4) = rate_u;
psi_ur_los(:,:,4)  = psi;

rate_u_ur_los_alg(:,:,:,4) = rate_u_alg;
psi_ur_los_alg(:,:,:,4)  = psi_alg;

load('../results/uplink/rate_uplink_mf_ur-los_M_64_K_18_L_13_SNR_0_dB_MC_10000.mat');

rate_u_ur_los(:,:,5) = rate_u;
psi_ur_los(:,:,5)  = psi;

rate_u_ur_los_alg(:,:,:,5) = rate_u_alg;
psi_ur_los_alg(:,:,:,5)  = psi_alg;

load('../results/uplink/rate_uplink_mf_ur-los_M_64_K_18_L_13_SNR_5_dB_MC_10000.mat');

rate_u_ur_los(:,:,6) = rate_u;
psi_ur_los(:,:,6)  = psi;

rate_u_ur_los_alg(:,:,:,6) = rate_u_alg;
psi_ur_los_alg(:,:,:,6)  = psi_alg;

load('../results/uplink/rate_uplink_mf_ur-los_M_64_K_18_L_13_SNR_10_dB_MC_10000.mat');

rate_u_ur_los(:,:,7) = rate_u;
psi_ur_los(:,:,7)  = psi;

rate_u_ur_los_alg(:,:,:,7) = rate_u_alg;
psi_ur_los_alg(:,:,:,7)  = psi_alg;

load('../results/uplink/rate_uplink_mf_sparse_M_64_K_18_L_13_SNR_-20_dB_MC_10000.mat');

rate_u_sparse(:,:,1) = rate_u;
psi_sparse(:,:,1)  = psi;

rate_u_sparse_alg(:,:,:,1) = rate_u_alg;
psi_sparse_alg(:,:,:,1)  = psi_alg;

load('../results/uplink/rate_uplink_mf_sparse_M_64_K_18_L_13_SNR_-15_dB_MC_10000.mat');

rate_u_sparse(:,:,2) = rate_u;
psi_sparse(:,:,2)  = psi;

rate_u_sparse_alg(:,:,:,2) = rate_u_alg;
psi_sparse_alg(:,:,:,2)  = psi_alg;

load('../results/uplink/rate_uplink_mf_sparse_M_64_K_18_L_13_SNR_-10_dB_MC_10000.mat');

rate_u_sparse(:,:,3) = rate_u;
psi_sparse(:,:,3)  = psi;

rate_u_sparse_alg(:,:,:,3) = rate_u_alg;
psi_sparse_alg(:,:,:,3)  = psi_alg;

load('../results/uplink/rate_uplink_mf_sparse_M_64_K_18_L_13_SNR_-5_dB_MC_10000.mat');

rate_u_sparse(:,:,4) = rate_u;
psi_sparse(:,:,4)  = psi;

rate_u_sparse_alg(:,:,:,4) = rate_u_alg;
psi_sparse_alg(:,:,:,4)  = psi_alg;

load('../results/uplink/rate_uplink_mf_sparse_M_64_K_18_L_13_SNR_0_dB_MC_10000.mat');

rate_u_sparse(:,:,5) = rate_u;
psi_sparse(:,:,5)  = psi;

rate_u_sparse_alg(:,:,:,5) = rate_u_alg;
psi_sparse_alg(:,:,:,5)  = psi_alg;

load('../results/uplink/rate_uplink_mf_sparse_M_64_K_18_L_13_SNR_5_dB_MC_10000.mat');

rate_u_sparse(:,:,6) = rate_u;
psi_sparse(:,:,6)  = psi;

rate_u_sparse_alg(:,:,:,6) = rate_u_alg;
psi_sparse_alg(:,:,:,6)  = psi_alg;

load('../results/uplink/rate_uplink_mf_sparse_M_64_K_18_L_13_SNR_10_dB_MC_10000.mat');

rate_u_sparse(:,:,7) = rate_u;
psi_sparse(:,:,7)  = psi;

rate_u_sparse_alg(:,:,:,7) = rate_u_alg;
psi_sparse_alg(:,:,:,7)  = psi_alg;

load('../results/uplink/rate_uplink_mf_rayleigh_M_64_K_18_L_13_SNR_-20_dB_MC_10000.mat');

rate_u_rayleigh(:,:,1) = rate_u;
psi_rayleigh(:,:,1)  = psi;

rate_u_rayleigh_alg(:,:,:,1) = rate_u_alg;
psi_rayleigh_alg(:,:,:,1)  = psi_alg;

load('../results/uplink/rate_uplink_mf_rayleigh_M_64_K_18_L_13_SNR_-15_dB_MC_10000.mat');

rate_u_rayleigh(:,:,2) = rate_u;
psi_rayleigh(:,:,2)  = psi;

rate_u_rayleigh_alg(:,:,:,2) = rate_u_alg;
psi_rayleigh_alg(:,:,:,2)  = psi_alg;

load('../results/uplink/rate_uplink_mf_rayleigh_M_64_K_18_L_13_SNR_-10_dB_MC_10000.mat');

rate_u_rayleigh(:,:,3) = rate_u;
psi_rayleigh(:,:,3)  = psi;

rate_u_rayleigh_alg(:,:,:,3) = rate_u_alg;
psi_rayleigh_alg(:,:,:,3)  = psi_alg;

load('../results/uplink/rate_uplink_mf_rayleigh_M_64_K_18_L_13_SNR_-5_dB_MC_10000.mat');

rate_u_rayleigh(:,:,4) = rate_u;
psi_rayleigh(:,:,4)  = psi;

rate_u_rayleigh_alg(:,:,:,4) = rate_u_alg;
psi_rayleigh_alg(:,:,:,4)  = psi_alg;

load('../results/uplink/rate_uplink_mf_rayleigh_M_64_K_18_L_13_SNR_0_dB_MC_10000.mat');

rate_u_rayleigh(:,:,5) = rate_u;
psi_rayleigh(:,:,5)  = psi;

rate_u_rayleigh_alg(:,:,:,5) = rate_u_alg;
psi_rayleigh_alg(:,:,:,5)  = psi_alg;

load('../results/uplink/rate_uplink_mf_rayleigh_M_64_K_18_L_13_SNR_5_dB_MC_10000.mat');

rate_u_rayleigh(:,:,6) = rate_u;
psi_rayleigh(:,:,6)  = psi;

rate_u_rayleigh_alg(:,:,:,6) = rate_u_alg;
psi_rayleigh_alg(:,:,:,6)  = psi_alg;

load('../results/uplink/rate_uplink_mf_rayleigh_M_64_K_18_L_13_SNR_10_dB_MC_10000.mat');

rate_u_rayleigh(:,:,7) = rate_u;
psi_rayleigh(:,:,7)  = psi;

rate_u_rayleigh_alg(:,:,:,7) = rate_u_alg;
psi_rayleigh_alg(:,:,:,7)  = psi_alg;

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

BIN_WIDTH_CDF  = 0.0005;

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
root_out_prob = '../figures/rate/out_prob_';

savefig = 0;

% p_o = 0.05;
% 
% R_sum_u_ns    = 27.37;
% R_sum_u_rs    = 22.43;
% R_sum_u_sos   = 22.43;
% R_sum_u_icibs = 23.39;
%  
% R_u_ns    = 5.470;
% R_u_rs    = 5.605;
% R_u_sos   = 5.605;
% R_u_icibs = 5.845;

colours = get(gca,'colororder');
close;

p_u_ur_los   = zeros(10,K,N_SNR);
p_u_sparse   = zeros(10,K,N_SNR);
p_u_rayleigh = zeros(10,K,N_SNR);

psi_ur_los_min = 0;
psi_ur_los_max = 0.45;

psi_sparse_min = 0;
psi_sparse_max = 0.7;

psi_rayleigh_min = 0;
psi_rayleigh_max = 0.2;

psi_step = 0.01;

psi_range_ur_los   = (psi_ur_los_min:psi_step:psi_ur_los_max);
psi_range_sparse   = (psi_sparse_min:psi_step:psi_sparse_max);
psi_range_rayleigh = (psi_rayleigh_min:psi_step:psi_rayleigh_max);

f_u_ur_los   = zeros(K,length(psi_range_ur_los),N_SNR);
f_u_sparse   = zeros(K,length(psi_range_sparse),N_SNR);
f_u_rayleigh = zeros(K,length(psi_range_rayleigh),N_SNR);

for snr_idx = 1:N_SNR
    for k = 1:1
        curvefit_ur_los   = fit(psi_ur_los(k,:,snr_idx)',rate_u_ur_los(k,:,snr_idx)','exp2');
        curvefit_sparse   = fit(psi_sparse(k,:,snr_idx)',rate_u_sparse(k,:,snr_idx)','exp2');
        curvefit_rayleigh = fit(psi_rayleigh(k,:,snr_idx)',rate_u_rayleigh(k,:,snr_idx)','poly3');
        
        % p_u_ur_los(1,k,snr_idx)  = curvefit_ur_los.p1;
        % p_u_ur_los(2,k,snr_idx)  = curvefit_ur_los.p2;
        % p_u_ur_los(3,k,snr_idx)  = curvefit_ur_los.p3;
        % p_u_ur_los(4,k,snr_idx)  = curvefit_ur_los.p4;
        % p_u_ur_los(5,k,snr_idx)  = curvefit_ur_los.p5;
        % p_u_ur_los(6,k,snr_idx)  = curvefit_ur_los.p6;
        % p_u_ur_los(7,k,snr_idx)  = curvefit_ur_los.p7;
        % p_u_ur_los(8,k,snr_idx)  = curvefit_ur_los.p8;
        % p_u_ur_los(9,k,snr_idx)  = curvefit_ur_los.p9;
        % p_u_ur_los(10,k,snr_idx) = curvefit_ur_los.p10;
                                                        
        p_u_ur_los(1,k,snr_idx) = curvefit_ur_los.a;
        p_u_ur_los(2,k,snr_idx) = curvefit_ur_los.b;
        p_u_ur_los(3,k,snr_idx) = curvefit_ur_los.c;
        p_u_ur_los(4,k,snr_idx) = curvefit_ur_los.d;
        
        % p_u_sparse(1,k,snr_idx) = curvefit_sparse.p1;
        % p_u_sparse(2,k,snr_idx) = curvefit_sparse.p2;
        % p_u_sparse(3,k,snr_idx) = curvefit_sparse.p3;
        % p_u_sparse(4,k,snr_idx) = curvefit_sparse.p4;
        
        p_u_sparse(1,k,snr_idx) = curvefit_sparse.a;
        p_u_sparse(2,k,snr_idx) = curvefit_sparse.b;
        p_u_sparse(3,k,snr_idx) = curvefit_sparse.c;
        p_u_sparse(4,k,snr_idx) = curvefit_sparse.d;
        
        p_u_rayleigh(1,k,snr_idx) = curvefit_rayleigh.p1;
        p_u_rayleigh(2,k,snr_idx) = curvefit_rayleigh.p2;
        p_u_rayleigh(3,k,snr_idx) = curvefit_rayleigh.p3;
        p_u_rayleigh(4,k,snr_idx) = curvefit_rayleigh.p4;
        
        % f_u_ur_los(k,:,snr_idx)   = p_u_ur_los(1,k,snr_idx)*psi_range_ur_los.^9 + p_u_ur_los(2,k,snr_idx)*psi_range_ur_los.^8 + ...
        %                             p_u_ur_los(3,k,snr_idx)*psi_range_ur_los.^7 + p_u_ur_los(4,k,snr_idx)*psi_range_ur_los.^6 + ...
        %                             p_u_ur_los(5,k,snr_idx)*psi_range_ur_los.^5 + p_u_ur_los(6,k,snr_idx)*psi_range_ur_los.^4 + ...
        %                             p_u_ur_los(7,k,snr_idx)*psi_range_ur_los.^3 + p_u_ur_los(8,k,snr_idx)*psi_range_ur_los.^2 + ...
        %                             p_u_ur_los(9,k,snr_idx)*psi_range_ur_los + p_u_ur_los(10,k,snr_idx);
        f_u_ur_los(k,:,snr_idx) = p_u_ur_los(1,k,snr_idx)*exp(p_u_ur_los(2,k,snr_idx)*psi_range_ur_los) + ...
                                  p_u_ur_los(3,k,snr_idx)*exp(p_u_ur_los(4,k,snr_idx)*psi_range_ur_los);
        % f_u_sparse(k,:,snr_idx)   = p_u_sparse(1,k,snr_idx)*psi_range_sparse.^3 + p_u_sparse(2,k,snr_idx)*psi_range_sparse.^2 + ...
        %                             p_u_sparse(3,k,snr_idx)*psi_range_sparse + p_u_sparse(4,k,snr_idx);
        f_u_sparse(k,:,snr_idx)   = p_u_sparse(1,k,snr_idx)*exp(p_u_sparse(2,k,snr_idx)*psi_range_sparse) + ...
                                    p_u_sparse(3,k,snr_idx)*exp(p_u_sparse(4,k,snr_idx)*psi_range_sparse);
        f_u_rayleigh(k,:,snr_idx) = p_u_rayleigh(1,k,snr_idx)*psi_range_rayleigh.^3 + p_u_rayleigh(2,k,snr_idx)*psi_range_rayleigh.^2 + ...
                                    p_u_rayleigh(3,k,snr_idx)*psi_range_rayleigh + p_u_rayleigh(4,k,snr_idx);
        
        figure;
        
        set(gcf,'position',[500 250 1200 600]);
        
        subplot(1,3,1);
        
        plot(psi_ur_los(k,:,snr_idx),rate_u_ur_los(k,:,snr_idx),'.','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(psi_range_ur_los,f_u_ur_los(k,:,snr_idx),'-','color',colours(2,:),'linewidth',linewidth);
        
        xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
        ylabel('Uplink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        
        legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([psi_ur_los_min psi_ur_los_max]);
        
        subplot(1,3,2);
        
        plot(psi_sparse(k,:,snr_idx),rate_u_sparse(k,:,snr_idx),'.','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(psi_range_sparse,f_u_sparse(k,:,snr_idx),'-','color',colours(2,:),'linewidth',linewidth);
        
        xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
        ylabel('Uplink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        
        legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([psi_sparse_min psi_sparse_max]);
        
        subplot(1,3,3)
        
        plot(psi_rayleigh(k,:,snr_idx),rate_u_rayleigh(k,:,snr_idx),'.','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(psi_range_rayleigh,f_u_rayleigh(k,:,snr_idx),'-','color',colours(2,:),'linewidth',linewidth);
        
        xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
        ylabel('Uplink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        
        legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([psi_rayleigh_min psi_rayleigh_max]);
        
        if(savefig == 1)
            saveas(gcf,[root_rate_fit 'M_' num2str(M) '_' num2str(k) '_user'],'fig');
            saveas(gcf,[root_rate_fit 'M_' num2str(M) '_' num2str(k) '_user'],'png');
            saveas(gcf,[root_rate_fit 'M_' num2str(M) '_' num2str(k) '_user'],'epsc2');
        end
    end
end

% for snr_idx = 1:N_SNR  
%     [values_u_ur_los,edges_u_ur_los]             = histcounts(mean(rate_u_ur_los(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     
%     [values_rs_u_ur_los,edges_rs_u_ur_los]       = histcounts(mean(rate_u_ur_los_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_sos_u_ur_los,edges_sos_u_ur_los]     = histcounts(mean(rate_u_ur_los_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_cbs_u_ur_los,edges_cbs_u_ur_los]     = histcounts(mean(rate_u_ur_los_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_icibs_u_ur_los,edges_icibs_u_ur_los] = histcounts(mean(rate_u_ur_los_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     
%     [values_u_sparse,edges_u_sparse]             = histcounts(mean(rate_u_sparse(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     
%     [values_rs_u_sparse,edges_rs_u_sparse]       = histcounts(mean(rate_u_sparse_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_sos_u_sparse,edges_sos_u_sparse]     = histcounts(mean(rate_u_sparse_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_cbs_u_sparse,edges_cbs_u_sparse]     = histcounts(mean(rate_u_sparse_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_icibs_u_sparse,edges_icibs_u_sparse] = histcounts(mean(rate_u_sparse_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     
%     [values_u_rayleigh,edges_u_rayleigh]             = histcounts(mean(rate_u_rayleigh(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     
%     [values_rs_u_rayleigh,edges_rs_u_rayleigh]       = histcounts(mean(rate_u_rayleigh_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_sos_u_rayleigh,edges_sos_u_rayleigh]     = histcounts(mean(rate_u_rayleigh_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_cbs_u_rayleigh,edges_cbs_u_rayleigh]     = histcounts(mean(rate_u_rayleigh_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_icibs_u_rayleigh,edges_icibs_u_rayleigh] = histcounts(mean(rate_u_rayleigh_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     
%     figure;
%     
%     set(gcf,'position',[0 0 800 600]);
%     
%     plot([edges_u_ur_los],[values_u_ur_los 1],'-','color',colours(1,:),'linewidth',linewidth);
%     hold on;
%     plot([edges_rs_u_ur_los],[values_rs_u_ur_los 1],'-','color',colours(2,:),'linewidth',linewidth);
%     plot([edges_sos_u_ur_los],[values_sos_u_ur_los 1],'-','color',colours(3,:),'linewidth',linewidth);
%     plot([edges_cbs_u_ur_los],[values_cbs_u_ur_los 1],'-','color',colours(4,:),'linewidth',linewidth);
%     plot([edges_icibs_u_ur_los],[values_icibs_u_ur_los 1],'-','color',colours(5,:),'linewidth',linewidth);
%     plot([edges_u_sparse],[values_u_sparse 1],'--','color',colours(1,:),'linewidth',linewidth);
%     plot([edges_rs_u_sparse],[values_rs_u_sparse 1],'--','color',colours(2,:),'linewidth',linewidth);
%     plot([edges_sos_u_sparse],[values_sos_u_sparse 1],'--','color',colours(3,:),'linewidth',linewidth);
%     plot([edges_cbs_u_sparse],[values_cbs_u_sparse 1],'--','color',colours(4,:),'linewidth',linewidth);
%     plot([edges_icibs_u_sparse],[values_icibs_u_sparse 1],'--','color',colours(5,:),'linewidth',linewidth);
%     plot([edges_u_rayleigh],[values_u_rayleigh 1],':','color',colours(1,:),'linewidth',linewidth);
%     plot([edges_rs_u_rayleigh],[values_rs_u_rayleigh 1],':','color',colours(2,:),'linewidth',linewidth);
%     plot([edges_sos_u_rayleigh],[values_sos_u_rayleigh 1],':','color',colours(3,:),'linewidth',linewidth);
%     plot([edges_cbs_u_rayleigh],[values_cbs_u_rayleigh 1],':','color',colours(4,:),'linewidth',linewidth);
%     plot([edges_icibs_u_rayleigh],[values_icibs_u_rayleigh 1],':','color',colours(5,:),'linewidth',linewidth);
%     
%     title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);    
%     xlabel('Uplink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
%     ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
%     
%     legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
%     
%     set(gca,'fontname',fontname,'fontsize',fontsize);
%     
%     % xlim([5 7]);
%     ylim([0 1]);
%     
%     if (savefig == 1)
%         saveas(gcf,[root_out_prob 'uplink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'fig');
%         saveas(gcf,[root_out_prob 'uplink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'png');
%         saveas(gcf,[root_out_prob 'uplink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'epsc2');
%     end
% end
% 
% for snr_idx = 1:N_SNR  
%     [values_u_ur_los,edges_u_ur_los]             = histcounts(mean(rate_u_ur_los(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     
%     [values_rs_u_ur_los,edges_rs_u_ur_los]       = histcounts(mean(rate_u_ur_los_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_sos_u_ur_los,edges_sos_u_ur_los]     = histcounts(mean(rate_u_ur_los_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_cbs_u_ur_los,edges_cbs_u_ur_los]     = histcounts(mean(rate_u_ur_los_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_icibs_u_ur_los,edges_icibs_u_ur_los] = histcounts(mean(rate_u_ur_los_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%         
%     figure;
%     
%     set(gcf,'position',[0 0 800 600]);
%     
%     plot([edges_u_ur_los],[values_u_ur_los 1],'-','color',colours(1,:),'linewidth',linewidth);
%     hold on;
%     plot([edges_rs_u_ur_los],[values_rs_u_ur_los 1],'-','color',colours(2,:),'linewidth',linewidth);
%     plot([edges_sos_u_ur_los],[values_sos_u_ur_los 1],'-','color',colours(3,:),'linewidth',linewidth);
%     plot([edges_cbs_u_ur_los],[values_cbs_u_ur_los 1],'-','color',colours(4,:),'linewidth',linewidth);
%     plot([edges_icibs_u_ur_los],[values_icibs_u_ur_los 1],'-','color',colours(5,:),'linewidth',linewidth);
%     
%     title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);    
%     xlabel('Uplink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
%     ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
%     
%     legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
%     
%     set(gca,'fontname',fontname,'fontsize',fontsize);
%     
%     % xlim([5 7]);
%     ylim([0 1]);
%     
%     if (savefig == 1)
%         saveas(gcf,[root_out_prob 'uplink_ur_los_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'fig');
%         saveas(gcf,[root_out_prob 'uplink_ur_los_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'png');
%         saveas(gcf,[root_out_prob 'uplink_ur_los_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'epsc2');
%     end
% end
% 
% for snr_idx = 1:N_SNR      
%     [values_u_sparse,edges_u_sparse]             = histcounts(mean(rate_u_sparse(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     
%     [values_rs_u_sparse,edges_rs_u_sparse]       = histcounts(mean(rate_u_sparse_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_sos_u_sparse,edges_sos_u_sparse]     = histcounts(mean(rate_u_sparse_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_cbs_u_sparse,edges_cbs_u_sparse]     = histcounts(mean(rate_u_sparse_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_icibs_u_sparse,edges_icibs_u_sparse] = histcounts(mean(rate_u_sparse_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%         
%     figure;
%     
%     set(gcf,'position',[0 0 800 600]);
%     
%     plot([edges_u_sparse],[values_u_sparse 1],'-','color',colours(1,:),'linewidth',linewidth);
%     hold on;
%     plot([edges_rs_u_sparse],[values_rs_u_sparse 1],'-','color',colours(2,:),'linewidth',linewidth);
%     plot([edges_sos_u_sparse],[values_sos_u_sparse 1],'-','color',colours(3,:),'linewidth',linewidth);
%     plot([edges_cbs_u_sparse],[values_cbs_u_sparse 1],'-','color',colours(4,:),'linewidth',linewidth);
%     plot([edges_icibs_u_sparse],[values_icibs_u_sparse 1],'-','color',colours(5,:),'linewidth',linewidth);
%     
%     title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);    
%     xlabel('Uplink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
%     ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
%     
%     legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
%     
%     set(gca,'fontname',fontname,'fontsize',fontsize);
%     
%     % xlim([5 7]);
%     ylim([0 1]);
%     
%     if (savefig == 1)
%         saveas(gcf,[root_out_prob 'uplink_sparse_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'fig');
%         saveas(gcf,[root_out_prob 'uplink_sparse_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'png');
%         saveas(gcf,[root_out_prob 'uplink_sparse_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'epsc2');
%     end
% end
% 
% for snr_idx = 1:N_SNR     
%     [values_u_rayleigh,edges_u_rayleigh]             = histcounts(mean(rate_u_rayleigh(:,:,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     
%     [values_rs_u_rayleigh,edges_rs_u_rayleigh]       = histcounts(mean(rate_u_rayleigh_alg(:,:,1,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_sos_u_rayleigh,edges_sos_u_rayleigh]     = histcounts(mean(rate_u_rayleigh_alg(:,:,2,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_cbs_u_rayleigh,edges_cbs_u_rayleigh]     = histcounts(mean(rate_u_rayleigh_alg(:,:,3,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     [values_icibs_u_rayleigh,edges_icibs_u_rayleigh] = histcounts(mean(rate_u_rayleigh_alg(:,:,4,snr_idx)),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
%     
%     figure;
%     
%     set(gcf,'position',[0 0 800 600]);
%     
%     plot([edges_u_rayleigh],[values_u_rayleigh 1],'-','color',colours(1,:),'linewidth',linewidth);
%     hold on;
%     plot([edges_rs_u_rayleigh],[values_rs_u_rayleigh 1],'-','color',colours(2,:),'linewidth',linewidth);
%     plot([edges_sos_u_rayleigh],[values_sos_u_rayleigh 1],'-','color',colours(3,:),'linewidth',linewidth);
%     plot([edges_cbs_u_rayleigh],[values_cbs_u_rayleigh 1],'-','color',colours(4,:),'linewidth',linewidth);
%     plot([edges_icibs_u_rayleigh],[values_icibs_u_rayleigh 1],'-','color',colours(5,:),'linewidth',linewidth);
%     
%     title(['SNR = ' num2str(snr(snr_idx)) 'dB'],'fontname',fontname,'fontsize',fontsize);    
%     xlabel('Uplink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
%     ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
%     
%     legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
%     
%     set(gca,'fontname',fontname,'fontsize',fontsize);
%     
%     % xlim([5 7]);
%     ylim([0 1]);
%     
%     if (savefig == 1)
%         saveas(gcf,[root_out_prob 'uplink_rayleigh_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'fig');
%         saveas(gcf,[root_out_prob 'uplink_rayleigh_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'png');
%         saveas(gcf,[root_out_prob 'uplink_rayleigh_avg_rate_ter_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_idx)],'epsc2');
%     end
% end