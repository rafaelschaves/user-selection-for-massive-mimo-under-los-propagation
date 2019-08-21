clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

r_k = 1.25;
r_l = 0.25;

M = 64;                                                                    % Number of antennas at base station
K = r_k*M;
L = r_l*K;

snr_db_ul = -7;
snr_db_dl = 10;

snr_ul = 10^(snr_db_ul/10);
snr_dl = 10^(snr_db_dl/10);

% Roots

root_load_ul = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Uplink/';
root_load_dl = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Downlink/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Selection/';

% Loading data

se = zeros(K,MC,2);                                                        % Rate using all K users

sum_se_1   = zeros(MC,2);
sum_se_2   = zeros(MC,2);
sum_se_upp = zeros(MC,2);
sum_se_low = zeros(MC,2);

load([root_load_ul 'spectral_efficiency_mf_rayleigh_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_db_ul) '_dB_MC_' num2str(MC) '.mat']);
load([root_load_dl 'spectral_efficiency_mf_rayleigh_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_db_dl) '_dB_MC_' num2str(MC) '.mat']);

se(:,:,1) = se_u;
se(:,:,2) = se_d;

psi       = reshape(sum(betarnd(1,M-1,K-1,K,MC),1),K,MC);
max_psi   = max(psi);
min_psi   = min(psi);
range_psi = (min_psi + max_psi)/2;
% range_psi = mean(psi);

sum_se_1(:,1)   = 1/2*sum(log2(1 + snr_ul*M./(1 + snr_ul*M*psi)));
sum_se_1(:,2)   = 1/2*sum(log2(1 + snr_dl*M./(K + snr_dl*M*psi)));

sum_se_2(:,1)   = 1/2*K*log2(1 + snr_ul*M./(1 + snr_ul*M*range_psi));
sum_se_2(:,2)   = 1/2*K*log2(1 + snr_dl*M./(K + snr_dl*M*range_psi));

sum_se_low(:,1) = 1/2*K*log2(1 + snr_ul*M./(1 + snr_ul*M*max_psi));
sum_se_low(:,2) = 1/2*K*log2(1 + snr_dl*M./(K + snr_dl*M*max_psi));

sum_se_upp(:,1) = 1/2*K*log2(1 + snr_ul*M./(1 + snr_ul*M*min_psi));
sum_se_upp(:,2) = 1/2*K*log2(1 + snr_dl*M./(K + snr_dl*M*min_psi));

% Post processing - Calculating the CDF

nbins = 20;

cdf_sum_se = zeros(nbins,2);
edg_sum_se = zeros(nbins+1,2);

cdf_sum_se_1 = zeros(nbins,2);
edg_sum_se_1 = zeros(nbins+1,2);

cdf_sum_se_2 = zeros(nbins,2);
edg_sum_se_2 = zeros(nbins+1,2);

cdf_sum_se_upp = zeros(nbins,2);
edg_sum_se_upp = zeros(nbins+1,2);

cdf_sum_se_low = zeros(nbins,2);
edg_sum_se_low = zeros(nbins+1,2);

[cdf_sum_se(:,1),edg_sum_se(:,1)] = histcounts(sum(se(:,:,1)),nbins,'normalization','cdf');
[cdf_sum_se(:,2),edg_sum_se(:,2)] = histcounts(sum(se(:,:,2)),nbins,'normalization','cdf');

[cdf_sum_se_1(:,1),edg_sum_se_1(:,1)] = histcounts(sum_se_1(:,1),nbins,'normalization','cdf');
[cdf_sum_se_1(:,2),edg_sum_se_1(:,2)] = histcounts(sum_se_1(:,2),nbins,'normalization','cdf');

[cdf_sum_se_2(:,1),edg_sum_se_2(:,1)] = histcounts(sum_se_2(:,1),nbins,'normalization','cdf');
[cdf_sum_se_2(:,2),edg_sum_se_2(:,2)] = histcounts(sum_se_2(:,2),nbins,'normalization','cdf');

[cdf_sum_se_upp(:,1),edg_sum_se_upp(:,1)] = histcounts(sum_se_upp(:,1),nbins,'normalization','cdf');
[cdf_sum_se_upp(:,2),edg_sum_se_upp(:,2)] = histcounts(sum_se_upp(:,2),nbins,'normalization','cdf');

[cdf_sum_se_low(:,1),edg_sum_se_low(:,1)] = histcounts(sum_se_low(:,1),nbins,'normalization','cdf');
[cdf_sum_se_low(:,2),edg_sum_se_low(:,2)] = histcounts(sum_se_low(:,2),nbins,'normalization','cdf');

% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 30;

legend_text = {'Simulated','Teoric','Approximated','Upper','Lower'};

savefig = 0;

colours = [0.0000 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840];

figure;

set(gcf,'position',[0 0 800 600]);

plot(edg_sum_se(:,1)    ,[cdf_sum_se(:,1); 1]    ,'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
hold on;
plot(edg_sum_se_1(:,1)  ,[cdf_sum_se_1(:,1); 1]  ,'o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(edg_sum_se_2(:,1)  ,[cdf_sum_se_2(:,1); 1]  ,'s','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(edg_sum_se_upp(:,1),[cdf_sum_se_upp(:,1); 1],'^','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(edg_sum_se_low(:,1),[cdf_sum_se_low(:,1); 1],'v','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(edg_sum_se(:,2)    ,[cdf_sum_se(:,2); 1]    ,'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
hold on;
plot(edg_sum_se_1(:,2)  ,[cdf_sum_se_1(:,2); 1]  ,'o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(edg_sum_se_2(:,2)  ,[cdf_sum_se_2(:,2); 1]  ,'s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(edg_sum_se_upp(:,2),[cdf_sum_se_upp(:,2); 1],'^','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(edg_sum_se_low(:,2),[cdf_sum_se_low(:,2); 1],'v','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);

xlabel('Sum-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);

legend(legend_text,'fontname',fontname,'fontsize',fontsize,'location','southeast');
legend box off;

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 1]);
xlim([10 40]);

%     if (savefig == 1)
%         saveas(gcf,[root_save 'cdf_sum_se_rayleigh_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'fig');
%         saveas(gcf,[root_save 'cdf_sum_se_rayleigh_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'png');
%         saveas(gcf,[root_save 'cdf_sum_se_rayleigh_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'epsc2');
%     end