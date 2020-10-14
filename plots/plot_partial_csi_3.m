clear;
close all;
clc;

% Macros

MC   = 500;                                                              % Size of the monte-carlo ensemble
N_MC = 10;

M = 100;                                                                   % Number of antennas at base station
K = 100;                                                                   % Number of users at the cell 

% M = 50  & K = [50  75]
% M = 100 & K = [100 150]

snr = -5;

err = [0 pi/(6*M) pi/(3*M) pi/(2*M)];

N_ALG = 3;                                                                % Number of algorithms for perform user scheduling
N_PRE = 3;
N_ERR = length(err);

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Downlink/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Selection/Downlink/';

zero_pad_1 = '%03d';
zero_pad_2 = '%02d';
zero_pad_3 = '%04d';

chn_type = 'ur_los';

if K > M
    L_max = M;
else
    L_max = K-1;
end

prob_sel   = zeros(L_max,N_ALG,N_ERR);
avg_sum_se = zeros(L_max,N_PRE,N_ERR);    
S_set_mc   = zeros(K,L_max,N_PRE,N_ERR,MC*N_MC);
se_all_mc  = zeros(K,N_PRE,N_ERR,MC*N_MC);
    
for n_mc = 1:N_MC
    load([root_load 'spectral_efficiency_all_L_' chn_type '_partial_csi_2_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' sprintf(zero_pad_3,MC) '_' sprintf(zero_pad_2,n_mc) '.mat']);
    
    idx_ini = (n_mc - 1)*MC + 1;
    idx_end = n_mc*MC;
    
    se_all_mc(:,:,:,idx_ini:idx_end)  = se;
    S_set_mc(:,:,:,:,idx_ini:idx_end) = S_set;
    
    clear se se_s_all_L S_set;
end

aux = repmat((1:L_max)'*N_MC*MC,1,N_ALG);

% aux_SOS_e0   = S_set_mc(:,:,1,1,:);
% aux_CBS_e0   = S_set_mc(:,:,2,1,:);
% aux_ICIBS_e0 = S_set_mc(:,:,3,1,:);
% 
% aux_SOS_e1   = S_set_mc(:,:,1,2,:);
% aux_CBS_e1   = S_set_mc(:,:,2,2,:);
% aux_ICIBS_e1 = S_set_mc(:,:,3,2,:);
% 
% aux_SOS_e2   = S_set_mc(:,:,1,3,:);
% aux_CBS_e2   = S_set_mc(:,:,2,3,:);
% aux_ICIBS_e2 = S_set_mc(:,:,3,3,:);
% 
% aux_SOS_e3   = S_set_mc(:,:,1,4,:);
% aux_CBS_e3   = S_set_mc(:,:,2,4,:);
% aux_ICIBS_e3 = S_set_mc(:,:,3,4,:);
% 
% for mc = 1:N_MC*MC
%     aux_prob_sos_1(:,mc) = sum((aux_SOS_e0(:,:,mc) - aux_SOS_e1(:,:,mc)<0),1)';
%     aux_prob_sos_2(:,mc) = sum((aux_SOS_e0(:,:,mc) - aux_SOS_e2(:,:,mc)<0),1)';
%     aux_prob_sos_3(:,mc) = sum((aux_SOS_e0(:,:,mc) - aux_SOS_e3(:,:,mc)<0),1)';
%     
%     aux_prob_cbs_1(:,mc) = sum((aux_CBS_e0(:,:,mc) - aux_CBS_e1(:,:,mc)<0),1)';
%     aux_prob_cbs_2(:,mc) = sum((aux_CBS_e0(:,:,mc) - aux_CBS_e2(:,:,mc)<0),1)';
%     aux_prob_cbs_3(:,mc) = sum((aux_CBS_e0(:,:,mc) - aux_CBS_e3(:,:,mc)<0),1)';
%     
%     aux_prob_icibs_1(:,mc) = sum((aux_ICIBS_e0(:,:,mc) - aux_ICIBS_e1(:,:,mc)<0),1)';
%     aux_prob_icibs_2(:,mc) = sum((aux_ICIBS_e0(:,:,mc) - aux_ICIBS_e2(:,:,mc)<0),1)';
%     aux_prob_icibs_3(:,mc) = sum((aux_ICIBS_e0(:,:,mc) - aux_ICIBS_e3(:,:,mc)<0),1)';
% end
% 
% prob_sel(:,1,2) = sum(aux_prob_sos_1,2)./aux(:,1);
% prob_sel(:,1,3) = sum(aux_prob_sos_2,2)./aux(:,1);
% prob_sel(:,1,4) = sum(aux_prob_sos_3,2)./aux(:,1);
% 
% prob_sel(:,2,2) = sum(aux_prob_cbs_1,2)./aux(:,1);
% prob_sel(:,2,3) = sum(aux_prob_cbs_2,2)./aux(:,1);
% prob_sel(:,2,4) = sum(aux_prob_cbs_3,2)./aux(:,1);
% 
% prob_sel(:,3,2) = sum(aux_prob_icibs_1,2)./aux(:,1);
% prob_sel(:,3,3) = sum(aux_prob_icibs_2,2)./aux(:,1);
% prob_sel(:,3,4) = sum(aux_prob_icibs_3,2)./aux(:,1);

prob_sel(:,:,2) = reshape(sum(sum((S_set_mc(:,:,:,1,:) - S_set_mc(:,:,:,2,:))<0,5),1),L_max,N_PRE)./aux;
prob_sel(:,:,3) = reshape(sum(sum((S_set_mc(:,:,:,1,:) - S_set_mc(:,:,:,3,:))<0,5),1),L_max,N_PRE)./aux;
prob_sel(:,:,4) = reshape(sum(sum((S_set_mc(:,:,:,1,:) - S_set_mc(:,:,:,4,:))<0,5),1),L_max,N_PRE)./aux;

%sum_se = reshape(sum(se_all_mc,1),N_PRE,N_ERR,MC*N_MC);

%avg_sum_se(:,:,k) = mean(sum_se,3);

% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 30;

marker = {'o','s','^'};

linestyle = {'-','--',':'};

savefig = 1;

% NS - No selection
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_prec   = {'MRT','ZF','MMSE'};
legend_algo   = {'NS','SOS','CBS','ICIBS'};
legend_algo_2 = {'SOS','CBS','ICIBS'};

location_1 = 'northwest';
location_2 = 'northeast';
location_3 = 'southwest';
location_4 = 'southeast';

colours = [0.0000 0.4470 0.7410;
           0.8500 0.3250 0.0980;
           0.9290 0.6940 0.1250;
           0.4940 0.1840 0.5560;
           0.4660 0.6740 0.1880;
           0.3010 0.7450 0.9330;
           0.6350 0.0780 0.1840];

figure;
       
set(gcf,'position',[0 0 800 600]);

if K > M
    plot(1:L_max,prob_sel(:,1,2),linestyle{1},'color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(1:L_max,prob_sel(:,2,2),linestyle{1},'color',colours(2,:),'linewidth',linewidth);
    plot(1:L_max,prob_sel(:,3,2),linestyle{1},'color',colours(3,:),'linewidth',linewidth);
    plot(1:L_max,prob_sel(:,1,3),linestyle{2},'color',colours(1,:),'linewidth',linewidth);
    plot(1:L_max,prob_sel(:,2,3),linestyle{2},'color',colours(2,:),'linewidth',linewidth);
    plot(1:L_max,prob_sel(:,3,3),linestyle{2},'color',colours(3,:),'linewidth',linewidth);
    plot(1:L_max,prob_sel(:,1,4),linestyle{3},'color',colours(1,:),'linewidth',linewidth);
    plot(1:L_max,prob_sel(:,2,4),linestyle{3},'color',colours(2,:),'linewidth',linewidth);
    plot(1:L_max,prob_sel(:,3,4),linestyle{3},'color',colours(3,:),'linewidth',linewidth);
else
    plot(1:K,[prob_sel(:,1,2)' 0],linestyle{1},'color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(1:K,[prob_sel(:,2,2)' 0],linestyle{1},'color',colours(2,:),'linewidth',linewidth);
    plot(1:K,[prob_sel(:,3,2)' 0],linestyle{1},'color',colours(3,:),'linewidth',linewidth);
    plot(1:K,[prob_sel(:,1,3)' 0],linestyle{2},'color',colours(1,:),'linewidth',linewidth);
    plot(1:K,[prob_sel(:,2,3)' 0],linestyle{2},'color',colours(2,:),'linewidth',linewidth);
    plot(1:K,[prob_sel(:,3,3)' 0],linestyle{2},'color',colours(3,:),'linewidth',linewidth);
    plot(1:K,[prob_sel(:,1,4)' 0],linestyle{3},'color',colours(1,:),'linewidth',linewidth);
    plot(1:K,[prob_sel(:,2,4)' 0],linestyle{3},'color',colours(2,:),'linewidth',linewidth);
    plot(1:K,[prob_sel(:,3,4)' 0],linestyle{3},'color',colours(3,:),'linewidth',linewidth);
end

xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
ylabel('$\textrm{Pr}\{\cal{S_{\varepsilon}} \neq \cal{S}\}$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');

legend(legend_algo_2,'fontname',fontname,'fontsize',fontsize,'location',location_2,'interpreter','latex');
legend box off;

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([1 K]);
ylim([0 1]);

if (savefig == 1)
    saveas(gcf,[root_save 'prob_user_selection_partial_csi_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
    saveas(gcf,[root_save 'prob_user_selection_partial_csi_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
    saveas(gcf,[root_save 'prob_user_selection_partial_csi_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
end