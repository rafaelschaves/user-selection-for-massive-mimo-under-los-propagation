clear;
close all;
clc;

% Macros

MC = 1000;                                                                 % Size of the monte-carlo ensemble

M = 200;                                                                   % Number of antennas at base station
K = 150;                                                                    % Number of users at the cell 

% M = 50  & K = [10 30 50 70]
% M = 100 & K = [10 30 50 70 100 130] 
% M = 200 & K = [10 30 50 70 100 130 150 200 220]

if K > M
    L_max = M-1;
else
    L_max = K-1;
end

snr = -5;

N_ALG = 3;                                                                 % Number of algorithms for perform user scheduling
N_PRE = 2;

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Downlink/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Selection/';

chn_type = 'ur_los';

% Loading data

sum_se_s = zeros(L_max,N_PRE,N_ALG,MC);

load([root_load 'spectral_efficiency_all_L_' chn_type '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC) '.mat']);

sum_se = reshape(sum(se,1),2,MC);

for L = 1:L_max
    sum_se_s(L,:,:,:) = sum(se_s_all_L(:,L,:,:,:),1);
end

avg_sum_se = mean(sum_se,2);
avg_sum_se_s = mean(sum_se_s,4);

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

legend_algo = {'NS','SOS','CBS','ICIBS'};

location_1 = 'northwest';
location_2 = 'northeast';

colours = [0.0000 0.4470 0.7410;
           0.8500 0.3250 0.0980;
           0.9290 0.6940 0.1250;
           0.4940 0.1840 0.5560;
           0.4660 0.6740 0.1880;
           0.3010 0.7450 0.9330;
           0.6350 0.0780 0.1840];

figure;
       
set(gcf,'position',[0 0 800 600]);

plot(1:L_max,avg_sum_se(1)*ones(L_max,1),'-k','linewidth',linewidth);
hold on;
plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
plot(1:L_max,avg_sum_se(2)*ones(L_max,1),'--k','linewidth',linewidth);
plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);

xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
ylabel('Sum-spectral efficiency','fontname',fontname,'fontsize',fontsize);

legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_2);
legend box off;

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([1 L_max]);

%         if chn_idx == 1
%             if m == 1
%                 dim = [0.16 0.18 0.035 0.65];
%             else
%                 dim = [0.16 0.18 0.035 0.725];
%             end
%             annotation('rectangle',dim,'linewidth',linewidth);
%
%             axes('position',[.6 .275 .25 .25]);
%             box on;
%
%             plot(pdf_psi_x(1,:,m,chn_idx),pdf_psi_y(1,:,m,chn_idx),'-','color',colours(1,:),'linewidth',linewidth);
%             hold on;
%             plot(pdf_psi_x(2,:,m,chn_idx),pdf_psi_y(2,:,m,chn_idx),'-','color',colours(2,:),'linewidth',linewidth);
%             plot(pdf_psi_x(3,:,m,chn_idx),pdf_psi_y(3,:,m,chn_idx),'-','color',colours(3,:),'linewidth',linewidth);
%             plot(pdf_psi_x(4,:,m,chn_idx),pdf_psi_y(4,:,m,chn_idx),'-','color',colours(4,:),'linewidth',linewidth);
%             plot(pdf_psi_x(5,:,m,chn_idx),pdf_psi_y(5,:,m,chn_idx),'-','color',colours(5,:),'linewidth',linewidth);
%
%             if m == 1
%                 xlim([0 0.04]);
%             elseif m == 2
%                 xlim([0 0.02]);
%             else
%                 xlim([0 0.015]);
%             end
%
%             set(gca,'fontname',fontname,'fontsize',fontsize);
%         end

if (savefig == 1)
    saveas(gcf,[root_save 'sum_se_all_L' chn_type '_M_' num2str(M) '_K_' num2str(K)],'fig');
    saveas(gcf,[root_save 'sum_se_all_L' chn_type '_M_' num2str(M) '_K_' num2str(K)],'png');
    saveas(gcf,[root_save 'sum_se_all_L' chn_type '_M_' num2str(M) '_K_' num2str(K)],'epsc2');
end

%         figure;
%
%         set(gcf,'position',[0 0 800 600]);
%
%         plot(1:L_max,avg_sum_se(2)*ones(L_max,1),'-k','linewidth',linewidth);
%         hold on;
%         plot(1:L_max,avg_sum_se_s(:,2,1),'-','color',colours(1,:),'linewidth',linewidth);
%         plot(1:L_max,avg_sum_se_s(:,2,2),'-','color',colours(2,:),'linewidth',linewidth);
%         plot(1:L_max,avg_sum_se_s(:,2,3),'-','color',colours(3,:),'linewidth',linewidth);
%
%         xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
%         ylabel('Sum-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
%
%         legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_2);
%         legend box off;
%
%         xlim([1 L_max]);
%
%         set(gca,'fontname',fontname,'fontsize',fontsize);