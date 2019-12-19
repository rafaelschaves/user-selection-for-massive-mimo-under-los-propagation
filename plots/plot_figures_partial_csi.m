clear;
close all;
clc;

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Partial CSI/spectral_efficiency_mf_ur_los_M_';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Partial CSI/';

% Macros

MC_1 = 1000;                                                                 % Size of the monte-carlo ensemble
N_XI = 11;

M = 64;                                                                    % Number of antennas at base station
K = 72;                                                                    % Number of mobile users
L = 8;                                                                     % Number of selected users

N_ALG = 4;                                                                 % Number of algorithms for perform user scheduling

BS_RADIUS = 100;
BS_POWER  = 1;

load([root_load num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_radius_' num2str(BS_RADIUS) '_m_BS_POWER_' num2str(BS_POWER) '_W_MC_' num2str(MC_1) '.mat'],'se_ep','se_max','se_ep_sel','se_max_sel');

se_ep_avg  = reshape(mean(se_ep,2),[],N_XI,MC_1);
se_max_avg = reshape(mean(se_max,2),[],N_XI,MC_1);

se_ep_sel_avg  = reshape(mean(se_ep_sel,2),[],N_XI,MC_1,N_ALG);
se_max_sel_avg = reshape(mean(se_max_sel,2),[],N_XI,MC_1,N_ALG);

N_BINS = 10;

cdf_se_ep = zeros(N_BINS,N_XI,N_ALG+1);
x_se_ep   = zeros(N_BINS+1,N_XI,N_ALG+1);

cdf_se_max = zeros(N_BINS,N_XI,N_ALG+1);
x_se_max   = zeros(N_BINS+1,N_XI,N_ALG+1);

for n_xi = 1:N_XI
    [cdf_se_ep(:,n_xi,1),x_se_ep(:,n_xi,1)]   = histcounts(reshape(se_ep_avg(:,n_xi,:),1,[]),N_BINS,'normalization','cdf');
    [cdf_se_max(:,n_xi,1),x_se_max(:,n_xi,1)] = histcounts(reshape(se_max_avg(:,n_xi,:),1,[]),N_BINS,'normalization','cdf');
    
    for n_alg = 1:N_ALG
        [cdf_se_ep(:,n_xi,n_alg+1),x_se_ep(:,n_xi,n_alg+1)]   = histcounts(reshape(se_ep_sel_avg(:,n_xi,:,n_alg),1,[]),N_BINS,'normalization','cdf');
        [cdf_se_max(:,n_xi,n_alg+1),x_se_max(:,n_xi,n_alg+1)] = histcounts(reshape(se_max_sel_avg(:,n_xi,:,n_alg),1,[]),N_BINS,'normalization','cdf');
    end
end

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

colours = [0.0000 0.4470 0.7410;
           0.8500 0.3250 0.0980;
           0.9290 0.6940 0.1250;
           0.4940 0.1840 0.5560;
           0.4660 0.6740 0.1880;
           0.3010 0.7450 0.9330;
           0.6350 0.0780 0.1840];

% figure;
% 
% set(gcf,'position',[0 0 800 600]);
%  
% plot(x_se_ep(:,end,1),[cdf_se_ep(:,end,1); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% hold on;
% plot(x_se_ep(:,end,2),[cdf_se_ep(:,end,2); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,end,3),[cdf_se_ep(:,end,3); 1],'-','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,end,4),[cdf_se_ep(:,end,4); 1],'-','color',colours(4,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,end,5),[cdf_se_ep(:,end,5); 1],'-','color',colours(5,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,end-1,1),[cdf_se_ep(:,end-1,1); 1],'--','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,end-1,2),[cdf_se_ep(:,end-1,2); 1],'--','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,end-1,3),[cdf_se_ep(:,end-1,3); 1],'--','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,end-1,4),[cdf_se_ep(:,end-1,4); 1],'--','color',colours(4,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,end-1,5),[cdf_se_ep(:,end-1,5); 1],'--','color',colours(5,:),'linewidth',linewidth,'markersize',markersize);
%  
% xlabel('Spectral efficiency per user (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
% 
% % legend({'CBS','ICIBS'},'fontname',fontname,'fontsize',fontsize,'location','southeast');
% % legend box off;
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 1]);
% 
% figure;
% 
% set(gcf,'position',[0 0 800 600]);
%  
% plot(x_se_max(:,end,1),[cdf_se_max(:,end,1); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% hold on;
% plot(x_se_max(:,end,2),[cdf_se_max(:,end,2); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,end,3),[cdf_se_max(:,end,3); 1],'-','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,end,4),[cdf_se_max(:,end,4); 1],'-','color',colours(4,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,end,5),[cdf_se_max(:,end,5); 1],'-','color',colours(5,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,end-1,1),[cdf_se_max(:,end-1,1); 1],'--','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,end-1,2),[cdf_se_max(:,end-1,2); 1],'--','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,end-1,3),[cdf_se_max(:,end-1,3); 1],'--','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,end-1,4),[cdf_se_max(:,end-1,4); 1],'--','color',colours(4,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,end-1,5),[cdf_se_max(:,end-1,5); 1],'--','color',colours(5,:),'linewidth',linewidth,'markersize',markersize);
%  
% xlabel('Spectral efficiency per user (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
% 
% % legend({'CBS','ICIBS'},'fontname',fontname,'fontsize',fontsize,'location','southeast');
% % legend box off;
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 1]);

figure;

set(gcf,'position',[0 0 800 600]);
 
plot(x_se_ep(:,11,4),[cdf_se_ep(:,11,4); 1],'o-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
hold on;
plot(x_se_ep(:,11,5),[cdf_se_ep(:,11,5); 1],'o-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_ep(:,10,4),[cdf_se_ep(:,10,4); 1],'s-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_ep(:,10,5),[cdf_se_ep(:,10,5); 1],'s-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_ep(:,9,4),[cdf_se_ep(:,9,4); 1],'^-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_ep(:,9,5),[cdf_se_ep(:,9,5); 1],'^-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_ep(:,8,4),[cdf_se_ep(:,8,4); 1],'v-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_ep(:,8,5),[cdf_se_ep(:,8,5); 1],'v-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,7,4),[cdf_se_ep(:,7,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,7,5),[cdf_se_ep(:,7,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,6,4),[cdf_se_ep(:,6,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,6,5),[cdf_se_ep(:,6,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,5,4),[cdf_se_ep(:,5,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,5,5),[cdf_se_ep(:,5,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,4,4),[cdf_se_ep(:,4,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,4,5),[cdf_se_ep(:,4,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,3,4),[cdf_se_ep(:,3,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,3,5),[cdf_se_ep(:,3,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,2,4),[cdf_se_ep(:,2,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,2,5),[cdf_se_ep(:,2,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,1,4),[cdf_se_ep(:,1,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_ep(:,1,5),[cdf_se_ep(:,1,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);

xlabel('Spectral efficiency per user (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);

legend({'CBS','ICIBS'},'fontname',fontname,'fontsize',fontsize,'location','southeast');
legend box off;

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 1]);
xlim([1.5 4]);

figure;

set(gcf,'position',[0 0 800 600]);
 
plot(x_se_max(:,11,4),[cdf_se_max(:,11,4); 1],'o-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
hold on;
plot(x_se_max(:,11,5),[cdf_se_max(:,11,5); 1],'o-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_max(:,10,4),[cdf_se_max(:,10,4); 1],'s-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_max(:,10,5),[cdf_se_max(:,10,5); 1],'s-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_max(:,9,4),[cdf_se_max(:,9,4); 1],'^-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_max(:,9,5),[cdf_se_max(:,9,5); 1],'^-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_max(:,8,4),[cdf_se_max(:,8,4); 1],'v-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(x_se_max(:,8,5),[cdf_se_max(:,8,5); 1],'v-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,7,4),[cdf_se_max(:,7,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,7,5),[cdf_se_max(:,7,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,6,4),[cdf_se_max(:,6,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,6,5),[cdf_se_max(:,6,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,5,4),[cdf_se_max(:,5,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,5,5),[cdf_se_max(:,5,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,4,4),[cdf_se_max(:,4,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,4,5),[cdf_se_max(:,4,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,3,4),[cdf_se_max(:,3,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,3,5),[cdf_se_max(:,3,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,2,4),[cdf_se_max(:,2,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,2,5),[cdf_se_max(:,2,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,1,4),[cdf_se_max(:,1,4); 1],'-','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
% plot(x_se_max(:,1,5),[cdf_se_max(:,1,5); 1],'-','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);

xlabel('Spectral efficiency per user (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);

legend({'CBS','ICIBS'},'fontname',fontname,'fontsize',fontsize,'location','southeast');
legend box off;

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 1]);
xlim([1.5 4]);

% [cdf_se_dl_user(:,1,1),se_dl_user(:,1,1)] = ecdf(reshape(se_dl_avg(:,:,end  ,:,1),1,[]));
% [cdf_se_dl_user(:,2,1),se_dl_user(:,2,1)] = ecdf(reshape(se_dl_avg(:,:,end-1,:,1),1,[]));
% [cdf_se_dl_user(:,3,1),se_dl_user(:,3,1)] = ecdf(reshape(se_dl_avg(:,:,end-2,:,1),1,[]));
% 
% [cdf_se_dl_user(:,1,2),se_dl_user(:,1,2)] = ecdf(reshape(se_dl_avg(:,:,end  ,:,2),1,[]));
% [cdf_se_dl_user(:,2,2),se_dl_user(:,2,2)] = ecdf(reshape(se_dl_avg(:,:,end-1,:,2),1,[]));
% [cdf_se_dl_user(:,3,2),se_dl_user(:,3,2)] = ecdf(reshape(se_dl_avg(:,:,end-2,:,2),1,[]));