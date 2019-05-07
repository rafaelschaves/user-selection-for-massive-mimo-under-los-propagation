clear;
close all;
clc;

M = 64;
K = 18;
L = 13;

snr = [-20 -15 -10 -5 0 5 10]';

likely_capacity_ur_los = [7.235 14.460 23.405 32.625 39.515 43.945 45.625;
                          5.080 10.550 17.685 25.365 31.530 35.245 37.460;
                          6.765 13.400 21.270 29.080 34.615 38.155 39.950;
                          5.445 12.265 22.840 35.530 46.330 53.975 57.520;
                          5.445 12.320 22.815 35.750 46.715 54.505 57.780];
                      
likely_capacity_sparse = [10.415 10.765 10.845 10.890 10.900 10.865 10.825;
                           9.045  9.410  9.420  9.530  9.570  9.550  9.530;
                           9.115  9.430  9.555  9.545  9.580  9.635  9.440;
                          14.515 15.005 15.155 15.265 15.195 15.210 15.185;
                          14.795 15.350 15.535 15.700 15.625 15.520 15.540];
                      
likely_capacity_rayleigh = [10.885 21.010 30.150 35.275 37.360 38.150 38.385;
                             8.080 16.205 24.250 29.170 31.225 32.060 32.310;
                             8.210 16.380 24.395 29.275 31.340 32.145 32.265;
                             8.210 16.815 25.765 31.510 34.005 35.010 35.285;
                             8.245 17.000 26.290 32.395 35.160 36.195 36.560];
                         
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

root_rate_lik = '../figures/rate/likely_';

savefig = 1;

colours = get(gca,'colororder');
close;                        

figure;
    
set(gcf,'position',[0 0 800 600]);
    
plot(snr,likely_capacity_ur_los(1,:),'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot(snr,likely_capacity_ur_los(2,:),'-','color',colours(2,:),'linewidth',linewidth);
plot(snr,likely_capacity_ur_los(3,:),'-','color',colours(3,:),'linewidth',linewidth);
plot(snr,likely_capacity_ur_los(4,:),'-','color',colours(4,:),'linewidth',linewidth);
plot(snr,likely_capacity_ur_los(5,:),'-','color',colours(5,:),'linewidth',linewidth);

xlabel('SNR (dB)','fontname',fontname,'fontsize',fontsize);
ylabel('95% likely sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    
legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
set(gca,'fontname',fontname,'fontsize',fontsize);

% dim = [0.65 0.75 0.2 0.15];
dim = [0.2 0.25 0.2 0.15];

annotation('ellipse',dim,'linewidth',linewidth);

axes('position',[.6 .20 .27 .27])
box on;

plot(snr(1:3),likely_capacity_ur_los(1,1:3),'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot(snr(1:3),likely_capacity_ur_los(2,1:3),'-','color',colours(2,:),'linewidth',linewidth);
plot(snr(1:3),likely_capacity_ur_los(3,1:3),'-','color',colours(3,:),'linewidth',linewidth);
plot(snr(1:3),likely_capacity_ur_los(4,1:3),'-','color',colours(4,:),'linewidth',linewidth);
plot(snr(1:3),likely_capacity_ur_los(5,1:3),'-','color',colours(5,:),'linewidth',linewidth);

set(gca,'fontname',fontname,'fontsize',fontsize);

% ylim([45 58]);

if (savefig == 1)
    saveas(gcf,[root_rate_lik 'uplink_ur_los_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'fig');
    saveas(gcf,[root_rate_lik 'uplink_ur_los_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'png');
    saveas(gcf,[root_rate_lik 'uplink_ur_los_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'epsc2');
end

figure;
    
set(gcf,'position',[0 0 800 600]);
    
plot(snr,likely_capacity_sparse(1,:),'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot(snr,likely_capacity_sparse(2,:),'-','color',colours(2,:),'linewidth',linewidth);
plot(snr,likely_capacity_sparse(3,:),'-','color',colours(3,:),'linewidth',linewidth);
plot(snr,likely_capacity_sparse(4,:),'-','color',colours(4,:),'linewidth',linewidth);
plot(snr,likely_capacity_sparse(5,:),'-','color',colours(5,:),'linewidth',linewidth);

xlabel('SNR (dB)','fontname',fontname,'fontsize',fontsize);
ylabel('95% likely sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    
legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
set(gca,'fontname',fontname,'fontsize',fontsize);
       
if (savefig == 1)
    saveas(gcf,[root_rate_lik 'uplink_sparse_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'fig');
    saveas(gcf,[root_rate_lik 'uplink_sparse_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'png');
    saveas(gcf,[root_rate_lik 'uplink_sparse_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'epsc2');
end

figure;
    
set(gcf,'position',[0 0 800 600]);
    
plot(snr,likely_capacity_rayleigh(1,:),'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot(snr,likely_capacity_rayleigh(2,:),'-','color',colours(2,:),'linewidth',linewidth);
plot(snr,likely_capacity_rayleigh(3,:),'-','color',colours(3,:),'linewidth',linewidth);
plot(snr,likely_capacity_rayleigh(4,:),'-','color',colours(4,:),'linewidth',linewidth);
plot(snr,likely_capacity_rayleigh(5,:),'-','color',colours(5,:),'linewidth',linewidth);

xlabel('SNR (dB)','fontname',fontname,'fontsize',fontsize);
ylabel('95% likely sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    
legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
set(gca,'fontname',fontname,'fontsize',fontsize);
       
if (savefig == 1)
    saveas(gcf,[root_rate_lik 'uplink_rayleigh_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'fig');
    saveas(gcf,[root_rate_lik 'uplink_rayleigh_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'png');
    saveas(gcf,[root_rate_lik 'uplink_rayleigh_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'epsc2');
end

% 
% bar_sum = [sum(mean(rate_u,2)) sum(mean(rate_d,2));
%            sum(mean(rate_rs_u,2))  sum(mean(rate_rs_d,2));
%            sum(mean(rate_sos_u,2)) sum(mean(rate_sos_d,2)); 
%            sum(mean(rate_icibs_u,2)) sum(mean(rate_icibs_d,2))];
% 
% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% bar(cat,bar_sum,BAR_SIZE);
% 
% xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
% ylabel('Average sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_link,'fontname',fontname,'fontsize',fontsize);
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 40]);
% 
% saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'fig');
% saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'png');
% saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'epsc2');
% 
% bar_sum_5 = [R_sum_u_ns R_sum_d_ns; 
%              R_sum_u_rs R_sum_d_rs; 
%              R_sum_u_sos R_sum_d_sos;
%              R_sum_u_icibs R_sum_d_icibs];
%          
% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% bar(cat,bar_sum_5,BAR_SIZE);
% 
% xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
% ylabel('95% likely sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_link,'fontname',fontname,'fontsize',fontsize);
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 40]);
% 
% saveas(gcf,[root_erg_rate '95_sum_rate_M_' num2str(M) '_K_' num2str(K)],'fig');
% saveas(gcf,[root_erg_rate '95_sum_rate_M_' num2str(M) '_K_' num2str(K)],'png');
% saveas(gcf,[root_erg_rate '95_sum_rate_M_' num2str(M) '_K_' num2str(K)],'epsc2');
% 
% 
% bar_u = [mean(mean(rate_u,2)) mean(mean(rate_d,2));
%          mean(mean(rate_rs_u,2)) mean(mean(rate_rs_d,2));
%          mean(mean(rate_sos_u,2)) mean(mean(rate_sos_d,2));
%          mean(mean(rate_icibs_u,2)) mean(mean(rate_icibs_d,2))];
% 
% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% bar(cat,bar_u,BAR_SIZE);
% 
% xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
% ylabel('Average rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_link,'fontname',fontname,'fontsize',fontsize);
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 10]);
% 
% saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'fig');
% saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'png');
% saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'epsc2');
% 
% bar_u_5 = [R_u_ns R_d_ns;
%            R_u_rs R_d_rs; 
%            R_u_sos R_d_sos;
%            R_u_icibs R_d_icibs];
% 
% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% bar(cat,bar_u_5,BAR_SIZE);
% 
% xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
% ylabel('95% likely average rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_link,'fontname',fontname,'fontsize',fontsize);
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 10]);
% 
% saveas(gcf,[root_erg_rate '95_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'fig');
% saveas(gcf,[root_erg_rate '95_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'png');
% saveas(gcf,[root_erg_rate '95_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'epsc2');