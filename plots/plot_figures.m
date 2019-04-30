clear;
close all;
clc;

load('../results/rate_mf_M_500_K_30_L_25_SNR_-7_dB_MC_10000.mat');

M = 500;
K = 5;

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
% ICIBS - ICI-based selection

legend_algo = {'NS','RS','SOS','ICIBS'};
legend_link = {'Uplink','Downlink'};

location = 'northwest';

cat = categorical(legend_algo);
cat = reordercats(cat,legend_algo);

root_rate_fit = '../figures/rate/fit_';
root_erg_rate = '../figures/rate/erg_cap_';
root_out_prob = '../figures/rate/out_prob_';

p_o = 0.05;

R_sum_u_ns    = 27.37;
R_sum_u_rs    = 22.43;
R_sum_u_sos   = 22.43;
R_sum_u_icibs = 23.39;

R_sum_d_ns    = 31.46;
R_sum_d_rs    = 26.16;
R_sum_d_sos   = 26.19;
R_sum_d_icibs = 28.16;

R_u_ns    = 5.470;
R_u_rs    = 5.605;
R_u_sos   = 5.605;
R_u_icibs = 5.845;

R_d_ns    = 6.285;
R_d_rs    = 6.535;
R_d_sos   = 6.540;
R_d_icibs = 7.035;

colours = get(gca,'colororder');
close;

p_u = zeros(2,K);

psi_range = (0:0.01:0.4);

f_u = zeros(K,length(psi_range));

for k = 1:K
    curvefit = fit(psi(k,:)',rate_u(k,:)','poly1');
    
    p_u(1,k) = curvefit.p1;
    p_u(2,k) = curvefit.p2;
    
    f_u(k,:) = p_u(1,k)*psi_range + p_u(2,k);
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot(psi(k,:),rate_u(k,:),'.','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(psi_range,f_u(k,:),'-','color',colours(2,:),'linewidth',linewidth);
    
    xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
    ylabel('Uplink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    
    legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    xlim([0 0.4]);
    
    saveas(gcf,[root_rate_fit 'uplink_M_' num2str(M) '_' num2str(k) '_user'],'fig');
    saveas(gcf,[root_rate_fit 'uplink_M_' num2str(M) '_' num2str(k) '_user'],'png');
    saveas(gcf,[root_rate_fit 'uplink_M_' num2str(M) '_' num2str(k) '_user'],'epsc2');
end

p_d = zeros(3,K);

psi_range = (0:0.01:0.4);

f_d = zeros(K,length(psi_range));

for k = 1:K
    curvefit = fit(psi(k,:)',rate_d(k,:)','poly2');
    
    p_d(1,k) = curvefit.p1;
    p_d(2,k) = curvefit.p2;
    p_d(3,k) = curvefit.p3;
        
    f_d(k,:) = p_d(1,k)*psi_range.^2 + p_d(2,k)*psi_range + p_d(3,k);
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot(psi(k,:),rate_d(k,:),'.','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(psi_range,f_d(k,:),'-','color',colours(2,:),'linewidth',linewidth);
    
    xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
    ylabel('Downlink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    
    legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    xlim([0 0.4]);
    
    saveas(gcf,[root_rate_fit 'downlink_M_' num2str(M) '_' num2str(k) '_user'],'fig');
    saveas(gcf,[root_rate_fit 'downlink_M_' num2str(M) '_' num2str(k) '_user'],'png');
    saveas(gcf,[root_rate_fit 'downlink_M_' num2str(M) '_' num2str(k) '_user'],'epsc2');
end

[values_u,edges_u]             = histcounts(sum(rate_u),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_rs_u,edges_rs_u]       = histcounts(sum(rate_rs_u),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_sos_u,edges_sos_u]     = histcounts(sum(rate_sos_u),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_icibs_u,edges_icibs_u] = histcounts(sum(rate_icibs_u),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');

x_sum_u = 20:35;
y_sum_u = p_o*ones(length(x_sum_u),1);

figure;

set(gcf,'position',[0 0 800 600]);

plot([20 edges_u 35],[0 values_u 1 1],'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot([20 edges_rs_u 35],[0 values_rs_u 1 1],'-','color',colours(2,:),'linewidth',linewidth);
plot([20 edges_sos_u 35],[0 values_sos_u 1 1],'-','color',colours(3,:),'linewidth',linewidth);
plot([20 edges_icibs_u 35],[0 values_icibs_u 1 1],'-','color',colours(4,:),'linewidth',linewidth);
plot(x_sum_u,y_sum_u,'--k','linewidth',linewidth);
plot(R_sum_u_ns,p_o,'o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(R_sum_u_rs,p_o,'o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(R_sum_u_sos,p_o,'o','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
plot(R_sum_u_icibs,p_o,'o','color',colours(4,:),'linewidth',linewidth,'markersize',markersize);

xlabel('Uplink sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);

legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([20 35]);
ylim([0 1]);

saveas(gcf,[root_out_prob 'uplink_sum_rate_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_out_prob 'uplink_sum_rate_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_out_prob 'uplink_sum_rate_M_' num2str(M) '_K_' num2str(K)],'epsc2');

[values_d,edges_d]              = histcounts(sum(rate_d),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_rs_d,edges_rs_d]        = histcounts(sum(rate_rs_d),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_sos_d,edges_sos_d]      = histcounts(sum(rate_sos_d),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_icibs_d, edges_icibs_d] = histcounts(sum(rate_icibs_d),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');

x_sum_d = 20:45;
y_sum_d = p_o*ones(length(x_sum_d),1);

figure;

set(gcf,'position',[0 0 800 600]);

plot([20 edges_d 45],[0 values_d 1 1],'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot([20 edges_rs_d 45],[0 values_rs_d 1 1],'-','color',colours(2,:),'linewidth',linewidth);
plot([20 edges_sos_d 45],[0 values_sos_d 1 1],'-','color',colours(3,:),'linewidth',linewidth);
plot([20 edges_icibs_d 45],[0 values_icibs_d 1 1],'-','color',colours(4,:),'linewidth',linewidth);
plot(x_sum_d,y_sum_d,'--k','linewidth',linewidth);
plot(R_sum_d_ns,p_o,'o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(R_sum_d_rs,p_o,'o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(R_sum_d_sos,p_o,'o','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
plot(R_sum_d_icibs,p_o,'o','color',colours(4,:),'linewidth',linewidth,'markersize',markersize);

xlabel('Downlink sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);

legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([20 45]);
ylim([0 1]);

saveas(gcf,[root_out_prob 'downlink_sum_rate_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_out_prob 'downlink_sum_rate_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_out_prob 'downlink_sum_rate_M_' num2str(M) '_K_' num2str(K)],'epsc2');

bar_sum = [sum(mean(rate_u,2)) sum(mean(rate_d,2));
           sum(mean(rate_rs_u,2))  sum(mean(rate_rs_d,2));
           sum(mean(rate_sos_u,2)) sum(mean(rate_sos_d,2)); 
           sum(mean(rate_icibs_u,2)) sum(mean(rate_icibs_d,2))];

figure;

set(gcf,'position',[0 0 800 600]);

bar(cat,bar_sum,BAR_SIZE);

xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
ylabel('Average sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);

legend(legend_link,'fontname',fontname,'fontsize',fontsize);

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 40]);

saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'epsc2');

bar_sum_5 = [R_sum_u_ns R_sum_d_ns; 
             R_sum_u_rs R_sum_d_rs; 
             R_sum_u_sos R_sum_d_sos;
             R_sum_u_icibs R_sum_d_icibs];
         
figure;

set(gcf,'position',[0 0 800 600]);

bar(cat,bar_sum_5,BAR_SIZE);

xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
ylabel('95% likely sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);

legend(legend_link,'fontname',fontname,'fontsize',fontsize);

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 40]);

saveas(gcf,[root_erg_rate '95_sum_rate_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_erg_rate '95_sum_rate_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_erg_rate '95_sum_rate_M_' num2str(M) '_K_' num2str(K)],'epsc2');

[values_u,edges_u]             = histcounts(mean(rate_u),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_rs_u,edges_rs_u]       = histcounts(mean(rate_rs_u),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_sos_u,edges_sos_u]     = histcounts(mean(rate_sos_u),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_icibs_u,edges_icibs_u] = histcounts(mean(rate_icibs_u),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');

x_u = 5:7;
y_u = p_o*ones(length(x_u),1);

figure;

set(gcf,'position',[0 0 800 600]);

plot([5 edges_u 7],[0 values_u 1 1],'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot([5 edges_rs_u 7],[0 values_rs_u 1 1],'-','color',colours(2,:),'linewidth',linewidth);
plot([5 edges_sos_u 7],[0 values_sos_u 1 1],'-','color',colours(3,:),'linewidth',linewidth);
plot([5 edges_icibs_u 7],[0 values_icibs_u 1 1],'-','color',colours(4,:),'linewidth',linewidth);
plot(x_u,y_u,'--k','linewidth',linewidth);
plot(R_u_ns,p_o,'o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(R_u_rs,p_o,'o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(R_u_sos,p_o,'o','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
plot(R_u_icibs,p_o,'o','color',colours(4,:),'linewidth',linewidth,'markersize',markersize);

xlabel('Uplink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);

legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([5 7]);
ylim([0 1]);

saveas(gcf,[root_out_prob 'uplink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_out_prob 'uplink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_out_prob 'uplink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'epsc2');

[values_d, edges_d]             = histcounts(mean(rate_d),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_rs_d, edges_rs_d]       = histcounts(mean(rate_rs_d),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_sos_d, edges_sos_d]       = histcounts(mean(rate_sos_d),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_icibs_d, edges_icibs_d] = histcounts(mean(rate_icibs_d),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');

x_d = 5:10;
y_d = p_o*ones(length(x_d),1);

figure;

set(gcf,'position',[0 0 800 600]);

plot([5 edges_d 9.5],[0 values_d 1 1],'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot([5 edges_rs_d 9.5],[0 values_rs_d 1 1],'-','color',colours(2,:),'linewidth',linewidth);
plot([5 edges_sos_d 9.5],[0 values_sos_d 1 1],'-','color',colours(3,:),'linewidth',linewidth);
plot([5 edges_icibs_d 9.5],[0 values_icibs_d 1 1],'-','color',colours(4,:),'linewidth',linewidth);
plot(x_d,y_d,'--k','linewidth',linewidth);
plot(R_d_ns,p_o,'o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(R_d_rs,p_o,'o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(R_d_sos,p_o,'o','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
plot(R_d_icibs,p_o,'o','color',colours(4,:),'linewidth',linewidth,'markersize',markersize);

xlabel('Downlink rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);

legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([5 9.5]);
ylim([0 1]);

saveas(gcf,[root_out_prob 'downlink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_out_prob 'downlink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_out_prob 'downlink_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'epsc2');

bar_u = [mean(mean(rate_u,2)) mean(mean(rate_d,2));
         mean(mean(rate_rs_u,2)) mean(mean(rate_rs_d,2));
         mean(mean(rate_sos_u,2)) mean(mean(rate_sos_d,2));
         mean(mean(rate_icibs_u,2)) mean(mean(rate_icibs_d,2))];

figure;

set(gcf,'position',[0 0 800 600]);

bar(cat,bar_u,BAR_SIZE);

xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
ylabel('Average rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);

legend(legend_link,'fontname',fontname,'fontsize',fontsize);

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 10]);

saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'epsc2');

bar_u_5 = [R_u_ns R_d_ns;
           R_u_rs R_d_rs; 
           R_u_sos R_d_sos;
           R_u_icibs R_d_icibs];

figure;

set(gcf,'position',[0 0 800 600]);

bar(cat,bar_u_5,BAR_SIZE);

xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
ylabel('95% likely average rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);

legend(legend_link,'fontname',fontname,'fontsize',fontsize);

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 10]);

saveas(gcf,[root_erg_rate '95_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_erg_rate '95_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_erg_rate '95_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'epsc2');