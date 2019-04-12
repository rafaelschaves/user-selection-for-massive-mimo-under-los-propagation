clear;
close all;
clc;

load('./results/example_rate_mf_M_500_K_4_MC_10000.mat');

rate_4 = rate;

load('./results/example_rate_mf_M_500_K_5_MC_10000.mat');

rate_5 = rate;

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

BIN_WIDTH_CDF  = 0.005;

legend_text = {'$K = 4$','$K = 5$'};
location = 'northwest';

% cat = categorical(legend_text);
% cat = reordercats(cat,legend_text);

root_example = './figures/rate/example_';

colours = get(gca,'colororder');
close;

p_o = 0.05;

R_4 = 6.425;
R_5 = 6.185;

R_sum_4 = 25.7;
R_sum_5 = 30.92;

x = 20:40;
y = p_o*ones(length(x),1);

[values_4, edges_4] = histcounts(sum(rate_4),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_5, edges_5] = histcounts(sum(rate_5),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');

figure;

set(gcf,'position',[0 0 800 600]);

plot([20 edges_4 40],[0 values_4 1 1],'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot([20 edges_5],[0 values_5 1],'-','color',colours(2,:),'linewidth',linewidth);
plot(x,y,'--k','linewidth',linewidth);
plot(R_sum_4,p_o,'o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(R_sum_5,p_o,'o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);

xlabel('Sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Outage Probability','fontname',fontname,'fontsize',fontsize);

legend(legend_text,'fontname',fontname,'fontsize',fontsize,'location',location,'interpreter','latex');

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 1]);

saveas(gcf,[root_example 'sum_rate'],'fig');
saveas(gcf,[root_example 'sum_rate'],'png');
saveas(gcf,[root_example 'sum_rate'],'epsc2');

x = 5:9;
y = 0.05*ones(length(x),1);

[values_4, edges_4] = histcounts(mean(rate_4),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_5, edges_5] = histcounts(mean(rate_5),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');

figure;

set(gcf,'position',[0 0 800 600]);

plot([5 edges_4 9],[0 values_4 1 1],'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot([5 edges_5 9],[0 values_5 1 1],'-','color',colours(2,:),'linewidth',linewidth);
plot(x,y,'--k','linewidth',linewidth);
plot(R_4,p_o,'o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(R_5,p_o,'o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);

xlabel('Rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Outage Probability','fontname',fontname,'fontsize',fontsize);

legend(legend_text,'fontname',fontname,'fontsize',fontsize,'location',location,'interpreter','latex');

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 1]);

saveas(gcf,[root_example 'avg_rate_ter'],'fig');
saveas(gcf,[root_example 'avg_rate_ter'],'png');
saveas(gcf,[root_example 'avg_rate_ter'],'epsc2');