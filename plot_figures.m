clear;
close all;
clc;

load('./results/rate_mf_M_500_K_5_MC_10000.mat');

M = 500;
K = 5;

% Ploting Figures

linewidth = 2;
fontname  = 'Times New Roman';
fontsize  = 20;

BIN_WIDTH_CDF  = 0.005;

legend_text = {'No selection', 'ICI-based selection'};

root_rate_fit = './figures/rate/fit';
root_erg_rate = './figures/rate/erg_cap';
root_out_prob = './figures/rate/out_prob';

colours = get(gca,'colororder');
close;

p = zeros(3,K);

psi_range = (0:0.01:0.4);

f = zeros(K,length(psi_range));

for k = 1:K
    curvefit = fit(psi(k,:)',rate(k,:)','poly2');
    
    p(1,k) = curvefit.p1;
    p(2,k) = curvefit.p2;
    p(3,k) = curvefit.p3;
    
    f(k,:) = p(1,k)*psi_range.^2 + p(2,k)*psi_range + p(3,k);
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot(psi(k,:),rate(k,:),'.','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(psi_range,f(k,:),'-','color',colours(2,:),'linewidth',linewidth);
    
    xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
    ylabel('Rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);

    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    xlim([0 0.4]);
    
    saveas(gcf,[root_rate_fit '_M_' num2str(M) '_' num2str(k) '_user'],'fig');
    saveas(gcf,[root_rate_fit '_M_' num2str(M) '_' num2str(k) '_user'],'png');
    saveas(gcf,[root_rate_fit '_M_' num2str(M) '_' num2str(k) '_user'],'epsc2');
end

[values, edges] = histcounts(sum(rate),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_prime, edges_prime] = histcounts(sum(rate_prime),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');

figure;

set(gcf,'position',[0 0 800 600]);

plot(edges,[values 1],'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot(edges_prime,[values_prime 1],'-','color',colours(2,:),'linewidth',linewidth);

xlabel('Sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Outage Probability','fontname',fontname,'fontsize',fontsize);

legend(legend_text,'fontname',fontname,'fontsize',fontsize,'location','southeast');

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 1]);

saveas(gcf,[root_out_prob 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_out_prob 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_out_prob 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'epsc2');

figure;

set(gcf,'position',[0 0 800 600]);

cat = categorical({'No selection','ICI-based selection'});

bar(cat,[sum(mean(rate,2)) sum(mean(rate_prime,2))],0.4);

ylabel('Average sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);

set(gca,'fontname',fontname,'fontsize',fontsize);

saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'epsc2');

[values, edges] = histcounts(mean(rate),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');
[values_prime, edges_prime] = histcounts(mean(rate_prime),'binwidth',BIN_WIDTH_CDF,'normalization','cdf');

figure;

set(gcf,'position',[0 0 800 600]);

plot(edges,[values 1],'-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot(edges_prime,[values_prime 1],'-','color',colours(2,:),'linewidth',linewidth);

xlabel('Rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
ylabel('Outage Probability','fontname',fontname,'fontsize',fontsize);

legend(legend_text,'fontname',fontname,'fontsize',fontsize,'location','southeast');

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 1]);

saveas(gcf,[root_out_prob 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_out_prob 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_out_prob 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'epsc2');

figure;

set(gcf,'position',[0 0 800 600]);

cat = categorical({'No selection','ICI-based selection'});

bar(cat,[mean(mean(rate,2)) mean(mean(rate_prime,2))],0.4);

ylabel('Average rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);

set(gca,'fontname',fontname,'fontsize',fontsize);

saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'fig');
saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'png');
saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'epsc2');
