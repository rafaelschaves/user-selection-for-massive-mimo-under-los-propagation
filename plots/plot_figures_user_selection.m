clear;
close all;
clc;

% Macros

MC   = 1000;                                                              % Size of the monte-carlo ensemble
N_MC = 1;

M = 50;                                                                  % Number of antennas at base station
K = 75;                                                                   % Number of users at the cell 

if K > M
    L_max = M;
else
    L_max = K-1;
end

N_ALG = 3;                                                                 % Number of algorithms for perform user scheduling                                                                           % Number of algorithms for perform user scheduling
N_PRE = 3;

snr = -5;

bandwidth   = 20e6;
dl_ul_ratio = 0.5;

% Roots

root_load = 'G:\My Drive\UFRJ\PhD\Codes\user-scheduling-massive-mimo\Results\Selection\Downlink\';
root_save = 'G:\My Drive\UFRJ\PhD\Codes\user-scheduling-massive-mimo\Results\Figures\Downlink\';

zero_pad_1 = '%03d';
zero_pad_2 = '%02d';

% Loading data

se_all_mc     = zeros(K,N_PRE,MC*N_MC);    
se_s_L_all_mc = zeros(L_max,L_max,N_PRE,N_ALG,MC*N_MC);
sum_se_s      = zeros(L_max,N_PRE,N_ALG,MC*N_MC);

for n_mc = 1:N_MC
    load([root_load 'spectral_efficiency_all_L_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC) '_' sprintf(zero_pad_2,n_mc) '.mat']);
     
    idx_ini = (n_mc - 1)*MC + 1;
    idx_end = n_mc*MC;
        
    se_all_mc(:,:,idx_ini:idx_end)         = bandwidth*dl_ul_ratio*se;
    se_s_L_all_mc(:,:,:,:,idx_ini:idx_end) = bandwidth*dl_ul_ratio*se_s_all_L;
    
    clear se se_s_all_L;
end

sum_se = reshape(sum(se_all_mc,1),N_PRE,MC*N_MC);

for L = 1:L_max
    sum_se_s(L,:,:,:) = sum(se_s_L_all_mc(:,L,:,:,:),1);
end

avg_sum_se   = mean(sum_se,2);
avg_sum_se_s = mean(sum_se_s,4);

[max_sum_se_s,L_star] = max(avg_sum_se_s,[],1);

max_sum_se_s = reshape(max_sum_se_s,N_PRE,N_ALG)';
L_star       = reshape(L_star,N_PRE,N_ALG)';

L = ceil(K/2);

for n_alg = 1:N_ALG
    for n_pre = 1:N_PRE
        sum_se_s_star(:,n_pre,n_alg) = sum_se_s(L_star(n_alg,n_pre),n_pre,n_alg,:);
    end    
end

N_BIN = 100;

cdf_sum_se = cell(N_ALG+1,N_PRE);
edg_sum_se = cell(N_ALG+1,N_PRE);

for n_pre = 1:N_PRE
    [cdf_sum_se{1,n_pre},edg_sum_se{1,n_pre}] = histcounts(sum_se(n_pre,:),N_BIN,'normalization','cdf');

    for n_alg = 1:N_ALG
        [cdf_sum_se{n_alg+1,n_pre},edg_sum_se{n_alg+1,n_pre}] = histcounts(sum_se_s_star(:,n_pre,n_alg),N_BIN,'normalization','cdf');
    end
end

% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 30;

marker    = {'o','s','^'};
linestyle = {'-','--',':'};

savefig = 0;
plotse  = 0;

if M == 50
    OM = 1e-6;
elseif M == 100
    OM = 1e-6;
end

% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_algo           = {'SOS','CBS','ICIBS'};
legend_algo_plus_prec = {'SOS','CBS','ICIBS','MRT','ZF','MMSE'};

um = {'(kbps)','(Mbps)','(Gbps)'};

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
           0.6350 0.0780 0.1840;
           0.0000 0.0000 0.0000];

figure;
       
set(gcf,'position',[0 0 800 600]);

if K <= M
    % Plot for user selection algorithm legends
    plot(1:K,OM*[avg_sum_se_s(:,1,1); avg_sum_se(1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(1:K,OM*[avg_sum_se_s(:,1,2); avg_sum_se(1)],'-' ,'color',colours(2,:),'linewidth',linewidth);
    plot(1:K,OM*[avg_sum_se_s(:,1,3); avg_sum_se(1)],'-' ,'color',colours(3,:),'linewidth',linewidth);
    % Plot for precoding algorithm legends
    plot(1:K,OM*[avg_sum_se_s(:,1,1); avg_sum_se(1)],'-' ,'color',colours(8,:),'linewidth',linewidth);
    plot(1:K,OM*[avg_sum_se_s(:,2,1); avg_sum_se(2)],'--','color',colours(8,:),'linewidth',linewidth);
    plot(1:K,OM*[avg_sum_se_s(:,3,1); avg_sum_se(3)],':' ,'color',colours(8,:),'linewidth',linewidth);
    % Plot for results
    plot(1:K,OM*[avg_sum_se_s(:,1,1); avg_sum_se(1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
    plot(1:K,OM*[avg_sum_se_s(:,1,2); avg_sum_se(1)],'-' ,'color',colours(2,:),'linewidth',linewidth);
    plot(1:K,OM*[avg_sum_se_s(:,1,3); avg_sum_se(1)],'-' ,'color',colours(3,:),'linewidth',linewidth);
    plot(1:K,OM*[avg_sum_se_s(:,2,1); avg_sum_se(2)],'--','color',colours(1,:),'linewidth',linewidth);
    plot(1:K,OM*[avg_sum_se_s(:,2,2); avg_sum_se(2)],'--','color',colours(2,:),'linewidth',linewidth);
    plot(1:K,OM*[avg_sum_se_s(:,2,3); avg_sum_se(2)],'--','color',colours(3,:),'linewidth',linewidth);
    plot(1:K,OM*[avg_sum_se_s(:,3,1); avg_sum_se(3)],':' ,'color',colours(1,:),'linewidth',linewidth);
    plot(1:K,OM*[avg_sum_se_s(:,3,2); avg_sum_se(3)],':' ,'color',colours(2,:),'linewidth',linewidth);
    plot(1:K,OM*[avg_sum_se_s(:,3,3); avg_sum_se(3)],':' ,'color',colours(3,:),'linewidth',linewidth);
else
    % Plot for user selection algorithm legends
    plot(1:L_max,OM*avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(1:L_max,OM*avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
    plot(1:L_max,OM*avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
    % Plot for precoding algorithm legends
    plot(1:L_max,OM*avg_sum_se_s(:,1,1),'-' ,'color',colours(8,:),'linewidth',linewidth);
    plot(1:L_max,OM*avg_sum_se_s(:,2,1),'--','color',colours(8,:),'linewidth',linewidth);
    plot(1:L_max,OM*avg_sum_se_s(:,3,1),':' ,'color',colours(8,:),'linewidth',linewidth);
    % Plot for results
    plot(1:L_max,OM*avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
    plot(1:L_max,OM*avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
    plot(1:L_max,OM*avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
    plot(1:L_max,OM*avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
    plot(1:L_max,OM*avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
    plot(1:L_max,OM*avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
    plot(1:L_max,OM*avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
    plot(1:L_max,OM*avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
    plot(1:L_max,OM*avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
end

xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);

if plotse == 1
    ylabel('Sum-spectral efficiency','fontname',fontname,'fontsize',fontsize);
else
    ylabel(['Throughput ' um{abs(log10(OM))/3}],'fontname',fontname,'fontsize',fontsize);
end

if M == 50 && K == 10
    legend(legend_algo_plus_prec,'fontname',fontname,'fontsize',fontsize,'location',location_4,'numcolumns',2);
    legend box off;
end

set(gca,'fontname',fontname,'fontsize',fontsize);

if K <= M
    xlim([1 K]);
else
    xlim([1 L_max]);
end

switch M
    case 50
        if plotse == 1
            ylim([0 20]);
        else
            ylim([0 OM*bandwidth*dl_ul_ratio*20]);
        end
        
        switch K
            case 10
                dim = [0.7 0.6 0.2 0.07];
                annotation('ellipse',dim,'linewidth',linewidth);
                
                axes('position',[.24 .585 .3 .3]);
                box on;
                
                plot(1:K,OM*[avg_sum_se_s(:,1,1); avg_sum_se(1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:K,OM*[avg_sum_se_s(:,1,2); avg_sum_se(1)],'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,1,3); avg_sum_se(1)],'-' ,'color',colours(3,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,2,1); avg_sum_se(2)],'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,2,2); avg_sum_se(2)],'--','color',colours(2,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,2,3); avg_sum_se(2)],'--','color',colours(3,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,3,1); avg_sum_se(3)],':' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,3,2); avg_sum_se(3)],':' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,3,3); avg_sum_se(3)],':' ,'color',colours(3,:),'linewidth',linewidth);
                
                set(gca,'fontname',fontname,'fontsize',fontsize);
                
                xlim([7 10]);
                
                if plotse == 1
                    ylim([12 12.5]);
                end
            case 25
                dim = [0.56 0.725 0.2 0.07];
                annotation('ellipse',dim,'linewidth',linewidth);
                
                axes('position',[.4 .275 .3 .3]);
                box on;
                
                plot(1:L_max,OM*avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:L_max,OM*avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
                
                set(gca,'fontname',fontname,'fontsize',fontsize);
                
                xlim([15 20]);
                
                if plotse == 1
                    ylim([15.5 16]);
                end
            case 50
                dim = [0.425 0.8 0.175 0.07];
                annotation('ellipse',dim,'linewidth',linewidth);
                
                axes('position',[.325 .275 .3 .3]);
                box on;
                
                plot(1:L_max,OM*avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:L_max,OM*avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
                
                set(gca,'fontname',fontname,'fontsize',fontsize);
                
                xlim([20 30]);
                
                if plotse == 1
                    ylim([17 18]);
                end
            case 75
                dim = [0.5 0.825 0.2 0.07];
                annotation('ellipse',dim,'linewidth',linewidth);
                
                axes('position',[.3 .275 .3 .3]);
                box on;
                
                plot(1:L_max,OM*avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:L_max,OM*avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
                
                set(gca,'fontname',fontname,'fontsize',fontsize);
                
                xlim([25 36]);
                
                if plotse == 1
                    ylim([17.5 18.55]);
                end
        end
    case 100
        if plotse == 1
            ylim([0 40]);
        else
            ylim([0 OM*bandwidth*dl_ul_ratio*40]);
        end
        
        switch K
            case 10
                dim = [0.75 0.49 0.15 0.075];
                annotation('ellipse',dim,'linewidth',linewidth);
                
                axes('position',[.25 .575 .3 .3]);
                box on;
                
                plot(1:K,OM*[avg_sum_se_s(:,1,1); avg_sum_se(1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:K,OM*[avg_sum_se_s(:,1,2); avg_sum_se(1)],'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,1,3); avg_sum_se(1)],'-' ,'color',colours(3,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,2,1); avg_sum_se(2)],'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,2,2); avg_sum_se(2)],'--','color',colours(2,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,2,3); avg_sum_se(2)],'--','color',colours(3,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,3,1); avg_sum_se(3)],':' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,3,2); avg_sum_se(3)],':' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_se_s(:,3,3); avg_sum_se(3)],':' ,'color',colours(3,:),'linewidth',linewidth);
                
                set(gca,'fontname',fontname,'fontsize',fontsize);
                
                xlim([7 10]);
                
                if plotse == 1
                    ylim([18 19.5]);
                end
            case 25
                dim = [0.7 0.65 0.15 0.05];
                annotation('ellipse',dim,'linewidth',linewidth);
                
                axes('position',[.575 .275 .3 .3]);
                box on;
                
                plot(1:L_max,OM*avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:L_max,OM*avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
                
                set(gca,'fontname',fontname,'fontsize',fontsize);
                
                xlim([18 22]);
                
                if plotse == 1
                    ylim([26 27]);
                end
            case 50
                dim = [0.6 0.75 0.15 0.05];
                annotation('ellipse',dim,'linewidth',linewidth);
                
                axes('position',[.45 .275 .3 .3]);
                box on;
                
                plot(1:L_max,OM*avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:L_max,OM*avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
                
                set(gca,'fontname',fontname,'fontsize',fontsize);
                
                xlim([31 37]);
                
                if plotse == 1
                    ylim([31 32]);
                end
            case 75
                dim = [0.525 0.775 0.15 0.05];
                annotation('ellipse',dim,'linewidth',linewidth);
                
                axes('position',[.35 .275 .3 .3]);
                box on;
                
                plot(1:L_max,OM*avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:L_max,OM*avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
                
                set(gca,'fontname',fontname,'fontsize',fontsize);
                
                xlim([35 50]);
                
                if plotse == 1
                    ylim([33 34.3]);
                end
            case 100
                dim = [0.45 0.8 0.15 0.05];
                annotation('ellipse',dim,'linewidth',linewidth);
                
                axes('position',[.325 .275 .3 .3]);
                box on;
                
                plot(1:L_max,OM*avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:L_max,OM*avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
                
                set(gca,'fontname',fontname,'fontsize',fontsize);
                
                xlim([45 60]);
                
                if plotse == 1
                    ylim([34.2 35.6]);
                end
            case 150
                dim = [0.5 0.825 0.225 0.05];
                annotation('ellipse',dim,'linewidth',linewidth);
                
                axes('position',[.35 .3 .3 .3]);
                box on;
                
                plot(1:L_max,OM*avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:L_max,OM*avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
                
                set(gca,'fontname',fontname,'fontsize',fontsize);
                
                xlim([55 76]);
                
                if plotse == 1
                    ylim([35 37.1]);
                end
        end
end
 
if savefig == 1
    if plotse == 1
        saveas(gcf,[root_save 'sum_se_all_L_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
        saveas(gcf,[root_save 'sum_se_all_L_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
        saveas(gcf,[root_save 'sum_se_all_L_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
    else
        saveas(gcf,[root_save 'throughput_all_L_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
        saveas(gcf,[root_save 'throughput_all_L_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
        saveas(gcf,[root_save 'throughput_all_L_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
    end
end

figure;

set(gcf,'position',[0 0 800 600]);

% Plots for user selection algorihtm legends
plot(OM*edg_sum_se{1,1},[cdf_sum_se{1,1} 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
hold on;
plot(OM*edg_sum_se{1,2},[cdf_sum_se{1,2} 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
plot(OM*edg_sum_se{1,3},[cdf_sum_se{1,3} 1],'-' ,'color',colours(3,:),'linewidth',linewidth);
% Plots for precoding algorithm legends
plot(OM*edg_sum_se{1,1},[cdf_sum_se{1,1} 1],'-' ,'color',colours(8,:),'linewidth',linewidth);
plot(OM*edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'--','color',colours(8,:),'linewidth',linewidth);
plot(OM*edg_sum_se{3,1},[cdf_sum_se{3,1} 1],':' ,'color',colours(8,:),'linewidth',linewidth);
% Plots for results
plot(OM*edg_sum_se{1,1},[cdf_sum_se{1,1} 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
plot(OM*edg_sum_se{1,2},[cdf_sum_se{1,2} 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
plot(OM*edg_sum_se{1,3},[cdf_sum_se{1,3} 1],'-' ,'color',colours(3,:),'linewidth',linewidth);
plot(OM*edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'--','color',colours(1,:),'linewidth',linewidth);
plot(OM*edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
plot(OM*edg_sum_se{2,3},[cdf_sum_se{2,3} 1],'--','color',colours(3,:),'linewidth',linewidth);
plot(OM*edg_sum_se{3,1},[cdf_sum_se{3,1} 1],':' ,'color',colours(1,:),'linewidth',linewidth);
plot(OM*edg_sum_se{3,2},[cdf_sum_se{3,2} 1],':' ,'color',colours(2,:),'linewidth',linewidth);
plot(OM*edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':' ,'color',colours(3,:),'linewidth',linewidth);

if plotse == 1
    xlabel('Sum-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
else
    xlabel(['Throughput ' um{abs(log10(OM))/3}],'fontname',fontname,'fontsize',fontsize);
end

ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);

%if K == 10
    legend(legend_algo_plus_prec,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',2);
    legend box off;
%end

set(gca,'fontname',fontname,'fontsize',fontsize);

ylim([0 1]);

switch M
    case 50
        switch K
            case 10
                switch L
                    case ceil(K/5)
                        xlim([63.8 63.92]);
                    case ceil(K/2)
                        xlim([103 104.75]);
                        % xlim([103 105]);
                end
            case 25
                switch L
                    case ceil(K/5)
                        xlim([104.2 104.65]);
                    case ceil(K/2)
                        xlim([148 152.5]);
                end
            case 50
                switch L
                    case ceil(K/5)
                        xlim([138.5 139.7]);
                    case ceil(K/2)
                        xlim([150 182]);
                end
            case 75
                switch L
                    case ceil(K/5)
                        xlim([157 159.4]);
                    case ceil(K/2)
                        xlim([50 207]);
                end
        end
    case 100
        switch K
            case 10
                switch L
                    case ceil(K/5)
                        xlim([82.2 82.3]);
                    case ceil(K/2)
                        xlim([144.5 145.7]);
                end
            case 25
                switch L
                    case ceil(K/5)
                        xlim([145.4 145.63]);
                    case ceil(K/2)
                        xlim([232 235.5]);
                end
            case 50
                switch L
                    case ceil(K/5)
                        xlim([208.5 209.25]);
                    case ceil(K/2)
                        xlim([295 301])
                end
            case 75
                switch L
                    case ceil(K/5)
                        xlim([248.5 250.1]);
                    case ceil(K/2)
                        xlim([320 340]);
                end
            case 100
                switch L
                    case ceil(K/5)
                        xlim([277 279.3]);
                    case ceil(K/2)
                        xlim([300 359]);
                end
            case 150
                switch L
                    case ceil(K/5)
                        xlim([315 318.5]);
                    case ceil(K/2)
                        xlim([100 385])
                end
        end        
end

if savefig == 1
    if plotse == 1
        saveas(gcf,[root_save 'cdf_sum_se_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
        saveas(gcf,[root_save 'cdf_sum_se_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
        saveas(gcf,[root_save 'cdf_sum_se_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
    else
        saveas(gcf,[root_save 'cdf_throughput_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
        saveas(gcf,[root_save 'cdf_throughput_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
        saveas(gcf,[root_save 'cdf_throughput_ur_los_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
    end
end