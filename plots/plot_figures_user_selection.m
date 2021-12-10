clear;
close all;
clc;

% Macros

<<<<<<< HEAD
MC   = 1000;                                                                                                                                          % Size of the monte-carlo ensemble
N_MC = 2;

M = 100;                                                                                                                                              % Number of antennas at base station
K = 10;                                                                                                                                              % Number of users at the cell 
=======
MC   = 1000;                                                              % Size of the monte-carlo ensemble
N_MC = 1;

M = 50;                                                                  % Number of antennas at base station
K = 75;                                                                   % Number of users at the cell 
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0

% M = 50  & K = [10 25 50 75]
% M = 100 & K = [10 25 50 75 100 150] 
% M = 200 & K = [10 25 50 75 100 150 200 250]

if K > M
    L_max = M;
else
    L_max = K-1;
end

<<<<<<< HEAD
N_ALG = 3;                                                                                                                                            % Number of algorithms for perform user scheduling
N_PRE = 3;

snr = -5;

bandwidth   = 20e6;
dl_ul_ratio = 0.5;

% Roots

%root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Downlink/';
root_load = 'D:\PhD\user-selection\Perfect CSI\';
=======
snr = -5;

N_ALG = 3;                                                                 % Number of algorithms for perform user scheduling
N_PRE = 3;

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Downlink/';
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Selection/Downlink/';

zero_pad_1 = '%03d';
zero_pad_2 = '%02d';

chn_type = 'ur_los';

% Loading data

<<<<<<< HEAD
se_all_mc     = zeros(K,N_PRE,MC*N_MC);    
se_s_L_all_mc = zeros(L_max,L_max,N_PRE,N_ALG,MC*N_MC);
sum_se_s      = zeros(L_max,N_PRE,N_ALG,MC*N_MC);

for n_mc = 1:N_MC
    load([root_load 'spectral_efficiency_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' ...
          num2str(MC) '_' sprintf(zero_pad_2,n_mc) '.mat']);
=======
se_all_mc   = zeros(K,N_PRE,MC*N_MC);    
se_s_L_all_mc = zeros(L_max,L_max,N_PRE,N_ALG,MC*N_MC);

sum_se_s      = zeros(L_max,N_PRE,N_ALG,MC*N_MC);

for n_mc = 1:N_MC
    load([root_load 'spectral_efficiency_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC) '_' sprintf(zero_pad_2,n_mc) '.mat']);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
     
    idx_ini = (n_mc - 1)*MC + 1;
    idx_end = n_mc*MC;
        
<<<<<<< HEAD
    se_all_mc(:,:,idx_ini:idx_end)         = bandwidth*dl_ul_ratio*se;
    se_s_L_all_mc(:,:,:,:,idx_ini:idx_end) = bandwidth*dl_ul_ratio*se_s_all_L;
=======
    se_all_mc(:,:,idx_ini:idx_end)         = se;
    se_s_L_all_mc(:,:,:,:,idx_ini:idx_end) = se_s_all_L;
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
    
    clear se se_s_all_L;
end

sum_se = reshape(sum(se_all_mc,1),N_PRE,MC*N_MC);

<<<<<<< HEAD
for l = 1:L_max
    sum_se_s(l,:,:,:) = sum(se_s_L_all_mc(:,l,:,:,:),1);
=======
for L = 1:L_max
    sum_se_s(L,:,:,:) = sum(se_s_L_all_mc(:,L,:,:,:),1);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
end

avg_sum_se   = mean(sum_se,2);
avg_sum_se_s = mean(sum_se_s,4);

<<<<<<< HEAD
[max_sum_se_s,L_star] = max(avg_sum_se_s,[],1); 
=======
[max_sum_se_s,L_star] = max(avg_sum_se_s,[],1);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0

max_sum_se_s = reshape(max_sum_se_s,N_PRE,N_ALG)';
L_star       = reshape(L_star,N_PRE,N_ALG)';

<<<<<<< HEAD
L = ceil(K/2);

N_BIN = 100;

cdf_sum_se = cell(N_PRE,N_ALG);
edg_sum_se = cell(N_PRE,N_ALG);

for n_alg = 1:N_ALG
    for n_pre = 1:N_PRE
        % [cdf_sum_se{n_pre,n_alg},edg_sum_se{n_pre,n_alg}] = histcounts(sum_se_s(L_star(n_alg,n_pre),n_pre,n_alg,:),N_BIN,'normalization','cdf');
        [cdf_sum_se{n_pre,n_alg},edg_sum_se{n_pre,n_alg}] = histcounts(sum_se_s(L,n_pre,n_alg,:),N_BIN,'normalization','cdf');
=======
% L_star = 80*ones(3,3);

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
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
    end
end

% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 30;

<<<<<<< HEAD
marker    = {'o','s','^'};
linestyle = {'-','--',':'};

savefig = 1;
plotse  = 0;

if M == 50
    OM = 1e-6;
elseif M == 100
    OM = 1e-6;
end

=======
marker = {'o','s','^'};

linestyle = {'-','--',':'};

savefig = 0;

% NS - No selection
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

<<<<<<< HEAD
legend_algo           = {'SOS','CBS','ICIBS'};
legend_algo_plus_prec = {'SOS','CBS','ICIBS','MRT','ZF','MMSE'};

um = {'(kbps)','(Mbps)','(Gbps)'};
=======
legend_algo   = {'NS','SOS','CBS','ICIBS'};
legend_algo_2 = {'SOS','CBS','ICIBS'};
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0

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
<<<<<<< HEAD
           0.6350 0.0780 0.1840;
           0.0000 0.0000 0.0000];
=======
           0.6350 0.0780 0.1840];
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0

figure;
       
set(gcf,'position',[0 0 800 600]);

% if K <= M
%     plot(1:L_max,avg_sum_se(1)*ones(L_max,1),'-k','linewidth',linewidth);
%     hold on;
%     plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se(2)*ones(L_max,1),'--k','linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se(3)*ones(L_max,1),':k','linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
% else
%     plot(1:L_max,avg_sum_se(1)*ones(L_max,1),'-k','linewidth',linewidth);
%     hold on;
%     plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
%     plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
% end

if K <= M
<<<<<<< HEAD
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
=======
    plot(1:K,[avg_sum_se_s(:,1,1); avg_sum_se(1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(1:K,[avg_sum_se_s(:,1,2); avg_sum_se(1)],'-' ,'color',colours(2,:),'linewidth',linewidth);
    plot(1:K,[avg_sum_se_s(:,1,3); avg_sum_se(1)],'-' ,'color',colours(3,:),'linewidth',linewidth);
    plot(1:K,[avg_sum_se_s(:,2,1); avg_sum_se(2)],'--','color',colours(1,:),'linewidth',linewidth);
    plot(1:K,[avg_sum_se_s(:,2,2); avg_sum_se(2)],'--','color',colours(2,:),'linewidth',linewidth);
    plot(1:K,[avg_sum_se_s(:,2,3); avg_sum_se(2)],'--','color',colours(3,:),'linewidth',linewidth);
    plot(1:K,[avg_sum_se_s(:,3,1); avg_sum_se(3)],':' ,'color',colours(1,:),'linewidth',linewidth);
    plot(1:K,[avg_sum_se_s(:,3,2); avg_sum_se(3)],':' ,'color',colours(2,:),'linewidth',linewidth);
    plot(1:K,[avg_sum_se_s(:,3,3); avg_sum_se(3)],':' ,'color',colours(3,:),'linewidth',linewidth);
else
    plot(1:L_max,avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(1:L_max,avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
    plot(1:L_max,avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
    plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
    plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
    plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
    plot(1:L_max,avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
    plot(1:L_max,avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
    plot(1:L_max,avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
end

xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
ylabel('Sum-spectral efficiency','fontname',fontname,'fontsize',fontsize);

if K == 10
    legend(legend_algo_2,'fontname',fontname,'fontsize',fontsize,'location',location_4);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
    legend box off;
end

set(gca,'fontname',fontname,'fontsize',fontsize);

if K <= M
    xlim([1 K]);
else
    xlim([1 L_max]);
end

switch chn_type
    case 'rayleigh'
    case 'ur_los'
        switch M
            case 50
<<<<<<< HEAD
                if plotse == 1
                    ylim([0 20]);
                else
                    ylim([0 OM*bandwidth*dl_ul_ratio*20]);
                end
=======
                ylim([0 20]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                
                switch K
                    case 10                        
                        dim = [0.7 0.6 0.2 0.07];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
                        axes('position',[.24 .585 .3 .3]);
                        box on;
                            
<<<<<<< HEAD
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
=======
                        plot(1:K,[avg_sum_se_s(:,1,1); avg_sum_se(1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:K,[avg_sum_se_s(:,1,2); avg_sum_se(1)],'-' ,'color',colours(2,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,1,3); avg_sum_se(1)],'-' ,'color',colours(3,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,2,1); avg_sum_se(2)],'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,2,2); avg_sum_se(2)],'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,2,3); avg_sum_se(2)],'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,3,1); avg_sum_se(3)],':' ,'color',colours(1,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,3,2); avg_sum_se(3)],':' ,'color',colours(2,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,3,3); avg_sum_se(3)],':' ,'color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        % xticks([7 10]);
                        % xticklabels({'0','1.25','2.5'})
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([7 10]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([12 12.5]);
                        end
=======
                        ylim([12 12.5]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 25
                        dim = [0.56 0.725 0.2 0.07];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
<<<<<<< HEAD
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
=======
                        axes('position',[.35 .275 .3 .3]);
                        box on;
                        
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([15 20]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([15.5 16]);
                        end
=======
                        ylim([15.5 16]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 50
                        dim = [0.425 0.8 0.175 0.07];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
<<<<<<< HEAD
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
=======
                        axes('position',[.275 .275 .3 .3]);
                        box on;
                        
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':' ,'color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':' ,'color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':' ,'color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([20 30]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([17 18]);
                        end
=======
                        ylim([17 18]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 75
                        dim = [0.5 0.825 0.2 0.07];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
                        axes('position',[.3 .275 .3 .3]);
                        box on;
                        
<<<<<<< HEAD
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
=======
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([25 36]);
<<<<<<< HEAD
                        
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
=======
                        ylim([17.5 18.55]);
                end
            case 100
                ylim([0 40]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                
                switch K
                    case 10                        
                        dim = [0.75 0.49 0.15 0.075];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
                        axes('position',[.25 .575 .3 .3]);
                        box on;
                        
<<<<<<< HEAD
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
=======
                        plot(1:K,[avg_sum_se_s(:,1,1); avg_sum_se(1)],'-' ,'color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:K,[avg_sum_se_s(:,1,2); avg_sum_se(1)],'-' ,'color',colours(2,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,1,3); avg_sum_se(1)],'-' ,'color',colours(3,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,2,1); avg_sum_se(2)],'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,2,2); avg_sum_se(2)],'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,2,3); avg_sum_se(2)],'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,3,1); avg_sum_se(3)],':' ,'color',colours(1,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,3,2); avg_sum_se(3)],':' ,'color',colours(2,:),'linewidth',linewidth);
                        plot(1:K,[avg_sum_se_s(:,3,3); avg_sum_se(3)],':' ,'color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([7 10]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([18 19.5]);
                        end
=======
                        ylim([18 19.5]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 25                 
                        dim = [0.7 0.65 0.15 0.05];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
<<<<<<< HEAD
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
=======
                        axes('position',[.55 .275 .3 .3]);
                        box on;
                        
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([18 22]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([26 27]);
                        end
=======
                        ylim([26 27]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 50   
                        dim = [0.6 0.75 0.15 0.05];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
                        axes('position',[.45 .275 .3 .3]);
                        box on;
                        
<<<<<<< HEAD
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
=======
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([31 37]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([31 32]);
                        end
=======
                        ylim([31 32]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 75
                        dim = [0.525 0.775 0.15 0.05];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
<<<<<<< HEAD
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
=======
                        axes('position',[.275 .275 .3 .3]);
                        box on;
                        
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([35 50]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([33 34.3]);
                        end
=======
                        ylim([33 34.3]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 100
                        dim = [0.45 0.8 0.15 0.05];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
                        axes('position',[.325 .275 .3 .3]);
                        box on;
                        
<<<<<<< HEAD
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
=======
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([45 60]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([34.2 35.6]);
                        end
=======
                        ylim([34.2 35.6]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 150
                        dim = [0.5 0.825 0.225 0.05];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
<<<<<<< HEAD
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
=======
                        axes('position',[.3 .3 .3 .3]);
                        box on;
                        
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([55 76]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([35 37.1]);
                        end
=======
                        ylim([35 37.1]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                end
            case 200
                switch K
                    case 10
                    case 25
                        ylim([5 43]);
                        
                        dim = [0.75 0.875 0.15 0.05];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
                        axes('position',[.475 .3 .3 .3]);
                        box on;
                        
<<<<<<< HEAD
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
=======
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([20 23]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([41.5 42.6]);
                        end
=======
                        ylim([41.5 42.6]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 50
                        ylim([5 55]);
                        
                        dim = [0.7 0.875 0.1 0.05];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
                        axes('position',[.45 .275 .3 .3]);
                        box on;
                        
<<<<<<< HEAD
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
=======
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([36 42]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([52.5 54]);
                        end
=======
                        ylim([52.5 54]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 75
                        ylim([5 60]);
                        
                        dim = [0.655 0.875 0.125 0.05];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
                        axes('position',[.45 .275 .3 .3]);
                        box on;
                        
<<<<<<< HEAD
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
=======
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([52 58]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([59.5 60.1]);
                        end
=======
                        ylim([59.5 60.1]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 100
                        ylim([5 64]);
                        
                        dim = [0.625 0.875 0.125 0.05];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
                        axes('position',[.45 .275 .3 .3]);
                        box on;
                        
<<<<<<< HEAD
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
=======
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([64 72]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([62 64]);
                        end
=======
                        ylim([62 64]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 150
                        ylim([5 69]);
                        
                        dim = [0.55 0.85 0.125 0.075];
                        annotation('ellipse',dim,'linewidth',linewidth);
                        
                        axes('position',[.4 .3 .3 .3]);
                        box on;
                        
<<<<<<< HEAD
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
=======
                        plot(1:L_max,avg_sum_se_s(:,1,1),'-','color',colours(1,:),'linewidth',linewidth);
                        hold on;
                        plot(1:L_max,avg_sum_se_s(:,1,2),'-','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,1,3),'-','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,1),'--','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,2),'--','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,2,3),'--','color',colours(3,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,1),':','color',colours(1,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,2),':','color',colours(2,:),'linewidth',linewidth);
                        plot(1:L_max,avg_sum_se_s(:,3,3),':','color',colours(3,:),'linewidth',linewidth);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                        
                        set(gca,'fontname',fontname,'fontsize',fontsize);
                        
                        xlim([80 95]);
<<<<<<< HEAD
                        
                        if plotse == 1
                            ylim([66 69]);
                        end
=======
                        ylim([66 69]);
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
                    case 200
                    case 250
                end
        end
    otherwise
end

<<<<<<< HEAD
if savefig == 1
    if plotse == 1
        saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
        saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
        saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
    else
        saveas(gcf,[root_save 'throughput_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
        saveas(gcf,[root_save 'throughput_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
        saveas(gcf,[root_save 'throughput_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
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

switch chn_type
    case 'ur_los'
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
            case 200
                switch K
                    case 10
                    case 25
                    case 50
                    case 75
                    case 100
                    case 150
                    case 200
                    case 250
                end
        end
    otherwise
end

=======
if (savefig == 1)
    saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
    saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
    saveas(gcf,[root_save 'sum_se_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
end

% figure;

% set(gcf,'position',[0 0 800 600]);
% 
% plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
% hold on;
% plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
% plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
% plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
% plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
% plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
% plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
% plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
% plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
% 
% xlabel('Sum-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_algo_2,'fontname',fontname,'fontsize',fontsize,'location',location_1);
% legend box off;
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
% switch chn_type
%     case 'rayleigh'
%     case 'ur_los'
%         switch M
%             case 50
%                 switch K
%                     case 10
%                         xlim([9 13])
%                         ylim([0 1]);
%                         
%                         dim = [0.2 0.2 0.4 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([9.5 11]);
%                         ylim([0.01 0.05]);
%                     case 25
%                         xlim([13 16.5]);
%                         ylim([0 1]);
%                         
%                         dim = [0.225 0.2 0.35 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([13.5 15]);
%                         ylim([0.03 0.05]);
%                     case 50
%                         xlim([15 18.2]);
%                         ylim([0 1]);
%                         
%                         dim = [0.375 0.2 0.3 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.275 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([16 17.3]);
%                         ylim([0.03 0.05]);
%                     case 75
%                         xlim([15.5 19]);
%                         ylim([0 1]);
%                         
%                         dim = [0.425 0.2 0.275 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([17.2 18]);
%                         ylim([0.03 0.05]);
%                 end
%             case 100
%                 switch K
%                     case 10
%                         xlim([14 20])
%                         ylim([0 1]);
%                         
%                         dim = [0.3 0.2 0.2 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([15 16.4]);
%                         ylim([0.03 0.05]);
%                     case 25
%                         xlim([22 28]);
%                         ylim([0 1]);
%                         
%                         dim = [0.225 0.2 0.35 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([23 24.6]);
%                         ylim([0.03 0.05]);
%                     case 50
%                         xlim([28 33]);
%                         ylim([0 1]);
%                         
%                         dim = [0.275 0.2 0.275 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([29 30.4]);
%                         ylim([0.03 0.05]);
%                     case 75
%                         xlim([30 35]);
%                         ylim([0 1]);
%                         
%                         dim = [0.425 0.2 0.2 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([31.8 33.2]);
%                         ylim([0.03 0.05]);
%                     case 100
%                         xlim([31.5 36.3]);
%                         ylim([0 1]);
%                         
%                         dim = [0.425 0.2 0.225 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([33.3 34.7]);
%                         ylim([0.03 0.05]);
%                     case 150
%                         xlim([32 37.7]);
%                         ylim([0 1]);
%                         
%                         dim = [0.475 0.2 0.275 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([34.6 36.4]);
%                         ylim([0.03 0.05]);
%                 end
%             case 200
%                 switch K
%                     case 10
%                         xlim([22 27.5]);
%                         ylim([0 1]);
%                         
%                         dim = [0.275 0.2 0.35 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([23 25.5]);
%                         ylim([0.015 0.05]);
%                     case 25
%                         xlim([35 44]);
%                         ylim([0 1]);
%                         
%                         dim = [0.225 0.2 0.2 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .375 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([36 38]);
%                         ylim([0.015 0.05]);
%                     case 50
%                         xlim([47 55.5]);
%                         ylim([0 1]);
%                         
%                         dim = [0.225 0.175 0.225 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([48 51.5]);
%                         ylim([0.02 0.05]);
%                     case 75
%                         xlim([55 61.8]);
%                         ylim([0 1]);
%                         
%                         dim = [0.15 0.175 0.225 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([55 57.5]);
%                         ylim([0.02 0.05]);
%                     case 100
%                         xlim([59 65.5]);
%                         ylim([0 1]);
%                         
%                         dim = [0.25 0.175 0.225 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .35 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([60 61.8]);
%                         ylim([0.02 0.05]);
%                     case 150
%                         xlim([62 70]);
%                         ylim([0 1]);
%                         
%                         dim = [0.3 0.175 0.3 0.035];
%                         annotation('ellipse',dim,'linewidth',linewidth);
%                         
%                         axes('position',[.25 .375 .3 .3]);
%                         box on;
%                         
%                         plot(edg_sum_se{2,1},[cdf_sum_se{2,1} 1],'-','color',colours(1,:),'linewidth',linewidth);
%                         hold on;
%                         plot(edg_sum_se{3,1},[cdf_sum_se{3,1} 1],'-','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,1},[cdf_sum_se{4,1} 1],'-','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,2},[cdf_sum_se{2,2} 1],'--','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,2},[cdf_sum_se{3,2} 1],'--','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,2},[cdf_sum_se{4,2} 1],'--','color',colours(3,:),'linewidth',linewidth);
%                         plot(edg_sum_se{2,3},[cdf_sum_se{2,3} 1],':','color',colours(1,:),'linewidth',linewidth);
%                         plot(edg_sum_se{3,3},[cdf_sum_se{3,3} 1],':','color',colours(2,:),'linewidth',linewidth);
%                         plot(edg_sum_se{4,3},[cdf_sum_se{4,3} 1],':','color',colours(3,:),'linewidth',linewidth);
%                         
%                         set(gca,'fontname',fontname,'fontsize',fontsize);
%                         
%                         xlim([64 67]);
%                         ylim([0.02 0.05]);
%                     case 200
%                     case 250
%                 end
%         end
%     otherwise
% end
<<<<<<< HEAD

if savefig == 1
    if plotse == 1
        saveas(gcf,[root_save 'cdf_sum_se_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
        saveas(gcf,[root_save 'cdf_sum_se_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
        saveas(gcf,[root_save 'cdf_sum_se_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
    else
        saveas(gcf,[root_save 'cdf_throughput_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
        saveas(gcf,[root_save 'cdf_throughput_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
        saveas(gcf,[root_save 'cdf_throughput_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_L_' sprintf(zero_pad_2,L) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
    end
end
=======
% 
% if (savefig == 1)
%     saveas(gcf,[root_save 'cdf_sum_se_star_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'fig');
%     saveas(gcf,[root_save 'cdf_sum_se_star_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'png');
%     saveas(gcf,[root_save 'cdf_sum_se_star_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(N_MC*MC)],'epsc2');
% end
>>>>>>> f3cbfc255096b547f744ae3295841f92e75e42f0
