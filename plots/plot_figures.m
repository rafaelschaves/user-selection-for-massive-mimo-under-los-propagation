clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = 128;                                                                    % Number of antennas at base station
K = 160;                                                                    % Number of mobile users
L = 80;                                                                     % Number of selected users

snr_ul = -7;
snr_dl = 10;

M_SIZ = length(M);                                                         % Size of the antennas set
N_ALG = 4;                                                                 % Number of algorithms for perform user scheduling
N_CHN = 2;                                                                 % Number of channel models simulated

% Roots

root_load_ul = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Uplink/';
root_load_dl = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Selection/Downlink/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Selection/';

chn_type = {'ur_los','rayleigh'};

% Loading data

se       = zeros(K,MC,M_SIZ,N_CHN,2);                                        % Rate using all K users
se_sel   = zeros(L,MC,N_ALG,M_SIZ,N_CHN,2);                                  % Rate using only L users
se_bound = zeros(MC,M_SIZ,N_CHN,2); 

for chn_idx = 1:N_CHN
    for m = 1:M_SIZ
        load([root_load_ul 'spectral_efficiency_mf_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_ul) '_dB_MC_' num2str(MC) '.mat']);
        load([root_load_dl 'spectral_efficiency_mf_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_dl) '_dB_MC_' num2str(MC) '.mat']);
        
        se(:,:,m,chn_idx,1) = se_u;
        se(:,:,m,chn_idx,2) = se_d;
        
        se_sel(:,:,:,m,chn_idx,1) = se_u_sel;
        se_sel(:,:,:,m,chn_idx,2) = se_d_sel;
        
        se_bound(:,m,chn_idx,1) = 0.5*sum(log2(1 + (10^(snr_ul/10)*M(m))./(1 + 10^(snr_ul/10)*(K-1)^(2)*M(m)*psi.^2)));
        se_bound(:,m,chn_idx,2) = 0.5*sum(log2(1 + (10^(snr_dl/10)*M(m)/K)./(1 + 10^(snr_dl/10)*(K-1)^(2)*M(m)*psi.^2/K)));
    end
end

% Post processing - Calculating the CDF

bin_width = 0.0005;

cdf_sum_se = cell(N_ALG+1,M_SIZ,N_CHN,2);
edg_sum_se = cell(N_ALG+1,M_SIZ,N_CHN,2); 

cdf_se_user = cell(N_ALG+1,M_SIZ,N_CHN,2);
edg_se_user = cell(N_ALG+1,M_SIZ,N_CHN,2);

cdf_se_bound = cell(M_SIZ,N_CHN,2);
edg_se_bound = cell(M_SIZ,N_CHN,2);

for chn_idx = 1:N_CHN
    for m = 1:M_SIZ
        [cdf_sum_se{1,m,chn_idx,1},edg_sum_se{1,m,chn_idx,1}] = histcounts(sum(se(:,:,m,chn_idx,1)),'binwidth',bin_width,'normalization','cdf');
        [cdf_sum_se{1,m,chn_idx,2},edg_sum_se{1,m,chn_idx,2}] = histcounts(sum(se(:,:,m,chn_idx,2)),'binwidth',bin_width,'normalization','cdf');
        
        [cdf_se_user{1,m,chn_idx,1},edg_se_user{1,m,chn_idx,1}] = histcounts(se(:,:,m,chn_idx,1),'binwidth',bin_width,'normalization','cdf');
        [cdf_se_user{1,m,chn_idx,2},edg_se_user{1,m,chn_idx,2}] = histcounts(se(:,:,m,chn_idx,2),'binwidth',bin_width,'normalization','cdf');
        
        [cdf_se_bound{m,chn_idx,1},edg_se_bound{m,chn_idx,1}] = histcounts(se_bound(:,m,chn_idx,1),'binwidth',bin_width,'normalization','cdf');
        [cdf_se_bound{m,chn_idx,2},edg_se_bound{m,chn_idx,2}] = histcounts(se_bound(:,m,chn_idx,2),'binwidth',bin_width,'normalization','cdf');
        
        cdf_sum_se{1,m,chn_idx,1} = [cdf_sum_se{1,m,chn_idx,1} 1];
        cdf_sum_se{1,m,chn_idx,2} = [cdf_sum_se{1,m,chn_idx,2} 1];
        
        edg_sum_se{1,m,chn_idx,1} = edg_sum_se{1,m,chn_idx,1} + bin_width/2;
        edg_sum_se{1,m,chn_idx,2} = edg_sum_se{1,m,chn_idx,2} + bin_width/2;
        
        cdf_se_user{1,m,chn_idx,1} = [cdf_se_user{1,m,chn_idx,1} 1];
        cdf_se_user{1,m,chn_idx,2} = [cdf_se_user{1,m,chn_idx,2} 1];
        
        edg_se_user{1,m,chn_idx,1} = edg_se_user{1,m,chn_idx,1} + bin_width/2;
        edg_se_user{1,m,chn_idx,2} = edg_se_user{1,m,chn_idx,2} + bin_width/2;

        cdf_se_bound{m,chn_idx,1} = [cdf_se_bound{m,chn_idx,1} 1];
        cdf_se_bound{m,chn_idx,2} = [cdf_se_bound{m,chn_idx,2} 1];
        
        edg_se_bound{m,chn_idx,1} = edg_se_bound{m,chn_idx,1} + bin_width/2;
        edg_se_bound{m,chn_idx,2} = edg_se_bound{m,chn_idx,2} + bin_width/2;
        
        for alg_idx = 1:N_ALG
            [cdf_sum_se{alg_idx+1,m,chn_idx,1},edg_sum_se{alg_idx+1,m,chn_idx,1}] = histcounts(sum(se_sel(:,:,alg_idx,m,chn_idx,1)),'binwidth',bin_width,'normalization','cdf');
            [cdf_sum_se{alg_idx+1,m,chn_idx,2},edg_sum_se{alg_idx+1,m,chn_idx,2}] = histcounts(sum(se_sel(:,:,alg_idx,m,chn_idx,2)),'binwidth',bin_width,'normalization','cdf');
            
            [cdf_se_user{alg_idx+1,m,chn_idx,1},edg_se_user{alg_idx+1,m,chn_idx,1}] = histcounts(se_sel(:,:,alg_idx,m,chn_idx,1),'binwidth',bin_width,'normalization','cdf');
            [cdf_se_user{alg_idx+1,m,chn_idx,2},edg_se_user{alg_idx+1,m,chn_idx,2}] = histcounts(se_sel(:,:,alg_idx,m,chn_idx,2),'binwidth',bin_width,'normalization','cdf');
            
            cdf_sum_se{alg_idx+1,m,chn_idx,1} = [cdf_sum_se{alg_idx+1,m,chn_idx,1} 1];
            cdf_sum_se{alg_idx+1,m,chn_idx,2} = [cdf_sum_se{alg_idx+1,m,chn_idx,2} 1];
                        
            edg_sum_se{alg_idx+1,m,chn_idx,1} = edg_sum_se{alg_idx+1,m,chn_idx,1} + bin_width/2;
            edg_sum_se{alg_idx+1,m,chn_idx,2} = edg_sum_se{alg_idx+1,m,chn_idx,2} + bin_width/2;
            
            cdf_se_user{alg_idx+1,m,chn_idx,1} = [cdf_se_user{alg_idx+1,m,chn_idx,1} 1];
            cdf_se_user{alg_idx+1,m,chn_idx,2} = [cdf_se_user{alg_idx+1,m,chn_idx,2} 1];
                        
            edg_se_user{alg_idx+1,m,chn_idx,1} = edg_se_user{alg_idx+1,m,chn_idx,1} + bin_width/2;
            edg_se_user{alg_idx+1,m,chn_idx,2} = edg_se_user{alg_idx+1,m,chn_idx,2} + bin_width/2;
        end
    end
end

% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 35;

savefig = 0;

% NS - No selection
% RS - Random selection
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_boun = {'Simulated','Bound'};
legend_algo = {'NS','RS','SOS','CBS','ICIBS'};

location_1 = 'northwest';
location_2 = 'northeast';

colours = [0.0000 0.4470 0.7410;
           0.8500 0.3250 0.0980;
           0.9290 0.6940 0.1250;
           0.4940 0.1840 0.5560;
           0.4660 0.6740 0.1880;
           0.3010 0.7450 0.9330;
           0.6350 0.0780 0.1840];

for chn_idx = 1:N_CHN
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot(edg_sum_se{1,1,chn_idx,1},cdf_sum_se{1,1,chn_idx,1},'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(edg_se_bound{1,chn_idx,1},cdf_se_bound{1,chn_idx,1},'o','color',colours(2,:),'linewidth',linewidth);
    %hold on;
    plot(edg_sum_se{1,1,chn_idx,2},cdf_sum_se{1,1,chn_idx,2},'--','color',colours(1,:),'linewidth',linewidth);
    plot(edg_se_bound{1,chn_idx,2},cdf_se_bound{1,chn_idx,2},'o','color',colours(3,:),'linewidth',linewidth);
    
    xlabel('Sum-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_boun,'fontname',fontname,'fontsize',fontsize,'location',location_1);
    legend box off;
        
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    ylim([0 1]);

    if (savefig == 1)
        saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'fig');
        saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'png');
        saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'epsc2');
    end
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot(edg_sum_se{1,1,chn_idx,1},cdf_sum_se{1,1,chn_idx,1},'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(edg_sum_se{2,1,chn_idx,1},cdf_sum_se{2,1,chn_idx,1},'-','color',colours(2,:),'linewidth',linewidth);
    plot(edg_sum_se{3,1,chn_idx,1},cdf_sum_se{3,1,chn_idx,1},'-','color',colours(3,:),'linewidth',linewidth);
    plot(edg_sum_se{4,1,chn_idx,1},cdf_sum_se{4,1,chn_idx,1},'-','color',colours(4,:),'linewidth',linewidth);
    plot(edg_sum_se{5,1,chn_idx,1},cdf_sum_se{5,1,chn_idx,1},'-','color',colours(5,:),'linewidth',linewidth);
    plot(edg_sum_se{1,1,chn_idx,2},cdf_sum_se{1,1,chn_idx,2},'--','color',colours(1,:),'linewidth',linewidth);
    plot(edg_sum_se{2,1,chn_idx,2},cdf_sum_se{2,1,chn_idx,2},'--','color',colours(2,:),'linewidth',linewidth);
    plot(edg_sum_se{3,1,chn_idx,2},cdf_sum_se{3,1,chn_idx,2},'--','color',colours(3,:),'linewidth',linewidth);
    plot(edg_sum_se{4,1,chn_idx,2},cdf_sum_se{4,1,chn_idx,2},'--','color',colours(4,:),'linewidth',linewidth);
    plot(edg_sum_se{5,1,chn_idx,2},cdf_sum_se{5,1,chn_idx,2},'--','color',colours(5,:),'linewidth',linewidth);
    
    xlabel('Sum-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_1);
    legend box off;
    
    if chn_idx == 1
        xticks(0:10:60);
    end
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    if chn_idx == 2
        xlim([10 inf])
    end
    
    ylim([0 1]);

    if (savefig == 1)
        saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'fig');
        saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'png');
        saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'epsc2');
    end
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot(edg_se_user{1,1,chn_idx,1},cdf_se_user{1,1,chn_idx,1},'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(edg_se_user{2,1,chn_idx,1},cdf_se_user{2,1,chn_idx,1},'-','color',colours(2,:),'linewidth',linewidth);
    plot(edg_se_user{3,1,chn_idx,1},cdf_se_user{3,1,chn_idx,1},'-','color',colours(3,:),'linewidth',linewidth);
    plot(edg_se_user{4,1,chn_idx,1},cdf_se_user{4,1,chn_idx,1},'-','color',colours(4,:),'linewidth',linewidth);
    plot(edg_se_user{5,1,chn_idx,1},cdf_se_user{5,1,chn_idx,1},'-','color',colours(5,:),'linewidth',linewidth);
    plot(edg_se_user{1,1,chn_idx,2},cdf_se_user{1,1,chn_idx,2},'--','color',colours(1,:),'linewidth',linewidth);
    plot(edg_se_user{2,1,chn_idx,2},cdf_se_user{2,1,chn_idx,2},'--','color',colours(2,:),'linewidth',linewidth);
    plot(edg_se_user{3,1,chn_idx,2},cdf_se_user{3,1,chn_idx,2},'--','color',colours(3,:),'linewidth',linewidth);
    plot(edg_se_user{4,1,chn_idx,2},cdf_se_user{4,1,chn_idx,2},'--','color',colours(4,:),'linewidth',linewidth);
    plot(edg_se_user{5,1,chn_idx,2},cdf_se_user{5,1,chn_idx,2},'--','color',colours(5,:),'linewidth',linewidth);
    
    xlabel('Spectral efficiency per user (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_2);
    legend box off;
    
    if chn_idx == 1
        xticks(0:0.5:4);
    end

    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    if chn_idx == 1
        xlim([0 4]);
    else
        xlim([0 2.5]);
    end
    ylim([0 1]);
    
    if (savefig == 1)
        saveas(gcf,[root_save 'cdf_se_per_user_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'fig');
        saveas(gcf,[root_save 'cdf_se_per_user_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'png');
        saveas(gcf,[root_save 'cdf_se_per_user_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'epsc2');
    end
end