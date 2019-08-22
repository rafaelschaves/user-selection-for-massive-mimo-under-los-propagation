clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = [64 128 256];                                                          % Number of antennas at base station
r_k = 1.25;
r_l = 0.25;

% K = 320;                                                                 % Number of mobile users
% L = 80;                                                                  % Number of selected users

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

se_cell  = cell(M_SIZ,2);
psi_cell = cell(M_SIZ,2);

for m = 1:M_SIZ
    K = r_k*M(m);
    L = r_l*K;
    
    se     = zeros(K,MC,N_CHN,2);                                        % Rate using all K users
    se_sel = zeros(L,MC,N_ALG,N_CHN,2);                                  % Rate using only L users
    
    psi_     = zeros(K,MC,N_CHN);
    psi_sel_ = zeros(L,MC,N_ALG,N_CHN);
    
    for chn_idx = 1:N_CHN
        load([root_load_ul 'spectral_efficiency_mf_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_ul) '_dB_MC_' num2str(MC) '.mat']);
        load([root_load_dl 'spectral_efficiency_mf_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_dl) '_dB_MC_' num2str(MC) '.mat']);
        
        se(:,:,chn_idx,1) = se_u;
        se(:,:,chn_idx,2) = se_d;
                
        psi_(:,:,chn_idx) = psi;
        
        se_sel(:,:,:,chn_idx,1) = se_u_sel;
        se_sel(:,:,:,chn_idx,2) = se_d_sel;
        
        psi_sel_(:,:,:,chn_idx) = psi_sel;
    end
    
    se_cell{m,1} = se;
    se_cell{m,2} = se_sel;
    
    psi_cell{m,1} = psi_;
    psi_cell{m,2} = psi_sel_;
end

% Post processing - Calculating the CDF

N_BIN = 50;
N_PDF = 100;

cdf_sum_se = cell(N_ALG+1,M_SIZ,N_CHN,2);
edg_sum_se = cell(N_ALG+1,M_SIZ,N_CHN,2);

cdf_se_user = cell(N_ALG+1,M_SIZ,N_CHN,2);
edg_se_user = cell(N_ALG+1,M_SIZ,N_CHN,2);

pdf_psi_y = zeros(N_ALG+1,N_PDF,M_SIZ,N_CHN);
pdf_psi_x = zeros(N_ALG+1,N_PDF,M_SIZ,N_CHN);

for m = 1:M_SIZ
    for chn_idx = 1:N_CHN
        % [pdf_psi{1,m,chn_idx},edg_psi{1,m,chn_idx}] = histcounts(psi_cell{m,1}(1,:,chn_idx),'binwidth',binwid,'normalization','pdf');
        
        [pdf_psi_y(1,:,m,chn_idx), pdf_psi_x(1,:,m,chn_idx)] = ksdensity(psi_cell{m,1}(1,:,chn_idx),'support','positive');
        
        [cdf_sum_se{1,m,chn_idx,1},edg_sum_se{1,m,chn_idx,1}] = histcounts(sum(se_cell{m,1}(:,:,chn_idx,1)),N_BIN,'normalization','cdf');
        [cdf_sum_se{1,m,chn_idx,2},edg_sum_se{1,m,chn_idx,2}] = histcounts(sum(se_cell{m,1}(:,:,chn_idx,2)),N_BIN,'normalization','cdf');
        
        [cdf_se_user{1,m,chn_idx,1},edg_se_user{1,m,chn_idx,1}] = histcounts(se_cell{m,1}(:,:,chn_idx,1),N_BIN,'normalization','cdf');
        [cdf_se_user{1,m,chn_idx,2},edg_se_user{1,m,chn_idx,2}] = histcounts(se_cell{m,1}(:,:,chn_idx,2),N_BIN,'normalization','cdf');
                
        cdf_sum_se{1,m,chn_idx,1} = [cdf_sum_se{1,m,chn_idx,1} 1];
        cdf_sum_se{1,m,chn_idx,2} = [cdf_sum_se{1,m,chn_idx,2} 1];
        
        cdf_se_user{1,m,chn_idx,1} = [cdf_se_user{1,m,chn_idx,1} 1];
        cdf_se_user{1,m,chn_idx,2} = [cdf_se_user{1,m,chn_idx,2} 1];
        
        for alg_idx = 1:N_ALG
            % [pdf_psi{alg_idx+1,m,chn_idx},edg_psi{alg_idx+1,m,chn_idx}] = histcounts(psi_cell{m,2}(1,:,alg_idx,chn_idx),N_BIN,'normalization','pdf');
            
            [pdf_psi_y(alg_idx+1,:,m,chn_idx), pdf_psi_x(alg_idx+1,:,m,chn_idx)] = ksdensity(psi_cell{m,2}(1,:,alg_idx,chn_idx),'support','positive');

            [cdf_sum_se{alg_idx+1,m,chn_idx,1},edg_sum_se{alg_idx+1,m,chn_idx,1}] = histcounts(sum(se_cell{m,2}(:,:,alg_idx,chn_idx,1)),N_BIN,'normalization','cdf');
            [cdf_sum_se{alg_idx+1,m,chn_idx,2},edg_sum_se{alg_idx+1,m,chn_idx,2}] = histcounts(sum(se_cell{m,2}(:,:,alg_idx,chn_idx,2)),N_BIN,'normalization','cdf');
            
            [cdf_se_user{alg_idx+1,m,chn_idx,1},edg_se_user{alg_idx+1,m,chn_idx,1}] = histcounts(se_cell{m,2}(:,:,alg_idx,chn_idx,1),N_BIN,'normalization','cdf');
            [cdf_se_user{alg_idx+1,m,chn_idx,2},edg_se_user{alg_idx+1,m,chn_idx,2}] = histcounts(se_cell{m,2}(:,:,alg_idx,chn_idx,2),N_BIN,'normalization','cdf');
            
            cdf_sum_se{alg_idx+1,m,chn_idx,1} = [cdf_sum_se{alg_idx+1,m,chn_idx,1} 1];
            cdf_sum_se{alg_idx+1,m,chn_idx,2} = [cdf_sum_se{alg_idx+1,m,chn_idx,2} 1];
            
            cdf_se_user{alg_idx+1,m,chn_idx,1} = [cdf_se_user{alg_idx+1,m,chn_idx,1} 1];
            cdf_se_user{alg_idx+1,m,chn_idx,2} = [cdf_se_user{alg_idx+1,m,chn_idx,2} 1];
        end
    end
end

% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 30;

marker = {'o','s','^'};

linestyle = {'-','--',':'};

savefig = 1;

% NS - No selection
% RS - Random selection
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

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
    for m = 1:M_SIZ
        K = r_k*M(m);
        L = r_l*K;
        
        if chn_idx == 1
            figure;
            
            set(gcf,'position',[0 0 800 600]);
            
            plot(pdf_psi_x(1,:,m,chn_idx),pdf_psi_y(1,:,m,chn_idx),'-','color',colours(1,:),'linewidth',linewidth);
            hold on;
            plot(pdf_psi_x(2,:,m,chn_idx),pdf_psi_y(2,:,m,chn_idx),'-','color',colours(2,:),'linewidth',linewidth);
            plot(pdf_psi_x(3,:,m,chn_idx),pdf_psi_y(3,:,m,chn_idx),'-','color',colours(3,:),'linewidth',linewidth);
            plot(pdf_psi_x(4,:,m,chn_idx),pdf_psi_y(4,:,m,chn_idx),'-','color',colours(4,:),'linewidth',linewidth);
            plot(pdf_psi_x(5,:,m,chn_idx),pdf_psi_y(5,:,m,chn_idx),'-','color',colours(5,:),'linewidth',linewidth);
        else
            if m == 1
                figure;
                
                set(gcf,'position',[0 0 800 600]);
                
                plot(pdf_psi_x(1,:,m,chn_idx),pdf_psi_y(1,:,m,chn_idx),linestyle{m},'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(pdf_psi_x(2,:,m,chn_idx),pdf_psi_y(2,:,m,chn_idx),linestyle{m},'color',colours(2,:),'linewidth',linewidth);
                plot(pdf_psi_x(3,:,m,chn_idx),pdf_psi_y(3,:,m,chn_idx),linestyle{m},'color',colours(3,:),'linewidth',linewidth);
                plot(pdf_psi_x(4,:,m,chn_idx),pdf_psi_y(4,:,m,chn_idx),linestyle{m},'color',colours(4,:),'linewidth',linewidth);
                plot(pdf_psi_x(5,:,m,chn_idx),pdf_psi_y(5,:,m,chn_idx),linestyle{m},'color',colours(5,:),'linewidth',linewidth);
            end
            
            plot(pdf_psi_x(1,:,m,chn_idx),pdf_psi_y(1,:,m,chn_idx),linestyle{m},'color',colours(1,:),'linewidth',linewidth);
            plot(pdf_psi_x(2,:,m,chn_idx),pdf_psi_y(2,:,m,chn_idx),linestyle{m},'color',colours(2,:),'linewidth',linewidth);
            plot(pdf_psi_x(3,:,m,chn_idx),pdf_psi_y(3,:,m,chn_idx),linestyle{m},'color',colours(3,:),'linewidth',linewidth);
            plot(pdf_psi_x(4,:,m,chn_idx),pdf_psi_y(4,:,m,chn_idx),linestyle{m},'color',colours(4,:),'linewidth',linewidth);
            plot(pdf_psi_x(5,:,m,chn_idx),pdf_psi_y(5,:,m,chn_idx),linestyle{m},'color',colours(5,:),'linewidth',linewidth);
        end
        
        xlabel('Inter-channel interference','fontname',fontname,'fontsize',fontsize);
        ylabel('Probability distribution','fontname',fontname,'fontsize',fontsize);
        
        legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_2);
        legend box off;
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
    end
end

% for m = 1:M_SIZ
%     K = r_k*M(m);
%     L = r_l*K;
%     
%     for chn_idx = 1:N_CHN
%         figure;
%         
%         set(gcf,'position',[0 0 800 600]);
%         
%         plot(pdf_psi_x(1,:,m,chn_idx),pdf_psi_y(1,:,m,chn_idx),'-','color',colours(1,:),'linewidth',linewidth);
%         hold on;
%         plot(pdf_psi_x(2,:,m,chn_idx),pdf_psi_y(2,:,m,chn_idx),'-','color',colours(2,:),'linewidth',linewidth);
%         plot(pdf_psi_x(3,:,m,chn_idx),pdf_psi_y(3,:,m,chn_idx),'-','color',colours(3,:),'linewidth',linewidth);
%         plot(pdf_psi_x(4,:,m,chn_idx),pdf_psi_y(4,:,m,chn_idx),'-','color',colours(4,:),'linewidth',linewidth);
%         plot(pdf_psi_x(5,:,m,chn_idx),pdf_psi_y(5,:,m,chn_idx),'-','color',colours(5,:),'linewidth',linewidth);
%         
%         xlabel('Inter-channel interference','fontname',fontname,'fontsize',fontsize);
%         ylabel('Probability distribution','fontname',fontname,'fontsize',fontsize);
%         
%         legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_2);
%         legend box off;
%         
%         set(gca,'fontname',fontname,'fontsize',fontsize);
%         
% %         if chn_idx == 1
% %             dim = [0.185 0.18 0.075 0.7];
% %         
% %             annotation('rectangle',dim,'linewidth',linewidth);
% %        
% %             axes('position',[.6 .275 .25 .25]);
% %             box on;
% %             
% %             histogram(psi_sel_(1,:,3,chn_idx),'normalization','pdf','facecolor',colours(4,:));
% %             hold on;
% %             histogram(psi_sel_(1,:,4,chn_idx),'normalization','pdf','facecolor',colours(5,:));
% %                         
% %             set(gca,'fontname',fontname,'fontsize',fontsize);
% %         end
%     end
%     
% %     for chn_idx = 1:N_CHN
% %         figure;
% %         
% %         set(gcf,'position',[0 0 800 600]);
% %         
% %         plot(edg_sum_se{1,m,chn_idx,1},cdf_sum_se{1,m,chn_idx,1},'-','color',colours(1,:),'linewidth',linewidth);
% %         hold on;
% %         plot(edg_sum_se{2,m,chn_idx,1},cdf_sum_se{2,m,chn_idx,1},'-','color',colours(2,:),'linewidth',linewidth);
% %         plot(edg_sum_se{3,m,chn_idx,1},cdf_sum_se{3,m,chn_idx,1},'-','color',colours(3,:),'linewidth',linewidth);
% %         plot(edg_sum_se{4,m,chn_idx,1},cdf_sum_se{4,m,chn_idx,1},'-','color',colours(4,:),'linewidth',linewidth);
% %         plot(edg_sum_se{5,m,chn_idx,1},cdf_sum_se{5,m,chn_idx,1},'-','color',colours(5,:),'linewidth',linewidth);
% %         plot(edg_sum_se{1,m,chn_idx,2},cdf_sum_se{1,m,chn_idx,2},'--','color',colours(1,:),'linewidth',linewidth);
% %         plot(edg_sum_se{2,m,chn_idx,2},cdf_sum_se{2,m,chn_idx,2},'--','color',colours(2,:),'linewidth',linewidth);
% %         plot(edg_sum_se{3,m,chn_idx,2},cdf_sum_se{3,m,chn_idx,2},'--','color',colours(3,:),'linewidth',linewidth);
% %         plot(edg_sum_se{4,m,chn_idx,2},cdf_sum_se{4,m,chn_idx,2},'--','color',colours(4,:),'linewidth',linewidth);
% %         plot(edg_sum_se{5,m,chn_idx,2},cdf_sum_se{5,m,chn_idx,2},'--','color',colours(5,:),'linewidth',linewidth);
% %         
% %         xlabel('Sum-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% %         ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
% %         
% %         legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_1);
% %         legend box off;
% %         
% %         set(gca,'fontname',fontname,'fontsize',fontsize);
% %         
% %         if m == 1 && chn_idx == 1
% %             xlim([0 50]);
% %         elseif m == 1 && chn_idx == 2
% %             xlim([10 35]);
% %         elseif m == 2 && chn_idx == 1
% %             xlim([20 110]);
% %         elseif m == 2 && chn_idx == 2
% %             xlim([25 70]);
% %         elseif m == 3 && chn_idx == 1
% %             xlim([40 225]);
% %         elseif m == 3 && chn_idx == 2
% %             xlim([50 140]);
% %         end
% %         
% %         ylim([0 1]);
% %         
% %         if (savefig == 1)
% %             saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'fig');
% %             saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'png');
% %             saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'epsc2');
% %         end
% %         
% %         figure;
% %         
% %         set(gcf,'position',[0 0 800 600]);
% %         
% %         plot(edg_se_user{1,m,chn_idx,1},cdf_se_user{1,m,chn_idx,1},'-','color',colours(1,:),'linewidth',linewidth);
% %         hold on;
% %         plot(edg_se_user{2,m,chn_idx,1},cdf_se_user{2,m,chn_idx,1},'-','color',colours(2,:),'linewidth',linewidth);
% %         plot(edg_se_user{3,m,chn_idx,1},cdf_se_user{3,m,chn_idx,1},'-','color',colours(3,:),'linewidth',linewidth);
% %         plot(edg_se_user{4,m,chn_idx,1},cdf_se_user{4,m,chn_idx,1},'-','color',colours(4,:),'linewidth',linewidth);
% %         plot(edg_se_user{5,m,chn_idx,1},cdf_se_user{5,m,chn_idx,1},'-','color',colours(5,:),'linewidth',linewidth);
% %         plot(edg_se_user{1,m,chn_idx,2},cdf_se_user{1,m,chn_idx,2},'--','color',colours(1,:),'linewidth',linewidth);
%         plot(edg_se_user{2,m,chn_idx,2},cdf_se_user{2,m,chn_idx,2},'--','color',colours(2,:),'linewidth',linewidth);
%         plot(edg_se_user{3,m,chn_idx,2},cdf_se_user{3,m,chn_idx,2},'--','color',colours(3,:),'linewidth',linewidth);
%         plot(edg_se_user{4,m,chn_idx,2},cdf_se_user{4,m,chn_idx,2},'--','color',colours(4,:),'linewidth',linewidth);
%         plot(edg_se_user{5,m,chn_idx,2},cdf_se_user{5,m,chn_idx,2},'--','color',colours(5,:),'linewidth',linewidth);
%         
%         xlabel('Spectral efficiency per user (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
%         ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
%         
%         legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_2);
%         legend box off;
%         
%         set(gca,'fontname',fontname,'fontsize',fontsize);
%         
%         if chn_idx == 1
%             xlim([0 4]);
%         else
%             xlim([0 2]);
%         end
%         
%         ylim([0 1]);
%         
%         if (savefig == 1)
%             saveas(gcf,[root_save 'cdf_se_per_user_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'fig');
%             saveas(gcf,[root_save 'cdf_se_per_user_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'png');
%             saveas(gcf,[root_save 'cdf_se_per_user_' chn_type{chn_idx} '_M_' num2str(M(m)) '_K_' num2str(K) '_L_' num2str(L)],'epsc2');
%         end
%     end
%end

% for chn_idx = 1:N_CHN
%     figure;
%     
%     set(gcf,'position',[0 0 800 600]);
%     
%     for m = 1:M_SIZ
%         plot(edg_sum_se{1,m,chn_idx,1},cdf_sum_se{1,m,chn_idx,1},['-' marker{m}],'color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%         hold on;
%         plot(edg_sum_se{2,m,chn_idx,1},cdf_sum_se{2,m,chn_idx,1},['-' marker{m}],'color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_sum_se{3,m,chn_idx,1},cdf_sum_se{3,m,chn_idx,1},['-' marker{m}],'color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_sum_se{4,m,chn_idx,1},cdf_sum_se{4,m,chn_idx,1},['-' marker{m}],'color',colours(4,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_sum_se{5,m,chn_idx,1},cdf_sum_se{5,m,chn_idx,1},['-' marker{m}],'color',colours(5,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_sum_se{1,m,chn_idx,2},cdf_sum_se{1,m,chn_idx,2},['--' marker{m}],'color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_sum_se{2,m,chn_idx,2},cdf_sum_se{2,m,chn_idx,2},['--' marker{m}],'color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_sum_se{3,m,chn_idx,2},cdf_sum_se{3,m,chn_idx,2},['--' marker{m}],'color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_sum_se{4,m,chn_idx,2},cdf_sum_se{4,m,chn_idx,2},['--' marker{m}],'color',colours(4,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_sum_se{5,m,chn_idx,2},cdf_sum_se{5,m,chn_idx,2},['--' marker{m}],'color',colours(5,:),'linewidth',linewidth,'markersize',markersize);
%     end
%     
%     xlabel('Sum-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
%     ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
%     
%     legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
%     set(gca,'fontname',fontname,'fontsize',fontsize);
%     
%     ylim([0 1]);
%     
%     if (savefig == 1)
%         saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_r_k_' num2str(r_k) '_r_l_' num2str(r_l)],'fig');
%         saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_r_k_' num2str(r_k) '_r_l_' num2str(r_l)],'png');
%         saveas(gcf,[root_save 'cdf_sum_se_' chn_type{chn_idx} '_r_k_' num2str(r_k) '_r_l_' num2str(r_l)],'epsc2');
%     end
%     
%     figure;
%     
%     set(gcf,'position',[0 0 800 600]);
%     
%     for m = 1:M_SIZ
%         plot(edg_se_user{1,m,chn_idx,1},cdf_se_user{1,m,chn_idx,1},['-' marker{m}],'color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%         hold on;
%         plot(edg_se_user{2,m,chn_idx,1},cdf_se_user{2,m,chn_idx,1},['-' marker{m}],'color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_se_user{3,m,chn_idx,1},cdf_se_user{3,m,chn_idx,1},['-' marker{m}],'color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_se_user{4,m,chn_idx,1},cdf_se_user{4,m,chn_idx,1},['-' marker{m}],'color',colours(4,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_se_user{5,m,chn_idx,1},cdf_se_user{5,m,chn_idx,1},['-' marker{m}],'color',colours(5,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_se_user{1,m,chn_idx,2},cdf_se_user{1,m,chn_idx,2},['--' marker{m}],'color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_se_user{2,m,chn_idx,2},cdf_se_user{2,m,chn_idx,2},['--' marker{m}],'color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_se_user{3,m,chn_idx,2},cdf_se_user{3,m,chn_idx,2},['--' marker{m}],'color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_se_user{4,m,chn_idx,2},cdf_se_user{4,m,chn_idx,2},['--' marker{m}],'color',colours(4,:),'linewidth',linewidth,'markersize',markersize);
%         plot(edg_se_user{5,m,chn_idx,2},cdf_se_user{5,m,chn_idx,2},['--' marker{m}],'color',colours(5,:),'linewidth',linewidth,'markersize',markersize);
%     end
%     
%     xlabel('Sum-spectral efficiency per user (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
%     ylabel('Outage probability','fontname',fontname,'fontsize',fontsize);
%     
%     legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
%     set(gca,'fontname',fontname,'fontsize',fontsize);
%     
%     ylim([0 1]);
%     
%     if (savefig == 1)
%         saveas(gcf,[root_save 'cdf_se_per_user' chn_type{chn_idx} '_r_k_' num2str(r_k) '_r_l_' num2str(r_l)],'fig');
%         saveas(gcf,[root_save 'cdf_se_per_user' chn_type{chn_idx} '_r_k_' num2str(r_k) '_r_l_' num2str(r_l)],'png');
%         saveas(gcf,[root_save 'cdf_se_per_user' chn_type{chn_idx} '_r_k_' num2str(r_k) '_r_l_' num2str(r_l)],'epsc2');
%     end
% end