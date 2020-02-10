clear;
close all;
clc;

% Macros

MC = 10000;                                                                 % Size of the monte-carlo ensemble

M = [64 128 256 512];                                                      % Number of antennas at base station
K = [8 16 32 64];                                                          % Number of mobile users

bs_power = [1 10];
radius   = [100 500 1000 2000];

M_SIZ = length(M);                                                         % Size of the antennas set
K_SIZ = length(K);
P_SIZ = length(bs_power);
R_SIZ = length(radius);

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Power Allocation/Downlink/';
root_save_eps = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Power Allocation/eps/';
root_save_fig = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Power Allocation/fig/';
root_save_png = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Power Allocation/png/';

% chn_type = {'ur_los','rayleigh'};

% Loading data

error_prob_all = zeros(M_SIZ,K_SIZ,R_SIZ,P_SIZ);

for p = 1:P_SIZ
    for r = 1:R_SIZ
        for m = 1:M_SIZ
            for k = 1:K_SIZ
                load([root_load 'results_error_prob_ur_los_M_' num2str(M(m)) ...
                      '_K_' num2str(K(k)) '_cell_radius_' num2str(radius(r)) ...
                      '_m_BS_power_' num2str(bs_power(p)) '_W_MC_' num2str(MC) '.mat']);
                error_prob_all(m,k,r,p) = sum(abs(gamma(:,1) - gamma(:,2)) > 1e-6)/MC;
            end
        end
    end
end

% Ploting Figures

linewidth  = 3;
markersize = 15;
fontname   = 'Times New Roman';
fontsize   = 30;

marker = {'o','s','^'};

linestyle = {'-','--',':'};

savefig = 1;

legend_M = {'$M = 64$','$M = 128$','$M = 256$','$M = 512$'};
legend_K = {'$K = 8$' ,'$K = 16$' ,'$K = 32$' ,'$K = 64$'};

location_1 = 'northwest';
location_2 = 'northeast';
location_3 = 'southeast';
location_4 = 'southwest';

colours = [0.0000 0.4470 0.7410;
           0.8500 0.3250 0.0980;
           0.9290 0.6940 0.1250;
           0.4940 0.1840 0.5560;
           0.4660 0.6740 0.1880;
           0.3010 0.7450 0.9330;
           0.6350 0.0780 0.1840];

for r = 1:R_SIZ
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot(K,error_prob_all(1,:,r,1),'-ok','linewidth',linewidth,'markersize',markersize);
    hold on;
    plot(K,error_prob_all(2,:,r,1),'-sk','linewidth',linewidth,'markersize',markersize);
    plot(K,error_prob_all(3,:,r,1),'-^k','linewidth',linewidth,'markersize',markersize);
    plot(K,error_prob_all(4,:,r,1),'-vk','linewidth',linewidth,'markersize',markersize);
    plot(K,error_prob_all(1,:,r,1),'-o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
    plot(K,error_prob_all(2,:,r,1),'-s','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
    plot(K,error_prob_all(3,:,r,1),'-^','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
    plot(K,error_prob_all(4,:,r,1),'-v','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
    plot(K,error_prob_all(1,:,r,2),'--o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
    plot(K,error_prob_all(2,:,r,2),'--s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
    plot(K,error_prob_all(3,:,r,2),'--^','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
    plot(K,error_prob_all(4,:,r,2),'--v','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
    
    xlabel('Number of users','fontname',fontname,'fontsize',fontsize);
    ylabel('Probability of failure','fontname',fontname,'fontsize',fontsize);
    
    if r == 1 
        legend(legend_M,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_3);
        legend box off;
    end
    
    % title(['$R = ' num2str(radius(r)) '$ m'],'fontname',fontname,'fontsize',fontsize,'interpreter','latex');
    
    % xticks([8 16 24 32 40 48 56 64])
    
    xlim([K(1) K(end)]);
    ylim([0 1]);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    if (savefig == 1)
        saveas(gcf,[root_save_fig 'prob_failure_R_' num2str(radius(r))],'fig');
        saveas(gcf,[root_save_png 'prob_failure_R_' num2str(radius(r))],'png');
        saveas(gcf,[root_save_eps 'prob_failure_R_' num2str(radius(r))],'epsc2');
    end
end

% for r = 1:R_SIZ
%     figure;
%     
%     set(gcf,'position',[0 0 800 600]);
%     
%     plot(M,error_prob_all(:,1,r,1),'-ok','linewidth',linewidth,'markersize',markersize);
%     hold on;
%     plot(M,error_prob_all(:,2,r,1),'-sk','linewidth',linewidth,'markersize',markersize);
%     plot(M,error_prob_all(:,3,r,1),'-^k','linewidth',linewidth,'markersize',markersize);
%     plot(M,error_prob_all(:,4,r,1),'-vk','linewidth',linewidth,'markersize',markersize);
%     plot(M,error_prob_all(:,1,r,1),'-o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%     plot(M,error_prob_all(:,2,r,1),'-s','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%     plot(M,error_prob_all(:,3,r,1),'-^','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%     plot(M,error_prob_all(:,4,r,1),'-v','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%     plot(M,error_prob_all(:,1,r,2),'--o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%     plot(M,error_prob_all(:,2,r,2),'--s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%     plot(M,error_prob_all(:,3,r,2),'--^','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%     plot(M,error_prob_all(:,4,r,2),'--v','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%     
%     xlabel('Number of antennas','fontname',fontname,'fontsize',fontsize);
%     ylabel('Probability of failure','fontname',fontname,'fontsize',fontsize);
%     
%     legend(legend_K,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_3);
%     legend box off;
%     
%     % title(['$R = ' num2str(radius(r)) '$ m'],'fontname',fontname,'fontsize',fontsize,'interpreter','latex');
%     
%     % xticks([8 16 24 32 40 48 56 64])
%     
%     xlim([M(1) M(end)]);
%     ylim([0 1]);
%     
%     set(gca,'fontname',fontname,'fontsize',fontsize);
%     
%     if (savefig == 1)
%         saveas(gcf,[root_save 'number_iterations_ur_los'],'fig');
%         saveas(gcf,[root_save 'number_iterations_ur_los'],'png');
%         saveas(gcf,[root_save 'number_iterations_ur_los'],'epsc2');
%     end
% end