clear;
close all;
clc;

% Macros

MC = 1000;                                                                 % Size of the monte-carlo ensemble

M = [64 128 256 512];                                                      % Number of antennas at base station
K = [8 16 32 64];                                                          % Number of mobile users

bs_power = [1 10];
radius   = [100 500 1000 2000];

M_SIZ = length(M);                                                         % Size of the antennas set
K_SIZ = length(K);
P_SIZ = length(bs_power);
R_SIZ = length(radius);
N_ALG = 3;

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Power Allocation/Downlink/';
root_save_eps = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Power Allocation/eps/';
root_save_fig = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Power Allocation/fig/';
root_save_png = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Power Allocation/png/';

% Loading data

n_it_all    = zeros(M_SIZ,K_SIZ,R_SIZ,P_SIZ,MC,N_ALG);
time_all    = zeros(M_SIZ,K_SIZ,R_SIZ,P_SIZ,MC,N_ALG);

for p = 1:P_SIZ
    for r = 1:R_SIZ
        for m = 1:M_SIZ
            for k = 1:K_SIZ
                load([root_load 'results_iterations_ur_los_M_' num2str(M(m)) ...
                      '_K_' num2str(K(k)) '_cell_radius_' num2str(radius(r)) ...
                      '_m_BS_power_' num2str(bs_power(p)) '_W_MC_' num2str(MC) '.mat']);
                n_it_all(m,k,r,p,:,:) = n_it;
                time_all(m,k,r,p,:,:) = time;
            end
        end
    end
end

% Post Processing

n_it_avg = reshape(mean(n_it_all,5),M_SIZ,K_SIZ,R_SIZ,P_SIZ,N_ALG);
time_avg = reshape(mean(time_all,5),M_SIZ,K_SIZ,R_SIZ,P_SIZ,N_ALG);

% Ploting Figures

linewidth  = 3;
markersize = 15;
fontname   = 'Times New Roman';
fontsize   = 30;

marker = {'o','s','^'};

linestyle = {'-','--',':'};

savefig = 1;

legend_M = {'$M = 64$','$M = 128$','$M = 256$','$M = 512$'};

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

for p = 1:P_SIZ
    for r = 1:R_SIZ
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        plot(K,reshape(n_it_avg(1,:,r,p,3),1,[]),'-ok','linewidth',linewidth,'markersize',markersize);
        hold on;
        plot(K,reshape(n_it_avg(2,:,r,p,3),1,[]),'-sk','linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(3,:,r,p,3),1,[]),'-^k','linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(4,:,r,p,3),1,[]),'-vk','linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(1,:,r,p,3),1,[]),'-o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(2,:,r,p,3),1,[]),'-s','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(3,:,r,p,3),1,[]),'-^','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(4,:,r,p,3),1,[]),'-v','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(1,:,r,p,1),1,[]),'--o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(2,:,r,p,1),1,[]),'--s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(3,:,r,p,1),1,[]),'--^','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(4,:,r,p,1),1,[]),'--v','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(1,:,r,p,2),1,[]),':o','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(2,:,r,p,2),1,[]),':s','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(3,:,r,p,2),1,[]),':^','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
        plot(K,reshape(n_it_avg(4,:,r,p,2),1,[]),':v','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
        
        xlabel('Number of users','fontname',fontname,'fontsize',fontsize);
        ylabel('Number of iterations','fontname',fontname,'fontsize',fontsize);
                
        xlim([K(1) K(end)]);
        
        if p == 1
            if r == 1
                ylim([0 45]);
            elseif r == 4
                legend(legend_M,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_3);
                legend box off;
                
                ylim([0 35]);
            else
                ylim([0 40]);
            end
        end
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        if (savefig == 1)
            saveas(gcf,[root_save_fig 'number_iterations_ur_los_R_' num2str(radius(r)) '_P_' num2str(bs_power(p))],'fig');
            saveas(gcf,[root_save_png 'number_iterations_ur_los_R_' num2str(radius(r)) '_P_' num2str(bs_power(p))],'png');
            saveas(gcf,[root_save_eps 'number_iterations_ur_los_R_' num2str(radius(r)) '_P_' num2str(bs_power(p))],'epsc2');
        end
    end
end

% for p = 1:P_SIZ
%     for r = 1:R_SIZ
%         figure;
%         
%         set(gcf,'position',[0 0 800 600]);
%         
%         plot(K,1e3*reshape(time_avg(1,:,r,p,1),1,[]),'-ok','linewidth',linewidth,'markersize',markersize);
%         hold on;
%         plot(K,1e3*reshape(time_avg(2,:,r,p,1),1,[]),'-sk','linewidth',linewidth,'markersize',markersize);
%         plot(K,1e3*reshape(time_avg(3,:,r,p,1),1,[]),'-^k','linewidth',linewidth,'markersize',markersize);
%         plot(K,1e3*reshape(time_avg(4,:,r,p,1),1,[]),'-vk','linewidth',linewidth,'markersize',markersize);
%         plot(K,1e3*reshape(time_avg(1,:,r,p,1),1,[]),'-o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%         plot(K,1e3*reshape(time_avg(2,:,r,p,1),1,[]),'-s','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%         plot(K,1e3*reshape(time_avg(3,:,r,p,1),1,[]),'-^','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%         plot(K,1e3*reshape(time_avg(4,:,r,p,1),1,[]),'-v','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
%         plot(K,1e3*reshape(time_avg(1,:,r,p,2),1,[]),'--o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%         plot(K,1e3*reshape(time_avg(2,:,r,p,2),1,[]),'--s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%         plot(K,1e3*reshape(time_avg(3,:,r,p,2),1,[]),'--^','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%         plot(K,1e3*reshape(time_avg(4,:,r,p,2),1,[]),'--v','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
%         
%         ylabel('Elapsed time (ms)','fontname',fontname,'fontsize',fontsize);
%         xlabel('Number of users','fontname',fontname,'fontsize',fontsize);
%         
%         legend(legend_M,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_1);
%         legend box off;
%         
%         xlim([K(1) K(end)]);
%         ylim([0 3]);
%         
%         set(gca,'fontname',fontname,'fontsize',fontsize);
%         
%         if (savefig == 1)
%             saveas(gcf,[root_save 'time_ur_los'],'fig');
%             saveas(gcf,[root_save 'time_ur_los'],'png');
%             saveas(gcf,[root_save 'time_ur_los'],'epsc2');
%         end
%     end
% end