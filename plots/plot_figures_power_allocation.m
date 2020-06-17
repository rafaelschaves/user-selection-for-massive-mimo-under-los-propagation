clear;
close all;
clc;

% Macros

MC = 1000;                                                                 % Size of the monte-carlo ensemble

M = [64 128 256 512];                                                          % Number of antennas at base station
K = [8 16 32 64];                                                          % Number of mobile users

bs_power = 1;
radius   = 500;

M_SIZ = length(M);                                                         % Size of the antennas set
K_SIZ = length(K);
N_ALG = 3;                                                                 % Number of algorithms for perform user scheduling

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Power Allocation/Downlink/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Power Allocation/';

% chn_type = {'ur_los','rayleigh'};

% Loading data

gamma_all   = zeros(MC,N_ALG,K_SIZ,M_SIZ);
gamma_u_all = zeros(MC,N_ALG,K_SIZ,M_SIZ);
n_it_all    = zeros(MC,N_ALG,K_SIZ,M_SIZ);
time_all    = zeros(MC,N_ALG,K_SIZ,M_SIZ);

for m = 1:M_SIZ
    for k = 1:K_SIZ
        load([root_load 'results_ur_los_M_' num2str(M(m)) '_K_' num2str(K(k)) '_cell_radius_' num2str(radius) '_m_BS_power_' num2str(bs_power) '_W_MC_' num2str(MC) '.mat']);
        gamma_all(:,:,k,m)   = gamma;
        gamma_u_all(:,:,k,m) = gamma_u;
        n_it_all(:,:,k,m)    = n_it;
        time_all(:,:,k,m)    = time;
    end
end

% Post Processing

n_it_avg = mean(n_it_all,1);
time_avg = mean(time_all,1);

% Ploting Figures

linewidth  = 3;
markersize = 15;
fontname   = 'Times New Roman';
fontsize   = 30;

marker = {'o','s','^'};

linestyle = {'-','--',':'};

savefig = 0;

% legend_algo = {'Algorithm 1','Algorithm 2','Algorithm 3'};
% legend_algo = {'Algorithm 1','Algorithm 2'};
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

figure;

set(gcf,'position',[0 0 800 600]);

plot(K,reshape(n_it_avg(:,1,:,1),1,[]),'-ok','linewidth',linewidth,'markersize',markersize);
hold on;
plot(K,reshape(n_it_avg(:,1,:,2),1,[]),'-sk','linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,1,:,3),1,[]),'-^k','linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,1,:,4),1,[]),'-vk','linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,1,:,1),1,[]),'-o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,1,:,2),1,[]),'-s','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,1,:,3),1,[]),'-^','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,1,:,4),1,[]),'-v','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,2,:,1),1,[]),'--o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,2,:,2),1,[]),'--s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,2,:,3),1,[]),'--^','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,2,:,4),1,[]),'--v','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,3,:,1),1,[]),':o','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,3,:,2),1,[]),':s','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,3,:,3),1,[]),':^','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(n_it_avg(:,3,:,4),1,[]),':v','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);

xlabel('Number of users','fontname',fontname,'fontsize',fontsize);
ylabel('Number of iterations','fontname',fontname,'fontsize',fontsize);

legend(legend_M,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_3);
legend box off;

xticks([8 16 24 32 40 48 56 64])

xlim([K(1) K(end)]);
ylim([0 40]);

set(gca,'fontname',fontname,'fontsize',fontsize);

if (savefig == 1)
    saveas(gcf,[root_save 'number_iterations_ur_los'],'fig');
    saveas(gcf,[root_save 'number_iterations_ur_los'],'png');
    saveas(gcf,[root_save 'number_iterations_ur_los'],'epsc2');
end

figure;

set(gcf,'position',[0 0 800 600]);

yyaxis left;

plot(K,1e3*reshape(time_avg(:,1,:,1),1,[]),'-ok','linewidth',linewidth,'markersize',markersize);
hold on;
plot(K,1e3*reshape(time_avg(:,1,:,2),1,[]),'-sk','linewidth',linewidth,'markersize',markersize);
plot(K,1e3*reshape(time_avg(:,1,:,3),1,[]),'-^k','linewidth',linewidth,'markersize',markersize);
plot(K,1e3*reshape(time_avg(:,1,:,4),1,[]),'-vk','linewidth',linewidth,'markersize',markersize);
plot(K,1e3*reshape(time_avg(:,1,:,1),1,[]),'-o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,1e3*reshape(time_avg(:,1,:,2),1,[]),'-s','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,1e3*reshape(time_avg(:,1,:,3),1,[]),'-^','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,1e3*reshape(time_avg(:,1,:,4),1,[]),'-v','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,1e3*reshape(time_avg(:,2,:,1),1,[]),'--o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,1e3*reshape(time_avg(:,2,:,2),1,[]),'--s','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,1e3*reshape(time_avg(:,2,:,3),1,[]),'--^','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,1e3*reshape(time_avg(:,2,:,4),1,[]),'--v','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);

ylabel('Elapsed time (ms)','fontname',fontname,'fontsize',fontsize);

ylim([0 3]);

yyaxis right;

plot(K,reshape(time_avg(:,3,:,1),1,[]),':o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
hold on;
plot(K,reshape(time_avg(:,3,:,2),1,[]),':s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(time_avg(:,3,:,3),1,[]),':^','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(K,reshape(time_avg(:,3,:,4),1,[]),':v','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);

ylabel('Elapsed time (s)','fontname',fontname,'fontsize',fontsize);

ylim([0 20]);

xlabel('Number of users','fontname',fontname,'fontsize',fontsize);

legend(legend_M,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_1);
legend box off;

xticks([8 16 24 32 40 48 56 64])

xlim([K(1) K(end)]);

set(gca,'fontname',fontname,'fontsize',fontsize);

if (savefig == 1)
    saveas(gcf,[root_save 'time_ur_los'],'fig');
    saveas(gcf,[root_save 'time_ur_los'],'png');
    saveas(gcf,[root_save 'time_ur_los'],'epsc2');
end