clear;
close all;
clc;

% Macros

MC = 1000;                                                                 % Size of the monte-carlo ensemble

M = [64 128 256 512];                                                      % Number of antennas at base station
K = [8 16 32 64];                                                          % Number of mobile users

% M = 128;
% K = 16;

bs_power = [1 10];
radius   = [100 500 1000 2000];

M_SIZ = length(M);                                                         % Size of the antennas set
K_SIZ = length(K);
P_SIZ = length(bs_power);
R_SIZ = length(radius);

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Power Allocation/Downlink/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Power Allocation/';

% chn_type = {'ur_los','rayleigh'};

% Loading data

results = zeros(M_SIZ,K_SIZ,R_SIZ,P_SIZ);

for m = 1:M_SIZ
    for k = 1:K_SIZ
        for p = 1:P_SIZ
            for r = 1:R_SIZ
                load([root_load 'results_upper_bound_ur_los_M_' num2str(M(m)) '_K_' num2str(K(k)) '_cell_radius_' num2str(radius(r)) '_m_BS_power_' num2str(bs_power(p)) '_W_MC_' num2str(MC) '.mat']);
                
                results(m,k,r,p) = result;
            end
        end
    end
end

% Ploting Figures

linewidth  = 3;
markersize = 15;
fontname   = 'Times New Roman';
fontsize   = 30;

BAR_SIZE = 0.8;

marker = {'o','s','^'};

linestyle = {'-','--',':'};

savefig = 0;

legend_power  = {'$1$ W','$10$ W'};
legend_M      = {'$M = 64$','$M = 128$','$M = 256$','$M = 512$'};
legend_K      = {'$K = 8$','$K = 16$','$K = 32$','$K = 64$'};

cell_radius = [100 500 1000 2000];

cat = categorical({'100','500','1000','2000'});
cat = reordercats(cat,{'100','500','1000','2000'});

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

% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% bar(cat,results(1,1,:,:)/MC,BAR_SIZE)
% 
% xlabel('Cell radius (m)','fontname',fontname,'fontsize',fontsize);
% ylabel('$\mathrm{Pr}\{\gamma_{\mathrm{u}}^{(3)} < \gamma^{*}\}$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
% 
% legend(legend_power,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_2);
% legend box off;
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% if (savefig == 1)
%     saveas(gcf,[root_save 'number_iterations_ur_los'],'fig');
%     saveas(gcf,[root_save 'number_iterations_ur_los'],'png');
%     saveas(gcf,[root_save 'number_iterations_ur_los'],'epsc2');
% end

figure;

set(gcf,'position',[0 0 800 600]);

yyaxis left;

plot(cell_radius,reshape(results(1,1,:,1)/MC,1,[]),'-ok','linewidth',linewidth,'markersize',markersize);
hold on;
plot(cell_radius,reshape(results(2,1,:,1)/MC,1,[]),'-sk','linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(3,1,:,1)/MC,1,[]),'-^k','linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(4,1,:,1)/MC,1,[]),'-vk','linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(1,1,:,1)/MC,1,[]),'-o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(2,1,:,1)/MC,1,[]),'-s','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(3,1,:,1)/MC,1,[]),'-^','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(4,1,:,1)/MC,1,[]),'-v','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);

ylabel('$\mathrm{Pr}\{\gamma_{\mathrm{u}}^{(3)} < \gamma^{*}\}$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');

% ylim([0.5 3]);

yyaxis right;

plot(cell_radius,reshape(results(1,1,:,2)/MC,1,[]),'-o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
hold on;
plot(cell_radius,reshape(results(2,1,:,2)/MC,1,[]),'-s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(3,1,:,2)/MC,1,[]),'-^','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(4,1,:,2)/MC,1,[]),'-v','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);

ylabel('$\mathrm{Pr}\{\gamma_{\mathrm{u}}^{(3)} < \gamma^{*}\}$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');

% ylim([5 20]);

xlabel('Cell radius (m)','fontname',fontname,'fontsize',fontsize);

legend(legend_M,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_2);
legend box off;

% xticks([8 16 24 32 40 48 56 64])

% xlim([K(1) K(end)]);

set(gca,'fontname',fontname,'fontsize',fontsize);

if (savefig == 1)
    saveas(gcf,[root_save 'time_ur_los'],'fig');
    saveas(gcf,[root_save 'time_ur_los'],'png');
    saveas(gcf,[root_save 'time_ur_los'],'epsc2');
end

figure;

set(gcf,'position',[0 0 800 600]);

yyaxis left;

plot(cell_radius,reshape(results(1,1,:,1)/MC,1,[]),'-ok','linewidth',linewidth,'markersize',markersize);
hold on;
plot(cell_radius,reshape(results(1,2,:,1)/MC,1,[]),'-sk','linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(1,3,:,1)/MC,1,[]),'-^k','linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(1,4,:,1)/MC,1,[]),'-vk','linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(1,1,:,1)/MC,1,[]),'-o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(1,2,:,1)/MC,1,[]),'-s','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(1,3,:,1)/MC,1,[]),'-^','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(1,4,:,1)/MC,1,[]),'-v','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);

ylabel('$\mathrm{Pr}\{\gamma_{\mathrm{u}}^{(3)} < \gamma^{*}\}$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');

% ylim([0.5 3]);

yyaxis right;

plot(cell_radius,reshape(results(1,1,:,2)/MC,1,[]),'-o','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
hold on;
plot(cell_radius,reshape(results(1,2,:,2)/MC,1,[]),'-s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(1,3,:,2)/MC,1,[]),'-^','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(cell_radius,reshape(results(1,4,:,2)/MC,1,[]),'-v','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);

ylabel('$\mathrm{Pr}\{\gamma_{\mathrm{u}}^{(3)} < \gamma^{*}\}$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');

% ylim([5 20]);

xlabel('Cell radius (m)','fontname',fontname,'fontsize',fontsize);

legend(legend_K,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_2);
legend box off;

% xticks([8 16 24 32 40 48 56 64])

% xlim([K(1) K(end)]);

set(gca,'fontname',fontname,'fontsize',fontsize);

if (savefig == 1)
    saveas(gcf,[root_save 'time_ur_los'],'fig');
    saveas(gcf,[root_save 'time_ur_los'],'png');
    saveas(gcf,[root_save 'time_ur_los'],'epsc2');
end