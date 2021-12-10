clear;
close all;
clc;

% Roots

root_save = 'C:\Users\RafaelChaves\Google Drive\UFRJ\PhD\Codes\user-scheduling-massive-mimo\Figures\Selection\Downlink\';

% SOS

A_sos = @(M,K,L) (M^3 + M^2 - M).*(L - 1).*(2*K - L) + (2*M - 1)*(K*L - L.*(L - 1))/2 + 2*M^2.*L.*(L - 1);
M_sos = @(M,K,L) 2*M^2.*(L - 1).*(2*M*K - M*L + L) + M*(2*K*L - L.*(L - 1));
D_sos = @(M,K,L) 2*M.*L;
S_sos = @(M,K,L) K*L - L.*(L - 1)/2;

% S-SOS

A_ssos = @(M,K,L) (M^3 + M^2 - M).*(L - 1).*(2*K - L) + (2*M - 1)*(K*L - L.*(L - 1))/2 + 4*M^2.*(L - 1);
M_ssos = @(M,K,L) 2*M^2.*(L - 1).*(2*M*K - M*L + 2) + M*(2*K*L - L.*(L - 1));
D_ssos = @(M,K,L) 2*M.*L;
S_ssos = @(M,K,L) K*L - L.*(L - 1)/2;

% CBS

A_cbs = @(M,K,L) (4*K^2*M + K*(K - 1)/2).*ones(1,length(L));
M_cbs = @(M,K,L) (4*K^2*M + K^2 - K).*ones(1,length(L));
D_cbs = @(M,K,L) 2*M.*K.*ones(1,length(L));
S_cbs = @(M,K,L) K*(K - 1)/2.*ones(1,length(L));

% ICIBS

A_icibs = @(M,K,L) 4*K^2*M + K*(K - 1)*(2*K + 5)/6 - L.*(L + 1).*(L - 1)/3;
M_icibs = @(M,K,L) (4*K^2*M + K^2 - K).*ones(1,length(L));
D_icibs = @(M,K,L) 2*M*K + K*(K + 1)/2 - L.*(L + 1)/2;
S_icibs = @(M,K,L) K*(K - 1)/2.*ones(1,length(L));

% Fast-ICIBS

A_sicibs = @(M,K,L) 4*K^2*M + K^2 - L.*(L + 1)/2;
M_sicibs = @(M,K,L) 4*K^2*M + K*(3*K - 1) - L.*(L + 1)/2;
D_sicibs = @(M,K,L) 2*M*K + K*(K + 1)/2 - L.*(L + 1)/2;
S_sicibs = @(M,K,L) K*(K - 1)/2.*ones(1,length(L));

M_1 = 50;
K_1 = 100;
L_1 = 1:50;

M_2 = 100;
K_2 = 100;
L_2 = 1:100;

flop_sos_1    = A_sos(M_1,K_1,L_1)    + M_sos(M_1,K_1,L_1)    + D_sos(M_1,K_1,L_1)    + S_sos(M_1,K_1,L_1);
flop_ssos_1   = A_ssos(M_1,K_1,L_1)   + M_ssos(M_1,K_1,L_1)   + D_ssos(M_1,K_1,L_1)   + S_ssos(M_1,K_1,L_1);
flop_cbs_1    = A_cbs(M_1,K_1,L_1)    + M_cbs(M_1,K_1,L_1)    + D_cbs(M_1,K_1,L_1)    + S_cbs(M_1,K_1,L_1);
flop_icibs_1  = A_icibs(M_1,K_1,L_1)  + M_icibs(M_1,K_1,L_1)  + D_icibs(M_1,K_1,L_1)  + S_icibs(M_1,K_1,L_1);
flop_sicibs_1 = A_sicibs(M_1,K_1,L_1) + M_sicibs(M_1,K_1,L_1) + D_sicibs(M_1,K_1,L_1) + S_sicibs(M_1,K_1,L_1);

flop_sos_2    = A_sos(M_2,K_2,L_2)    + M_sos(M_2,K_2,L_2)    + D_sos(M_2,K_2,L_2)    + S_sos(M_2,K_2,L_2);
flop_ssos_2   = A_ssos(M_2,K_2,L_2)   + M_ssos(M_2,K_2,L_2)   + D_ssos(M_2,K_2,L_2)   + S_ssos(M_2,K_2,L_2);
flop_cbs_2    = A_cbs(M_2,K_2,L_2)    + M_cbs(M_2,K_2,L_2)    + D_cbs(M_2,K_2,L_2)    + S_cbs(M_2,K_2,L_2);
flop_icibs_2  = A_icibs(M_2,K_2,L_2)  + M_icibs(M_2,K_2,L_2)  + D_icibs(M_2,K_2,L_2)  + S_icibs(M_2,K_2,L_2);
flop_sicibs_2 = A_sicibs(M_2,K_2,L_2) + M_sicibs(M_2,K_2,L_2) + D_sicibs(M_2,K_2,L_2) + S_sicibs(M_2,K_2,L_2);

% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 30;

marker    = {'o','s','^','v','x'};
linestyle = {'-','--',':'};

savefig = 1;
 
legend_algo = {'SOS','S-SOS','CBS','ICIBS','S-ICIBS'};

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

plot(L_1,flop_sos_1*1e-9   ,'-' ,'color',colours(1,:),'linewidth',linewidth);
hold on;
plot(L_1,flop_ssos_1*1e-9  ,'-' ,'color',colours(5,:),'linewidth',linewidth);
plot(L_1,flop_cbs_1*1e-9   ,'-' ,'color',colours(2,:),'linewidth',linewidth);
plot(L_1,flop_icibs_1*1e-9 ,'-' ,'color',colours(3,:),'linewidth',linewidth);
plot(L_1,flop_sicibs_1*1e-9,'-' ,'color',colours(4,:),'linewidth',linewidth);
plot(L_2,flop_sos_2*1e-9   ,'--','color',colours(1,:),'linewidth',linewidth);
plot(L_2,flop_ssos_2*1e-9  ,'--','color',colours(5,:),'linewidth',linewidth);
plot(L_2,flop_cbs_2*1e-9   ,'--','color',colours(2,:),'linewidth',linewidth);
plot(L_2,flop_icibs_2*1e-9 ,'--','color',colours(3,:),'linewidth',linewidth);
plot(L_2,flop_sicibs_2*1e-9,'--','color',colours(4,:),'linewidth',linewidth);

ylabel('GFlops','fontname',fontname,'fontsize',fontsize);
xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);

legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',2);
legend box off;

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([1 max(L_2)]);

axes('position',[.5 .3 .35 .35]);
box on;

plot(L_1,flop_cbs_1*1e-6   ,'-' ,'color',colours(2,:),'linewidth',linewidth);
hold on;
plot(L_1,flop_icibs_1*1e-6 ,'-' ,'color',colours(3,:),'linewidth',linewidth);
plot(L_1,flop_sicibs_1*1e-6,'-' ,'color',colours(4,:),'linewidth',linewidth);
plot(L_2,flop_cbs_2*1e-6   ,'--','color',colours(2,:),'linewidth',linewidth);
plot(L_2,flop_icibs_2*1e-6 ,'--','color',colours(3,:),'linewidth',linewidth);
plot(L_2,flop_sicibs_2*1e-6,'--','color',colours(4,:),'linewidth',linewidth);

ylabel('MFlops','fontname',fontname,'fontsize',fontsize);

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([1 max(L_2)]);
ylim([3 9]);

if savefig == 1
    saveas(gcf,[root_save 'computational_complexity'],'fig');
    saveas(gcf,[root_save 'computational_complexity'],'png');
    saveas(gcf,[root_save 'computational_complexity'],'epsc2');
end

% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% yyaxis left
% 
% plot(L_1,flop_sos_1*1e-9 ,'-' ,'linewidth',linewidth);
% hold on;
% plot(L_1,flop_ssos_1*1e-9,'--','linewidth',linewidth);
% 
% ylabel('GFlops Count','fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 2]);
% 
% yyaxis right
% 
% plot(L_1,flop_cbs_1*1e-6   ,'-' ,'linewidth',linewidth);
% hold on;
% plot(L_1,flop_icibs_1*1e-6 ,'--','linewidth',linewidth);
% plot(L_1,flop_sicibs_1*1e-6,':' ,'linewidth',linewidth);
% 
% % ylabel('MFlops Count','fontname',fontname,'fontsize',fontsize);
% 
% ylim([2.2 2.5]);
% 
% xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',3);
% legend box off;
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% xlim([1 max(L_1)]);
% 
% if savefig == 1
%     saveas(gcf,[root_save 'computational_complexity_M_' num2str(M_1) '_K_' num2str(K_1)],'fig');
%     saveas(gcf,[root_save 'computational_complexity_M_' num2str(M_1) '_K_' num2str(K_1)],'png');
%     saveas(gcf,[root_save 'computational_complexity_M_' num2str(M_1) '_K_' num2str(K_1)],'epsc2');
% end
% 
% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% yyaxis left
% 
% plot(L_2,flop_sos_2*1e-9 ,'-' ,'linewidth',linewidth);
% hold on;
% plot(L_2,flop_ssos_2*1e-9,'--','linewidth',linewidth);
% 
% % ylabel('GFlops Count','fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 60]);
% 
% yyaxis right
% 
% plot(L_2,flop_cbs_2*1e-6   ,'-' ,'linewidth',linewidth);
% hold on;
% plot(L_2,flop_icibs_2*1e-6 ,'--','linewidth',linewidth);
% plot(L_2,flop_sicibs_2*1e-6,':' ,'linewidth',linewidth);
% 
% ylabel('MFlops Count','fontname',fontname,'fontsize',fontsize);
% 
% ylim([18 20]);
% 
% xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
% 
% % legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',3);
% % legend box off;
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% xlim([1 max(L_2)]);
% 
% if savefig == 1
%     saveas(gcf,[root_save 'computational_complexity_M_' num2str(M_2) '_K_' num2str(K_2)],'fig');
%     saveas(gcf,[root_save 'computational_complexity_M_' num2str(M_2) '_K_' num2str(K_2)],'png');
%     saveas(gcf,[root_save 'computational_complexity_M_' num2str(M_2) '_K_' num2str(K_2)],'epsc2');
% end