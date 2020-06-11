clear;
close all;
clc;

% Macros

M = [50 100 200]';
K = [10 25 50 75 100 150]';

M_SIZ = length(M);
K_SIZ = length(K);
N_ALG = 3;                                                                 % Number of algorithms for perform user scheduling
N_PRE = 3;

% Roots

root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Figures/Selection/Downlink/';

zero_pad_1 = '%03d';
zero_pad_2 = '%02d';

chn_type = 'ur_los';

sum_se = zeros(K_SIZ,M_SIZ,N_PRE);
L_star = zeros(K_SIZ,N_ALG,M_SIZ,N_PRE);

sum_se(:,1,1) = [11.8593 14.2438 15.2947 15.6819 0 0]';
sum_se(:,1,2) = [10.5841 7.6095  2.5325  0.5366  0 0]';
sum_se(:,1,3) = [11.6354 11.9242 10.7754 10.4741 0 0]';

sum_se(:,2,1) = [18.6069 24.8489 28.2778 29.6633 30.4046 31.1983]';
sum_se(:,2,2) = [18.3592 20.2795 15.1036 9.8819 5.7731 0.9448]';
sum_se(:,2,3) = [18.8707 22.9825 22.1466 20.3252 19.0472 18.2715]';

sum_se(:,3,1) = [26.8831 40.1241 49.3859 53.7560 56.3416 59.1804]';
sum_se(:,3,2) = [27.3405 38.2939 39.8475 35.6565 30.4775 20.5494]';
sum_se(:,3,3) = [27.5728 39.6931 44.4186 43.6460 41.5768 37.1698]';

L_star(:,:,1,1) = [8 8 8; 17 17 17; 27 28 27; 33 35 33; 0 0 0; 0 0 0];
L_star(:,:,1,2) = [8 8 8; 16 16 16; 24 24 24; 28 28 29; 0 0 0; 0 0 0];
L_star(:,:,1,3) = [8 8 8; 16 16 16; 24 24 25; 29 29 29; 0 0 0; 0 0 0];

L_star(:,:,2,1) = [9 9 9; 20 20 20; 35 35 34; 46 46 45; 54 55 54; 67 69 65];
L_star(:,:,2,2) = [9 9 9; 19 19 19; 32 32 32; 41 41 42; 48 48 49; 57 57 58];
L_star(:,:,2,3) = [9 9 9; 20 20 20; 32 32 33; 42 42 42; 48 49 49; 58 58 59];

L_star(:,:,3,1) = [9 9 9; 22 22 22; 40 40 40; 56 56 55; 69 70 69; 91 93 91];
L_star(:,:,3,2) = [9 9 9; 22 22 22; 39 39 39; 53 53 53; 65 65 65; 84 84 84];
L_star(:,:,3,3) = [9 9 9; 22 22 23; 39 39 39; 54 53 54; 65 65 66; 84 84 84];

% Ploting Figures

linewidth  = 3;
markersize = 15;
fontname   = 'Times New Roman';
fontsize   = 30;

BAR_SIZE = 0.8;

marker = {'o','s','^'};

linestyle = {'-','--',':'};

savefig = 1;

legend_pre = {'MRT','ZF','MMSE'};
legend_alg = {'SOS','CBS', 'ICIBS'};
legend_M = {'$M = 50$','$M = 100$','$M = 200$'};

cat1 = categorical({'10','25','50','75'});
cat1 = reordercats(cat1,{'10','25','50','75'});

cat2 = categorical({'10','25','50','75','100','150'});
cat2 = reordercats(cat2,{'10','25','50','75','100','150'});

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

plot(K(1:4),sum_se(1:4,1,1),'-o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
hold on;
plot(K,sum_se(:,2,1),'-s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(K,sum_se(:,3,1),'-v','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
plot(K(1:4),sum_se(1:4,1,2),'--o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,sum_se(:,2,2),'--s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(K,sum_se(:,3,2),'--v','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);
plot(K(1:4),sum_se(1:4,1,3),':o','color',colours(1,:),'linewidth',linewidth,'markersize',markersize);
plot(K,sum_se(:,2,3),':s','color',colours(2,:),'linewidth',linewidth,'markersize',markersize);
plot(K,sum_se(:,3,3),':v','color',colours(3,:),'linewidth',linewidth,'markersize',markersize);

xlabel('Number of users','fontname',fontname,'fontsize',fontsize);
ylabel('Sum-spectral efficiency','fontname',fontname,'fontsize',fontsize,'interpreter','latex');

legend(legend_M,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_1);
legend box off;

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([10 150]);
%ylim([0 40]);

if (savefig == 1)
    saveas(gcf,[root_save 'NS_sum_se'],'fig');
    saveas(gcf,[root_save 'NS_sum_se'],'png');
    saveas(gcf,[root_save 'NS_sum_se'],'epsc2');
end

for m = 1:M_SIZ
    for n_pre = 1:N_PRE
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        if m == 1
            bar(cat1,L_star(1:4,:,m,n_pre),BAR_SIZE);
        else
            bar(cat2,L_star(:,:,m,n_pre),BAR_SIZE);
        end
        
        %         xtipsl = b(1).XData - BAR_SIZE;
        %         ytips1 = b(1).YData;
        %         labels1 = string(b(1).YData);
        %         text(xtipsl,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom');
        %
        %         xtips2 = b(2).XData;
        %         ytips2 = b(2).YData;
        %         labels2 = string(b(2).YData);
        %         text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom');
        %
        %         xtips3 = b(3).XData + BAR_SIZE;
        %         ytips3 = b(3).YData;
        %         labels3 = string(b(3).YData);
        %         text(xtips3,ytips3,labels3,'HorizontalAlignment','center','VerticalAlignment','bottom');
        
        xlabel('Number of users','fontname',fontname,'fontsize',fontsize);
        ylabel('$L^{\star}$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
        
        legend(legend_alg,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_1);
        legend box off;
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        if (savefig == 1)
            saveas(gcf,[root_save 'L_star_M_' num2str(M(m)) '_' legend_pre{n_pre}],'fig');
            saveas(gcf,[root_save 'L_star_M_' num2str(M(m)) '_' legend_pre{n_pre}],'png');
            saveas(gcf,[root_save 'L_star_M_' num2str(M(m)) '_' legend_pre{n_pre}],'epsc2');
        end
    end
end