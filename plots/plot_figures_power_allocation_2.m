clear;
close all;
clc;

% Macros

MC = 1000;                                                                 % Size of the monte-carlo ensemble

M = [64 128 256 512];                                                      % Number of antennas at base station
K = [8 16 32 64];                                                          % Number of mobile users

% M = 64;
% K = 8;

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

savefig = 1;

legend_power  = {'$P = 1$ W','$P = 10$ W'};
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

for m = 1:M_SIZ
    for k = 1:K_SIZ
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        bar(cat,reshape(results(m,k,:,:),[],2)/MC,BAR_SIZE)
        
        xlabel('Cell radius (m)','fontname',fontname,'fontsize',fontsize);
        ylabel('$\mathrm{Pr}\{\gamma_{\mathrm{r}}^{(3)} < \gamma^{\star}\}$','fontname',fontname,'fontsize',fontsize,'interpreter','latex');
        
        legend(legend_power,'fontname',fontname,'fontsize',fontsize,'interpreter','latex','location',location_2);
        legend box off;
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        if (savefig == 1)
            saveas(gcf,[root_save_fig 'bound_test_M_' num2str(M(m)) '_K_' num2str(K(k))],'fig');
            saveas(gcf,[root_save_png 'bound_test_M_' num2str(M(m)) '_K_' num2str(K(k))],'png');
            saveas(gcf,[root_save_eps 'bound_test_M_' num2str(M(m)) '_K_' num2str(K(k))],'epsc2');
        end
    end
end