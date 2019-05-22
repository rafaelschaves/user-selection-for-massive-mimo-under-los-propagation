clear;
close all;
clc;

root_load = '../results/comp_time/comp_time_';
root_fig  = '../figures/comp_time/comp_time_';

% Macros

MC    = 5000;                                                              % Size of the monte-carlo ensemble
N_ALG = 4;                                                                 % Number of algorithms for perform user scheduling
M_SIZ = 2;                                                                 % Size of the antennas set
N_CHN = 3;                                                                 % Number of channel models simulated

M = [64 256];                                                              % Number of antennas at base station
K = 18;                                                                    % Number of mobile users
L = K - 1;                                                                 % Maximum number of selected users

% Loading data

all_comp_time = zeros(MC,L,N_ALG,N_CHN,M_SIZ);                             % Computational time

for m_idx = 1:M_SIZ
    load([root_load 'M_' num2str(M(m_idx)) '_K_' num2str(K) '_MC_' num2str(MC) '.mat']);
    
    all_comp_time(:,:,:,:,m_idx) = comp_time(:,:,:,:);
end

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

savefig = 1;

% RS - Random selection
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_algo = {'RS','SOS','CBS','ICIBS'};
channel_mod = {'ur_los','sparse','rayleigh'};

location_1 = 'northeast';
location_2 = 'northwest';

colours = get(gca,'colororder');
close;

ylim_1 = [0 5];
ylim_2 = [0 18];

for chn_idx = 1:N_CHN
    for m_idx = 1:M_SIZ
        figure;
        
        set(gcf,'position',[0 0 800 600]);
        
        plot(1:L,1e3*mean(all_comp_time(:,:,1,chn_idx,m_idx)),'-','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:L,1e3*mean(all_comp_time(:,:,2,chn_idx,m_idx)),'-','color',colours(2,:),'linewidth',linewidth);
        plot(1:L,1e3*mean(all_comp_time(:,:,3,chn_idx,m_idx)),'-','color',colours(3,:),'linewidth',linewidth);
        plot(1:L,1e3*mean(all_comp_time(:,:,4,chn_idx,m_idx)),'-','color',colours(4,:),'linewidth',linewidth);
        
        xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
        ylabel('Computational time (ms)','fontname',fontname,'fontsize',fontsize);
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        if(m_idx == 1)
            legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_2);

            xlim([1 L]);
            ylim(ylim_1);
        else
            legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_2);

            xlim([1 L]);
            ylim(ylim_2);
        end
        
        if(savefig == 1)
            saveas(gcf,[root_fig channel_mod{chn_idx} '_M_' num2str(M(m_idx)) '_K_' num2str(K) '_MC_' num2str(MC)],'fig');
            saveas(gcf,[root_fig channel_mod{chn_idx} '_M_' num2str(M(m_idx)) '_K_' num2str(K) '_MC_' num2str(MC)],'png');
            saveas(gcf,[root_fig channel_mod{chn_idx} '_M_' num2str(M(m_idx)) '_K_' num2str(K) '_MC_' num2str(MC)],'epsc2');
        end
    end
end