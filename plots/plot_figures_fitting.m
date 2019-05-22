clear;
close all;
clc;

% Macros

MC = 10000;                                                                % Size of the monte-carlo ensemble

M = [64 256];                                                              % Number of antennas at base station
K = 18;                                                                    % Number of mobile users
L = 13;                                                                    % Number of selected users

snr = (-20:5:10)';                                                         % SNR in dB

M_SIZ = length(M);                                                         % Size of the antennas set
N_SNR = length(snr);                                                       % Size of the SNR set 
N_CHN = 3;                                                                 % Number of channel models simulated

% Loading data

rate     = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                  % Rate using all K users
psi_all  = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                  % Interchannel interference for all users

for m = 1:length(M)
    for snr_idx = 1:N_SNR
        load(['../results/scheduling/uplink/rate_uplink_mf_ur-los_M_' num2str(M(m)) ...
              '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) ...
              '_dB_MC_' num2str(MC) '.mat']);
        
        rate(:,:,m,snr_idx,1)    = rate_u;
        psi_all(:,:,m,snr_idx,1) = psi;
                
        load(['../results/scheduling/uplink/rate_uplink_mf_sparse_M_' num2str(M(m)) ...
              '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) ...
              '_dB_MC_' num2str(MC) '.mat']);
        
        rate(:,:,m,snr_idx,2)    = rate_u;
        psi_all(:,:,m,snr_idx,2) = psi;
                
        load(['../results/scheduling/uplink/rate_uplink_mf_rayleigh_M_' num2str(M(m)) ...
              '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) ...
              '_dB_MC_' num2str(MC) '.mat']);
        
        rate(:,:,m,snr_idx,3)    = rate_u;
        psi_all(:,:,m,snr_idx,3) = psi;        
    end
end

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

savefig = 1;

location = 'northwest';

root = '../figures/fitting/fit_uplink_';

channel_type = {'ur_los','sparse','rayleigh'};

colours = get(gca,'colororder');
close;

p   = zeros(10,N_SNR,M_SIZ,N_CHN);

psi_min = 0;
psi_max = 1;

psi_step = 0.01;

psi_range = (psi_min:psi_step:psi_max);

f   = zeros(length(psi_range),N_SNR,M_SIZ,N_CHN);

for chn_idx = 1:N_CHN
    for m = 1:M_SIZ
        for snr_idx = 1:N_SNR
            if(chn_idx == 3)
                curvefit = fit(psi_all(1,:,m,snr_idx,chn_idx)',rate(1,:,m,snr_idx,chn_idx)','poly4');
            
                p(1,snr_idx,m,chn_idx) = curvefit.p1;
                p(2,snr_idx,m,chn_idx) = curvefit.p2;
                p(3,snr_idx,m,chn_idx) = curvefit.p3;
                p(4,snr_idx,m,chn_idx) = curvefit.p4;
                p(5,snr_idx,m,chn_idx) = curvefit.p5;
            
                f(:,snr_idx,m,chn_idx) = p(1,snr_idx,m,chn_idx)*psi_range.^4 + p(2,snr_idx,m,chn_idx)*psi_range.^3 + ...
                                         p(3,snr_idx,m,chn_idx)*psi_range.^2 + p(4,snr_idx,m,chn_idx)*psi_range + ...
                                         p(5,snr_idx,m,chn_idx);
            else
                curvefit = fit(psi_all(1,:,m,snr_idx,chn_idx)',rate(1,:,m,snr_idx,chn_idx)','poly8');
            
                p(1,snr_idx,m,chn_idx) = curvefit.p1;
                p(2,snr_idx,m,chn_idx) = curvefit.p2;
                p(3,snr_idx,m,chn_idx) = curvefit.p3;
                p(4,snr_idx,m,chn_idx) = curvefit.p4;
                p(5,snr_idx,m,chn_idx) = curvefit.p5;
                p(6,snr_idx,m,chn_idx) = curvefit.p6;
                p(7,snr_idx,m,chn_idx) = curvefit.p7;
                p(8,snr_idx,m,chn_idx) = curvefit.p8;
                p(9,snr_idx,m,chn_idx) = curvefit.p9;        
            
                f(:,snr_idx,m,chn_idx) = p(1,snr_idx,m,chn_idx)*psi_range.^8 + p(2,snr_idx,m,chn_idx)*psi_range.^7 + ...
                                         p(3,snr_idx,m,chn_idx)*psi_range.^6 + p(4,snr_idx,m,chn_idx)*psi_range.^5 + ...
                                         p(5,snr_idx,m,chn_idx)*psi_range.^4 + p(6,snr_idx,m,chn_idx)*psi_range.^3 + ...
                                         p(7,snr_idx,m,chn_idx)*psi_range.^2 + p(8,snr_idx,m,chn_idx)*psi_range + ...
                                         p(9,snr_idx,m,chn_idx); 
            end
            figure;
            
            set(gcf,'position',[0 0 800 600]);
            
            plot(psi_all(1,:,m,snr_idx,chn_idx),rate(1,:,m,snr_idx,chn_idx),'.','color',colours(1,:),'linewidth',linewidth);
            hold on;
            plot(psi_range,f(:,snr_idx,m,chn_idx),'-','color',colours(2,:),'linewidth',linewidth);
            
            xlabel('Interchannel interference','fontname',fontname,'fontsize',fontsize);
            ylabel('Uplink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
            
            legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
            
            set(gca,'fontname',fontname,'fontsize',fontsize);
            
            if(chn_idx == 1)
                if(m == 1)
                    xlim([0 0.4]);
                elseif(m == 2)
                    xlim([0 0.2]);
                end
            elseif(chn_idx == 2)
                if(m == 1)
                    xlim([0 0.6]);
                    ylim([0 4]);
                elseif(m == 2)
                    xlim([0 0.4]);
                    ylim([0 5]);
                end
            elseif(chn_idx == 3)
                if(m == 1)
                    xlim([0.05 0.15]);
                elseif(m==2)
                    xlim([0.03 0.08]);
                end
            end
            
            if(savefig == 1)
                saveas(gcf,[root channel_type{chn_idx} '_M_' num2str(M(m)) '_SNR_' num2str(snr(snr_idx))],'fig');
                saveas(gcf,[root channel_type{chn_idx} '_M_' num2str(M(m)) '_SNR_' num2str(snr(snr_idx))],'png');
                saveas(gcf,[root channel_type{chn_idx} '_M_' num2str(M(m)) '_SNR_' num2str(snr(snr_idx))],'epsc2');
            end           
         end
    end
end