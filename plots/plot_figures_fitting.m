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
N_ALG = 4;                                                                 % Number of algorithms for perform user scheduling
N_SNR = length(snr);                                                       % Size of the SNR set 
N_CHN = 3;                                                                 % Number of channel models simulated

% Loading data

rate = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                      % Rate using all K users
psi_  = zeros(K,MC,M_SIZ,N_SNR,N_CHN);                                      % Interchannel interference for all users

for m_idx = 1:length(M)
    for snr_idx = 1:N_SNR
        load(['../results/uplink/rate_uplink_mf_ur-los_M_' num2str(M(m_idx)) ...
            '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) ...
            '_dB_MC_' num2str(MC) '.mat']);
        
        rate(:,:,m_idx,snr_idx,1) = rate_u;
        psi_(:,:,m_idx,snr_idx,1)  = psi;
                
        load(['../results/uplink/rate_uplink_mf_sparse_M_' num2str(M(m_idx)) ...
            '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) ...
            '_dB_MC_' num2str(MC) '.mat']);
        
        rate(:,:,m_idx,snr_idx,2) = rate_u;
        psi_(:,:,m_idx,snr_idx,2)  = psi;
                
        load(['../results/uplink/rate_uplink_mf_rayleigh_M_' num2str(M(m_idx)) ...
            '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr(snr_idx)) ...
            '_dB_MC_' num2str(MC) '.mat']);
        
        rate(:,:,m_idx,snr_idx,3) = rate_u;
        psi_(:,:,m_idx,snr_idx,3)  = psi;        
    end
end

% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

savefig = 0;

% NS - No selection
% RS - Random selection
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_algo = {'NS','RS','SOS','CBS','ICIBS'};
channel_mod = {'ur_los','sparse','rayleigh'};

location = 'northwest';

root_out_prob = '../figures/rate/out_prob_';

colours = get(gca,'colororder');
close;

p_u_ur_los   = zeros(10,K,N_SNR);
p_u_sparse   = zeros(10,K,N_SNR);
p_u_rayleigh = zeros(10,K,N_SNR);

psi_ur_los_min = 0;
psi_ur_los_max = 0.45;

psi_sparse_min = 0;
psi_sparse_max = 0.7;

psi_rayleigh_min = 0;
psi_rayleigh_max = 0.2;

psi_step = 0.01;

psi_range_ur_los   = (psi_ur_los_min:psi_step:psi_ur_los_max);
psi_range_sparse   = (psi_sparse_min:psi_step:psi_sparse_max);
psi_range_rayleigh = (psi_rayleigh_min:psi_step:psi_rayleigh_max);

f_u_ur_los   = zeros(K,length(psi_range_ur_los),N_SNR);
f_u_sparse   = zeros(K,length(psi_range_sparse),N_SNR);
f_u_rayleigh = zeros(K,length(psi_range_rayleigh),N_SNR);

for snr_idx = 1:N_SNR
    for k = 1:1
        curvefit_ur_los   = fit(psi_(k,:,1,snr_idx,1)',rate(k,:,1,snr_idx,1)','exp2');
        curvefit_sparse   = fit(psi_(k,:,1,snr_idx,2)',rate(k,:,1,snr_idx,2)','exp2');
        curvefit_rayleigh = fit(psi_(k,:,1,snr_idx,3)',rate(k,:,1,snr_idx,3)','poly3');
        
        % p_u_ur_los(1,k,snr_idx)  = curvefit_ur_los.p1;
        % p_u_ur_los(2,k,snr_idx)  = curvefit_ur_los.p2;
        % p_u_ur_los(3,k,snr_idx)  = curvefit_ur_los.p3;
        % p_u_ur_los(4,k,snr_idx)  = curvefit_ur_los.p4;
        % p_u_ur_los(5,k,snr_idx)  = curvefit_ur_los.p5;
        % p_u_ur_los(6,k,snr_idx)  = curvefit_ur_los.p6;
        % p_u_ur_los(7,k,snr_idx)  = curvefit_ur_los.p7;
        % p_u_ur_los(8,k,snr_idx)  = curvefit_ur_los.p8;
        % p_u_ur_los(9,k,snr_idx)  = curvefit_ur_los.p9;
        % p_u_ur_los(10,k,snr_idx) = curvefit_ur_los.p10;
                                                        
        p_u_ur_los(1,k,snr_idx) = curvefit_ur_los.a;
        p_u_ur_los(2,k,snr_idx) = curvefit_ur_los.b;
        p_u_ur_los(3,k,snr_idx) = curvefit_ur_los.c;
        p_u_ur_los(4,k,snr_idx) = curvefit_ur_los.d;
        
        % p_u_sparse(1,k,snr_idx) = curvefit_sparse.p1;
        % p_u_sparse(2,k,snr_idx) = curvefit_sparse.p2;
        % p_u_sparse(3,k,snr_idx) = curvefit_sparse.p3;
        % p_u_sparse(4,k,snr_idx) = curvefit_sparse.p4;
        
        p_u_sparse(1,k,snr_idx) = curvefit_sparse.a;
        p_u_sparse(2,k,snr_idx) = curvefit_sparse.b;
        p_u_sparse(3,k,snr_idx) = curvefit_sparse.c;
        p_u_sparse(4,k,snr_idx) = curvefit_sparse.d;
        
        p_u_rayleigh(1,k,snr_idx) = curvefit_rayleigh.p1;
        p_u_rayleigh(2,k,snr_idx) = curvefit_rayleigh.p2;
        p_u_rayleigh(3,k,snr_idx) = curvefit_rayleigh.p3;
        p_u_rayleigh(4,k,snr_idx) = curvefit_rayleigh.p4;
        
        % f_u_ur_los(k,:,snr_idx)   = p_u_ur_los(1,k,snr_idx)*psi_range_ur_los.^9 + p_u_ur_los(2,k,snr_idx)*psi_range_ur_los.^8 + ...
        %                             p_u_ur_los(3,k,snr_idx)*psi_range_ur_los.^7 + p_u_ur_los(4,k,snr_idx)*psi_range_ur_los.^6 + ...
        %                             p_u_ur_los(5,k,snr_idx)*psi_range_ur_los.^5 + p_u_ur_los(6,k,snr_idx)*psi_range_ur_los.^4 + ...
        %                             p_u_ur_los(7,k,snr_idx)*psi_range_ur_los.^3 + p_u_ur_los(8,k,snr_idx)*psi_range_ur_los.^2 + ...
        %                             p_u_ur_los(9,k,snr_idx)*psi_range_ur_los + p_u_ur_los(10,k,snr_idx);
        f_u_ur_los(k,:,snr_idx) = p_u_ur_los(1,k,snr_idx)*exp(p_u_ur_los(2,k,snr_idx)*psi_range_ur_los) + ...
                                  p_u_ur_los(3,k,snr_idx)*exp(p_u_ur_los(4,k,snr_idx)*psi_range_ur_los);
        % f_u_sparse(k,:,snr_idx)   = p_u_sparse(1,k,snr_idx)*psi_range_sparse.^3 + p_u_sparse(2,k,snr_idx)*psi_range_sparse.^2 + ...
        %                             p_u_sparse(3,k,snr_idx)*psi_range_sparse + p_u_sparse(4,k,snr_idx);
        f_u_sparse(k,:,snr_idx)   = p_u_sparse(1,k,snr_idx)*exp(p_u_sparse(2,k,snr_idx)*psi_range_sparse) + ...
                                    p_u_sparse(3,k,snr_idx)*exp(p_u_sparse(4,k,snr_idx)*psi_range_sparse);
        f_u_rayleigh(k,:,snr_idx) = p_u_rayleigh(1,k,snr_idx)*psi_range_rayleigh.^3 + p_u_rayleigh(2,k,snr_idx)*psi_range_rayleigh.^2 + ...
                                    p_u_rayleigh(3,k,snr_idx)*psi_range_rayleigh + p_u_rayleigh(4,k,snr_idx);
        
        figure;
        
        set(gcf,'position',[500 250 1200 600]);
        
        subplot(1,3,1);
        
        plot(psi_(k,:,1,snr_idx,1),rate(k,:,1,snr_idx,1),'.','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(psi_range_ur_los,f_u_ur_los(k,:,snr_idx),'-','color',colours(2,:),'linewidth',linewidth);
        
        xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
        ylabel('Uplink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        
        legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([psi_ur_los_min psi_ur_los_max]);
        
        subplot(1,3,2);
        
        plot(psi_(k,:,1,snr_idx,2),rate(k,:,1,snr_idx,2),'.','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(psi_range_sparse,f_u_sparse(k,:,snr_idx),'-','color',colours(2,:),'linewidth',linewidth);
        
        xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
        ylabel('Uplink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        
        legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([psi_sparse_min psi_sparse_max]);
        
        subplot(1,3,3)
        
        plot(psi_(k,:,1,snr_idx,3),rate(k,:,1,snr_idx,3),'.','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(psi_range_rayleigh,f_u_rayleigh(k,:,snr_idx),'-','color',colours(2,:),'linewidth',linewidth);
        
        xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
        ylabel('Uplink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        
        legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([psi_rayleigh_min psi_rayleigh_max]);
        
        if(savefig == 1)
            saveas(gcf,[root_rate_fit 'M_' num2str(M) '_' num2str(k) '_user'],'fig');
            saveas(gcf,[root_rate_fit 'M_' num2str(M) '_' num2str(k) '_user'],'png');
            saveas(gcf,[root_rate_fit 'M_' num2str(M) '_' num2str(k) '_user'],'epsc2');
        end
    end
end