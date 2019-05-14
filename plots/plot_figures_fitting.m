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
        curvefit_ur_los   = fit(psi_ur_los(k,:,snr_idx)',rate_u_ur_los(k,:,snr_idx)','exp2');
        curvefit_sparse   = fit(psi_sparse(k,:,snr_idx)',rate_u_sparse(k,:,snr_idx)','exp2');
        curvefit_rayleigh = fit(psi_rayleigh(k,:,snr_idx)',rate_u_rayleigh(k,:,snr_idx)','poly3');
        
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
        
        plot(psi_ur_los(k,:,snr_idx),rate_u_ur_los(k,:,snr_idx),'.','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(psi_range_ur_los,f_u_ur_los(k,:,snr_idx),'-','color',colours(2,:),'linewidth',linewidth);
        
        xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
        ylabel('Uplink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        
        legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([psi_ur_los_min psi_ur_los_max]);
        
        subplot(1,3,2);
        
        plot(psi_sparse(k,:,snr_idx),rate_u_sparse(k,:,snr_idx),'.','color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(psi_range_sparse,f_u_sparse(k,:,snr_idx),'-','color',colours(2,:),'linewidth',linewidth);
        
        xlabel('Interchannel inteference','fontname',fontname,'fontsize',fontsize);
        ylabel('Uplink rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
        
        legend({'Data','Fitted curve'},'fontname',fontname,'fontsize',fontsize);
        
        set(gca,'fontname',fontname,'fontsize',fontsize);
        
        xlim([psi_sparse_min psi_sparse_max]);
        
        subplot(1,3,3)
        
        plot(psi_rayleigh(k,:,snr_idx),rate_u_rayleigh(k,:,snr_idx),'.','color',colours(1,:),'linewidth',linewidth);
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
