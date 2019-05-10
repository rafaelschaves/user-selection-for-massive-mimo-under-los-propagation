clear;
close all;
clc;

M = [64 256];
K = 18;
L = 13;

N_CHN = 3;

snr = [-20 -15 -10 -5 0 5 10]';

likely_capacity(:,:,1,1) = [0.54625 1.02980 1.58030 2.06970 2.42570 2.61270 2.72580;
                            0.56375 1.08370 1.68720 2.25370 2.69870 2.94720 3.11330;
                            0.59725 1.17830 1.88730 2.55770 3.07680 3.37570 3.58220;
                            0.68025 1.44830 2.44570 3.42930 4.19080 4.64480 4.89530;
                            0.68225 1.45380 2.44930 3.44930 4.20180 4.66230 4.89670];

likely_capacity(:,:,1,2) = [0.58625 0.59675 0.59425 0.60725 0.60225 0.60175 0.60725;
                            0.70525 0.72075 0.72775 0.72975 0.73575 0.72975 0.73575;
                            0.69875 0.71975 0.71875 0.72575 0.72725 0.71975 0.72325;
                            1.11730 1.15980 1.16880 1.17630 1.16780 1.17980 1.17730;
                            1.14030 1.18630 1.19780 1.20430 1.20430 1.20530 1.21080];
                        
likely_capacity(:,:,1,3) = [0.60575 1.16630 1.67780 1.96230 2.08020 2.11670 2.12920;
                            0.62275 1.24720 1.86630 2.24570 2.41120 2.46670 2.47620;
                            0.63225 1.25930 1.88130 2.25230 2.41720 2.47020 2.49070;
                            0.63275 1.29280 1.98430 2.42720 2.62470 2.68980 2.71470;
                            0.63525 1.30730 2.02770 2.49520 2.70770 2.78030 2.80930];
                        
likely_capacity(:,:,2,1) = [1.4848 2.4077 3.3257 4.1358 4.7928 5.2687 5.4833;
                            1.5123 2.4632 3.4347 4.2943 5.0498 5.5988 5.9278;
                            1.6113 2.6343 3.7203 4.6698 5.5428 6.1573 6.5443;
                            1.8228 3.1498 4.6193 6.0413 7.2293 8.1267 8.6542;
                            1.8228 3.1503 4.6223 6.0443 7.2538 8.1677 8.6932];
                        
likely_capacity(:,:,2,2) = [1.1793 1.2003 1.2023 1.2028 1.2007 1.2178 1.2018;
                            1.3618 1.3968 1.3843 1.3938 1.3993 1.4038 1.3998;
                            1.3682 1.4143 1.4093 1.4058 1.4088 1.4198 1.4173;
                            2.1942 2.2463 2.2483 2.2602 2.2462 2.2877 2.2612;
                            2.2348 2.2977 2.3022 2.3073 2.2992 2.3323 2.3132];

likely_capacity(:,:,2,3) = [1.6438 2.5903 3.2942 3.6547 3.7978 3.8492 3.8582;
                            1.6827 2.7132 3.5397 3.9973 4.1878 4.2542 4.2753;
                            1.6898 2.7212 3.5448 4.0003 4.1873 4.2553 4.2728;
                            1.7043 2.7862 3.6923 4.2168 4.4437 4.5258 4.5507;
                            1.7108 2.8102 3.7422 4.2988 4.5398 4.6298 4.6573];
% Ploting Figures

linewidth  = 2;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

BIN_WIDTH_CDF  = 0.005;

BAR_SIZE = 0.8;

% NS - No selection
% RS - Random selection
% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_algo = {'NS','RS','SOS','CBS','ICIBS'};
legend_link = {'Uplink','Downlink'};
channel_mod = {'ur_los','sparse','rayleigh'};

location = 'northwest';

cat = categorical(legend_algo);
cat = reordercats(cat,legend_algo);

root_rate_lik = '../figures/rate/likely_';

savefig = 1;

colours = get(gca,'colororder');
close;                        

% figure;
%     
% set(gcf,'position',[0 0 800 600]);
%     
% plot(snr,likely_capacity_ur_los(1,:),'-','color',colours(1,:),'linewidth',linewidth);
% hold on;
% plot(snr,likely_capacity_ur_los(2,:),'-','color',colours(2,:),'linewidth',linewidth);
% plot(snr,likely_capacity_ur_los(3,:),'-','color',colours(3,:),'linewidth',linewidth);
% plot(snr,likely_capacity_ur_los(4,:),'-','color',colours(4,:),'linewidth',linewidth);
% plot(snr,likely_capacity_ur_los(5,:),'-','color',colours(5,:),'linewidth',linewidth);
% 
% xlabel('SNR (dB)','fontname',fontname,'fontsize',fontsize);
% ylabel('95% likely sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
%     
% legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
%     
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% % dim = [0.65 0.75 0.2 0.15];
% dim = [0.2 0.25 0.2 0.15];
% 
% annotation('ellipse',dim,'linewidth',linewidth);
% 
% axes('position',[.6 .20 .27 .27])
% box on;
% 
% plot(snr(1:3),likely_capacity_ur_los(1,1:3),'-','color',colours(1,:),'linewidth',linewidth);
% hold on;
% plot(snr(1:3),likely_capacity_ur_los(2,1:3),'-','color',colours(2,:),'linewidth',linewidth);
% plot(snr(1:3),likely_capacity_ur_los(3,1:3),'-','color',colours(3,:),'linewidth',linewidth);
% plot(snr(1:3),likely_capacity_ur_los(4,1:3),'-','color',colours(4,:),'linewidth',linewidth);
% plot(snr(1:3),likely_capacity_ur_los(5,1:3),'-','color',colours(5,:),'linewidth',linewidth);
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% % ylim([45 58]);
% 
% if (savefig == 1)
%     saveas(gcf,[root_rate_lik 'uplink_ur_los_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'fig');
%     saveas(gcf,[root_rate_lik 'uplink_ur_los_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'png');
%     saveas(gcf,[root_rate_lik 'uplink_ur_los_sum_rate_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L)],'epsc2');
% end

for chn_idx = 1:N_CHN
    
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    plot(snr,likely_capacity(1,:,1,chn_idx),'-','color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(snr,likely_capacity(2,:,1,chn_idx),'-','color',colours(2,:),'linewidth',linewidth);
    plot(snr,likely_capacity(3,:,1,chn_idx),'-','color',colours(3,:),'linewidth',linewidth);
    plot(snr,likely_capacity(4,:,1,chn_idx),'-','color',colours(4,:),'linewidth',linewidth);
    plot(snr,likely_capacity(5,:,1,chn_idx),'-','color',colours(5,:),'linewidth',linewidth);
    plot(snr,likely_capacity(1,:,2,chn_idx),'--','color',colours(1,:),'linewidth',linewidth);
    plot(snr,likely_capacity(2,:,2,chn_idx),'--','color',colours(2,:),'linewidth',linewidth);
    plot(snr,likely_capacity(3,:,2,chn_idx),'--','color',colours(3,:),'linewidth',linewidth);
    plot(snr,likely_capacity(4,:,2,chn_idx),'--','color',colours(4,:),'linewidth',linewidth);
    plot(snr,likely_capacity(5,:,2,chn_idx),'--','color',colours(5,:),'linewidth',linewidth);

    
    
    xlabel('SNR (dB)','fontname',fontname,'fontsize',fontsize);
    ylabel('95% likely sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    
    legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location);
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    if (savefig == 1)
        saveas(gcf,[root_rate_lik 'uplink' channel_mod{chn_idx} '_avg_rate_ter_K_' num2str(K) '_L_' num2str(L)],'fig');
        saveas(gcf,[root_rate_lik 'uplink' channel_mod{chn_idx} '_avg_rate_ter_K_' num2str(K) '_L_' num2str(L)],'png');
        saveas(gcf,[root_rate_lik 'uplink' channel_mod{chn_idx} '_avg_rate_ter_K_' num2str(K) '_L_' num2str(L)],'epsc2');
    end
    
end

% 
% bar_sum = [sum(mean(rate_u,2)) sum(mean(rate_d,2));
%            sum(mean(rate_rs_u,2))  sum(mean(rate_rs_d,2));
%            sum(mean(rate_sos_u,2)) sum(mean(rate_sos_d,2)); 
%            sum(mean(rate_icibs_u,2)) sum(mean(rate_icibs_d,2))];
% 
% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% bar(cat,bar_sum,BAR_SIZE);
% 
% xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
% ylabel('Average sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_link,'fontname',fontname,'fontsize',fontsize);
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 40]);
% 
% saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'fig');
% saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'png');
% saveas(gcf,[root_erg_rate 'sum_rate_M_' num2str(M) '_K_' num2str(K)],'epsc2');
% 
% bar_sum_5 = [R_sum_u_ns R_sum_d_ns; 
%              R_sum_u_rs R_sum_d_rs; 
%              R_sum_u_sos R_sum_d_sos;
%              R_sum_u_icibs R_sum_d_icibs];
%          
% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% bar(cat,bar_sum_5,BAR_SIZE);
% 
% xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
% ylabel('95% likely sum-rate (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_link,'fontname',fontname,'fontsize',fontsize);
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 40]);
% 
% saveas(gcf,[root_erg_rate '95_sum_rate_M_' num2str(M) '_K_' num2str(K)],'fig');
% saveas(gcf,[root_erg_rate '95_sum_rate_M_' num2str(M) '_K_' num2str(K)],'png');
% saveas(gcf,[root_erg_rate '95_sum_rate_M_' num2str(M) '_K_' num2str(K)],'epsc2');
% 
% 
% bar_u = [mean(mean(rate_u,2)) mean(mean(rate_d,2));
%          mean(mean(rate_rs_u,2)) mean(mean(rate_rs_d,2));
%          mean(mean(rate_sos_u,2)) mean(mean(rate_sos_d,2));
%          mean(mean(rate_icibs_u,2)) mean(mean(rate_icibs_d,2))];
% 
% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% bar(cat,bar_u,BAR_SIZE);
% 
% xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
% ylabel('Average rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_link,'fontname',fontname,'fontsize',fontsize);
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 10]);
% 
% saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'fig');
% saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'png');
% saveas(gcf,[root_erg_rate 'avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'epsc2');
% 
% bar_u_5 = [R_u_ns R_d_ns;
%            R_u_rs R_d_rs; 
%            R_u_sos R_d_sos;
%            R_u_icibs R_d_icibs];
% 
% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% bar(cat,bar_u_5,BAR_SIZE);
% 
% xlabel('Algorithms','fontname',fontname,'fontsize',fontsize);
% ylabel('95% likely average rate per terminal (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_link,'fontname',fontname,'fontsize',fontsize);
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% ylim([0 10]);
% 
% saveas(gcf,[root_erg_rate '95_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'fig');
% saveas(gcf,[root_erg_rate '95_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'png');
% saveas(gcf,[root_erg_rate '95_avg_rate_ter_M_' num2str(M) '_K_' num2str(K)],'epsc2');