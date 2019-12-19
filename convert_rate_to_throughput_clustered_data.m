clear;
close all;
clc;

% Clustered folder

root_uplink   = '../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Clustered/Uplink/';
root_downlink = '../../../Google Drive/UFRJ/PhD/Codes/user-scheduling-massive-mimo/Results/Clustered/Downlink/';

L = 13;

R = 0.2;                                                                   % Radius in km
r = sqrt(3)/2*R;

coherence_time = 15000;
% pilot_time     = [18 36 72];
ratio          = 0.5;
bandwidth      = 20e6;
% cell_area      = 3*sqrt(3)*R^2/2;

for m = [64 256]
    for k = [18 36 72]
        for theta_mid  = [0 45 90]
            for theta_step = [1 10 20 30 90]
                % cell_area = 1/2*r^2*theta_step*pi/180;
                cell_area = 1;
                
                load([root_uplink 'rate_mf_ur_los_M_' num2str(m) '_K_' ...
                      num2str(k) '_L_' num2str(L) '_theta_mid_' ...
                      num2str(theta_mid) '_theta_step_' ...
                      num2str(theta_step) '_SNR_10_dB_MC_10000.mat']);
                
                thrput_u     = bandwidth*ratio*(1 - k/coherence_time)*rate_u/cell_area;
                thrput_u_sel = bandwidth*ratio*(1 - k/coherence_time)*rate_u_sel/cell_area;
                
                load([root_downlink 'rate_mf_ur_los_M_' num2str(m) '_K_' ...
                      num2str(k) '_L_' num2str(L) '_theta_mid_' ...
                      num2str(theta_mid) '_theta_step_' ...
                      num2str(theta_step) '_SNR_10_dB_MC_10000.mat']);
              
                thrput_d     = bandwidth*ratio*(1 - k/coherence_time)*rate_d/cell_area;
                thrput_d_sel = bandwidth*ratio*(1 - k/coherence_time)*rate_d_sel/cell_area;
                
                save([root_uplink 'throughput_outdoors_pedestrian_mf_ur_los_M_' ...
                      num2str(m) '_K_' num2str(k) '_L_' num2str(L) ...
                      '_theta_mid_' num2str(theta_mid) '_theta_step_' ...
                      num2str(theta_step) '_SNR_10_dB_MC_10000.mat'], ...
                      'thrput_u','psi','thrput_u_sel','psi_sel');

                save([root_downlink 'throughput_outdoors_pedestrian_mf_ur_los_M_' ...
                      num2str(m) '_K_' num2str(k) '_L_' num2str(L) ...
                      '_theta_mid_' num2str(theta_mid) '_theta_step_' ...
                      num2str(theta_step) '_SNR_10_dB_MC_10000.mat'], ...
                      'thrput_d','psi','thrput_d_sel','psi_sel');
            end
        end
    end
end