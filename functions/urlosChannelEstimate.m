function [G_hat] = urlosChannelEstimate(commcell,theta,max_err)

c = 3e8;                                                                   % Light speed

n_antenna = commcell.nAntennas;                                            % Number of transmit antennas at base station
n_user    = commcell.nUsers;                                               % Number of user terminals
f_c       = commcell.frequency;                                            % Carrier frequency of the transmitted signal

antenna_spacing = c/(2*f_c);                                               % Antenna spacing of transmitt array

theta_err = -max_err + 2*max_err*rand(n_user,1);
theta_user = theta + theta_err;

A = steeringVector(n_antenna, theta_user, antenna_spacing, c/f_c);
                
G_hat = A;                                                                 % Estimate of the small-scale fading coefficient matrix
end

function steering_vector = steeringVector(M,theta,d,lambda)

m   = (0:M-1)';
tau = d*sin(theta)/lambda;
Tau = repmat(tau',M,1);

steering_vector = exp(-2*pi*1i*Tau.*repmat(m,1,size(theta,1)));
end