function [rx_signal,chann_mtx,large_scale_vec] = channel(tx_signal, ...
                                                         power_tx, ...
                                                         snr, ...
                                                         commcell, ...
                                                         fading_type)
                                                     
n_block   = size(tx_signal,2);
n_antenna = commcell.nAntennas;
                                 
[chann_mtx, large_scale_vec] = massiveMIMOChannel(commcell,fading_type);

chann_mtx = chann_mtx*sqrt(diag(1./large_scale_vec));
  
noise       = randn(n_antenna,n_block) + 1i*randn(n_antenna,n_block);      % AWGN noise
power_noise = norm(noise(:),2)^2/(n_antenna*n_block);                      % Noise's power

power_noise_bar = (power_tx/power_noise)/snr;                              % Normalized noise's power in function of SNR
noise_bar       = sqrt(power_noise_bar)*noise;                             % Normalized AWGN noise

rx_signal = chann_mtx*tx_signal + noise_bar;                               % Received signal

end

