function [rx_signal] = channel(tx_signal,pow_tx,chnl_mtx,snr,type)
                                                     
type = upper(type);

n_block   = size(tx_signal,2);

switch type
    case 'UPLINK'
        n_antenna = size(chnl_mtx,1);
        
        noise     = randn(n_antenna,n_block) + 1i*randn(n_antenna,n_block);% AWGN noise
        pow_noise = norm(noise(:),2)^2/(n_antenna*n_block);                % Noise's power
    case 'DOWNLINK'
        n_user = size(chnl_mtx,1);
        
        noise       = randn(n_user,n_block) + 1i*randn(n_user,n_block);    % AWGN noise
        pow_noise = norm(noise(:),2)^2/(n_user*n_block);                   % Noise's power
    otherwise
        error('Invalid type of channel');
end
                                   
pow_noise_bar = (pow_tx/pow_noise)/snr;                                    % Normalized noise's power in function of SNR
noise_bar     = sqrt(pow_noise_bar)*noise;                                 % Normalized AWGN noise

rx_signal = chnl_mtx*tx_signal + noise_bar;                                % Received signal

end