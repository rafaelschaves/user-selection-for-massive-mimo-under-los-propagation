function [signal_hat,varargout] = decoder(rx_signal,chnl_mtx,decpar)

% Initialization

decoder   = upper(decpar.decoder);
ref_power = decpar.power;

n_antenna = size(chnl_mtx,1);                                              % Number of antennas at base station
n_user    = size(chnl_mtx,2);                                              % Number of users in the cell
n_block   = size(rx_signal,2);                                             % Number of transmitted blocks

switch decoder
    
    case 'MF'
        
        decod_mtx  = chnl_mtx'/n_antenna;                                  % MF decoding matrix
        signal_hat = decod_mtx*rx_signal;                                  % MF decoded signal
        
        power_signal_hat = norm(signal_hat(:),2)^2/(n_user*n_block);       % Received signal power calculation
        
        signal_hat = sqrt(ref_power/power_signal_hat)*signal_hat;          % Received signal power normalization
        
        varargout{1} = decod_mtx;
        varargout{2} = power_signal_hat;
        
    case 'ZF'
        
        decod_mtx  = (chnl_mtx'*chnl_mtx)\chnl_mtx';                       % ZF decoding matrix
        signal_hat = decod_mtx*rx_signal;                                  % ZF decoded signal
        
        power_signal_hat = norm(signal_hat(:),2)^2/(n_user*n_block);       % Received signal power calculation
        
        signal_hat = sqrt(ref_power/power_signal_hat)*signal_hat;          % Received signal power normalization
        
        varargout{1} = decod_mtx;
        varargout{2} = power_signal_hat;
        
    case 'MMSE'
        
        I_M = eye(n_antenna);
        snr = decpar.snr;
        
        decod_mtx  = (chnl_mtx'*chnl_mtx + I_M/snr)\chnl_mtx';             % MMSE decoding matrix
        signal_hat = decod_mtx*rx_signal;                                  % MMSE decoded signal
        
        power_signal_hat = norm(signal_hat(:),2)^2/(n_user*n_block);       % Received signal power calculation
        
        signal_hat = sqrt(ref_power/power_signal_hat)*signal_hat;          % Received signal power normalization

        varargout{1} = decod_mtx;
        varargout{2} = power_signal_hat;
              
    otherwise
        
        error('Invalid precoder.');
        
end

end