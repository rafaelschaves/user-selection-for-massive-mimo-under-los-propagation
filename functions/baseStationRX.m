function [tx_bit_hat,varargout] = baseStationRX(rx_signal, ...
                                                chnl_mtx, ...
                                                decpar, ...
                                                n_bit)

n_user  = size(chnl_mtx,2);
n_block = size(rx_signal,2);

% Initialization

tx_bit_hat = zeros(n_bit*n_block,n_user);                                  % Estimated message in bits

[tx_signal_hat,decod_mtx] = decoder(rx_signal,chnl_mtx,decpar);            % Decoding received signal

% Signal decodification for each user

for k = 1:n_user
    tx_bit_hat(:,k) = qamdemod(tx_signal_hat(k,:).',2^n_bit,'OutputType','bit');
end

varargout{1} = tx_signal_hat;
varargout{2} = decod_mtx;

end

