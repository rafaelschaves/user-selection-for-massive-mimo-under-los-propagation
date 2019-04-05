function [tx_bit_hat,varargout] = baseStationRX(rx_signal, ...
                                                chnl_mtx, ...
                                                decpar, ...
                                                n_bit)
                                            

const_size = 2^n_bit;                                                      % Constellation size

% Receiver processing

[tx_signal_hat,decod_mtx] = decoder(rx_signal,chnl_mtx,decpar);            % Decoding received signal

tx_bit_hat = demapper(tx_signal_hat,const_size);                           % Signal demapping for each user

varargout{1} = tx_signal_hat;
varargout{2} = decod_mtx;

end

