function [tx_signal,pow_tx_signal,varargout] = userTX(n_user,n_block,n_bit)

% Initialization

tx_bit    = zeros(n_bit*n_block,n_user);                                   % Message in bits
tx_signal = zeros(n_user,n_block);                                         % Message modulated in 2^B-QAM

const_size = 2^n_bit;                                                      % Constellation size

for k = 1:n_user
    tx_bit(:,k) = randi([0 1],n_bit*n_block,1);                            % Message in bits of the kth user
    tx_signal(k,:) = qammod(tx_bit(:,k),const_size,'InputType','bit').';   % Message modulated in 2^B-QAM for the kth user
end
    
pow_tx_signal = norm(tx_signal(:),2)^2/(n_user*n_block);                   % Transmitted signal's power

varargout{1} = tx_bit;

end