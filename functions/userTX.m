function [tx_signal,pow_tx_signal,varargout] = userTX(n_user, ...
                                                      n_block, ...
                                                      n_bit, ...
                                                      varargin)
                                                  
% MACROS

N_ARGIN  = 4;                                                              % Number of input arguments
N_ARGOUT = 3;                                                              % Number of output arguments

% Initialization

tx_signal = zeros(n_user,n_block);                                         % Transmitted modulated message

const_size = 2^n_bit;                                                      % Constellation size

switch nargin
    case N_ARGIN-1
        tx_bit = zeros(n_bit*n_block,n_user);                              % Transmitted message in bit

        for k = 1:n_user
            tx_bit(:,k)    = randi([0 1],n_bit*n_block,1);                 % Message in bits of the kth user
            tx_signal(k,:) = qammod(tx_bit(:,k), ...                       % Genarating message to be transmitted
                                    const_size, ...
                                    'InputType','bit').';                  % Message modulated in 2^B-QAM for the kth user
        end
    case N_ARGIN
        tx_bit = varargin{1};                                              % Transmitted message in bit
        
        for k = 1:n_user
            tx_signal(k,:) = qammod(tx_bit(:,k), ...                       % Genarating message to be transmitted
                                    const_size, ...
                                    'InputType','bit').';                  % Message modulated in 2^B-QAM for the kth user
        end
    otherwise
end

pow_tx_signal = norm(tx_signal(:),2)^2/(n_user*n_block);                   % Transmitted signal's power

varargout{1} = tx_bit;

end