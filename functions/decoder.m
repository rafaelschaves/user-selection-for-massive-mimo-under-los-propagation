function [s_hat,varargout] = decoder(y,H,decoder_parameters)

% Initialization

decoding = upper(decoder_parameters.precoder);

M = size(H,1);                                                             % Number of antennas at base station
% K = size(H,2);                                                             % Number of users

switch decoding
    
    case 'MF'
        
        Q     = H'/M;                                                      % MF decoding matrix
        s_hat = Q*y;                                                       % MF decoded signal
        
        varargout{1} = Q;
        % varargout{2} = norm(s_hat,2)/M;
        
    % case 'MF-SIC'
    %    
    %     for k = 1:K
    %        
    %     end
        
    case 'ZF'
        
        Q = (H'*H)\H';                                                     % ZF decoding matrix
        s_hat = Q*y;                                                       % ZF decoded signal
        
        varargout{1} = Q;
        % varargout{2} = norm(s_hat,2)/M;
        
    case 'MMSE'
        
        Q = (H'*H + 1/snr*eye(M))\H';                                      % MMSE decoding matrix
        s_hat = Q*y;                                                       % MMSE decoded signal
        
        varargout{1} = Q;
        % varargout{2} = norm(s_hat,2)/M;
              
    otherwise
        
        error('Invalid precoder.');
        
end

end