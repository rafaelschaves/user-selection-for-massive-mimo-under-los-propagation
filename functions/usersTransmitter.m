function [message_sym,power_sym,varargout] = usersTransmitter(n_user, ...
                                                              n_block, ...
                                                              n_bit)

constellation_size = 2^n_bit;                                              % Constellation size

message_bit = zeros(n_bit*n_block,n_user);                                 % Message in bits
message_sym = zeros(n_user,n_block);                                       % Message modulated in 2^B-QAM

for k = 1:n_user
    message_bit(:,k) = randi([0 1],n_bit*n_block,1);                       % Message in bits of the kth user
    message_sym(k,:) = qammod(message_bit(:,k),constellation_size, ...
                              'InputType','bit').';                        % Message modulated in 2^B-QAM for the kth user
end
    
power_sym = norm(message_sym(:),2)^2/(n_user*n_block);                     % Transmitted signal's power

varargout{1} = message_bit;

end

