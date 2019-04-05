function [demapped_signal] = demapper(input_signal,const_size)

% Initialization

n_user  = size(input_signal,1);
n_block = size(input_signal,2);
n_bit   = log2(const_size); 

demapped_signal = zeros(n_block*n_bit,n_user);

for k = 1:n_user
    demapped_signal(:,k) = qamdemod(input_signal(k,:).', ...
                                    const_size, ...
                                    'OutputType','bit');
end

end

