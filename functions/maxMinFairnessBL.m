function [gamma,eta,varargout] = maxMinFairnessBL(channel_mtx,large_scale,snr,user_index)

n_antennas = size(channel_mtx,1);
n_users    = size(channel_mtx,2);
n_blocks   = size(user_index,1);

I_k = eye(n_users);

tolerance = 1e-6;

n_iterations = 1;

channel_norm = vecnorm(channel_mtx,2);
Channel_norm = repmat(channel_norm,n_antennas,1);

channel_mtx_norm = channel_mtx./Channel_norm;

R = abs(channel_mtx_norm'*channel_mtx_norm).^2 - I_k;
b = 1./(snr*large_scale.*(channel_norm.^2)');

R_bar = zeros(n_users,n_users);
N     = zeros(n_blocks,n_users);

for n = 1:n_blocks
    R_bar(user_index{n},user_index{n}) = R(user_index{n},user_index{n});
    N(n,user_index{n}) = 1;
end

if R_bar == 0
    gamma = 1/b;
    eta   = 1;
else
    gamma_l = 0;
    gamma_u = 1/max(eig(R_bar));
    % gamma_u = n_antennas*snr*max(large_scale);

    varargout{1} = gamma_u;
    
    
    while (gamma_u - gamma_l) > tolerance
        gamma = (gamma_l + gamma_u)/2;
        
        cvx_begin quiet
            variable eta(n_users) nonnegative;
            subject to
            (I_k/gamma - R_bar)*eta == b;
            N*eta <= ones(n_blocks,1);
        cvx_end
        
        if strcmp(cvx_status,'Solved')
            gamma_l = gamma;
        else
            gamma_u = gamma;
        end
        n_iterations = n_iterations + 1;
    end
end

varargout{2} = n_iterations;
varargout{3} = R;
varargout{4} = b;
end