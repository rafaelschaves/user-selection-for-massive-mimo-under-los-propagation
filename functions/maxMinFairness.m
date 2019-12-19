function [gamma,eta,varargout] = maxMinFairness(channel_mtx,large_scale,snr,algorithm)

algorithm = upper(algorithm);

n_antennas = size(channel_mtx,1);
n_users    = size(channel_mtx,2);

I_k = eye(n_users);

tolerance = 1e-6;

n_iterations = 1;

channel_norm = vecnorm(channel_mtx,2);
Channel_norm = repmat(channel_norm,n_antennas,1);

channel_mtx_norm = channel_mtx./Channel_norm;

R = abs(channel_mtx_norm'*channel_mtx_norm).^2 - I_k;
b = 1./(snr*large_scale.*(channel_norm.^2)');

if R == 0
    gamma = 1/b;
    eta   = 1;
else
    switch algorithm
        case 'ALGORITHM 1'
            gamma_l = 0;
            gamma_u = 1/max(eig(R));
            
            varargout{1} = gamma_u;
            
            while (gamma_u - gamma_l) > tolerance
                gamma = (gamma_l + gamma_u)/2;
                eta   = (I_k/gamma - R)\b;
                
                if norm(eta,1) <= 1
                    gamma_l = gamma;
                else
                    gamma_u = gamma;
                end
                
                n_iterations = n_iterations + 1;
            end
            
        case 'ALGORITHM 2'
            gamma = 1/max(sum(R,2));
            eta   = (I_k/gamma - R)\b;
            
            if norm(eta,1) <= 1
                gamma_l = gamma;
                gamma_u = 1/max(eig(R));
            else
                gamma_l = 0;
                gamma_u = gamma;
            end
            
            varargout{1} = gamma_u;
            
            while (gamma_u - gamma_l) > tolerance
                gamma = (gamma_l + gamma_u)/2;
                eta   = (I_k/gamma - R)\b;
                
                if norm(eta,1) <= 1
                    gamma_l = gamma;
                else
                    gamma_u = gamma;
                end
                
                n_iterations = n_iterations + 1;
            end
                        
        case 'ALGORITHM 3'
            gamma_l = 0;
            gamma_u = n_antennas*snr*max(large_scale);
            
            varargout{1} = gamma_u;
            
            while (gamma_u - gamma_l) > tolerance
                gamma = (gamma_l + gamma_u)/2;
                
                cvx_begin quiet
                    variable eta(n_users) nonnegative;
                    subject to
                        (I_k/gamma - R)*eta == b;
                cvx_end
                
                if norm(eta,1) <= 1
                    gamma_l = gamma;
                else
                    gamma_u = gamma;
                end
                
                n_iterations = n_iterations + 1;
            end
            
        case 'ALGORITHM 4'
            gamma_l = 0;
            gamma_u = n_antennas*snr*max(large_scale);
            
            varargout{1} = gamma_u;
            
            while (gamma_u - gamma_l) > tolerance
                gamma = (gamma_l + gamma_u)/2;
                eta   = (I_k/gamma - R)\b;
                
                if norm(eta,1) <= 1
                    gamma_l = gamma;
                else
                    gamma_u = gamma;
                end
                
                n_iterations = n_iterations + 1;
            end
            
        otherwise
    end
end

varargout{2} = n_iterations;
varargout{3} = R;
varargout{4} = b;
end