function [gamma,eta,varargout] = maxMinFairness(channel_mtx,large_scale,snr)

n_antennas = size(channel_mtx,1);
n_users    = size(channel_mtx,2);

I_k = eye(n_users);

tolerance = 1e-6;

channel_norm = vecnorm(channel_mtx);
Channel_norm = repmat(channel_norm,n_antennas,1);

channel_mtx_norm = channel_mtx./Channel_norm;

Lambda = abs(channel_mtx_norm'*channel_mtx_norm).^2 - I_k;
b      = 1./(large_scale.*(channel_norm.^2)');

gamma_l = 0;

if Lambda == 0
    gamma_u = 1000;
else
    gamma_u = 1/max(eig(Lambda));
end

% gamma_u = n_antennas*snr*max(large_scale);

solve = false;

while solve == false
    gamma = (gamma_l + gamma_u)/2;        
    
%    if n_users == 2
        eta = (I_k - gamma*Lambda)\(gamma*b/snr);
%     else
%         cvx_begin quiet
%             variable eta(n_users) nonnegative;
%             subject to
%                 (I_k - gamma*Lambda)*eta == gamma*b./snr;
%         cvx_end
%     end
    
    if sum(isnan(eta) ~= 0)
        gamma_u = gamma_u/2;
    elseif (gamma_u - gamma_l)/2 < tolerance
        solve = true;
    elseif norm(eta,1) <= 1
        gamma_l = gamma;
    elseif norm(eta,1) > 1
        gamma_u = gamma;
    end 
end

varargout{1} = Lambda;
varargout{2} = b;
end