function [psi] = ici(chnl_mtx)

n_antenna = size(chnl_mtx,1);
n_user    = size(chnl_mtx,2);

chnl_mtx_norm = zeros(n_antenna,n_user);

for k = 1:n_user
    norm_h_k = norm(chnl_mtx(:,k),2);
    
    if(norm_h_k == 0)
        chnl_mtx_norm(:,k) = 0;
    else
        chnl_mtx_norm(:,k) = chnl_mtx(:,k)/norm_h_k;
    end
end

interf_mtx_norm = chnl_mtx_norm'*chnl_mtx_norm;

psi = (sum(abs(interf_mtx_norm),2) - 1)/(n_user - 1);

psi(psi <= 0) = 0;



% chnl_vec_norm = vecnorm(chnl_mtx,2);
% Chnl_vec_norm = repmat(chnl_vec_norm,n_antenna,1);

% chnl_mtx_norm = chnl_mtx./Chnl_vec_norm;

end