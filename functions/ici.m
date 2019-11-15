function [psi] = ici(chnl_mtx)

n_antenna  = size(chnl_mtx,1);
n_user     = size(chnl_mtx,2);
n_user_aux = n_user - length(find(all(chnl_mtx == 0)));

norm_h_k = vecnorm(chnl_mtx,2);
idx_aux  = find(norm_h_k ~= 0);

chnl_mtx_norm            = zeros(n_antenna,n_user);
chnl_mtx_norm(:,idx_aux) = chnl_mtx(:,idx_aux)./norm_h_k(idx_aux);

psi = (sum(abs(chnl_mtx_norm'*chnl_mtx_norm),2) - 1)/(n_user_aux - 1);

end