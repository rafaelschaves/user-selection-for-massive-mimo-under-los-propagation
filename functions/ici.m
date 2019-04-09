function [psi] = ici(chnl_mtx)

n_antenna = size(chnl_mtx,1);

chnl_vec_norm = vecnorm(chnl_mtx,2);
Chnl_vec_norm = repmat(chnl_vec_norm,n_antenna,1);

chnl_mtx_norm = chnl_mtx./Chnl_vec_norm;

interf_mtx_norm = chnl_mtx_norm'*chnl_mtx_norm;

psi = sum(abs(interf_mtx_norm),2) - 1;

end
