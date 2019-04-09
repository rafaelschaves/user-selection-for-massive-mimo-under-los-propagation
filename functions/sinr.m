function [gamma] = sinr(chnl_mtx,snr)

n_antenna = size(chnl_mtx,1);

chnl_vec_norm = vecnorm(chnl_mtx,2)';
Chnl_vec_norm = repmat(chnl_vec_norm',n_antenna,1);

chnl_mtx_norm = chnl_mtx./Chnl_vec_norm;

interf_mtx = chnl_mtx_norm'*chnl_mtx;

pow_interf = sum(abs(interf_mtx).^2,2) - chnl_vec_norm.^2;

gamma = chnl_vec_norm.^2./(1./snr + pow_interf);

end
