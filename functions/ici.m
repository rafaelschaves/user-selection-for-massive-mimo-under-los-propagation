function [psi] = ici(chnl_mtx)

n_antenna  = size(chnl_mtx,1);
n_user     = size(chnl_mtx,2);
n_user_aux = size(chnl_mtx,2) - length(find(all(chnl_mtx == 0)));

chnl_mtx_norm = zeros(n_antenna,n_user);

norm_h_k = vecnorm(chnl_mtx,2);

for k = 1:n_user    
    if(norm_h_k(k) ~= 0)
        chnl_mtx_norm(:,k) = chnl_mtx(:,k)/norm_h_k;
    end
end

interf_mtx_norm = chnl_mtx_norm'*chnl_mtx_norm;

psi = (sum(abs(interf_mtx_norm),2) - 1)/(n_user_aux - 1);

end