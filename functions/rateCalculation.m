function [rate,varargout] = rateCalculation(chnl_mtx,snr,type)

gamma = sinr(chnl_mtx,snr,type); 
rate  = log2(1 + gamma);

varargout{1} = gamma;

end

function [gamma] = sinr(chnl_mtx,snr,type)

type = upper(type);

n_antenna = size(chnl_mtx,1);

switch type
    case 'UPLINK'
        chnl_vec_norm = vecnorm(chnl_mtx,2)';
        Chnl_vec_norm = repmat(chnl_vec_norm',n_antenna,1);

        chnl_mtx_norm = chnl_mtx./Chnl_vec_norm;

        interf_mtx = chnl_mtx_norm'*chnl_mtx;

        pow_interf = sum(abs(interf_mtx).^2,2) - chnl_vec_norm.^2;

        gamma = chnl_vec_norm.^2./(1./snr + pow_interf);
    case 'DOWNLINK'
        chnl_vec_norm = vecnorm(chnl_mtx,2)';
        Chnl_vec_norm = repmat(chnl_vec_norm',n_antenna,1);
        
        chnl_mtx_norm = chnl_mtx./Chnl_vec_norm;

        interf_mtx = chnl_mtx.'*conj(chnl_mtx_norm);
        
        pow_interf = sum(abs(interf_mtx).^2,2) - chnl_vec_norm.^2;

        gamma = chnl_vec_norm.^2./(1./snr + pow_interf);
    otherwise
        error('Invalid link');
end

end