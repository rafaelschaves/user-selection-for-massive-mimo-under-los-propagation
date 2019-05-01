function [rate,varargout] = rateCalculation(chnl_mtx,proc_mtx,pow_vec,snr,type)

gamma = sinr(chnl_mtx,proc_mtx,pow_vec,snr,type); 
rate  = log2(1 + gamma);

varargout{1} = gamma;

end

function [gamma] = sinr(chnl_mtx,proc_mtx,pow_vec,snr,type)

type = upper(type);

switch type
    case 'UPLINK'
        aux_mtx        = abs(proc_mtx'*chnl_mtx).^2;
        pow_signal_vec = pow_vec.*diag(aux_mtx);
        pow_interf_vec = sum(pow_vec.*aux_mtx,2) - pow_signal_vec;
        proc_mtx_norm  = vecnorm(proc_mtx,2).^2;
        
        gamma = (snr.*pow_signal_vec)./(proc_mtx_norm' + snr.*pow_interf_vec);
    case 'DOWNLINK'
        aux_mtx        = abs(chnl_mtx.'*proc_mtx).^2;
        pow_signal_vec = pow_vec.*diag(aux_mtx);
        pow_interf_vec = sum(pow_vec.*aux_mtx,2) - pow_signal_vec;
        
        gamma = (snr.*pow_signal_vec)./(1 + snr.*pow_interf_vec);
    otherwise
        error('Invalid link');
end

end