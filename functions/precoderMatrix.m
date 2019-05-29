function [prec_mtx,varargout] = precoderMatrix(chn_mtx,precoder)

precoder = upper(precoder);

n_antenna = size(chn_mtx,1);

switch precoder
    case 'MATCHED FILTER'
            chn_norm     = vecnorm(chn_mtx);
            chn_norm_mtx = repmat(chn_norm,n_antenna,1);
            
            chn_mtx_norm = chn_mtx./chn_norm_mtx;
            
            prec_mtx     = conj(chn_mtx_norm);
            varargout{1} = chn_mtx_norm;
    case 'MF'
            chn_norm     = vecnorm(chn_mtx);
            chn_norm_mtx = repmat(chn_norm,n_antenna,1);
            
            chn_mtx_norm = chn_mtx./chn_norm_mtx;
            
            prec_mtx     = conj(chn_mtx_norm);
            varargout{1} = chn_mtx_norm;
    otherwise
        error('Invalid precoder');
end

end