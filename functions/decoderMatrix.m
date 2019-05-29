function [dec_mtx,varargout] = decoderMatrix(chn_mtx,decoder)

decoder = upper(decoder);

n_antenna = size(chn_mtx,1);

switch decoder
    case 'MATCHED FILTER'
            chn_norm     = vecnorm(chn_mtx);
            chn_norm_mtx = repmat(chn_norm,n_antenna,1);
            
            chn_mtx_norm = chn_mtx./chn_norm_mtx;
            
            dec_mtx      = chn_mtx_norm;
            varargout{1} = conj(chn_mtx_norm);
    case 'MF'
            chn_norm     = vecnorm(chn_mtx);
            chn_norm_mtx = repmat(chn_norm,n_antenna,1);
            
            chn_mtx_norm = chn_mtx./chn_norm_mtx;
            
            dec_mtx      = chn_mtx_norm;
            varargout{1} = conj(chn_mtx_norm);
    otherwise
        error('Invalid precoder');
end

end