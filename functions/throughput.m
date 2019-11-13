function [thrput,varargout] = throughput(chnl_mtx,proc_mtx,pow_vec,link,snr,settings)

tau_c     = settings.coherenceTime;                                        % Coherence time in samples
tau_p     = settings.PilotTime;                                            % Pilot time in samples
ratio     = settings.uplinkDownlinkTimeRatio;                              % Ratio between the uplink and downlink payload time
bandwidth = settings.bandwidth;                                            % Sytem bandwidth in Hz
cell_area = settings.cellArea;                                             % Cell area in km^2

gamma  = sinr(chnl_mtx,proc_mtx,pow_vec,snr,link);                         % SINR
se     = ratio*(1 - tau_p/tau_c)*log2(1 + gamma);                          % Spectral efficiency in b/s/Hz/cell
thrput = bandwidth*se/cell_area;                                           % Throughput in b/s/km^2

varargout{1} = se;
varargout{2} = gamma;

end

function [gamma] = sinr(chnl_mtx,proc_mtx,pow_vec,snr,link)

link = upper(link);

switch link
    case 'UPLINK'
        aux_mtx        = abs(proc_mtx'*chnl_mtx).^2;
        pow_signal_vec = pow_vec.*diag(aux_mtx);
        pow_interf_vec = sum(snr.*pow_vec.*aux_mtx,2) - snr.*pow_signal_vec;
        proc_mtx_norm  = vecnorm(proc_mtx,2).^2;
        
        gamma = (snr.*pow_signal_vec)./(proc_mtx_norm' + pow_interf_vec);
    case 'DOWNLINK'
        aux_mtx        = abs(chnl_mtx.'*proc_mtx).^2;
        pow_signal_vec = pow_vec.*diag(aux_mtx);
        pow_interf_vec = sum(pow_vec.*aux_mtx,2) - pow_signal_vec;
        
        gamma = (snr.*pow_signal_vec)./(1 + snr.*pow_interf_vec);
    otherwise
        error('Invalid link');
end

end