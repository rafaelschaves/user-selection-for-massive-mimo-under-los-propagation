function [se_mr,varargout] = DLspectralEfficiency(H,beta,rho,eta,varargin)

N_ARGIN = 6;

if nargin == N_ARGIN-2
    tau_d = 1;
    tau_c = 1;
elseif nargin == N_ARGIN-1
    tau_d = varargin{1};
    tau_c = 200;
elseif nargin == N_ARGIN
    tau_d = varargin{1};
    tau_c = varargin{2};
end

% MR Processing

eta_mr  = eta(:,1);             % Power allocation vector for MR processing
H_mr    = H(:,:,1);
beta_mr = beta(:,1);

V_mr      = conj(H_mr);
v_mr_norm = vecnorm(V_mr);

W_mr = V_mr*diag(1./v_mr_norm);

[P_s,P_i] = computeSignalAndInterferencePower(H_mr,beta_mr,rho,eta_mr,W_mr);

gamma_mr = P_s./(P_i + 1);
se_mr    = (tau_d/tau_c)*log2(1 + gamma_mr);

% ZF Processing

if nargout > 1
    if size(H,3) == 1
        H_zf = H;
    elseif size(H,3) == 2
        H_zf = H(:,:,2);
    end
    
    if size(beta,2) == 1
        beta_zf = beta;
    elseif size(beta,2) == 2
        beta_zf = beta(:,2);
    end
    
    if size(eta,2) == 1
        eta_zf = eta;
    elseif size(eta,2) == 2
        eta_zf = eta(:,2);
    end
    
    V_zf      = conj(H_zf)/(H_zf.'*conj(H_zf));
    v_zf_norm = vecnorm(V_zf);
    
    W_zf = V_zf*diag(1./v_zf_norm);
    
    [P_s,P_i] = computeSignalAndInterferencePower(H_zf,beta_zf,rho,eta_zf,W_zf);
    
    gamma_zf = P_s./(P_i + 1);
    se_zf    = (tau_d/tau_c)*log2(1 + gamma_zf);
    
    varargout{1} = se_zf;
end

end

function [P_s,P_i] = computeSignalAndInterferencePower(H,beta,rho,eta,W)

R = (H.'*W);

P_s = rho*beta.*eta.*diag(abs(R).^2);                                      % Signal power
P_i = rho*beta.*(sum((abs(R).^2 - diag(diag(abs(R).^2)))*diag(eta),2));    % Interference signal

end
