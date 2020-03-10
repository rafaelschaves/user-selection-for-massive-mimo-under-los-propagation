function [H_s,varargout] = userSelector(H,beta,rho,alg,type,varargin)

N_ARGIN = 7;                                                  

alg = upper(alg);                                                          % Algorithm
type = upper(type);                                                        % Type of selection

if(nargin == N_ARGIN-2)
    L    = 5;
    tau  = 0.5;
elseif(nargin == N_ARGIN-1)
    L    = varargin{1};
    tau  = 0.5;
elseif(nargin == N_ARGIN)
    L    = varargin{1};
    tau  = varargin{2};
else
    error('Invalid number of input arguments');
end

switch alg
    case 'EXHAUSTIVE SEARCH SELECTION EP'
        [H_s,varargout{1}] = exhaustiveSearchSelectionEP(H,beta,rho,L);
    case 'EXHAUSTIVE SEARCH SELECTION MMF'
        [H_s,varargout{1}] = exhaustiveSearchSelectionMMF(H,beta,rho,L);
    case 'SEMI-ORTHOGONAL SELECTION'
        [H_s,varargout{1}] = semiOrthogonalSelection(H,beta,L);
    case 'CORRELATION-BASED SELECTION'
        [H_s,varargout{1}] = correlationBasedSelection(H,L,tau,type);
    case 'ICI-BASED SELECTION'
        [H_s,varargout{1}] = iciBasesSelection(H,L,tau,type);
    case 'FR-BASED SELECTION'
        [H_s,varargout{1}] = frBasedSelection(H,beta,L,tau,type);
    otherwise
        error('Invalid algorithm');
end

end

function [corr_mtx] = correlationMatrix(mtx)

n_row = size(mtx,1);
n_col = size(mtx,2);

norm_mtx_k = vecnorm(mtx,2);
idx_aux    = find(norm_mtx_k ~= 0);

mtx_norm            = zeros(n_row,n_col);
mtx_norm(:,idx_aux) = mtx(:,idx_aux)./norm_mtx_k(idx_aux);

corr_mtx = abs(mtx_norm'*mtx_norm);
              
end

function [H_s,varargout] = exhaustiveSearchSelectionEP(H,beta,rho,L)

K = size(H,2); % Number users

if size(beta,1) == 1
    beta = repmat(beta,K,1);
end

S_set_aux = nchoosek(1:K,L);
N_S_set   = size(S_set_aux,1);

se_mr = zeros(L,N_S_set);
se_zf = zeros(L,N_S_set);

for n = 1:N_S_set    
    [se_mr(:,n),se_zf(:,n)] = DLspectralEfficiency(H(:,S_set_aux(n,:)),beta(S_set_aux(n,:)),rho,1/L);
end

[~,idx_S_set_mr] = max(sum(se_mr,1));
[~,idx_S_set_zf] = max(sum(se_zf,1));

S_set_mr = S_set_aux(idx_S_set_mr,:);
S_set_zf = S_set_aux(idx_S_set_zf,:);

H_s_mr = H(:,S_set_aux(idx_S_set_mr,:));
H_s_zf = H(:,S_set_aux(idx_S_set_zf,:));

H_s(:,:,1) = H_s_mr;
H_s(:,:,2) = H_s_zf;

S_set(:,1) = S_set_mr';
S_set(:,2) = S_set_zf';

varargout{1} = S_set;

end

function [H_s,varargout] = exhaustiveSearchSelectionMMF(H,beta,rho,L)

K = size(H,2);                                              % Number users

if size(beta,1) == 1
    beta = repmat(beta,K,1);
end

S_set_aux = nchoosek(1:K,L);
N_S_set   = size(S_set_aux,1);

se_mr = zeros(L,N_S_set);
se_zf = zeros(L,N_S_set);

for n = 1:N_S_set
    [~,eta_mr] = maxMinFairness(H(:,S_set_aux(n,:)),beta(S_set_aux(n,:)),rho,'algorithm 2');
    eta_zf = (1/sum(1./beta(S_set_aux(n,:))))./beta(S_set_aux(n,:));
    
    eta = [eta_mr eta_zf];
    
    [se_mr(:,n),se_zf(:,n)] = DLspectralEfficiency(H(:,S_set_aux(n,:)),beta(S_set_aux(n,:)),rho,eta);
end

[~,idx_S_set_mr]  = max(min(se_mr,[],1));
[~,idx_S_set_zf]  = max(min(se_zf,[],1));

S_set_mr  = S_set_aux(idx_S_set_mr,:);
S_set_zf  = S_set_aux(idx_S_set_zf,:);

H_s_mr = H(:,S_set_aux(idx_S_set_mr,:));
H_s_zf = H(:,S_set_aux(idx_S_set_zf,:));

H_s(:,:,1) = H_s_mr;
H_s(:,:,2) = H_s_zf;

S_set(:,1) = S_set_mr';
S_set(:,2) = S_set_zf';

varargout{1} = S_set;

end

function [H_s,varargout] = semiOrthogonalSelection(H,beta,L)

M = size(H,1);                                              % Number of antennas at base station

S_set     = zeros(L,1);
g_ort_sel = zeros(M,L);

I_M = eye(M);

G_aux = H*diag(sqrt(beta));

for l = 1:L
    if(l == 1)
        g_ort = G_aux;
    else
        g_ort = (I_M - g_ort_sel*g_ort_sel')*G_aux;
    end
    
    [~,S_set(l)] = max(vecnorm(g_ort));
    
    g_ort_sel(:,l) = g_ort(:,S_set(l))/norm(g_ort(:,S_set(l)),2);
    
    G_aux(:,S_set(l)) = zeros(M,1);
end

S_set = sort(S_set);

H_s = H(:,S_set);

varargout{1} = S_set;

end

function [H_s,varargout] = correlationBasedSelection(H,L,tau,type)

M = size(H,1);                                              % Number of antennas at base station
K = size(H,2);                                              % Number users

I_K = eye(K);

H_aux = H;

S_set = (1:K)';

switch type
    case 'FIXED'
        S_set_drop = zeros(K-L,1);
        
        for i = 1:K-L
            R     = correlationMatrix(H_aux);
            R_aux = R - I_K;
            
            [~,idx_corr] = max(R_aux(:));
            [corr_i,corr_j] = ind2sub(size(R_aux),idx_corr);
            
            h_i = R_aux(:,corr_i);
            h_j = R_aux(:,corr_j);
            
            h_i(corr_j) = [];
            h_j(corr_i) = [];
            
            if(max(h_i) > max(h_j))
                S_set_drop(i) = corr_i;
            else
                S_set_drop(i) = corr_j;
            end
            
            H_aux(:,S_set_drop(i)) = zeros(M,1);
        end
        
        S_set(S_set_drop) = [];
        
        H_s = H(:,S_set);
                
        varargout{1} = S_set;
    case 'AUTOMATIC'
        selection = 1;
        idx_while = 1;
        
        if(tau == 1)
            idx_while = 1;
        else
            while(selection == 1)
                R = correlationMatrix(H_aux);
                R_aux = R - I_K;
                
                [max_corr,idx_corr] = max(R_aux(:));
                
                if(max_corr < tau || idx_while == K)
                    selection = 0;
                else
                    [corr_i,corr_j] = ind2sub(size(R_aux),idx_corr);
                    
                    h_i = R_aux(:,corr_i);
                    h_j = R_aux(:,corr_j);
                    
                    h_i(corr_j) = [];
                    h_j(corr_i) = [];
                    
                    if(max(h_i) > max(h_j))
                        S_set_drop(idx_while) = corr_i;
                    else
                        S_set_drop(idx_while) = corr_j;
                    end
                    
                    H_aux(:,S_set_drop(idx_while)) = zeros(M,1);
                    
                    idx_while = idx_while + 1;
                end
            end
        end
        
        L = K - idx_while + 1;
        
        if(L == K)
            H_s = H;
        else
            S_set(S_set_drop) = [];
            
            H_s = H(:,S_set);
        end
                
        varargout{1} = S_set;
    otherwise
        error('Invalid type of selection');
end

end

function [H_s,varargout] = iciBasesSelection(H,L,tau,type)

M = size(H,1);                                              % Number of antennas at base station
K = size(H,2);                                              % Number users

H_aux = H;

S_set = (1:K)';

switch type
    case 'FIXED'
        S_set_drop = zeros(K-L,1);
                
        for i = 1:K-L
            psi = ici(H_aux);
            
            [~,S_set_drop(i)] = max(psi);
            
            H_aux(:,S_set_drop(i)) = zeros(M,1);
        end
        
        S_set(S_set_drop) = [];
        
        H_s = H(:,S_set);
                
        varargout{1} = S_set;
    case 'AUTOMATIC'
        selection = 1;
        idx_while = 1;
        
        if(tau == 1)
            idx_while = 1;
        else
            while(selection == 1)
                psi = ici(H_aux);
                
                decision = (psi > tau);
                
                if(sum(decision) == 0 || idx_while == K)
                    selection = 0;
                else
                    idx_decision = S_set(decision);
                    
                    [~,idx_aux] = max(psi(decision));
                    
                    S_set_drop(idx_while) = idx_decision(idx_aux);
                    
                    H_aux(:,S_set_drop(idx_while)) = zeros(M,1);
                    
                    idx_while = idx_while + 1;
                end
            end
        end
        
        L = K - idx_while + 1;
        
        if(L == K)
            H_s = H;
        else
            S_set(S_set_drop) = [];
            
            H_s = H(:,S_set);
        end
        
        varargout{1} = S_set;
    otherwise
        error('Invalid type of selection');
end

end

function [H_s,varargout] = frBasedSelection(H,beta,L,tau,type)

M = size(H,1);                                              % Number of antennas at base station
K = size(H,2);                                              % Number users

H_aux = H;

S_set = (1:K)';

switch type
    case 'FIXED'
        S_set_drop = zeros(K-L,1);
                
        for i = 1:K-L
            psi = ici(H_aux);
            
            [~,S_set_drop(i)] = max(psi./beta);
            
            H_aux(:,S_set_drop(i)) = zeros(M,1);
        end
        
        S_set(S_set_drop) = [];
        
        H_s = H(:,S_set);
           
        varargout{1} = S_set;
    case 'AUTOMATIC'
        selection = 1;
        idx_while = 1;
        
        if(tau == 1)
            idx_while = 1;
        else
            while(selection == 1)
                psi = ici(H_aux);
                
                decision = (psi > tau);
                
                if(sum(decision) == 0 || idx_while == K)
                    selection = 0;
                else
                    idx_decision = S_set(decision);
                    
                    [~,idx_aux] = max(psi(decision));
                    
                    S_set_drop(idx_while) = idx_decision(idx_aux);
                    
                    H_aux(:,S_set_drop(idx_while)) = zeros(M,1);
                    
                    idx_while = idx_while + 1;
                end
            end
        end
        
        L = K - idx_while + 1;
        
        if(L == K)
            H_s = H;
        else
            S_set(S_set_drop) = [];
            
            H_s = H(:,S_set);
        end
        
        varargout{1} = S_set;
    otherwise
        error('Invalid type of selection');
end

end