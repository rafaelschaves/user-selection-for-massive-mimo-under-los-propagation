function [user_sel,sel_chnl_mtx,varargout] = userSelector(chnl_mtx, ...
                                                          lrg_scl, ...
                                                          algorithm, ...
                                                          type, ...
                                                          varargin)

N_ARGIN = 6;                                                  

algorithm = upper(algorithm);                                              % Algorithm
type = upper(type);                                                        % Type of selection

if(nargin == N_ARGIN-2)
    n_selected = 5;
    threshold  = 0.5;
elseif(nargin == N_ARGIN-1)
    n_selected = varargin{1};
    threshold  = 0.5;
elseif(nargin == N_ARGIN)
    n_selected = varargin{1};
    threshold  = varargin{2};
else
    error('Invalid number of input arguments');
end

n_antenna = size(chnl_mtx,1);                                              % Number of antennas at base station
n_user    = size(chnl_mtx,2);                                              % Number users

switch algorithm
     case 'RANDOM SELECTION'
        [user_sel,sel_chnl_mtx,varargout{1},varargout{2}] = randomSelection(chnl_mtx,n_selected);
    case 'SEMI-ORTHOGONAL SELECTION'
        [user_sel,sel_chnl_mtx,varargout{1},varargout{2}] = semiOrthogonalSelection(chnl_mtx,n_selected);
    case 'CORRELATION-BASED SELECTION'
        [user_sel,sel_chnl_mtx,varargout{1},varargout{2}] = correlationBasedSelection(chnl_mtx,n_selected,threshold,type);
    case 'ICI-BASED SELECTION'
        switch type
            case 'FIXED'
                user_drop = zeros(n_user - n_selected,1);
        
                chnl_mtx_aux = chnl_mtx;
        
                for i = 1:n_user-n_selected
                    psi = ici(chnl_mtx_aux);
            
                    [~,user_drop(i)] = max(psi);
 
                    chnl_mtx_aux(:,user_drop(i)) = zeros(n_antenna,1);
                end
        
                user_sel = (1:n_user)';
                user_sel(user_drop) = [];
        
                sel_chnl_mtx = chnl_mtx;
                sel_chnl_mtx(:,user_drop) = [];
                
                drop_chnl_mtx = chnl_mtx;
                drop_chnl_mtx(:,user_sel) = [];
        
                varargout{1} = user_drop;
                varargout{2} = drop_chnl_mtx;
            case 'AUTOMATIC'
                user_sel = (1:n_user)';
                
                chnl_mtx_aux = chnl_mtx;
                
                selection = 1;
                idx_while = 1;
                
                if(threshold == 1)
                    idx_while = 1;
                else
                    while(selection == 1)
                        psi = ici(chnl_mtx_aux);
                        
                        decision = (psi > threshold);
                        
                        if(sum(decision) == 0 || idx_while == n_user)
                            selection = 0;
                        else
                            idx_decision = user_sel(decision);
                            
                            [~,idx_aux] = max(psi(decision));
                            
                            user_drop(idx_while) = idx_decision(idx_aux);
                            
                            chnl_mtx_aux(:,user_drop(idx_while)) = zeros(n_antenna,1);
                            
                            idx_while = idx_while + 1;
                        end
                    end
                end
                
                n_selected = n_user - idx_while + 1;
                
                if(n_selected == n_user)
                    user_drop = [];
                    sel_chnl_mtx = chnl_mtx;
                else
                    user_sel(user_drop) = [];
        
                    sel_chnl_mtx = chnl_mtx;
                    sel_chnl_mtx(:,user_drop) = [];
                end
                
                drop_chnl_mtx = chnl_mtx;
                drop_chnl_mtx(:,user_sel) = [];
        
                varargout{1} = user_drop;
                varargout{2} = drop_chnl_mtx;
            otherwise
                error('Invalid type of selection');
        end
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

function [user_sel,sel_chnl_mtx,varargout] = randomSelection(chnl_mtx,n_selected)

n_user = size(chnl_mtx,2);                                                 % Number users

user_drop = randperm(n_user,n_user - n_selected)';                         % User that won't transmitt or receive data

user_sel = (1:n_user)';
user_sel(user_drop) = [];

sel_chnl_mtx = chnl_mtx;
sel_chnl_mtx(:,user_drop) = [];

drop_chnl_mtx = chnl_mtx;
drop_chnl_mtx(:,user_sel) = [];

varargout{1} = user_drop;
varargout{2} = drop_chnl_mtx;

end

function [user_sel,sel_chnl_mtx,varargout] = semiOrthogonalSelection(chnl_mtx,n_selected)

n_antenna = size(chnl_mtx,1);                                              % Number of antennas at base station
n_user    = size(chnl_mtx,2);                                              % Number users

user_sel = zeros(n_selected,1);
ort_proj_sel = zeros(n_antenna,n_selected);

eye_M = eye(n_antenna);

chnl_mtx_aux = chnl_mtx;

for i = 1:n_selected
    if(i == 1)
        ort_proj = chnl_mtx;
        g_norm = vecnorm(ort_proj,2);
        
        [~,user_sel(i)] = max(g_norm);
        
        ort_proj_sel(:,i) = ort_proj(:,user_sel(i))/norm(ort_proj(:,user_sel(i)),2);
        
        chnl_mtx_aux(:,user_sel(i)) = zeros(n_antenna,1);
    else
        proj_mtx = eye_M - ort_proj_sel*ort_proj_sel';
        ort_proj = proj_mtx*chnl_mtx_aux;
        
        g_norm = vecnorm(ort_proj,2);
        
        [~,user_sel(i)] = max(g_norm);
        
        ort_proj_sel(:,i) = ort_proj(:,user_sel(i))/norm(ort_proj(:,user_sel(i)),2);
        
        chnl_mtx_aux(:,user_sel(i)) = zeros(n_antenna,1);
    end
end

user_drop = (1:n_user)';
user_drop(user_sel) = [];

user_sel = sort(user_sel);

sel_chnl_mtx = chnl_mtx;
sel_chnl_mtx(:,user_drop) = [];

drop_chnl_mtx = chnl_mtx;
drop_chnl_mtx(:,user_sel) = [];

varargout{1} = user_drop;
varargout{2} = drop_chnl_mtx;

end

function [user_sel,sel_chnl_mtx,varargout] = correlationBasedSelection(chnl_mtx,n_selected,threshold,type)

n_antenna = size(chnl_mtx,1);                                              % Number of antennas at base station
n_user    = size(chnl_mtx,2);                                              % Number users

I_K = eye(n_user);

chnl_mtx_aux = chnl_mtx;

switch type
    case 'FIXED'
        user_drop = zeros(n_user - n_selected,1);
        
        for i = 1:n_user-n_selected
            corr_mtx     = correlationMatrix(chnl_mtx_aux);
            corr_mtx_aux = corr_mtx - I_K;
            
            [~,idx_corr] = max(corr_mtx_aux(:));
            [corr_i,corr_j] = ind2sub(size(corr_mtx_aux),idx_corr);
            
            h_i = corr_mtx_aux(:,corr_i);
            h_j = corr_mtx_aux(:,corr_j);
            
            h_i(corr_j) = [];
            h_j(corr_i) = [];
            
            if(max(h_i) > max(h_j))
                user_drop(i) = corr_i;
            else
                user_drop(i) = corr_j;
            end
            
            chnl_mtx_aux(:,user_drop(i)) = zeros(n_antenna,1);
        end
        
        user_sel = (1:n_user)';
        user_sel(user_drop) = [];
        
        sel_chnl_mtx = chnl_mtx;
        sel_chnl_mtx(:,user_drop) = [];
        
        drop_chnl_mtx = chnl_mtx;
        drop_chnl_mtx(:,user_sel) = [];
        
        varargout{1} = user_drop;
        varargout{2} = drop_chnl_mtx;
    case 'AUTOMATIC'
        selection = 1;
        idx_while = 1;
        
        if(threshold == 1)
            idx_while = 1;
        else
            while(selection == 1)
                corr_mtx = correlationMatrix(chnl_mtx_aux);
                corr_mtx_aux = corr_mtx - I_K;
                
                [max_corr,idx_corr] = max(corr_mtx_aux(:));
                
                if(max_corr < threshold || idx_while == n_user)
                    selection = 0;
                else
                    [corr_i,corr_j] = ind2sub(size(corr_mtx_aux),idx_corr);
                    
                    h_i = corr_mtx_aux(:,corr_i);
                    h_j = corr_mtx_aux(:,corr_j);
                    
                    h_i(corr_j) = [];
                    h_j(corr_i) = [];
                    
                    if(max(h_i) > max(h_j))
                        user_drop(idx_while) = corr_i;
                    else
                        user_drop(idx_while) = corr_j;
                    end
                    
                    chnl_mtx_aux(:,user_drop(idx_while)) = zeros(n_antenna,1);
                    
                    idx_while = idx_while + 1;
                end
            end
        end
        
        n_selected = n_user - idx_while + 1;
        
        if(n_selected == n_user)
            user_drop = [];
            user_sel = (1:n_user)';
            sel_chnl_mtx = chnl_mtx;
        else
            user_sel = (1:n_user)';
            user_sel(user_drop) = [];
            
            sel_chnl_mtx = chnl_mtx;
            sel_chnl_mtx(:,user_drop) = [];
        end
        
        drop_chnl_mtx = chnl_mtx;
        drop_chnl_mtx(:,user_sel) = [];
        
        varargout{1} = user_drop;
        varargout{2} = drop_chnl_mtx;
    otherwise
        error('Invalid type of selection');
end

end