function [user_set_sel,sel_chnl_mtx] = userScheduling(chnl_mtx, ...
                                                      algorithm, ...
                                                      type, ...
                                                      varargin)

N_ARGIN = 5;                                                  

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
        user_set_drop  = randperm(n_user,n_user - n_selected)';            % User that won't transmitt or receive data
        
        user_set_sel = (1:n_user)';
        user_set_sel(user_set_drop) = []; 
    
        sel_chnl_mtx = chnl_mtx;
        sel_chnl_mtx(:,user_set_drop) = [];
    case 'SEMI-ORTHOGONAL SELECTION'
        user_set_sel = zeros(n_selected,1);
        ort_proj_sel = zeros(n_antenna,n_selected);
        
        eye_M = eye(n_antenna);
        
        chnl_mtx_aux = chnl_mtx;
        
        for i = 1:n_selected
            if(i == 1)
                ort_proj = chnl_mtx;
                g_norm = vecnorm(ort_proj,2);
                
                [~,user_set_sel(i)] = max(g_norm);
                
                ort_proj_sel(:,i) = ort_proj(:,user_set_sel(i))/norm(ort_proj(:,user_set_sel(i)),2);
                
                chnl_mtx_aux(:,user_set_sel(i)) = zeros(n_antenna,1);
            else
                proj_mtx = eye_M - sum(ort_proj_sel*ort_proj_sel',2);
                ort_proj = proj_mtx*chnl_mtx_aux;
                
                g_norm = vecnorm(ort_proj,2);
                
                [~,user_set_sel(i)] = max(g_norm);
                
                ort_proj_sel(:,i) = ort_proj(:,user_set_sel(i))/norm(ort_proj(:,user_set_sel(i)),2);
                
                chnl_mtx_aux(:,user_set_sel(i)) = zeros(n_antenna,1);
            end
        end
        
        user_set_drop = (1:n_user)';
        user_set_drop(user_set_sel) = [];
        
        user_set_sel = sort(user_set_sel);
        
        sel_chnl_mtx = chnl_mtx;
        sel_chnl_mtx(:,user_set_drop) = [];
    case 'CORRELATION-BASED SELECTION'
        switch type
            case 'FIXED'
                user_set_drop = zeros(n_user - n_selected,1);
                
                eye_K = eye(n_user);
                
                chnl_mtx_aux = chnl_mtx;
                
                for i = 1:n_user-n_selected
                    corr_mtx     = correlationMatrix(chnl_mtx_aux);
                    corr_mtx_aux = corr_mtx - eye_K;
                    
                    [~,idx_corr] = max(corr_mtx_aux(:));
                    [corr_i,corr_j] = ind2sub(size(corr_mtx_aux),idx_corr);
                    
                    h_i = corr_mtx_aux(:,corr_i);
                    h_j = corr_mtx_aux(:,corr_j);
                    
                    h_i(corr_j) = [];
                    h_j(corr_i) = [];
                    
                    if(max(h_i) > max(h_j))
                        user_set_drop(i) = corr_i;
                    else
                        user_set_drop(i) = corr_j;
                    end
                    
                    chnl_mtx_aux(:,user_set_drop(i)) = zeros(n_antenna,1);
                end
                
                user_set_sel = (1:n_user)';
                user_set_sel(user_set_drop) = [];
                
                sel_chnl_mtx = chnl_mtx;
                sel_chnl_mtx(:,user_set_drop) = [];
            case 'AUTOMATIC'
                eye_K = eye(n_user);
                
                chnl_mtx_aux = chnl_mtx;
                
                selection = 1;
                idx_while = 1;
                
                while(selection == 1)
                    corr_mtx = correlationMatrix(chnl_mtx_aux);
                    corr_mtx_aux = corr_mtx - eye_K;
                    
                    [max_corr,idx_corr] = max(corr_mtx_aux(:));
                    
                    if(max_corr < threshold)
                        selection = 0;
                    else
                        [corr_i,corr_j] = ind2sub(size(corr_mtx_aux),idx_corr);
                    
                        h_i = corr_mtx_aux(:,corr_i);
                        h_j = corr_mtx_aux(:,corr_j);
                    
                        h_i(corr_j) = [];
                        h_j(corr_i) = [];
                    
                        if(max(h_i) > max(h_j))
                            user_set_drop(idx_while) = corr_i;
                        else
                            user_set_drop(idx_while) = corr_j;
                        end
                    
                        chnl_mtx_aux(:,user_set_drop(idx_while)) = zeros(n_antenna,1);
                        
                        idx_while = idx_while + 1;
                    end
                end
                
                n_selected = n_user - idx_while + 1;
                
                if(n_selected == n_user)
                    user_set_sel = (1:n_user)';
                    sel_chnl_mtx = chnl_mtx;
                else
                    user_set_sel = (1:n_user)';
                    user_set_sel(user_set_drop) = [];
        
                    sel_chnl_mtx = chnl_mtx;
                    sel_chnl_mtx(:,user_set_drop) = [];
                end
            otherwise
                error('Invalid type of selection');
        end
    case 'ICI-BASED SELECTION'
        switch type
            case 'FIXED'
                user_set_drop = zeros(n_user - n_selected,1);
        
                chnl_mtx_aux = chnl_mtx;
        
                for i = 1:n_user-n_selected
                    psi = ici(chnl_mtx_aux);
            
                    [~,user_set_drop(i)] = max(psi);
 
                    chnl_mtx_aux(:,user_set_drop(i)) = zeros(n_antenna,1);
                end
        
                user_set_sel = (1:n_user)';
                user_set_sel(user_set_drop) = [];
        
                sel_chnl_mtx = chnl_mtx;
                sel_chnl_mtx(:,user_set_drop) = [];
            case 'AUTOMATIC'
                user_set_sel = (1:n_user)';
                
                chnl_mtx_aux = chnl_mtx;
                
                selection = 1;
                idx_while = 1;
                
                while(selection == 1)
                    psi = ici(chnl_mtx_aux);
                    
                    decision = (psi > threshold);
                    
                    if(sum(decision) == 0)
                        selection = 0;
                    else
                        idx_decision = user_set_sel(decision);
                    
                        [~,idx_aux] = max(psi(decision));
                    
                        user_set_drop(idx_while) = idx_decision(idx_aux);
                    
                        chnl_mtx_aux(:,user_set_drop(idx_while)) = zeros(n_antenna,1);
                        
                        idx_while = idx_while + 1;
                    end
                end
                
                n_selected = n_user - idx_while + 1;
                
                if(n_selected == n_user)
                    sel_chnl_mtx = chnl_mtx;
                else
                    user_set_sel(user_set_drop) = [];
        
                    sel_chnl_mtx = chnl_mtx;
                    sel_chnl_mtx(:,user_set_drop) = [];
                end
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

mtx_norm = zeros(n_row,n_col);

for k = 1:n_col
    norm_mtx_k = norm(mtx(:,k),2);
    
    if(norm_mtx_k == 0)
        mtx_norm(:,k) = 0;
    else
        mtx_norm(:,k) = mtx(:,k)/norm_mtx_k;
    end
end

corr_mtx     = abs(mtx_norm'*mtx_norm);
              
end