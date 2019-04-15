function [user_set_sel,sel_chnl_mtx] = userScheduling(chnl_mtx, ...
                                                      n_selected, ...
                                                      algorithm)

algorithm = upper(algorithm);                                              % Algorithm

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
    case 'ICI-BASED SELECTION'
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
    otherwise
        error('Invalid algorithm');
end

end