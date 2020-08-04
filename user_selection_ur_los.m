addpath('./functions/')

% Cheking deirectory

dir_save  = './results/scheduling/downlink/';
root_save = [dir_save 'spectral_efficiency_all_L_'];

if ~exist(dir_save,'dir')
    mkdir(dir_save);
end

% Checking variables

if ~exist('MC','var')
    MC = 5;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
end

if ~exist('M','var')
    M = 200;                                                               % Number of antennas at the base station
end

if ~exist('K','var')
    K = 20;                                                                % Number of users at the cell
end

err = [0 pi/(6*M) pi/(3*M) pi/(2*M)];

N_ALG = 3;
N_PRE = 3;
N_ERR = length(err);

commcell.nAntennas       = M;                                              % Number of Antennas
commcell.nUsers          = K;                                              % Number of Users
commcell.radius          = 500;                                            % Cell's raidus (circumradius) in meters
commcell.bsHeight        = 32;                                             % Height of base station in meters
commcell.userHeight      = [1 2];                                          % Height of user terminals in meters ([min max])
commcell.frequency       = 1.9e9;                                          % Carrier frequency in Hz
commcell.meanShadowFad   = 0;                                              % Shadow fading mean in dB
commcell.stdDevShadowFad = 8;                                              % Shadow fading standard deviation in dB
commcell.city            = 'large';                                        % Type of city

linkprop.bsPower         = 10;                                             % in Watts
linkprop.userPower       = 0.2;                                            % in Watts
linkprop.AntennaGainBS   = 0;                                              % in dBi
linkprop.AntennaGainUser = 0;                                              % in dBi
linkprop.noiseFigureBS   = 9;                                              % in dB
linkprop.noiseFigureUser = 9 ;                                             % in dB
linkprop.bandwidth       = 20e6;                                           % in Hz

channel_type = 'ur-los';

[~,snr_db] = linkBudgetCalculation(linkprop);                              % SNR in dB
snr        = 10.^(snr_db/10);

beta_db = -148 - 37.6*log10(commcell.radius/1000);
beta    = 10.^(beta_db/10);

snr_eff = round(snr_db + beta_db);

% Initialization

algorithm_type = {'semi-orthogonal selection','correlation-based selection','ici-based selection'};

if K > M
    L_max = M;
else
    L_max = K-1;
end


se            = zeros(K,N_PRE,N_ERR,MC);
se_s_all_L    = zeros(L_max,L_max,N_PRE,N_ALG,N_ERR,MC);
S_set         = zeros(K,L_max,N_ALG,N_ERR,MC);
pos_and_theta = zeros(K,3);

for mc = 1:MC
    mc
    
    [H,~,pos_and_theta] = massiveMIMOChannel(commcell,channel_type);
    
    for err_idx = 1:N_ERR
        H_hat = urlosChannelEstimate(commcell,pos_and_theta(:,3),err(err_idx));
        
        [se(:,1,err_idx,mc),se(:,2,err_idx,mc),se(:,3,err_idx,mc)] = DLspectralEfficiency(H,beta,snr,1/K,H_hat);                                      % No Selection
        
        for L = 1:L_max                                                                                                                               % Number of selected users
            for alg_idx = 1:N_ALG
                [H_s, S_set_aux] = userSelector(H_hat,beta,snr,algorithm_type{alg_idx},'fixed',L,[]);
                
                %         if alg_idx == 1
                %             beta_s = [beta(S_set(:,1)) beta(S_set(:,2))];
                %         else
                %             beta_s = beta(S_set);
                %         end
                
                S_set(S_set_aux,L,alg_idx,err_idx,mc) = 1;
                
                [se_s_mf,se_s_zf,se_s_mmse] = DLspectralEfficiency(H_s,beta,snr,1/L);
                
                se_s_all_L(:,L,1,alg_idx,err_idx,mc) = [se_s_mf; zeros(L_max-L,1)];
                se_s_all_L(:,L,2,alg_idx,err_idx,mc) = [se_s_zf; zeros(L_max-L,1)];
                se_s_all_L(:,L,3,alg_idx,err_idx,mc) = [se_s_mmse; zeros(L_max-L,1)];
            end
        end
    end
end

save([root_save strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr_eff) '_dB_MC_' num2str(MC) '.mat'],'se','se_s_all_L','S_set');