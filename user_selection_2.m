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

if ~exist('channel_type','var')
    channel_type = 'ur-los';
end

N_ALG = 3;
N_PRE = 3;

commcell.nAntennas       = M;                                              % Number of Antennas
commcell.nUsers          = K;                                              % Number of Users
commcell.radius          = 500;                                           % Cell's raidus (circumradius) in meters
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

[~,snr_db] = linkBudgetCalculation(linkprop);                              % SNR in dB
snr        = 10.^(snr_db/10);

beta_db = -148 - 37.6*log10(commcell.radius/1000);
beta    = 10.^(beta_db/10);

snr_eff = round(snr_db + beta_db);

% Initialization

algorithm_type = {'semi-orthogonal selection', ...
                  'correlation-based selection', ...
                  'ici-based selection'};

if K > M
    L_max = M;
else
    L_max = K-1;
end

se  = zeros(K,N_PRE,MC);
% psi = zeros(K,MC);
    
se_s_all_L = zeros(L_max,L_max,N_PRE,N_ALG,MC);
% psi_all_L = cell(L_max,MC);
              
for mc = 1:MC
    mc
    
    [H,~] = massiveMIMOChannel(commcell,channel_type);
    
    [se(:,1,mc),se(:,2,mc),se(:,3,mc)] = DLspectralEfficiency(H,beta,snr,1/K);        % No Selection
        
    for L = 1:L_max                                                      % Number of selected users
        L
        
        for alg_idx = 1:N_ALG
            H_s = userSelector(H,beta,snr,algorithm_type{alg_idx},'fixed',L,[]);
            
            %         if alg_idx == 1
            %             beta_s = [beta(S_set(:,1)) beta(S_set(:,2))];
            %         else
            %             beta_s = beta(S_set);
            %         end
            
            [se_s_mf,se_s_zf,se_s_mmse] = DLspectralEfficiency(H_s,beta,snr,1/L);
            
            se_s_all_L(:,L,1,alg_idx,mc) = [se_s_mf; zeros(L_max-L,1)];
            se_s_all_L(:,L,2,alg_idx,mc) = [se_s_zf; zeros(L_max-L,1)];
            se_s_all_L(:,L,3,alg_idx,mc) = [se_s_mmse; zeros(L_max-L,1)];
        end
    end
end

save([root_save strrep(channel_type,'-','_') '_M_' num2str(M) ...
      '_K_' num2str(K) '_SNR_' num2str(snr_eff) '_dB_MC_' ...
      num2str(MC) '.mat'],'se','se_s_all_L');