addpath('./functions/')

% Cheking deirectory

dir_save_dl  = './results/scheduling/downlink/';
root_save_dl = [dir_save_dl 'spectral_efficiency_mf_'];

if ~exist(dir_save_dl,'dir')
    mkdir(dir_save_dl);
end

% Checking variables

if ~exist('MC','var')
    MC = 1;                                                            % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
end

if ~exist('M','var')
    M = 200;                                                               % Number of antennas at the base station
end

if ~exist('K','var')
    K = 20;                                                                % Number of users at the cell
end

if ~exist('L','var')
    L = 10;                                                                % Number of selected users
end

if ~exist('channel_type','var')
    channel_type = 'ur-los';
end

N_ALG = 3;
N_PRE = 2;

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

se   = zeros(K,N_PRE,MC);
se_s = zeros(L,N_PRE,MC,N_ALG);

psi   = zeros(K,MC);
psi_s = zeros(L,MC,N_ALG);

% algorithm_type = {'exhaustive search selection ep', ...
%                   'semi-orthogonal selection', ...
%                   'correlation-based selection', ...
%                   'ici-based selection', ...
%                   'fr-based selection'};

algorithm_type = {'semi-orthogonal selection', ...
                  'correlation-based selection', ...
                  'ici-based selection'};

for mc = 1:MC
    mc
    
    [H,~] = massiveMIMOChannel(commcell,channel_type);
    
    [se(:,1,mc),se(:,2,mc)] = DLspectralEfficiency(H,beta,snr,1/K);      % No Selection
          
    for alg_idx = 1:N_ALG
        H_s = userSelector(H,beta,snr,algorithm_type{alg_idx},'fixed',L,[]);
        
%         if alg_idx == 1
%             beta_s = [beta(S_set(:,1)) beta(S_set(:,2))];
%         else
%             beta_s = beta(S_set);
%         end
        
        [se_s(:,1,mc,alg_idx),se_s(:,2,mc,alg_idx)] = DLspectralEfficiency(H_s,beta,snr,1/L);
    end
end

% save([root_save_dl strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
%       num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_eff) '_dB_MC_' ...
%       num2str(MC) '.mat'],'se','se_s','psi','psi_s');