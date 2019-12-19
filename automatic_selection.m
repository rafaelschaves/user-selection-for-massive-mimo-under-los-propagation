addpath('./functions/')

% Cheking deirectory

dir_save_dow = './results/auto_scheduling/downlink/';
dir_save_upl = './results/auto_scheduling/uplink/';

root_save_dow = [dir_save_dow 'rate_mf_'];
root_save_upl = [dir_save_upl 'rate_mf_'];

if ~exist(dir_save_dow,'dir')
    mkdir(dir_save_dow);
end

if ~exist(dir_save_upl,'dir')
    mkdir(dir_save_upl);
end

% Checking variables

if ~exist('MC','var')
    MC = 100;                                                            % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
end

if ~exist('M','var')
    M = 64;                                                                % Number of antennas at the base station
end

if ~exist('K','var')
    K = 72;                                                                % Number of users at the cell
end

if ~exist('snr_db','var')
    snr_db = 10;                                                           % SNR in dB
end

if ~exist('channel_type','var')
    channel_type = 'ur-los';
end

commcell.nAntennas       = M;                                              % Number of Antennas
commcell.nUsers          = K;                                              % Number of Users
commcell.radius          = 200;                                            % Cell's raidus (circumradius) in meters
commcell.bsHeight        = 32;                                             % Height of base station in meters
commcell.userHeight      = [1 2];                                          % Height of user terminals in meters ([min max])
commcell.nPaths          = 30;                                             % Number of Multipaths
commcell.frequency       = 1.9e9;                                          % Carrier frequency in Hz
commcell.meanShadowFad   = 0;                                              % Shadow fading mean in dB
commcell.stdDevShadowFad = 8;                                              % Shadow fading standard deviation in dB
commcell.city            = 'large';                                        % Type of city

settings.coherenceTime           = 1;                                      % Coherence time in samples
settings.PilotTime               = 0;                                      % Pilot time in samples
settings.uplinkDownlinkTimeRatio = 0.5;                                    % Ratio between the uplink and downlink payload time
settings.bandwidth               = 20e6;                                   % Sytem bandwidth in Hz
settings.cellArea                = 1;                                      % Cell area in km^2

N_ALG = 2;

snr = 10.^((snr_db)/10);                                                   % SNR

% Normalized threshold

tau_norm_min = 0;
tau_norm_max = 1;

tau_norm_step = 0.01;

tau_norm = tau_norm_min:tau_norm_step:tau_norm_max;

N_TAU = length(tau_norm);

tau_icibs_max = 0.25;

tau(:,1) = tau_norm';                                                      % Threshold - CBS
tau(:,2) = tau_icibs_max*tau_norm';                                        % Threshold - ICIBS

% Initialization

L = zeros(MC,N_TAU,N_ALG);

user_sel = cell(MC,N_TAU,N_ALG);
se_ul    = cell(MC,N_TAU,N_ALG);
se_dl    = cell(MC,N_TAU,N_ALG);
psi      = cell(MC,N_TAU,N_ALG);

algorithm_type = {'correlation-based selection','ici-based selection'};

for mc = 1:MC
    mc
    
    [G,~] = massiveMIMOChannel(commcell,channel_type);
        
    for alg_idx = 1:N_ALG
        for tau_idx = 1:N_TAU
            [user_sel{mc,tau_idx,alg_idx},H_sel] = userSelector(G,algorithm_type{alg_idx},'automatic',[],tau(tau_idx,alg_idx));
            
            L(mc,tau_idx,alg_idx) = size(H_sel,2);
            
            [Q,W] = decoderMatrix(H_sel,'mf');
                        
            pow_ul = ones(L(mc,tau_idx,alg_idx),1);
            pow_dl = ones(L(mc,tau_idx,alg_idx),1)/L(mc,tau_idx,alg_idx);
            
            [~,se_ul{mc,tau_idx,alg_idx}] = throughput(H_sel,Q,pow_ul,'uplink',snr,settings);
            [~,se_dl{mc,tau_idx,alg_idx}] = throughput(H_sel,W,pow_dl,'downlink',snr,settings);
            
            psi{mc,tau_idx,alg_idx} = ici(H_sel);
        end
    end
end

save([root_save_dow strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_db) '_dB_MC_' ...
      num2str(MC) '.mat'],'user_sel','se_dl','psi','L');

save([root_save_upl strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_tau_' num2str(N_TAU) '_SNR_' num2str(snr_db) '_dB_MC_' ...
      num2str(MC) '.mat'],'user_sel','se_ul','psi','L');