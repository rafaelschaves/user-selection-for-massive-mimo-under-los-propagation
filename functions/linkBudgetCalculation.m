function [snr_u_effe,snr_d_effe] = linkBudgetCalculation(linkprop,varargin)

% MACROS

N_ARGIN  = 3;
N_ARGOUT = 2;

% Constants

k_B = 1.38e-23;                                                            % Boltzman constant in Joule/Kelvin

% Link Properties

% n_user            = linkprop.nUsers;
bs_power          = linkprop.bsPower;              % in Watts
user_power        = linkprop.userPower;            % in Watts
antenna_gain_bs   = linkprop.AntennaGainBS;        % in dBi
antenna_gain_user = linkprop.AntennaGainUser;      % in dBi
noise_figure_bs   = linkprop.noiseFigureBS;        % in dB
noise_figure_user = linkprop.noiseFigureUser;      % in dB
bandwidth         = linkprop.bandwidth;            % in Hz

% Testing for errors

% Erros for wrong numbers of input and output arguments

if (nargin > N_ARGIN)
    error('Wrong number of input arguments');
elseif (nargout > N_ARGOUT)
    error('Wrong number of output arguments');
end

% % Erros for wrong values in numeric variables
% 
% if (n_antenna <= 0)
%     error('Number of antennas must be a positive integer number');
% elseif (n_user <= 0)
%     error('Number of user terminals must be a positive integer number');
% elseif (R <=0)
%     error('Cell radius must be a positive real number');
% elseif (bs_height <= 0)
%     error('Base station height must be a positive real number');
% end

if (nargin == N_ARGIN - 2)
    large_scale_fading = 1;
    T_noise = 300;
elseif (nargin == N_ARGIN - 1)
    large_scale_fading = varargin{1};
    T_noise = 300;
elseif (nargin == N_ARGIN)
    large_scale_fading = varargin{1};
    T_noise = varargin{2};
end

bs_power_dbm   = 10*log10(1000*bs_power);
user_power_dbm = 10*log10(1000*user_power);

BN_0     = bandwidth*k_B*T_noise;
BN_0_dbm = 10*log10(1000*BN_0);

effective_noise_u = BN_0_dbm + noise_figure_bs;
effective_noise_d = BN_0_dbm + noise_figure_user;

snr_u = user_power_dbm + antenna_gain_bs + antenna_gain_user - effective_noise_u;
snr_d = bs_power_dbm + antenna_gain_bs + antenna_gain_user - effective_noise_d;

snr_u_effe = snr_u*large_scale_fading;
snr_d_effe = snr_d*large_scale_fading;

end