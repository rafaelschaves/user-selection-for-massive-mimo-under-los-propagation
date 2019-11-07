% Federal University of Rio de Janeiro - UFRJ
% Electrical Engineering Program - COPPE
% Signals, Multimedia, and Telecommunications Laboratory - SMT
%
% Author: Rafael da Silva Chaves
% email: rafael.chves@smt.ufrj.br
%
% Abstract: This function generate a multipath channel for multi-user
% multiple-input multiple-output (MU-MIMO) transmission, following a time
% division duplex (TDD) transmission. This function generates three
% different types of fading channel. The first one is independent Rayleigh
% fading channel (rich scattering fading channel), which models an
% isotropic rich scattering fading. The second is the uniformly random 
% line-of-sight (UR-LoS), which models an environment with only LoS.
%
% References:
%
% [1] - T. L. Marzetta, E. G. Larsson, H. Yang, et al., "Fundamentals of
% Massive MIMO", 1 ed., Cambridge, Cambridge University Press, 2016
%
% [2] - A. M. Sayeed, "Deconstructing Multiantenna Fading Channels", IEEE
% Transaction on Signal Processing, vol. 50, no. 10, pp. 2563-2579, Oct.
% 2006
%
% Help:
%
% [G, beta, varargout] = multipathMUMIMOChannel(cell,   ...
%                                               fading, ...
%                                               varargin)
%
% Inputs:
%
%         -- 'commcell' is a structure with the parameters of the 
%         communication cell.
%
%           * 'commcell.nAntennas' is the field with the number of antennas
%           in the array, it must be a positive integer number.
%
%           * 'commcell.nUsers' is the field with the number of users in
%           the cell, it must be a postive integer number.
%
%           * 'commcell.radius' is the field with the cell radius, it must 
%           be a positive real number.
%
%           * 'commcell.bsHeight' is the filed with the heigth of BS
%           antenna array, it must be a postive real number.
%
%           * 'commcell.userHeight' is the field with the maximum and
%           minimum user heights, it must be a positve real vector with 
%           size 2.
%
%           * 'commcell.frequency' is the field with the carrier frequency
%           of the transmitted signal, it must be a positive real number.
%
%           * 'commcell.meanShdowFad' is the field with the mean of shadow
%           fading distribution in dB, it must be a positve real number.
%
%           * 'commcell.stdDevShadowFad' is the field with the standard
%           deviation of the shadow fading distribution in dB, it must be a
%           postive real number.
%
%           * 'commcell.city' is the field with the type of the city, it 
%           must be a string. It can be 'LARGE', 'MEDIUM' or 'SMALL'.
%
%           * 'fading' is the type of desired fading channel, it must be a
%           string. It can be 'RAYLEIGH', 'UR-LOS' or 'SPARSE'.
%
%         -- 'varargin{1}' is a structure that contains predefined position 
%         for the users and the objects interfering in the transmission.
%
%           * 'varargin{1}.x_user' is the field with the x coordinate for 
%           all users, it must be a postive real vector with size 
%           'commcell.nUsers'.
%
%           * 'varargin{1}.y_user' is the field with the y coordinate for
%           all users, it mus be a positive real vector with size
%           'commcell.nUsers'.
%
%           * 'varargin{2}.x_object' is the field with the x coordinate for
%           all interfering objects in the transmission, it must be a
%           positve real matrix with size 'commcell.nUsers' x 
%           'commcell.nPaths'.
%
%           * 'varargin{2}.y_object' is the field with the y coordinate for
%           all interfering objects in the transmission, it must be a
%           positive real matrix with size 'commcell.nUsers' x
%           'commcell.nPaths'.
%
% Outputs:
%
%         -- 'G' is the small-scale channel matrix, it is a complex matrix 
%         with size 'commcell.nAntennas' x 'commcell.nUsers'.
%
%         -- 'beta' is the large-scale fading vector, it is a real vector
%         with size 'commcell.nUsers'.
%
%         -- 'varargout{1}' is a vector with the position and the angle of 
%         arrival for all users.
%
%         -- 'varargout{2}' is a vector with the position and the angle of
%         arrival for all interfering objects.

function [G, beta, varargout] = massiveMIMOChannel(commcell, ...
                                                   fading,   ...
                                                   varargin)
                                              
% MACROS

N_ARGIN  = 4;                                                              % Number of input arguments
N_ARGOUT = 4;                                                              % Number of output arguments

% Constants

c = 3e8;                                                                   % Light speed

% Cell parameters

n_antenna             = commcell.nAntennas;                                % Number of transmit antennas at base station
n_user                = commcell.nUsers;                                   % Number of user terminals
R                     = commcell.radius;                                   % Cell's raidus (circumradius) in meters
bs_height             = commcell.bsHeight;                                 % Height of base station in meters
user_height           = commcell.userHeight;                               % Minimum and maximum heights of user terminals in meters
f_c                   = commcell.frequency;                                % Carrier frequency of the transmitted signal
mean_shadow_fad_dB    = commcell.meanShadowFad;                            % Shadow fading mean in dB
std_dev_shadow_fad_dB = commcell.stdDevShadowFad;                          % Shadow fading standard deviation in dB
city                  = commcell.city;                                     % Type of the city

city   = upper(city);
fading = upper(fading);                                                    % Type of fading that occurs in the transmission

% Testing for errors

% Erros for wrong numbers of input and output arguments

if (nargin > N_ARGIN)
    error('Wrong number of input arguments');
elseif (nargout > N_ARGOUT)
    error('Wrong number of output arguments');
end

% Erros for wrong values in numeric variables

if (n_antenna <= 0)
    error('Number of antennas must be a positive integer number');
elseif (n_user <= 0)
    error('Number of user terminals must be a positive integer number');
elseif (R <=0)
    error('Cell radius must be a positive real number');
elseif (bs_height <= 0)
    error('Base station height must be a positive real number');
end

H_user = user_height(1) + (user_height(2) - user_height(1))*rand(n_user,1);% Height of user terminals in meters

antenna_spacing = c/(2*f_c);                                               % Antenna spacing of transmitt array

if (nargin == N_ARGIN - 2)
    [x_user,y_user] = userPositionGenerator(n_user,R);
elseif (nargin >= N_ARGIN-1)
    coordinate = varargin{1};
    
    x_user = coordinate.x_user;
    y_user = coordinate.y_user;
end

theta_user = atan2(y_user,x_user);                                         % Departure angle in rad

varargout{1} = [x_user y_user theta_user];

% Large-scale Fading Calculation
% The large-scale coefficient is independent of the small-scale effects and
% depende only the distances of BS and users

mean_shadow_fad    = 10^(mean_shadow_fad_dB/20);                           % Shadow fading mean
std_dev_shadow_fad = 10^(std_dev_shadow_fad_dB/20);                        % Shadow fading standard deviation

mu_shadow_fad      = log(mean_shadow_fad^2/sqrt(std_dev_shadow_fad^2 + ...
    mean_shadow_fad^2));
sigma_shadow_fad   = sqrt(log((std_dev_shadow_fad/mean_shadow_fad)^2 + 1));

d_bs_user = sqrt(x_user.^2 + y_user.^2);                                   % Distance between base station and users in meters
r_bs_user = sqrt(d_bs_user.^2 + (bs_height - H_user).^2);                  % Length of the path traveled by the signal in meters

z = lognrnd(mu_shadow_fad,sigma_shadow_fad,n_user,1);                      % Shadow fading                 

beta = z.*pathLoss(city,bs_height,H_user,r_bs_user,f_c);                   % Large-scale fading

switch fading
    case 'RAYLEIGH'
        G = (randn(n_antenna,n_user) + 1i*randn(n_antenna,n_user))/sqrt(2);% Small-scale fading coefficient matrix
    case 'UR-LOS'
        A = steeringVector(n_antenna, theta_user, antenna_spacing, c/f_c);
        
        phi   = -pi + 2*pi*rand(n_user,1);                                 % Phase shift
        phase = exp(1i*phi);                                               
        Phase = repmat(phase.',n_antenna,1);                                
        
        G = Phase.*A;                                                      % Small-scale fading coefficient matrix
    otherwise
        error('Invalid type of fading');
end

end

function [x,y] = userPositionGenerator(n_coord,R)

r = sqrt(3)/2*R;

aux_cord = rand(n_coord,1);

K_1 = sum(aux_cord < 1/3);
K_2 = sum(aux_cord < 2/3 & aux_cord > 1/3);

u = rand(n_coord,1);
v = rand(n_coord,1);

u_1 = u(1:K_1,1);
v_1 = v(1:K_1,1);

u_2 = u(K_1+1:K_1+K_2,1);
v_2 = v(K_1+1:K_1+K_2,1);

u_3 = u(K_1+K_2+1:n_coord,1);
v_3 = v(K_1+K_2+1:n_coord,1);

x_1 = -R/2*u_1 + R*v_1;
y_1 = r*u_1;

x_2 = -R/2*u_2 - R/2*v_2;
y_2 = -r*u_2 + r*v_2;

x_3 = R*u_3 - R/2*v_3;
y_3 = -r*v_3;

x = [x_1' x_2' x_3']';
y = [y_1' y_2' y_3']';

end

function L = pathLoss(city,h_bs,h_user,distance,frequency)

% COST-231 Hata Model

frequency = frequency/1e6;
distance = distance/1e3;

switch city
    case 'LARGE'
        if((150 <= frequency) && (frequency <= 200))
            C_h = 8.29*(log10(1.54*h_user)).^2 - 1.1;
        elseif((200 < frequency) && (frequency <= 2000))
            C_h = 3.2*(log10(11.75*h_user)).^2 - 4.97;
        else
        end
        C = 3;
    case 'MEDIUM'        
        C_h = 0.8 + (1.1*log10(frequency) - 0.7)*h_user - 1.56*log10(frequency);
        C = 0;
    case 'SMALL'
        C_h = 0.8 + (1.1*log10(frequency) - 0.7)*h_user - 1.56*log10(frequency);
        C = 0;
    otherwise
        error('Invalid type of city');
end

L_dB = -(46.30 + 33.90*log10(frequency) - 13.82*log10(h_bs) - C_h + ...
        (44.9 - 6.55*log10(h_bs))*log10(distance) + C);
    
L = 10.^(L_dB/10);

end

function steering_vector = steeringVector(M,theta,d,lambda)

m   = (0:M-1)';
tau = d*sin(theta)/lambda;
Tau = repmat(tau',M,1);

steering_vector = exp(-2*pi*1i*Tau.*repmat(m,1,size(theta,1)));

end