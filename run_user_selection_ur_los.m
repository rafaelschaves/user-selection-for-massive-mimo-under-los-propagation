clear;
close all;
clc;

rng('shuffle');                                                            % Necessary for different seeds in each run of this script

MC = 1000;                                                                 % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)

M = 50;                                                                    % Number of antennas at the base station

% M = 50  and K = [10 25 50 75]
% M = 100 and K = [10 25 50 75 100 150]
% M = 200 and K = [10 25 50 75 100 150 200 250]

for K = 75                                                            % Number of users at the cell
    K
    run user_selection_ur_los.m
end