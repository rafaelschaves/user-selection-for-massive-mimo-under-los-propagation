clear;
close all;
clc;

rng('shuffle');                                                          % Necessary for different seeds in each run of this script

MC = 5;                                                                  % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
MC_ERR = 100;

M = 50;                                                                  % Number of antennas at the base station

for K = 75                                                               % Number of users at the cell
    K
    run user_selection_ur_los.m
end