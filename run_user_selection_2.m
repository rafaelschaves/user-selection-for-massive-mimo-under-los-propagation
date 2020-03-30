clear;
close all;
clc;

MC = 1000;                                                                 % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)

channel_type = 'ur-los';

M = 100;                                                                   % Number of antennas at the base station

for K = [10 25 50 75 100 150]                                              % Number of users at the cell
    K
    run user_selection_2.m
end