clear;
close all;
clc;

rng('shuffle');    

MC = 500;                                                             % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)

for M = 100                                                           % Number of antennas at the base station
    for K = 150                                                       % Number of users at the cell
        for theta_mid = [0 pi/4]
            for theta_step = pi/180:pi/360:pi/36 
                run user_selection_clustered.m
            end
        end
    end
end