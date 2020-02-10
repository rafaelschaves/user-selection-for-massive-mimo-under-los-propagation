clear;
close all;
clc;

MC = 1000;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)

channel_type = 'ur-los';

for bs_power = [1 10]                                                      % Power available at BS
    for radius = [100 500 1000 2000]                                       % Cell's radius
        for M = [64 128 256 512]                                           % Number of antennas at the base station
            for K = [8 16 32 64]                                           % Number of users
                run max_min_error_prob.m
            end
        end
    end
end