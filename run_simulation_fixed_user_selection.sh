#!/bin/bash

MC = 10000
L = 13
channel_type = ('ur-los' 'rayleigh')

for idx_channel in 0 1
do
    for M in 64 256
    do
        for K in 18 36 72
        do
            for snr in -20 -15 -10 -5 0 5 10 
            do
                arguments = "MC = $MC; M = $M; K = $K; L = $L; snr_db = $snr; channel_type = ${channel_type[$idx_channel]};"
                nice -n 19 nohup matlab -nojvm -r $arguments <user_selection.m> output_${channel_type[$idx_channel]}_$M_$K_$L_$snr.log &
            done
        done
    done
done