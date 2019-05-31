#!/bin/bash

MC=10000
L=13
CHANNEL=('ur-los' 'rayleigh')

for iDX in 0 1
do
    for M in 64 256
    do
        for K in 18 36 72
        do
            for SNR in -20 -15 -10 -5 0 5 10 
            do
                ARGUMENTS="MC = $MC; M = $M; K = $K; L = $L; snr_db = $SNR; channel_type = ${CHANNEL[$IDX]};"
                nice -n 19 nohup matlab -nojvm -r $ARGUMENTS <user_selection.m> output_${CHANNEL[$IDX]}_$M_$K_$L_$SNR.log &
            done
        done
    done
done