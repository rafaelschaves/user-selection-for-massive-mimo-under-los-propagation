#!/bin/bash

MC=10000
M=64
K=18
L=13
CHANNEL=ur-los

for SNR in -20 -15 -10 -5 0 5 10 
do
    ARGUMENTS="MC=$MC; M=$M; K=$K; L=$L; snr_db=$SNR; channel_type='$CHANNEL';"
    nice -n 19 nohup matlab -nojvm -r $ARGUMENTS <user_selection.m> ./log/output_$SNR.log &
done