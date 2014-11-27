#!/bin/bash

fileout="$1"
fileinp="$1"

gnuplot << EOF
set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 1000, 700
set title "Morse Potential Vibration"
set xlabel "Displacement"
#unset key
set xrange [-2:8]
set yrange [-2:75]

set output "${fileout}.png"
plot "$1" u 1:2 w lines title 'Morse Potential', for [i=3:10] "$1" u 1:i w lines title 'eigenstate '.i
EOF
