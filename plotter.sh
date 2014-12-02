#!/bin/bash

if [[ "$1" == "all" ]]
    then
    for i in dynamic/data*.dat; do ./bin/plotter.sh $i;done
    exit 0
fi


col=$(cat $1 | tail -1 | wc -w)

fileout="$1"
fileinp="$1"

gnuplot << EOF
set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 1000, 700
set title "Morse Potential Vibration"
set xlabel "Displacement"
#unset key
set xrange [-1:3]
set yrange [-0.2:4]

set output "${fileout}.png"
plot "$1" u 1:2 w lines title 'Morse Potential', for [i=3:$col] "$1" u 1:i w lines  title 'eigenstate '.(i-3)
EOF
