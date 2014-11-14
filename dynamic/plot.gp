set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 1000, 700
set output 'plot.png'
set title "Morse Potential Vibration"
set xlabel "Displacement"
#set ylabel ""
unset key


set xrange [-2:8]
set yrange [-2:75]

plot 'data.dat' u 1:2 w lines title "Potential",\
    '' u 1:3 w lines title "Wavefunction",\
    '' u 1:4 w lines title "Probability Density",\
    '' u 1:5 w lines title "Energy Eigenfunction",\
    '' u 1:6 w lines
#plot for [i=2:20] "../usr/dynamic/data.dat" u 1:i w lines

set output 'time_dependent.png'
set title "Time excitation"
set xlabel "Time (ps)"
plot for [i=2:20] "../usr/dynamic/tdata.dat" u 1:i w lines
