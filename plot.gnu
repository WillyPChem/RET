#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set output 'Test.eps'

set xlabel 'Time'
set ylabel 'Population'

plot 'Out.txt' u 1:3 w l lw 3 title 'Chromophore 1', \
'Out.txt' u 1:4 w l lw 3 title 'Chromophore 2', \
'Out.txt' u 1:2 w l lw 3 title 'Ground State'
