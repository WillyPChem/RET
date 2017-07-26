#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 30
set output 'PopulationMG.eps'
set xlabel 'Time (femtoseconds)'
a = 1
set xrange [0:20]
set pointsize 2
#set yrange [0:($1*0.02418)]
set ylabel 'Population'
plot '../Emol_300nm_Pop_1.txt' u ($1*0.02418):3 w l lw 2 title '1', \
'../Emol_300nm_Pop_2.txt' u ($1*0.02418):3 w l lw 2 title '2', \
'../Emol_300nm_Pop_3.txt' u ($1*0.02418):3 w l lw 2 title '3', \
'../Emol_300nm_Pop_4.txt' u ($1*0.02418):3 w l lw 2 title '4', \
'../Emol_300nm_Pop_5.txt' u ($1*0.02418):3 w l lw 2 title '5', \
'../Emol_300nm_Pop_6.txt' u ($1*0.02418):3 w l lw 2 title '6', \
'../Emol_300nm_Pop_7.txt' u ($1*0.02418):3 w l lw 2 title '7', \
'../Emol_300nm_Pop_8.txt' u ($1*0.02418):3 w l lw 2 title '8', \
'../Emol_300nm_Pop_9.txt' u ($1*0.02418):3 w l lt 0 lw 2 title '9'
