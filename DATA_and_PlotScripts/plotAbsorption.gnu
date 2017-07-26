#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 30
set output 'AbsorptionMG.eps'
set xlabel 'Wavelength (nm)'
a = 1
set pointsize 2
#set yrange [0:3.2]
set xrange [200:800]
set ylabel 'Extinction Cross Section'
plot '../Emol_300nm_Abs_1.txt' u (1240/$1):($3/a) w l lw 2 title '1', \
'../Emol_300nm_Abs_2.txt' u (1240/$1):($3/a) w l lw 2 title '2', \
'../Emol_300nm_Abs_3.txt' u (1240/$1):($3/a) w l lw 2 title '3', \
'../Emol_300nm_Abs_4.txt' u (1240/$1):($3/a) w l lw 2 title '4', \
'../Emol_300nm_Abs_5.txt' u (1240/$1):($3/a) w l lw 2 title '5', \
'../Emol_300nm_Abs_6.txt' u (1240/$1):($3/a) w l lw 2 title '6', \
'../Emol_300nm_Abs_7.txt' u (1240/$1):($3/a) w l lw 2 title '7', \
'../Emol_300nm_Abs_8.txt' u (1240/$1):($3/a) w l lw 2 title '8', \
'../Emol_300nm_Abs_9.txt' u (1240/$1):($3/a) w l lw 2 dt 2 title '9'
