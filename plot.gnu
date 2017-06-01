#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set output 'Population.eps'

set xlabel 'Time (atomic units)'
set ylabel 'Population'

plot 'Output.txt' u 1:2 w l lw 3 title 'Ground State', \
'Output.txt' u 1:3 w l lw 3 title 'Excited State 1', \
'Output.txt' u 1:4 w l lw 3 title 'Excited State 2'


set output 'Dipole.eps'
set ylabel 'Dipole Moment (atomic units)'
plot 'DipoleMoment.txt' u 1:2 w l lw 3 title 'Dipole Moment'

a = 1.86e9
b = 0.0414
c = 9.075e8
set key top left
set output 'Absorption.eps'
set xlabel 'Energy (eV)'
set ylabel 'Absorption (arb. units)'
plot 'Absorption_Spectrum.txt' u 1:($2/c) w l lw 3 title 'Malachite Green', \
'Absorption_SpectrumAu.txt' u 1:($2/a) w l lw 3 title 'Gold', \
'spectrum.txt' u (1240/$1):($2/b) w l lw 3 title 'Mie Gold'
