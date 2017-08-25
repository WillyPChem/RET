#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25

set output 'Emission_AgWater.eps'
a=3.5e-58
set yrange [0:1]
set ylabel 'Emission (arb. units)'
set xlabel 'Energy (eV)'
set xrange [1:3.2]
set key top left
plot 'EmissionSpectru_AgWater_inf.dat' u 1:($3/a/0.2) w l lw 3 title 'Uncoupled', \
'EmissionSpectru_AgWater_80.dat' u 1:($3/a/0.2) w l lw 3 lc rgb 'blue' title '4.23 nm', \
'EmissionSpectru_AgWater_60.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'dark-green' title '3.17 nm', \
'EmissionSpectru_AgWater_50.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'red' title '2.64 nm'
#'EmissionSpectru_AgWater_50.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'title '2.64 nm'

f = 1
reset
set yrange [-0.2:0.2]
set xrange [0:2000]
set output 'Dipole_AgWater.eps'
set xlabel 'Time (atomic units)'
set ylabel 'Dipole Moment (atomic units)'
plot 'DipoleMomentMG_AgWater_inf.dat' u ($1*f):2 w l lw 4 title 'Uncoupled', \
'DipoleMomentMG_AgWater_100.dat' u ($1*f):2 w l lw 4 title '5.29 nm', \
'DipoleMomentMG_AgWater_80.dat' u ($1*f):2 w l lw 4 title '4.23 nm', \
'DipoleMomentMG_AgWater_60.dat' u ($1*f):2 w l lw 4 title '3.17 nm', \
'DipoleMomentMG_AgWater_50.dat' u ($1*f):2 w l lw 4 title '2.64 nm'

reset
set xrange [0:3000]
set xlabel 'Time (atomic units)'
set ylabel 'Excited-state 1 Population'
set output 'Population1_AgWater.eps'
plot 'PopulationMG_AgWater_inf.dat' u ($1*f):3 w l lw 4 title 'Uncoupled', \
'PopulationMG_AgWater_100.dat' u ($1*f):3 w l lw 4 title '5.29 nm', \
'PopulationMG_AgWater_80.dat' u ($1*f):3 w l lw 4 title '4.23 nm', \
'PopulationMG_AgWater_60.dat' u ($1*f):3 w l lw 4 title '3.17 nm', \
'PopulationMG_AgWater_50.dat' u ($1*f):3 w l lw 4 title '2.64 nm'

set ylabel 'Excited-state 2 Population'
set output 'Population2_AgWater.eps'
plot 'PopulationMG_AgWater_inf.dat' u ($1*f):4 w l lw 4 title 'Uncoupled', \
'PopulationMG_AgWater_100.dat' u ($1*f):4 w l lw 4 title '5.29 nm', \
'PopulationMG_AgWater_80.dat' u ($1*f):4 w l lw 4 title '4.23 nm', \
'PopulationMG_AgWater_60.dat' u ($1*f):4 w l lw 4 title '3.17 nm', \
'PopulationMG_AgWater_50.dat' u ($1*f):4 w l lw 4 title '2.64 nm'

