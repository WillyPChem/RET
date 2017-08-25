#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25

set output 'Emission_AgWater.eps'
a=0.25*3.5e-58
set ylabel 'Emission (arb. units)'
set xlabel 'Energy (eV)'
set xrange [1.5:3.5]
set key top left
plot 'EmissionSpectru_AgWater_inf.dat' u 1:($3/a/0.2) w l lw 3 title 'Uncoupled', \
'EmissionSpectru_AgWater_80.dat' u 1:($3/a/0.2) w l lw 3 lc rgb 'blue' title 'r = 4.23 nm', \
'EmissionSpectru_AgWater_60.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'dark-green' title 'r = 3.17 nm', \
'EmissionSpectru_AgWater_50.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'red' title 'r = 2.64 nm'
#'EmissionSpectru_AgWater_50.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'title '2.64 nm'

set key top right
set output 'Emission_AuWater.eps'
plot 'EmissionSpectru_AuWater_inf.dat' u 1:($3/a/0.2) w l lw 3 title 'Uncoupled', \
'EmissionSpectru_AuWater_80.dat' u 1:($3/a/0.2) w l lw 3 lc rgb 'blue' title 'r = 4.23 nm', \
'EmissionSpectru_AuWater_60.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'dark-green' title 'r = 3.17 nm', \
'EmissionSpectru_AuWater_50.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'red' title 'r = 2.64 nm'


set output 'Emission_PtWater.eps'
plot 'EmissionSpectru_PtWater_inf.dat' u 1:($3/a/0.2) w l lw 3 title 'Uncoupled', \
'EmissionSpectru_PtWater_80.dat' u 1:($3/a/0.2) w l lw 3 lc rgb 'blue' title 'r = 4.23 nm', \
'EmissionSpectru_PtWater_60.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'dark-green' title 'r = 3.17 nm', \
'EmissionSpectru_PtWater_50.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'red' title 'r = 2.64 nm'

set output 'Emission_PtAir.eps'
plot 'EmissionSpectru_PtAir_inf.dat' u 1:($3/a/0.2) w l lw 3 title 'Uncoupled', \
'EmissionSpectru_PtAir_80.dat' u 1:($3/a/0.2) w l lw 3 lc rgb 'blue' title 'r = 4.23 nm', \
'EmissionSpectru_PtAir_60.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'dark-green' title 'r = 3.17 nm', \
'EmissionSpectru_PtAir_50.dat' u 1:($3/a/0.2) w l lw 3 lc rgb  'red' title 'r = 2.64 nm'


