#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25

set output 'Pop1_AgWater.eps'
set xrange [0:2000]
set ylabel 'Excited State 1 Population'
set xlabel 'Time (atomic units)'
set key top left
plot 'PopulationMG_AgWater_inf.dat' u 1:($3) w l lw 3 title 'Uncoupled', \
'PopulationMG_AgWater_80.dat' u 1:($3) w l lw 3 lc rgb 'blue' title 'r = 4.23 nm', \
'PopulationMG_AgWater_60.dat' u 1:($3) w l lw 3 lc rgb  'dark-green' title 'r = 3.17 nm', \
'PopulationMG_AgWater_50.dat' u 1:($3) w l lw 3 lc rgb  'red' title 'r = 2.64 nm'
#'PopulationMG_AgWater_50.dat' u 1:($3) w l lw 3 lc rgb  'title '2.64 nm'

set yrange [0:0.0001]
set key top right
set output 'Pop1_AuWater.eps'
plot 'PopulationMG_AuWater_inf.dat' u 1:($3) w l lw 3 title 'Uncoupled', \
'PopulationMG_AuWater_80.dat' u 1:($3) w l lw 3 lc rgb 'blue' title 'r = 4.23 nm', \
'PopulationMG_AuWater_60.dat' u 1:($3) w l lw 3 lc rgb  'dark-green' title 'r = 3.17 nm', \
'PopulationMG_AuWater_50.dat' u 1:($3) w l lw 3 lc rgb  'red' title 'r = 2.64 nm'

set yrange [0:0.00012]
set output 'Pop1_PtWater.eps'
plot 'PopulationMG_PtWater_inf.dat' u 1:($3) w l lw 3 title 'Uncoupled', \
'PopulationMG_PtWater_80.dat' u 1:($3) w l lw 3 lc rgb 'blue' title 'r = 4.23 nm', \
'PopulationMG_PtWater_60.dat' u 1:($3) w l lw 3 lc rgb  'dark-green' title 'r = 3.17 nm', \
'PopulationMG_PtWater_50.dat' u 1:($3) w l lw 3 lc rgb  'red' title 'r = 2.64 nm'

set yrange [0:0.0003]
set output 'Pop1_PtAir.eps'
plot 'PopulationMG_PtAir_inf.dat' u 1:($3) w l lw 3 title 'Uncoupled', \
'PopulationMG_PtAir_80.dat' u 1:($3) w l lw 3 lc rgb 'blue' title 'r = 4.23 nm', \
'PopulationMG_PtAir_60.dat' u 1:($3) w l lw 3 lc rgb  'dark-green' title 'r = 3.17 nm', \
'PopulationMG_PtAir_50.dat' u 1:($3) w l lw 3 lc rgb  'red' title 'r = 2.64 nm'

set yrange [0:0.0001]
set key top right
set output 'Pop2_AuWater.eps'
plot 'PopulationMG_AuWater_inf.dat' u 1:($4) w l lw 3 title 'Uncoupled', \
'PopulationMG_AuWater_80.dat' u 1:($4) w l lw 3 lc rgb 'blue' title 'r = 4.23 nm', \
'PopulationMG_AuWater_60.dat' u 1:($4) w l lw 3 lc rgb  'dark-green' title 'r = 3.17 nm', \
'PopulationMG_AuWater_50.dat' u 1:($4) w l lw 3 lc rgb  'red' title 'r = 2.64 nm'

set yrange [0:0.00012]
set output 'Pop2_PtWater.eps'
plot 'PopulationMG_PtWater_inf.dat' u 1:($4) w l lw 3 title 'Uncoupled', \
'PopulationMG_PtWater_80.dat' u 1:($4) w l lw 3 lc rgb 'blue' title 'r = 4.23 nm', \
'PopulationMG_PtWater_60.dat' u 1:($4) w l lw 3 lc rgb  'dark-green' title 'r = 3.17 nm', \
'PopulationMG_PtWater_50.dat' u 1:($4) w l lw 3 lc rgb  'red' title 'r = 2.64 nm'

set yrange [0:0.0003]
set output 'Pop2_PtAir.eps'
plot 'PopulationMG_PtAir_inf.dat' u 1:($4) w l lw 3 title 'Uncoupled', \
'PopulationMG_PtAir_80.dat' u 1:($4) w l lw 3 lc rgb 'blue' title 'r = 4.23 nm', \
'PopulationMG_PtAir_60.dat' u 1:($4) w l lw 3 lc rgb  'dark-green' title 'r = 3.17 nm', \
'PopulationMG_PtAir_50.dat' u 1:($4) w l lw 3 lc rgb  'red' title 'r = 2.64 nm'


