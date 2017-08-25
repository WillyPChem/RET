#!/usr/bin/gnuplot

set palette rgb 21,22,23
set size square
set terminal postscript enhanced color 'Helvetica' 25
set xlabel 'Distance (atomic units)'
set ylabel 'Energy (eV)'
set output 'Ag_Scan.eps'
set pm3d map
set xrange [40:70]
set yrange [2:6]
splot 'SCAN_AG_WATER_FINER.txt' u 1:2:7 notitle 

set output 'Au_Scan.eps'
splot 'SCAN_AU_WATER_FINER.txt' u 1:2:7 notitle

set output 'SMA_Water_Scan.eps'
splot 'SCAN_SMA_WATER_FINER.txt' u 1:2:7 notitle

set output 'SMA_Air_Scan.eps'
splot 'SCAN_SMA_AIR_FINER.txt' u 1:2:7 notitle
