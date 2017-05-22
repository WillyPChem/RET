#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25

set output 'C2_Iron_Abs_vs_Energy.eps'
set xlabel 'Energy (eV)'
set ylabel 'Efficiency'

plot 'spectrum.txt' u (1240/$1):2 w l lw 3 title 'Iron Absorption C2'
