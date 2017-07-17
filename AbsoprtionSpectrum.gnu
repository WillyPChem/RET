#!/usr/bin/gnuplot

set terminal png
set output "Absorption_A.png"

set xlabel 'wavenumber'
set ylabel 'Absoprtion Intensity' 
plot 'AbosprtionSpectrum.txt' u 1:2 title 'NanoParticle', \ 
'AbsorptionSpectrum.txt' u 1:3 title 'Molecule'
