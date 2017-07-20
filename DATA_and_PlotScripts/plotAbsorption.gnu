#!/Users/jay/gnuplot/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 30
set output 'AbsorptionMG.eps'
set xlabel 'Wavelength (nm)'
a = 8.2e29
set pointsize 2
set yrange [0:3.2]
set xrange [400:800]
set ylabel 'Extinction Cross Section'
plot 'DATA/TRIALA/AbsorptionSpectrum.txt' u (1240/$1):($3/a) w l lw 7 title 'A', \
'DATA/TRIALB/AbsorptionSpectrum.txt' u (1240/$1):($3/a) w l lw 7 title 'B', \
'DATA/TRIALC/AbsorptionSpectrum.txt' u (1240/$1):($3/a) w l lw 7 title 'C', \
'DATA/TRIALD/AbsorptionSpectrum.txt' u (1240/$1):($3/a) w l lw 7 title 'D', \
'DATA/TRIALE/AbsorptionSpectrum.txt' u (1240/$1):($3/a) w l lw 7 title 'E', \
'DATA/TRIALF/AbsorptionSpectrum.txt' u (1240/$1):($3/a) w l lw 7 title 'F', \
'DATA/TRIALG/AbsorptionSpectrum.txt' u (1240/$1):($3/a) w l lw 7 title 'G', \
'DATA/TRIALH/AbsorptionSpectrum.txt' u (1240/$1):($3/a) w l lw 7 title 'H', \
'DATA/TRIALI/AbsorptionSpectrum.txt' u (1240/$1):($3/a) w lp lw 7 pointinterval 5 title 'I'
