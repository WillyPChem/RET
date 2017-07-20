#!/Users/jay/gnuplot/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 35
set output 'DipoleA.eps'
set xrange [0:20]
set yrange [-9:9]
set xlabel 'Time (femtoseconds)'
set ylabel 'Dipole Moment (atomic units)'
plot 'DATA/TRIALA/DipoleMoment.txt' u ($1*0.02418):2 w l lw 7 title 'Au', \
'DATA/TRIALA/DipoleMomentMG.txt' u ($1*0.02418):2 w l lw 7 title 'MG'

set output 'DipoleD.eps'
plot 'DATA/TRIALD/DipoleMoment.txt' u ($1*0.02418):2 w l lw 7 title 'Au', \
'DATA/TRIALD/DipoleMomentMG.txt' u ($1*0.02418):2 w l lw 7 title 'MG'

set output 'DipoleG.eps'
plot 'DATA/TRIALG/DipoleMoment.txt' u ($1*0.02418):2 w l lw 7 title 'Au', \
'DATA/TRIALG/DipoleMomentMG.txt' u ($1*0.02418):2 w l lw 7 title 'MG'

unset ytics
unset ylabel
set output 'DipoleB.eps'
plot 'DATA/TRIALB/DipoleMoment.txt' u ($1*0.02418):2 w l lw 7 title 'Au', \
'DATA/TRIALB/DipoleMomentMG.txt' u ($1*0.02418):2 w l lw 7 title 'MG'

set output 'DipoleC.eps'
plot 'DATA/TRIALC/DipoleMoment.txt' u ($1*0.02418):2 w l lw 7 title 'Au', \
'DATA/TRIALC/DipoleMomentMG.txt' u ($1*0.02418):2 w l lw 7 title 'MG'

set output 'DipoleE.eps'
plot 'DATA/TRIALE/DipoleMoment.txt' u ($1*0.02418):2 w l lw 7 title 'Au', \
'DATA/TRIALE/DipoleMomentMG.txt' u ($1*0.02418):2 w l lw 7 title 'MG'

set output 'DipoleF.eps'
plot 'DATA/TRIALF/DipoleMoment.txt' u ($1*0.02418):2 w l lw 7 title 'Au', \
'DATA/TRIALF/DipoleMomentMG.txt' u ($1*0.02418):2 w l lw 7 title 'MG'

set output 'DipoleH.eps'
plot 'DATA/TRIALH/DipoleMoment.txt' u ($1*0.02418):2 w l lw 7 title 'Au', \
'DATA/TRIALH/DipoleMomentMG.txt' u ($1*0.02418):2 w l lw 7 title 'MG'

set output 'DipoleI.eps'
plot 'DATA/TRIALI/DipoleMoment.txt' u ($1*0.02418):2 w l lw 7 title 'Au', \
'DATA/TRIALI/DipoleMomentMG.txt' u ($1*0.02418):2 w l lw 7 title 'MG'
