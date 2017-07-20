#!/Users/jay/gnuplot/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 30
set output 'PopulationMG.eps'
set xlabel 'Time (femtoseconds)'
a = 8e29
set xrange [0:20]
set pointsize 2
#set yrange [0:($1*0.02418)]
set ylabel 'Population'
plot 'DATA/TRIALA/PopulationMG.txt' u ($1*0.02418):3 w l lw 7 title 'A', \
'DATA/TRIALB/PopulationMG.txt' u ($1*0.02418):3 w l lw 7 title 'B', \
'DATA/TRIALC/PopulationMG.txt' u ($1*0.02418):3 w l lw 7 title 'C', \
'DATA/TRIALD/PopulationMG.txt' u ($1*0.02418):3 w l lw 7 title 'D', \
'DATA/TRIALE/PopulationMG.txt' u ($1*0.02418):3 w l lw 7 title 'E', \
'DATA/TRIALF/PopulationMG.txt' u ($1*0.02418):3 w l lw 7 title 'F', \
'DATA/TRIALG/PopulationMG.txt' u ($1*0.02418):3 w l lw 7 title 'G', \
'DATA/TRIALH/PopulationMG.txt' u ($1*0.02418):3 w l lw 7 title 'H', \
'DATA/TRIALI/PopulationMG.txt' u ($1*0.02418):3 w l lt 0 lw 7 title 'I'
