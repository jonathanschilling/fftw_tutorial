#!/usr/bin/gnuplot -persist

set grid
set xtics 6
set xlabel "grid index l"
set xrange [ -2.0 : 39.0 ]
set ylabel "R / m"
set yrange [ 5.2 : 6.3 ]

plot '../axis_R_symm.dat'          u          1 w lp title   "DFT"   lw 2 ps 2, \
     '../axis_R_halfGrid.dat' u ($0+0.5):1 w lp title "REDFT01" lw 2 ps 2
