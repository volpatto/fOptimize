set terminal postscript eps enhanced color 
set encoding utf8
set output 'f1.eps'

#set terminal latex

#set key left top
set key rmargin right vertical
set xlabel 'x'
set ylabel 'y'
#set format x "%.3g
#set tics font ", 10"

set autoscale
set xrange [-100:100]
set yrange [-100:100]

set grid
set view 30,60,1
set contour base
set cntrparam levels auto 20

splot x**2.0+y**2.0 + 2.0*x*y
