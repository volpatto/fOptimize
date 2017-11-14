set terminal postscript eps enhanced color 
set encoding utf8
set output 'f3sit.eps'

#set terminal latex

#set key left top
set key rmargin right vertical
set xlabel 'x'
set ylabel 'y'
#set format x "%.3g
#set tics font ", 10"

set autoscale
set xrange [-4:3]
#set yrange [-10:]

set grid

filename = "secant3.dat"
data1 = 2
data2 = 3

f1(x) = 3.0*x**4.0 - 4.0*x**3.0 + 1.0
df1(x) = 12.0*x**3.0 - 12.0*x**2.0
ddf1(x) = 36.0*x**2.0 - 24.0*x
f2(x) = (1.0/4.0)*x**4 - (5.0/3.0)*x**3.0 - 6.0*x**2.0 + 19.0*x - 7.0
df2(x) = x**3.0 - 5.0*x**2.0 - 12.0*x + 19.0
ddf2(x) = 3.0*x**2.0 - 10.0*x - 12.0

set style data points
plot f2(x) title "f(x)" lw 2, filename using data1:data2 title "Iterações" with linespoints pt 7
