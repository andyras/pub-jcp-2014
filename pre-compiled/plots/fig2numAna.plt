#!/usr/bin/env gnuplot

reset

set style data lines

set terminal pdfcairo dashed enhanced color lw 4 size 3,2 font 'Arial,14' rounded
set output 'fig2numAna.pdf'

set tics scale 0.5 out nomirror
set xtics 0.9
set ytics 0.4
set y2tics 0.01 offset -0.5,0
set lmargin 6
set rmargin 6

set xlabel 'Time (fraction of recurrence)' offset 0,1.5
set ylabel 'P_{r}(t)' norotate offset 4

set title ''

set border 11

# change x range just for fitting part of curve
set xr [0.3:0.5]

f(x) = a*exp(-b*x)
a = 1
b = 0.0001
# time constant is 2*pi*V**2/(\hbar*\Delta)
b = 2*3.14159265*0.0005**2/(0.0735/100)*8548.551438
set fit logfile '/dev/null' quiet
fit f(x) '../data/sysA_numana/tcprob.out' u ($1/8548.551438):2 via a

set xrange [0:0.95]
set yrange [0:0.40]
set y2range [0:0.01]

set y2label ''

x1 = 0.02; y1 = 1.02; x2 = 0.10; y2 = 1.02
set arrow 1 from graph x1,y1 to graph x2,y2 filled heads size 0.02,15
set label 1 'FET' center at graph (x1+x2)/2.,0.96
x1 = 0.125; x2 = 0.45
set arrow 2 from graph x1,y1 to graph x2,y2 filled heads size 0.02,15
set label 2 'BET' center at graph (x1+x2)/2.,0.96

set key spacing 1.1

plot '../data/sysA_numana/tcprob.out' u ($1/8548.551438):2 title 'Numerical', \
'../data/sysA_numana/ms_est.out' u ($1/8548.551438):2 lc rgb 'blue' title 'Analytical', \
f(x) lc rgb 'black' title 'Exponential fit', \
'../data/sysA_onestate/tcprob.out' u ($1/8548.551438):2 lc rgb 'slateblue1' title 'Eq. 14, E_{lr} = 0', \
'../data/sysA_incoherent/avg/tcprob_avg.out' u ($1/8548.551438):2 lc rgb 'dark-green' title "Average using Eq. B3\n(right y axis)" axes x1y2

reset
