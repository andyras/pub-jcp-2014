#!/usr/bin/env gnuplot
# plot of bridge vs. no bridge

reset

myfont = 'Arial,14'

set terminal pdf dashed enhanced size 6,3 font myfont lw 3 rounded
set output 'fig6bridgeNoBridge.pdf'

set linetype 99 lt 1
set linetype 1 lc rgb 'red'
set linetype 2 lc rgb 'blue'
set linetype 3 lc rgb 'black'
set linetype 4 lc rgb 'dark-chartreuse'
set linetype 5 lc rgb 'steelblue'
set linetype 6 lc rgb 'red' lt 2
set linetype 7 lc rgb 'blue' lt 3
set linetype 8 lc rgb 'black' lt 4
set linetype 9 lc rgb 'dark-chartreuse' lt 5
set linetype 10 lc rgb 'steelblue' lt 1
set linetype 11 lc rgb 'red' lt 3

set style data lines

set terminal pdf dashed enhanced size 3,2 font myfont lw 3

set xlabel 'Time (fs)' offset 0,2
set ylabel 'P_{QD}(t)' norotate offset 4
set border 3
set tics scale 0.5 nomirror out
set ytics 0.4 offset 0.6
set xtics 100 offset 0,0.4
set xr [0:100]
set yr [0:0.4]

# center line
xpos=0.28
set arrow 1 from graph xpos,0 to graph xpos,1 nohead lc rgb 'black' lw 0.5

set key font ',12'

plot '../data/sysB_m0.015_s0.001/tcprob.out' u ($1/41.3414):2 notitle, \
'../data/sysC_g0.00_b0.03_m0.015_s0.001/tcprob.out' u ($1/41.3414):2 t 'Model C: V_{lb} = V_{br}', \
'../data/sysC_g0.00_b0.03_m0.015_s0.001_Vkb0.001_Vbc0.0004/tcprob.out' u ($1/41.3414):2 t 'Model C: V_{lb} > V_{br}', \
'../data/sysC_g0.00_b0.03_m0.015_s0.001_Vkb0.001_Vbc0.004/tcprob.out' u ($1/41.3414):2 t 'Model C: V_{lb} < V_{br}', \
2 lt 1 t 'Model B'
