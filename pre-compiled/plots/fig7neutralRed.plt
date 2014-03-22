#!/usr/bin/env gnuplot
# ---- plot of vibrational states on bridge reduced vs. neutral ---- #

myfont = 'Arial,14'

set terminal pdf dashed enhanced size 6,3 font myfont lw 3 rounded

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

set output 'fig7neutralRed.pdf'

set style data lines

set size 0.45,0.5
unset key
set tics scale 0.5 nomirror out
set xtics 20
set tmargin 2
set bmargin 1.5
set autoscale y
set border 3

set multiplot

set origin 0.55,0.5
set xr [0:80]
set xtics 80 offset 0,0.4
set ytics 1 offset 0.6
set xlabel 'Time (fs)' offset 0,2
set title '{/Arial-Bold B)} P_{neutral,ν}(t)    (E_b - {/Symbol m}_E) = 0.41 eV' offset 0,-1
plotstr="plot '<paste ../data/sysC_g2.00_b0.03_m0.015_s0.001/vibprob_qd.out ../data/sysC_g2.00_b0.03_m0.015_s0.001/vibprob_bu.out' u ($1/41.3414):($2+$14) ls 1"
do for [ii=3:12] {
 plotstr=sprintf("%s, '' u ($1/41.3414):($%d+$%d) ls (%d-1)", plotstr, ii, ii+12, ii)
}
eval plotstr

set origin 0.125,0.5
set xr [0:40]
set xtics 40
set ytics 0.010
set title '{/Arial-Bold A)} P_{ionic,ν}(t)    (E_b - {/Symbol m}_E) = 0.41 eV' offset 0,-1
plot for [ii=2:12] '../data/sysC_g2.00_b0.03_m0.015_s0.001/vibprob_br.out' u ($1/41.3414):ii lt (ii-1)

set origin 0.55,0.0
set xr [0:80]
set xtics 80
set ytics 1
set title '{/Arial-Bold D)} P_{neutral,ν}(t)    ({/Symbol m}_E - E_b) = 0.41 eV' offset 0,-1
plotstr="plot '<paste ../data/sysC_g2.00_b0.015_m0.03_s0.001/vibprob_qd.out ../data/sysC_g2.00_b0.015_m0.03_s0.001/vibprob_bu.out' u ($1/41.3414):($2+$14) ls 1"
do for [ii=3:12] {
 plotstr=sprintf("%s, '' u ($1/41.3414):($%d+$%d) ls (%d-1)", plotstr, ii, ii+12, ii)
}
eval plotstr

y1 = 0.475
set arrow 1 from screen 0.16,y1 to screen 1,y1 nohead

set origin 0.125,0.0
set xr [0:40]
set xtics 40
set ytics 0.02
set title '{/Arial-Bold C)} P_{ionic,ν}(t)    ({/Symbol m}_E - E_b) = 0.41 eV' offset 0,-1
plot for [ii=2:12] '../data/sysC_g2.00_b0.015_m0.03_s0.001/vibprob_br.out' u ($1/41.3414):ii lt (ii-1)

# dummy plot for key
unset xlabel
unset title
unset border
unset tics
set size 1,1
set origin -0.1,0
#set key center left maxcols 1 title 'i' box font ',16'
set key center left box height 0.5 width 0.5
#set key maxcols 1 box font ',14' at graph 0.11,0.80 height 0.5 width 0.5
set yr [0:1]
plot for [ii=0:10] 2 t 'ν='.ii lt (ii+1)

unset multiplot
