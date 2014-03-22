#!/usr/bin/env gnuplot
# ---- plot varying reorganization, sigma, and normalized version ---- #

reset

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

set terminal pdf dashed enhanced rounded size 3,2 font 'Arial,10' lw 2.5
set style data lines

set output 'fig8reorg.pdf'

set xlabel 'Time (fs)' offset 0,1.9
set ylabel 'P_{r}(t)' offset 6 norotate

set border 3

set tics scale 0.5 nomirror out
set xtics 80 offset 0,0.5
set ytics offset 0.5

set xrange [:80]

set key spacing 0.45

set multiplot layout 2,2
# E_b > mu
set yrange [0:0.25]
set ytics 0.25
set title 'A)  E_b above {/Symbol m}_{E}' offset 0,-0.5
#set key out right center title '{/Symbol l} (eV)' font 'Arial,6' samplen 2.5
set key at graph 0.7,1 title '{/Symbol l} (eV)' font 'Arial,6' samplen 2.5
set tmargin 1.5
set bmargin 2
set lmargin 5
set rmargin 2
plot '../data/sysC_g0.00_b0.03_m0.015_s0.001/tcprob.out' u ($1/41.3414):2 t '0.00', \
     '../data/sysC_g2.00_b0.03_m0.015_s0.001/tcprob.out' u ($1/41.3414):2 t '0.40', \
     '../data/sysC_g2.50_b0.03_m0.015_s0.001/tcprob.out' u ($1/41.3414):2 t '0.63'

# E_b < mu
set yrange [0:0.25]
set ytics 0.25
set title 'B)  E_b below {/Symbol m}_{E}' offset 0,-0.5
set key at graph 0.7,0.4 title '{/Symbol l} (eV)' font 'Arial,6' samplen 2.5
plot '../data/sysC_g0.00_b0.015_m0.03_s0.001/tcprob.out' u ($1/41.3414):2 t '0.00', \
     '../data/sysC_g2.00_b0.015_m0.03_s0.001/tcprob.out' u ($1/41.3414):2 t '0.40', \
     '../data/sysC_g2.50_b0.015_m0.03_s0.001/tcprob.out' u ($1/41.3414):2 t '0.63'

# E_b > mu, normalized
stats '../data/sysC_g0.00_b0.015_m0.03_s0.001/tcprob.out' u ($1/41.3414):2 name 'g1' nooutput
stats '../data/sysC_g2.00_b0.015_m0.03_s0.001/tcprob.out' u ($1/41.3414):2 name 'g2' nooutput
stats '../data/sysC_g2.50_b0.015_m0.03_s0.001/tcprob.out' u ($1/41.3414):2 name 'g3' nooutput

set yrange [0:1]
set ytics 1
set ylabel 'P_{r}(t)/P_{r,max}' offset 3 rotate
set title 'C)  E_b below {/Symbol m}_{E}' offset 0,-0.5
set key at graph 0.7,0.4 title '{/Symbol l} (eV)' font 'Arial,6' samplen 2.5
plot '../data/sysC_g0.00_b0.015_m0.03_s0.001/tcprob.out' u ($1/41.3414):($2/g1_max_y) t '0.00', \
     '../data/sysC_g2.00_b0.015_m0.03_s0.001/tcprob.out' u ($1/41.3414):($2/g2_max_y) t '0.40', \
     '../data/sysC_g2.50_b0.015_m0.03_s0.001/tcprob.out' u ($1/41.3414):($2/g3_max_y) t '0.63'

# E_b > mu, change sigma
set yrange [0:0.25]
set ytics 0.25
set ylabel 'P_{r}(t)' offset 6 norotate
#set key out right title '{/Symbol s}_E (meV)' font 'Arial,6' samplen 2.5
set key at graph 0.7,0.5 title '{/Symbol s}_E (meV)' font 'Arial,6' samplen 2.5
set title 'D)  E_b below {/Symbol m}_E, changing {/Symbol s}_E' offset 0,-0.5
plot '../data/sysC_g2.00_b0.015_m0.03_s0.001/tcprob.out' u ($1/41.3414):2 t '27' lc rgb 'blue', \
     '../data/sysC_g2.00_b0.015_m0.03_s0.002/tcprob.out' u ($1/41.3414):2 t '54' lc rgb 'blue', \
     '../data/sysC_g2.00_b0.015_m0.03_s0.004/tcprob.out' u ($1/41.3414):2 t '108' lc rgb 'blue'
unset multiplot

reset
