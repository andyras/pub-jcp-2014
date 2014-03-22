#!/usr/bin/env gnuplot
myfont = 'Arial,12'
set terminal pdfcairo dashed enhanced color size 3,4 font myfont lw 4 dl 1 rounded
set output 'fig4sysAchangeV.pdf'

d1 = '../data/sysA_V.000020_sigma0.001/tcprob.out'
stats d1 name 'd1' nooutput
d2 = '../data/sysA_V.000040_sigma0.001/tcprob.out'
stats d2 name 'd2' nooutput
d3 = '../data/sysA_V.000060_sigma0.001/tcprob.out'
stats d3 name 'd3' nooutput
d4 = '../data/sysA_V.000080_sigma0.001/tcprob.out'
stats d4 name 'd4' nooutput
d5 = '../data/sysA_V.000100_sigma0.001/tcprob.out'
stats d5 name 'd5' nooutput
set linetype 99 lt 1
set linetype 1 lc rgb 'black' lt 3
set linetype 2 lc rgb 'blue'
set linetype 3 lc rgb 'red' lt 99
set linetype 4 lc rgb 'dark-chartreuse'
set linetype 5 lc rgb 'steelblue'

set samples 1000

## Distributions centered around 0
gau(x,x0,s) = exp(-((x-x0)/s)**2/2)
# Lorentzian broadening of a state depends on density of states, coupling.
Delta = 0.0735/1000*27.211
p = 1/Delta
# not normalized Lorentzian
#lor(x) = V**2/(x**2 + (pi*V**2*p)**2)
# normalized Lorentzian
lor(x,x0,V) = 1/(1 + ((x-x0)/(pi*V**2*p))**2)

set xlabel 'Time (fs)' offset 0,0.8
set ylabel "P_{SS}(t)" norotate
set xr [0:200]
set key right Left top font ',10' reverse samplen 3.2 spacing 0.75

set border 3

set tics out scale 0.5 nomirror
set xtics 200 offset 0,0.35
set ytics 0.2 offset 0.35
set yrange [0:0.25]


set multiplot title "Model {/Arial-Bold A} On Resonance\nChanging V for {/Symbol s}_E = 27.2 meV"


### first plot: raw data
set size 1,0.30
set origin 0,0.635

set lmargin 11
set rmargin 2

set label "A)" at screen 0.015,0.895 font "Arial-Bold,16"

#plot d1 u ($1/41.3414):($2) t 'V = 0.544 meV', \
#d2 u ($1/41.3414):($2) t 'V = 1.09 meV', \
#d3 u ($1/41.3414):($2) t 'V = 1.63 meV', \
#d4 u ($1/41.3414):($2) t 'V = 2.18 meV', \
#d5 u ($1/41.3414):($2) t 'V = 2.72 meV'
plot d1 u ($1/41.3414):($2) t 'V^* = 0.465 meV', \
d2 u ($1/41.3414):($2) t 'V^* = 1.87 meV', \
d3 u ($1/41.3414):($2) t 'V^* = 4.17 meV', \
d4 u ($1/41.3414):($2) t 'V^* = 7.47 meV', \
d5 u ($1/41.3414):($2) t 'V^* = 11.6 meV'

set arrow 1 from screen 0,0.625 to screen 1,0.625 nohead


### second plot: normalized data
set origin 0,0.315

unset key

set yr [0:1]

set ytics 1

set label "B)" at screen 0.015,0.59 font "Arial-Bold,16"
set ylabel "Normalized\nP_{SS}(t) (a.u.)" norotate offset 1.5,0.5

plot d1 u ($1/41.3414):($2/d1_max_y), \
d2 u ($1/41.3414):($2/d2_max_y), \
d3 u ($1/41.3414):($2/d3_max_y), \
d4 u ($1/41.3414):($2/d4_max_y), \
d5 u ($1/41.3414):($2/d5_max_y)

set arrow 1 from screen 0,0.305 to screen 1,0.305 nohead


### third plot: broadenings vs. wavepacket width
set origin 0,0.025

set border 1
set lmargin 2
set rmargin 2

set yr [0:1.25]
set xr [-0.020:0.020]
set xtics 0.02 scale 0

unset ytics
unset ylabel
set xlabel 'Energy Offset from Single State (eV)'
set title 'Wavepacket Breadth vs. Single State Broadening'

set label "C)" at screen 0.015,0.275 font "Arial-Bold,16"

set key at 0.02,1.4

do for [ii=2:30] {set arrow ii from 2e-3*(ii-15),0 to 2e-3*(ii-15),1 nohead lw 0.5 }
plot gau(x, 0, 27.2e-3) lt 3 lc rgb 'gray70' lw 3 t "Wavepacket\nWidth", \
lor(x, 0, 0.544e-3) lt 1 notitle, \
lor(x, 0, 1.09e-3) lt 2 notitle, \
lor(x, 0, 1.63e-3) lt 3 notitle, \
lor(x, 0, 2.18e-3) lt 4 notitle, \
lor(x, 0, 2.72e-3) lt 5 notitle

unset multiplot
