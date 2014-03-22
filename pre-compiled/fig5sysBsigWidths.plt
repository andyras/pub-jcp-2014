#!/usr/bin/env gnuplot

myfont = 'Arial,12'
set terminal pdfcairo dashed enhanced color size 3,4 font myfont lw 4 dl 1 rounded
set output 'fig5sysBsigWidths.pdf'

d1 = '../data/sysB_V.000140_sigma.0010/tcprob.out'
stats d1 name 'd1' nooutput
d2 = '../data/sysB_V.000140_sigma.0020/tcprob.out'
stats d2 name 'd2' nooutput
d3 = '../data/sysB_V.000140_sigma.0040/tcprob.out'
stats d3 name 'd3' nooutput
d4 = '../data/sysB_V.000140_sigma.0080/tcprob.out'
stats d4 name 'd4' nooutput

set linetype 99 lt 1
set linetype 1 lc rgb 'black' lt 3
set linetype 2 lc rgb 'blue'
set linetype 3 lc rgb 'red' lt 99
set linetype 4 lc rgb 'dark-chartreuse'

set xlabel 'Time (fs)' offset 0,0.8
set ylabel "P_{QD}(t)" norotate
set xr [0:80]
set key right Left top font ',10' reverse samplen 3.2

set border 3

set tics out scale 0.5 nomirror
set xtics 20 offset 0,0.35
set ytics 0.3 offset 0.35
set yrange [0:0.45]

set multiplot title "Model {/Arial-Bold B} On Resonance\nChanging {/Symbol s}_E for V = 3.81 meV"

### first plot: raw data
set origin 0,0.635
set size 1,0.3

set lmargin 11
set rmargin 2

set label "A)" at screen 0.015,0.895 font "Arial-Bold,16"

plot d1 u ($1/41.3414):2 t '{/Symbol s}_E = 0.596 V^*', \
d2 u ($1/41.3414):2 t '{/Symbol s}_E = 2.39 V^*', \
d3 u ($1/41.3414):2 t '{/Symbol s}_E = 4.78 V^*', \
d4 u ($1/41.3414):2 t '{/Symbol s}_E = 9.56 V^*'

set arrow 1 from screen 0,0.625 to screen 1,0.625 nohead

### second plot: normalized data
set origin 0,0.315
set yr [0:1]
set ytics 1
set ylabel "Normalized\nP_{QD}(t) (a.u.)" norotate offset 1.5,0.5

set label "B)" at screen 0.015,0.59 font "Arial-Bold,16"

#plot d1 u ($1/41.3414):($2/d1_max_y) t '{/Symbol s}_E = 27.2 meV', \
#d2 u ($1/41.3414):($2/d2_max_y) t '{/Symbol s}_E = 54.4 meV', \
#d3 u ($1/41.3414):($2/d3_max_y) t '{/Symbol s}_E = 109 meV', \
#d4 u ($1/41.3414):($2/d4_max_y) t '{/Symbol s}_E = 218 meV'
plot d1 u ($1/41.3414):($2/d1_max_y) t '{/Symbol s}_E = 0.596 V^*', \
d2 u ($1/41.3414):($2/d2_max_y) t '{/Symbol s}_E = 2.39 V^*', \
d3 u ($1/41.3414):($2/d3_max_y) t '{/Symbol s}_E = 4.78 V^*', \
d4 u ($1/41.3414):($2/d4_max_y) t '{/Symbol s}_E = 9.56 V^*'

set arrow 1 from screen 0,0.305 to screen 1,0.305 nohead

### third plot: widths and broadening
set label "C)" at screen 0.015,0.275 font "Arial-Bold,16"
set origin 0,0.025
set samples 1000

set lmargin
set bmargin 2

set xr [-0.3:1.1]
set yr [0:1]
set key top left
#set key at 0,1.5

set border 1

unset ytics
set xtics 1 scale 0 out nomirror

set xlabel 'Energy Offset from Bulk Band Edge (eV)'
unset ylabel


set title 'Wavepacket Breadth vs. Single State Broadening'

## Distributions centered around 0
gau(x,x0,s) = exp(-((x-x0)/s)**2/2)
# Lorentzian broadening of a state depends on density of states, coupling.
V = 3.81e-3
Delta = 0.0735/1000*27.211
p = 1/Delta
# normalized Lorentzian
lor(x,x0) = 1/(1 + ((x-x0)/(pi*V**2*p))**2)

plot lor(x, 27.211*0.0364926) + \
lor(x, 27.211*0.0358741) + \
lor(x, 27.211*0.0352556) + \
lor(x, 27.211*0.034637) + \
lor(x, 27.211*0.0340185) + \
lor(x, 27.211*0.0334) + \
lor(x, 27.211*0.0327815) + \
lor(x, 27.211*0.032163) + \
lor(x, 27.211*0.0290704) + \
lor(x, 27.211*0.0284519) + \
lor(x, 27.211*0.0278333) + \
lor(x, 27.211*0.0272148) + \
lor(x, 27.211*0.0265963) + \
lor(x, 27.211*0.0228852) + \
lor(x, 27.211*0.0216481) + \
lor(x, 27.211*0.0210296) + \
lor(x, 27.211*0.0185556) + \
lor(x, 27.211*0.0167) + \
lor(x, 27.211*0.0160815) + \
lor(x, 27.211*0.015463) + \
lor(x, 27.211*0.00865926) + \
lor(x, 27.211*0.00804074) + \
lor(x, 27.211*0.00742222) + \
lor(x,0) axis x1y2 t "Broadened\nQD States" lw 3 lt 3 lc rgb 'gray70', \
gau(x, 27.211*0.015, 2.72e-2) notitle lt 1, \
gau(x, 27.211*0.015, 5.44e-2) notitle lt 2, \
gau(x, 27.211*0.015, 1.09e-1) notitle lt 3, \
gau(x, 27.211*0.015, 2.18e-1) notitle lt 4

unset key
unset ylabel
unset xlabel
unset tics
unset title
set autoscale y


unset multiplot
