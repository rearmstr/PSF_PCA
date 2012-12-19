reset
#
PS=1
PDF=0
if(PS==1) set term post eps enh rounded color dashed dl 4  'Times-bold' 32
if(PS==1) set out 'histo.ps'
if(PS==1) set size 2
# if(PS==1) set size square
#
set xlabel 'rms of e (per component)'
set ylabel ' '
# set logscale x
# set logscale y
# set xr [1.0E-4:1.5E-2]
# set yr [1.0E-1:2.0E3]
#
set parametric
plot 'histo.dat' u 1:2 w histeps lw 3 notitle
#
if(PS==1) set terminal x11
if(PS==1) set size 1,1
