reset
#
PS=1
if(PS==1) set term post eps enh rounded color dashed dl 4  'Times' 28
if(PS==1) set out 'whisker.ps'
if(PS==1) set size 2
if(PS==1) set size square
#
set xlabel 'x'
set ylabel 'y'
set nologscale x
set nologscale y
# set xr [-19:19]
# set yr [-19:19]
set xr [-200:8800]
set yr [-200:8500]
#
plot 'whiskerEigen.dat' using 1:2:3:4 with vectors nohead lt 1 notitle
#
if(PS==1) set terminal x11
if(PS==1) set size 1,1
