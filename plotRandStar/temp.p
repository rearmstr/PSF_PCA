reset
#
PS=1
PDF=0
if(PS==1) set term post eps enh rounded color dashed dl 4  'Times' 28
if(PS==1) set out 'whisker.ps'
if(PS==1) set size 2,2.5
# if(PS==1) set size square
#
# set xlabel 'x'
# set ylabel 'y'
set nologscale x
set nologscale y
set xr [336.16:338.5]
set yr [-16:-14]
#
load "drawDESchipBound.p"

plot 'LL.dat' using 4:5 with p ps 0.01 lc 3 lw 3 notitle

if(PS==1) set terminal x11
if(PS==1) set size 1,1
