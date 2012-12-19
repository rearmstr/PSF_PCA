reset
#
PS=1
PDF=0
# if(PS==1) set term post eps enh rounded color dashed dl 4  'Times' 28
if(PS==1) set term png
if(PS==1) set out "/astro/tutti0/www/mzm/temp/whisker_400.png"
if(PS==1) set size 2,1.22
# if(PS==1) set size square
#
set xlabel 'RA (degree)'
# set xtics out
# set ylabel 'y'
set nologscale x
set nologscale y
set xr [336.16:338.5]
set yr [-16:-14]
set bmargin 4
set tmargin 1
#
set multiplot
set size 1,1.22
set origin 0,0
set lmargin 6
set rmargin 0
set xtics 336.5,0.5,338
set label "data" at 336.4,-14.1
set label "decam--27--41-i-12" at 338,-14.1
load "drawDESchipBound.p"

plot 'whiskerData.dat' using 1:2:3:4 with vectors nohead lt 1 notitle
#
set size 1,1.22
set origin 1,0
set lmargin 0
set rmargin 6
set format y ""
set y2label "DEC (degree)"
set xtics 336.5,0.5,338.5
unset label; set label "residual" at 338.1,-14.1
# load "labelDESchipCenter.p"
load "drawDESchipBound.p"

plot 'whiskerDiff.dat' using 1:2:3:4 with vectors nohead lt 1 lw 1 notitle
#
set nomultiplot
# if(PS==1) set terminal x11
# if(PS==1) set size 1,1
