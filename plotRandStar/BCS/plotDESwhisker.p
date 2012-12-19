reset
#
PS=1
PDF=0
# if(PS==1) set term post eps enh rounded color dashed dl 4  'Times' 28
if(PS==1) set term png size 1200,600
if(PS==1) set out "/astro/tutti0/www/mzm/temp_bcs_r_new/whisker_30.png"
#
set xlabel 'x'
# set xtics out
# set ylabel 'y'
set nologscale x
set nologscale y
set xr [-200:8800]
set yr [-200:8800]
# set bmargin 4
# set tmargin 1
#
set multiplot
set size 0.5,1.0
set origin 0,0
set lmargin 6
set rmargin 0
# set xtics 336.5,0.5,338
set label "data" at 0,8500
set label "BCS0511-5448Ar.061218_0755.112" at 2000,8500
load "draw_bcs_chipBound.p"

plot 'whiskerData.dat' using 1:2:3:4 with vectors nohead lt 1 notitle
#
set size 0.5,1.0
set origin 0.5,0
set lmargin 0
set rmargin 6
set format y ""
# set y2label "DEC (degree)"
# set xtics 336.5,0.5,338.5
unset label; set label "residual" at 0,8500
# load "labelDESchipCenter.p"
load "draw_bcs_chipBound.p"

plot 'whiskerDiff.dat' using 1:2:3:4 with vectors nohead lt 1 lw 1 notitle
#
set nomultiplot
# if(PS==1) set terminal x11
# if(PS==1) set size 1,1
