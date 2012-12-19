reset
#
PS=1
PDF=0
if(PS==1) set term post eps enh rounded color dashed dl 4  'Times' 28
if(PS==1) set out 'whisker.ps'
if(PS==1) set size 2,1.5
# if(PS==1) set size square
#
# set xlabel 'x'
# set ylabel 'y'
set nologscale x
set nologscale y
set xr [-0.32:0.32]
set yr [-0.32:0.32]
set xtics -0.3,0.1,0.2
set parametric
set tr [0.0:2*pi]
rin=0.12
rout=0.3
set bmargin 3
set tmargin 0.6
#
set multiplot
set size 1,1.5
set origin 0,0
set lmargin 4
set rmargin 0
load "drawChipBound.p"
set arrow 1 from -0.02,0 to 0,0 nohead lw 2 lc 1
set label at 0.01,0 "4% e"
# set label at -0.03,0 "0.3 X"
set label at 0.18,0.29 "exposure 16"

plot 'whiskerData.dat' using 1:2:3:4 with vectors nohead lt 1 notitle, \
     rout*sin(t),rout*cos(t) w l lt 1 lc 3 notitle, \
     rin*sin(t),rin*cos(t) w l lt 1 lc 3 notitle
#
set size 1,1.5
set origin 1,0
set lmargin 0
set rmargin 4
set format y ""
load "drawChipBound.p"
set arrow 2 from -0.02,0 to 0,0 nohead lw 2 lc 1
unset label; set label at 0.01,0 "1% e"
# unset label; set label at -0.02,0 "3 X"
set label at 0.18,0.29 "residual"
# set object 50 circle at 0,0 size 0.12
# set object 51 circle at 0,0 size 0.3

# plot "focalPlaneOuter.dat" w circles lc rgb "gold" fs transparent solid 0.1 noborder notitle
#      "focalPlaneInner.dat" w circles lc rgb "white" fs solid 1.0 noborder notitle
plot 'whiskerDiff.dat' using 1:2:3:4 with vectors nohead lt 1 lw 1 notitle, \
     rout*sin(t),rout*cos(t) w l lt 1 lc 3 notitle, \
     rin*sin(t),rin*cos(t) w l lt 1 lc 3 notitle
#     'whiskerDiff_BS-BSNN_16.dat' using 1:2:3:4 with vectors nohead lt 1 lc 3 notitle
# plot 'whiskerEigen.dat' using 1:2:3:4 with vectors nohead lt 1 notitle
#
set nomultiplot
if(PS==1) set terminal x11
if(PS==1) set size 1,1
