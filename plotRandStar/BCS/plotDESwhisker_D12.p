reset
#
# whisker and D1, D2 on the same page.
#
PS=1
PDF=0
# for D1, D2
NX=1; NY=2
DX=0.06;DY=0.01;SX=0.9;SY=0.22
D12_dy=SY*NY+DY*1.5            # D1, D2 plot y size
# if(PS==1) set term post eps enh rounded color dashed dl 4  'Times' 28
if(PS==1) set term png size 1100,1000
if(PS==1) set out "/astro/tutti0/www/mzm/11t_D12_20Eigen_tol-7/whisker_3.png"
# if(PS==1) set size 2,3.22+D12_dy
#
# whisker plot
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
set size 0.5,0.55
set origin 0,D12_dy
set lmargin 5
set rmargin 0
set xtics 336.5,0.5,338
set label "data" at 336.4,-14.1
set label "decam--18--38-i-12" at 338,-14.1
load "drawDESchipBound.p"

plot 'whiskerData.dat' using 1:2:3:4 with vectors nohead lt 1 notitle
#
set size 0.5,0.55
set origin 0.5,D12_dy
set lmargin 0
set rmargin 5
set format y ""
set y2label "DEC (degree)"
set xtics 336.5,0.5,338.5
unset label; set label "residual" at 338.1,-14.1
# load "labelDESchipCenter.p"
load "drawDESchipBound.p"

plot 'whiskerDiff.dat' using 1:2:3:4 with vectors nohead lt 1 lw 1 notitle
#
# D1, D2 plot
#
unset y2label
set rmargin 1
set lmargin 8
xmin=0
xmax=2.3
set xrange [xmin:xmax]
set xtics xmin,0.5,xmax
# set format y "%3.1t {/Symbol \264}10^{%L}"
set format y "%8.1e"
#
# set xlabel '{/Symbol q} [degree]'
set xlabel 'separation r [degree]'
set pointsize 0.5
#
set key top right
#
#
set size SX,SY
set origin DX,DY
set tmargin 0
set xlabel 'separation r [degree]'
set format x "%4.2f"
set yrange [-1.0E-6:1.0E-6]
set ytics -1.0E-6,5E-7,5E-7
set ylabel 'D_2'
set arrow from xmin,0 to xmax,0 nohead
plot "D1_D2.dat" u 1:6:7 with yerrorbars pt 6  lw 1 lt 3 notitle
#
set origin DX,SY+DY
set bmargin 0
set tmargin 1
set yrange [-1.0E-7:1.0E-7]
set ytics -5.0E-8,5E-8
set noxlabel; set format x ""
set ylabel 'D_1'
set arrow from xmin,0 to xmax,0 nohead
plot "D1_D2.dat" u 1:4:5 with yerrorbars pt 4  lw 1 lt 1 notitle
#
unset multiplot
# if(PS==1) set terminal x11
# if(PS==1) set size 1,1
