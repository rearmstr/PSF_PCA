#-----------------------------------------------------------------------
reset
#
PS=1
NX=1; NY=2
DX=0.06;DY=0.05;SX=2.0;SY=1.0
set rmargin 0
set lmargin 12
if(PS==1) set term post eps enh rounded color dashed dl 4  'Times' 30
if(PS==1) set out 'D1_D2.ps'
if(PS==1) set size SX*NX+DX*1.5,SY*NY+DY*1.5
#
xmin=0.02
xmax=0.6
f_len=22
# xmax=2.3
set xrange [xmin/f_len:xmax/f_len]
set logscale x
set format y "%3.1t {/Symbol \264}10^{%L}"
#
set xlabel '{/Symbol q} [degree]'
set ylabel 'D_1({/Symbol q}) and D_2({/Symbol q})'
set pointsize 0.5
#
set key top right
set samples 1025
#
#
set multiplot
set size SX,SY
set origin DX,DY
set tmargin 0
set yrange [-8.0E-6:5.0E-6]
# set ytics -1.0E-6,5E-7,5E-7
set xlabel '{/Symbol q} [degree]'
set ylabel 'D_2({/Symbol q})'
set key bottom right
set arrow from xmin/f_len,0 to xmax/f_len,0 nohead lc 3
plot "D1_D2_SNAP_2k_10x10_valSet.dat" u ($1/f_len):6:7 with yerrorbars pt 6  lw 3 lt 1 lc 1 t "validation set 2k", \
     "D1_D2_SNAP_2k_10x10.dat"        u ($1/f_len):6:7 with yerrorbars pt 4  lw 3 lt 1 lc 3 t "data 2k", \
     "D1_D2_SNAP_4k_10x10_valSet.dat" u ($1/f_len):6:7 with yerrorbars pt 4  lw 3 lt 1 lc 4 t "validation set 4k", \
     "D1_D2_SNAP_4k_10x10.dat"        u ($1/f_len):6:7 with yerrorbars pt 4  lw 3 lt 1 lc 0 t "data 4k"
#
set origin DX,SY+DY
set bmargin 0
set tmargin 1
set yrange [-5.0E-6:5.0E-6]
# set ytics -5.0E-8,5E-8
set noxlabel; set format x ""
set ylabel 'D_1({/Symbol q})'
set key top right
set arrow from xmin/f_len,0 to xmax/f_len,0 nohead lc 3
plot "D1_D2_SNAP_2k_10x10_valSet.dat" u ($1/f_len):4:5 with yerrorbars pt 4  lw 3 lt 1 lc 1 t "validation set 2k", \
     "D1_D2_SNAP_2k_10x10.dat"        u ($1/f_len):4:5 with yerrorbars pt 4  lw 3 lt 1 lc 3 t "data 2k", \
     "D1_D2_SNAP_2k_10x10_1-8thN_valSet.dat" u ($1/f_len):4:5 with yerrorbars pt 4  lw 3 lt 1 lc 4 t "validation set 4k", \
     "D1_D2_SNAP_2k_10x10_1-8thN.dat"        u ($1/f_len):4:5 with yerrorbars pt 4  lw 3 lt 1 lc 0 t "data 4k"
#
unset multiplot
if(PS==1) set terminal x11
if(PS==1) set size 1,1
