#-----------------------------------------------------------------------
reset
#
PS=1
if(PS==1) set term post eps enh rounded color dashed dl 4  'Times' 30
if(PS==1) set out '2pt.ps'
if(PS==1) set size 2
#
# set xrange [0:37]
# set yrange [-0.001:0.001]
# set logscale y
set format y "%3.1t {/Symbol \264}10^{%L}"
#
set xlabel 'r (in degree)'
set ylabel 'Ellipticity 2pt correlation function'
set pointsize 2.0
#
set key top right
set samples 1025
#
# plot "2ptFuncResidual.dat" t "residual" with yerrorbars pt 4
plot "2ptFuncResidual.dat" t "residual" with p pt 4, \
     "2ptFuncData.dat" with p t "data"
#     "2ptFuncRecon.dat" t "reconstruction" with p
#
# if(PS==1) set terminal x11
# if(PS==1) set size 1,1
