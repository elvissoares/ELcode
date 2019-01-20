#!/usr/bin/gnuplot
#
# Plotting a color map using the default gnuplot palette
#
# AUTHOR: Elvis Soares

reset

# epslatex
set terminal epslatex size 8cm,6cm color colortext 11 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"
#set terminal epslatex size 12cm,9cm color colortext header "\\newcommand{\\ft}[0]{\\footnotesize}"

set output "silicon-explosive.tex"

# color definitions
load 'palletes/parula.pal'

# Axes
set style line 101 lc rgb '#000000' lt 1
set border back ls 101
#set tics nomirror out scale 0.75
# Grid
set style line 102 lc rgb'#808080' dt 4 lw 1
#set grid back ls 102

unset key

set logscale 
set format '$10^{%T}$'

set xlabel 'time (s)'
set ylabel 'Mass Fraction'

set xr[1e-7:1e-1]
set xtics 1e2

set yr[1e-8:3]
set ytics 1e2

set label 1 '\ft ${}^{4}$He' at 1e-7,7 tc ls 11
set label 2 '\ft ${}^{12}$C' at 5e-7,7 tc ls 12
set label 3 '\ft ${}^{16}$O' at 2.5e-6,7 tc ls 13
set label 4 '\ft ${}^{20}$Ne' at 1.5e-5,7 tc ls 14
set label 5 '\ft ${}^{24}$Mg' at 1e-4,7 tc ls 15
set label 6 '\ft [${}^{28}$Si]' at 9e-4,7 tc ls 16
set label 7 '\ft [${}^{56}$Ni]' at 1e-2,7 tc ls 17

#set rmargin at screen 0.85

plot 'silicon-explosive.dat' w l ls 11 lw 4,\
     'silicon-explosive.dat' u 1:3  w l ls 12 lw 4,\
     'silicon-explosive.dat' u 1:4 w l ls 13 lw 4,\
     'silicon-explosive.dat' u 1:5 w l ls 14 lw 4,\
     'silicon-explosive.dat' u 1:6 w l ls 15 lw 4,\
     'silicon-explosive.dat' u 1:7 w l ls 16 lw 4,\
     'silicon-explosive.dat' u 1:8 w l ls 17 lw 4
     
set output "energy-silicon-explosive.tex"

set yr[1e10:1e23]
set ytics 1e4

set ylabel 'Nuclear Energy Rate (erg/g.s)'

set arrow from 2.5e-6,1e16 to 2.5e-6,1e22 nohead

set label 1 '\ft $\epsilon>0$' at 4e-6,1e17 tc rgb '#000000'
set label 2 '\ft $\epsilon<0$' at 3e-7,1e17 tc rgb '#000000'

f(x) = (x >= 0 ? x : abs(x))

plot 'silicon-explosive.dat' u 1:($9 >= 0 ? $9 : 1/0) w l lc rgb '#808080' lw 4,\
'silicon-explosive.dat' u 1:($9 < 0 ? -$9 : 1/0) w l lc rgb '#808080' dt 3 lw 4

exit
#Necessario para compilar o arquivo diretamente com pdflatex
