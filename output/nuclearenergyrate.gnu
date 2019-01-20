#!/usr/bin/gnuplot
#
# Plotting a color map using the default gnuplot palette
#
# AUTHOR: Elvis Soares

reset

# wxt
#set terminal wxt size 420,420 enhanced font 'Verdana,10' persist
# png
#set terminal pngcairo size 420,420 enhanced font 'Verdana,10'

# epslatex#
set terminal epslatex size 9cm,7cm color colortext standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"
#set terminal epslatex size 9cm,7cm color colortext header "\\newcommand{\\ft}[0]{\\footnotesize}"
#set terminal epslatex size 12cm,9cm color colortext header "\\newcommand{\\ft}[0]{\\footnotesize}"

set output "nuclearenergyrate.tex"

# color definitions
load 'palletes/mypalette.pal'

# Axes
set style line 101 lc rgb '#000000' lt 1
set border back ls 101
#set tics nomirror out scale 0.75
# Grid
set style line 102 lc rgb'#808080' dt 4 lw 1
#set grid back ls 102

unset key

set logscale y
set format y '$10^{%T}$'

set xlabel 'Time (s)'
set ylabel 'Nuclear Energy Rate (erg/g.s)'

#set xr[3.7:3.75]
#set xtics 0.001

# getting slope for text placing
#set label 2 '\ft  $u_e = 10^{28}$ erg.cm$^{-3}$'         at 9.9,3.8     rotate by  -90 center tc ls 11

plot 'nuclearenergyrate.dat' u 1:($2) w l ls 1 lw 3,\
     'nuclearenergyrate.dat' u 1:($11) w l ls 2 lw 3,\
     'nuclearenergyrate.dat' u 1:($21) w l ls 3 lw 3,\
     'nuclearenergyrate.dat' u 1:($31) w l ls 4 lw 3,\
     'nuclearenergyrate.dat' u 1:($41) w l ls 5 lw 3,\
     'nuclearenergyrate.dat' u 1:($51) w l ls 6 lw 3,\
     'nuclearenergyrate.dat' u 1:($61) w l ls 7 lw 3,\
     'nuclearenergyrate.dat' u 1:($71) w l ls 8 lw 3,\
     'nuclearenergyrate.dat' u 1:($81) w l ls 9 lw 3,\
     'nuclearenergyrate.dat' u 1:($91) w l ls 10 lw 3,\
     'nuclearenergyrate.dat' u 1:($101) w l ls 11 lw 3,\
     'nuclearenergyrate.dat' u 1:($111) w l ls 12 lw 3,\
     'nuclearenergyrate.dat' u 1:($121) w l ls 12 lw 3,\
     'nuclearenergyrate.dat' u 1:($131) w l ls 12 lw 3,\
     'nuclearenergyrate.dat' u 1:($141) w l ls 12 lw 3

exit
#Necessario para compilar o arquivo diretamente com pdflatex
