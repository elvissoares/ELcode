#!/usr/bin/gnuplot
#
# Plotting a color map using the default gnuplot palette
#
# AUTHOR: Elvis Soares

reset


# epslatex
#set terminal epslatex size 9cm,6.75cm color colortext 10 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"
set terminal epslatex size 12cm,9cm color colortext 12 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"
set output "shellinfo.tex"

# jpeg
#set terminal jpeg size 480,360 enhanced font 'Verdana,10'
#set output "shellinfo.jpg"

# color definitions
load '../palletes/parula.pal'

unset key

# Axes
set style line 101 lc rgb '#000000' lt 1
set border back ls 101
#set tics nomirror out scale 0.75
# Grid
set style line 102 lc rgb'#808080' dt 3 lw 1
set grid back ls 102

set mxtics
set mytics

set format '$%g$'

set ylabel 'Temperature (K)'
set xlabel 'Density (g/cm$^{3}$)'

set yrange [6.:8.7]
set xrange [4:9.5]

set xtics 1
set ytics 0.5

# getting slope for text placing
set label 3 '\ft  $\eta = -2.4$' at 6,9.75 rotate by 35 center tc ls 11

set label 1 '\ft  $\tau_{He} = 10^{6}\; yr$' at 4.8,8 rotate by 0 center tc ls 13
set label 2 '\ft  $\tau_{He} = \tau_{hd}$' at 5.5,8.5 rotate by 0 center tc ls 13

set label 3 '\ft  $\tau_{C} = 10^{6}\; yr$' at 8,8.52 rotate by 0 center tc ls 12
set label 4 '\ft  $\tau_{C} = \tau_{hd}$' at 8.5,9.05 rotate by 0 center tc ls 12

plot '../thermonuclear/reactionrateHe.dat' w l lc rgb '#edb120' dt 3 lw 3,\
     '../thermonuclear/reactionrateHe.dat' u 1:3 w l lc rgb '#edb120' dt 3 lw 3,\
     '../thermonuclear/reactionrateC.dat' w l lc rgb '#d95319' dt 3 lw 3,\
     '../thermonuclear/reactionrateC.dat' u 1:3 w l lc rgb '#d95319' dt 3 lw 3,\
     'shell-accreted.dat' u 2:3 w l ls 17 lw 3,\
     'shell-central.dat' u 2:3 w l ls 17 lw 3,\
     'WD-profile-M1.08-N108.dat' u (log10($1)):(log10($2)) w l lc rgb '#a2142f' dt 3 lw 3,\
     '../CaseA-shellhelium.dat' u (log10($2)):(log10($3)) w l ls 17 lw 3,\
     '/home/elvis/Dropbox/Doutorado/Delayed Thermalization/code/output/CaseA-profile-N116.dat' u (log10($2)):(log10($3)) w l lc rgb '#a2142f' dt 3 lw 3


exit
#Necessario para compilar o arquivo diretamente com pdflatex
