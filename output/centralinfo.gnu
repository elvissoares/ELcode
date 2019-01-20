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

# epslatex
set terminal epslatex size 12cm,9cm color colortext 12 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"
#set terminal epslatex size 9cm,7cm color colortext standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"

#set terminal epslatex size 12cm,9cm color colortext header "\\newcommand{\\ft}[0]{\\footnotesize}"

set output "centralinfo.tex"

# color definitions
load 'palletes/parula.pal'

unset key

# Axes
set style line 101 lc rgb '#000000' lt 1
set border back ls 101
#set tics nomirror out scale 0.75
# Grid
set style line 102 lc rgb'#808080' dt 4 lw 1
#set grid back ls 102

set mxtics
set mytics

set format '$10^{%g}$'

set ylabel 'Temperature (K)'
set xlabel 'Density (g/cm$^{3}$)'

set yrange [7:9.3]
set xrange [4:9.5]

set xtics 1
set ytics 0.5

# getting slope for text placing
set label 3 '\ft  $\eta = -2.4$' at 6,9.75 rotate by 35 center tc ls 11

set label 1 '\ft  $\tau_{He} = 10^{6}\; yr$' at 4.8,8 rotate by 0 center tc ls 13
set label 2 '\ft  $\tau_{He} = \tau_{hd}$' at 5.5,8.5 rotate by 0 center tc ls 13

set label 3 '\ft  $\tau_{C} = 10^{6}\; yr$' at 8,8.52 rotate by 0 center tc ls 12
set label 4 '\ft  $\tau_{C} = \tau_{hd}$' at 8.5,9.05 rotate by 0 center tc ls 12

plot '../../eos/degenerateline.dat' w l ls 11 lw 3,\
     '../../thermonuclear/reactionrateHe.dat' w l lc rgb '#edb120' dt 3 lw 3,\
     '../../thermonuclear/reactionrateHe.dat' u 1:3 w l lc rgb '#edb120' dt 3 lw 3,\
     '../../thermonuclear/reactionrateC.dat' w l lc rgb '#d95319' dt 3 lw 3,\
     '../../thermonuclear/reactionrateC.dat' u 1:3 w l lc rgb '#d95319' dt 3 lw 3,\
     '../../equilibrium/N100/initial-condition-M1-Ts6.14051e+07-N100.dat' u (log10($2)):(log10($3/3)) w l lc rgb '#a2142f' dt 2 lw 3,\
     'shell-accreted.dat' u 2:3 w l ls 17 lw 3,\
     'shell-central.dat' u 2:3 w l ls 17 lw 3,\
     'shell-accreted0.dat' u 2:3 w l ls 12 lw 3,\
     'shell-central0.dat' u 2:3 w l ls 12 lw 3


exit
#Necessario para compilar o arquivo diretamente com pdflatex
