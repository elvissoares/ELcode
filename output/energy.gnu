#!/usr/bin/gnuplot
#
# Plotting a color map using the default gnuplot palette
#
# AUTHOR: Elvis Soares

reset

# epslatex#
set terminal epslatex size 12cm,9cm color colortext 12 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"

set output "energy.tex"

# color definitions
load 'palletes/parula.pal'

# Axes
set style line 101 lc rgb '#000000' lt 1
set border back ls 101
#set tics nomirror out scale 0.75
# Grid
set style line 102 lc rgb'#808080' dt 4 lw 1
#set grid back ls 102

set key horiz
set key center top

#set logscale y
#set format x '$10^{%T}$'

set xlabel 'time (s)'
set ylabel 'Energy ($10^{51}$ erg)'

# getting slope for text placing
#set label 2 '\ft  $u_e = 10^{28}$ erg.cm$^{-3}$'         at 9.9,3.8     rotate by  -90 center tc ls 11

plot 'energy.dat' u 1:($2/1e51) ti '\ft $K$' w l ls 15 lw 3,\
     'energy.dat' u 1:($3/1e51) ti '\ft $V_{g}$' w l ls 11 lw 3,\
     'energy.dat' u 1:($4/1e51) ti '\ft $U_{i}$' w l ls 13 lw 3,\
     'energy.dat' u 1:($5/1e51) ti '\ft $Q_{n}$' w l ls 12 lw 3,\
     'energy.dat' u 1:($6/1e51) ti '\ft $E$' w l ls 14 lw 3  
     
set output "energy-total.tex"
unset key
plot 'energy.dat' u 1:($6/1e51) w l ls 14 lw 3  

exit
#Necessario para compilar o arquivo diretamente com pdflatex
