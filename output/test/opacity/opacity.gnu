#!/usr/bin/gnuplot
#
# Plotting a color map using the default gnuplot palette
#
# AUTHOR: Elvis Soares

reset

# epslatex#
set terminal epslatex size 12cm,9cm color colortext 11 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"

set output "opacity-carbon.tex"

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

set logscale y
set format y '$%T$'

set xlabel '$\log{\rho}$ (g cm${}^{-3}$)'
set ylabel '$\log{\kappa}$ (cm${}^{2}$ g${}^{-1}$)'

set xr[-5:12]
#set xtics 0.001

# getting slope for text placing
set label 1 '\ft  $\kappa_{cond}$' at 0.0,1e8 center
set label 2 '\ft  $\kappa_{rad}$' at 10.0,5e6 center
set label 3 '\ft  $\kappa_{tot}$' at 1.0,5 center

plot 'opacity-carbon-logT8.dat' u 1:2 w l lc rgb '#0072bd' dt 2 lw 3,\
     'opacity-carbon-logT8.dat' u 1:3 w l lc rgb  '#d95319' dt 2 lw 3,\
     'opacity-carbon-logT8.dat' u 1:4 w l ls 17 lw 4

exit
#Necessario para compilar o arquivo diretamente com pdflatex
