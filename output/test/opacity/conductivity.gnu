#!/usr/bin/gnuplot
#
# Plotting a color map using the default gnuplot palette
#
# AUTHOR: Elvis Soares

reset

# epslatex#
set terminal epslatex size 12cm,9cm color colortext 12 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"

set output "conductivity.tex"

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

set xlabel '$\log{\rho}$ (g cm${}^{-3}$)'
set ylabel '$\log{\lambda_e}$ (erg cm${}^{-1}$ s${}^{-1}$ K${}^{-1}$)'

set yr[11:19]
#set xtics 0.001

# getting slope for text placing
set label 1 '\ft  $T = 10^9$ K' at -2.0,16.6 center
set label 2 '\ft  $10^8$ K' at -2.0,14.3 center
set label 3 '\ft  $10^7$ K' at -2.0,12 center

plot 'thermalconductivity-carbon.dat' u 1:4 w l ls 11 lw 4,\
     'thermalconductivity-carbon.dat' u 1:5 w l ls 11 lw 4,\
     'thermalconductivity-carbon.dat' u 1:6 w l ls 11 lw 4,\

exit
#Necessario para compilar o arquivo diretamente com pdflatex
