#!/usr/bin/gnuplot
#
# Plotting a color map using the default gnuplot palette
#
# AUTHOR: Elvis Soares

reset


# epslatex
set terminal epslatex size 9cm,6.75cm color colortext 10 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"
#set terminal epslatex size 12cm,9cm color colortext 10 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"

set output "densityprofile.tex"

# color definitions
load '../palletes/parula.pal'

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

set format '$%g$'

set ylabel 'Density (g/cm${}^3$)'
set xlabel '$M_r$ ($M_\odot$)'

set yrange [5:8]
set xrange [0.0:1.2]

set xtics 0.2
set ytics 1

# getting slope for text placing
set label 3 '\ft  $\eta = -2.4$' at 6,9.75 rotate by 35 center tc ls 11

plot 'WD-profile-M1.08-N108.dat' u 1:(log10($2)) w l ls 12 lw 4,\
     '../CaseA-profile-N116.dat' u 1:(log10($2)) w l lc rgb '#a2142f' dt 3 lw 4


exit
#Necessario para compilar o arquivo diretamente com pdflatex
