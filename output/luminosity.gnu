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

set output "luminosity.tex"

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
set format y '$10^{%T}$'

set xlabel 'time (s)'
set ylabel 'Luminosity (erg/s)'

# getting slope for text placing
#set label 2 '\ft  $u_e = 10^{28}$ erg.cm$^{-3}$'         at 9.9,3.8     rotate by  -90 center tc ls 11

#set xr[0:4]

plot 'luminosity.dat' u 1:2  w l ls 11 lw 2


exit
#Necessario para compilar o arquivo diretamente com pdflatex
