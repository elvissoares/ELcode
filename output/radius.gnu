#!/usr/bin/gnuplot
#
# Plotting a color map using the default gnuplot palette
#
# AUTHOR: Elvis Soares

reset

# epslatex#
set terminal epslatex size 12cm,9cm color colortext 12 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"

set output "radius.tex"

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

set xlabel 'time (seconds)'
set ylabel 'Radius (earth radius)'

set yr[0.03:3]
set xr[0.02:4]
#set xtics 0.001

REarth = 6.3781e8

# getting slope for text placing
#set label 2 '\ft  $u_e = 10^{28}$ erg.cm$^{-3}$'         at 9.9,3.8     rotate by  -90 center tc ls 11

plot 'radius.dat' u 1:($2/REarth) w l ls 1 lw 3,\
     'radius.dat' u 1:($3/REarth) w l ls 1 lw 3,\
     'radius.dat' u 1:($4/REarth) w l ls 1 lw 3,\
     'radius.dat' u 1:($5/REarth) w l ls 1 lw 3,\
     'radius.dat' u 1:($6/REarth) w l ls 1 lw 3,\
     'radius.dat' u 1:($7/REarth) w l ls 1 lw 3,\
     'radius.dat' u 1:($8/REarth) w l ls 1 lw 3,\
     'radius.dat' u 1:($9/REarth) w l ls 1 lw 3,\
     'radius.dat' u 1:($10/REarth) w l ls 1 lw 3,\
     'radius.dat' u 1:($11/REarth) w l ls 2 lw 3,\
     'radius.dat' u 1:($12/REarth) w l ls 3 lw 3,\
     'radius.dat' u 1:($13/REarth) w l ls 3 lw 3,\
     'radius.dat' u 1:($14/REarth) w l ls 3 lw 3,\
     'radius.dat' u 1:($15/REarth) w l ls 3 lw 3,\
     'radius.dat' u 1:($16/REarth) w l ls 3 lw 3,\
     'radius.dat' u 1:($17/REarth) w l ls 3 lw 3,\
     'radius.dat' u 1:($18/REarth) w l ls 3 lw 3,\
     'radius.dat' u 1:($19/REarth) w l ls 3 lw 3,\
     'radius.dat' u 1:($20/REarth) w l ls 3 lw 3,\
     'radius.dat' u 1:($21/REarth) w l ls 3 lw 3,\
     'radius.dat' u 1:($31/REarth) w l ls 4 lw 3,\
     'radius.dat' u 1:($41/REarth) w l ls 5 lw 3,\
     'radius.dat' u 1:($51/REarth) w l ls 6 lw 3,\
     'radius.dat' u 1:($61/REarth) w l ls 7 lw 3,\
     'radius.dat' u 1:($71/REarth) w l ls 8 lw 3,\
     'radius.dat' u 1:($81/REarth) w l ls 9 lw 3,\
     'radius.dat' u 1:($91/REarth) w l ls 10 lw 3,\
     'radius.dat' u 1:($101/REarth) w l ls 11 lw 3,\
     'radius.dat' u 1:($109/REarth) w l ls 12 lw 3,\
     'radius.dat' u 1:($110/REarth) w l ls 12 lw 3,\
     'radius.dat' u 1:($111/REarth) w l ls 12 lw 3,\
     'radius.dat' u 1:($112/REarth) w l ls 12 lw 3,\
     'radius.dat' u 1:($113/REarth) w l ls 12 lw 3,\
     'radius.dat' u 1:($114/REarth) w l ls 12 lw 3,\
     'radius.dat' u 1:($115/REarth) w l ls 12 lw 3,\
     'radius.dat' u 1:($116/REarth) w l ls 12 lw 3,\
     'radius.dat' u 1:($117/REarth) w l ls 12 lw 3

exit
#Necessario para compilar o arquivo diretamente com pdflatex
