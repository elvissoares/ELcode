#!/usr/bin/gnuplot
#
# Plotting a color map using the default gnuplot palette
#
# AUTHOR: Elvis Soares

reset

# epslatex#
set terminal epslatex size 8cm,6cm color colortext 11 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"

set output "sod-density-profile.tex"

# color definitions
load 'parula.pal'

# Axes
set style line 101 lc rgb '#000000' lt 1
set border back ls 101
#set tics nomirror out scale 0.75
# Grid
set style line 102 lc rgb'#808080' dt 4 lw 1
#set grid back ls 102

unset key

#set logscale y
#set format y '$10^{%T}$'

set xlabel 'Radius (earth radius)'
set ylabel 'Density ($10^{7}$ g/cm$^3$)'

set yr[0.0:1.2]
set xr[0.0:0.8]
#set xtics 0.001

REarth = 6.3781e8

# getting slope for text placing
#set label 2 '\ft  $u_e = 10^{28}$ erg.cm$^{-3}$'         at 9.9,3.8     rotate by  -90 center tc ls 11

plot 'sod-profile-t=0.dat' u ($1/REarth):($2/1e7) w l ls 17 lw 4,\
'sod-profile-medium.dat' u ($1/REarth):($2/1e7) w l ls 12 lw 4,\
'sod-profile-final.dat' u ($1/REarth):($2/1e7) w l ls 13 lw 4

set output "sod-pressure-profile.tex"

set yr[0.0:10]

set ylabel 'Pressure ($10^{23}$ erg/cm$^3$)'

plot 'sod-profile-t=0.dat' u ($1/REarth):($3/1e23) w l ls 17 lw 4,\
     'sod-profile-medium.dat' u ($1/REarth):($3/1e23) w l ls 12 lw 4,\
     'sod-profile-final.dat' u ($1/REarth):($3/1e23) w l ls 13 lw 4

set output "sod-velocity-profile.tex"

set yr[-0.5:4]
set ytics 1
set ylabel 'Velocity ($10^{8}$ cm/s)'
 
plot 'sod-profile-t=0.dat' u ($1/REarth):($4/1e8) w l ls 17 lw 4,\
     'sod-profile-medium.dat' u ($1/REarth):($4/1e8) w l ls 12 lw 4,\
     'sod-profile-final.dat' u ($1/REarth):($4/1e8) w l ls 13 lw 4


exit
#Necessario para compilar o arquivo diretamente com pdflatex
