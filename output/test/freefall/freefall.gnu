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
set terminal epslatex size 8cm,6cm color colortext 11 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"
#set terminal epslatex size 12cm,9cm color colortext 11 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"

# color definitions
load 'parula.pal'

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

# getting slope for text placing
#set label 1 '\ft  $\rho_0 = 10^7$ g/cm$^3$' at 0.001,0.5 rotate by 0 center tc ls 101
#set label 2 '\ft  $R_0 = 10^9$ cm' at 0.001,0.7 rotate by 0 center tc ls 101
#set label 3 '\ft $N =100$' at 0.5,0.015 rotate by 0 center tc ls 101

set output "freefall.tex"

set logscale

set xr[0.0001:2]
set yr[0.01:1]

set ylabel '$R/R_0$'
set xlabel '$\pi/2-(t-t_0)(8\pi G \rho_0/3)^{1/2}$'

plot 'freefall.dat' u 3:1 w l ls 17 lw 3,\
     'freefall.dat' u 2:1 pt 7 ps 1 lc rgb '#0072bd'

set terminal epslatex size 8cm,6cm color colortext 11 standalone header "\\newcommand{\\ft}[0]{\\footnotesize}"

set output "freefall-rho.tex"

unset logscale
set logscale y

set xr[0.0:0.8]
set yr[1e7:1e13]

set format y '$10^{%T}$'

set ylabel '$\rho$ (g/cm$^3$)'
set xlabel '$t$ (s)'

plot 'freefall-rho.dat' w l ls 17 lw 3


set output "freefall-profile.tex"

set logscale

set xr[1:1e12]
set yr[1e6:2e14]

set format '$10^{%T}$'

set ylabel '$\rho$ (g/cm$^3$)'
set xlabel '$R$ (cm)'

set ytics 1e2
set xtics 1e3

unset label

set label 1 '\ft  $t=0$ s ' at 9e2,1.7e7 rotate by 0 center tc ls 101
set label 2 '\ft  $t=0.66410$ s ' at 1e3,4e13 rotate by 0 center tc ls 101
set label 3 '\ft  $t=0.55009$ s ' at 1e3,1.8e8 rotate by 0 center tc ls 101
set label 4 '\ft  $t=0.64505$ s ' at 1e3,4.5e9 rotate by 0 center tc ls 101
set label 5 '\ft  $t=0.66001$ s ' at 1e3,8e10 rotate by 0 center tc ls 101
set label 6 '\ft  $t=0.66350$ s ' at 1e3,2.7e12 rotate by 0 center tc ls 101
set label 7 '\ft  $\rho \propto R^{-3}$' at 5e9,1e10 rotate by 0 center tc ls 101


plot 'freefall-profile-t=0.dat' w l ls 17 lw 3,\
     'freefall-profile-t=0.55009.dat' w l ls 17 lw 3,\
     'freefall-profile-t=0.645052.dat' w l ls 17 lw 3,\
     'freefall-profile-t=0.660012.dat' w l ls 17 lw 3,\
     'freefall-profile-t=0.663501.dat' w l ls 17 lw 3,\
     'freefall-profile-t=0.6641.dat' w l ls 17 lw 3,\
     'freefall-rho.dat' u 3:2 w l lc rgb '#0072bd' dt 2 lw 3

exit
#Necessario para compilar o arquivo diretamente com pdflatex
