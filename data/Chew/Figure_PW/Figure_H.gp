reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure_H.tex'

set key top right

set xrange [-2.0:+1.0]
set yrange [0.0:+4.5]

set xtics 0.5
set ytics 1.0

# set mxtics 0.1
# set mytics 5

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

set arrow from 0.0, 0.0 to 0.0, 4.5 nohead ls 1 dt 1 lc 'black' lw 1
set arrow from -0.2, 0.0 to -0.2, 4.5 nohead ls 1 dt 1 lc 'black' lw 1
set arrow from -0.5, 0.0 to -0.5, 4.5 nohead ls 1 dt 1 lc 'black' lw 1
set arrow from -1.0, 0.0 to -1.0, 4.5 nohead ls 1 dt 1 lc 'black' lw 1
set arrow from -1.3, 0.0 to -1.3, 4.5 nohead ls 1 dt 1 lc 'black' lw 1
set arrow from -1.5, 0.0 to -1.5, 4.5 nohead ls 1 dt 1 lc 'black' lw 1

set xlabel '$z$ [m]'
set ylabel '$|H|$ [mA/m]'
set title ''

plot '../data_PW.txt' using 1:(1000*$5) with lines lw 4 dt 1 lt 6 title '$H_{x}$', \
     '../data_PW.txt' using 1:(1000*$6)  with lines lw 4 dt 1 lt 7 title '$H_{y}$', \
     '../data_PW.txt' using 1:(1000*$7)  with lines lw 4 dt 1 lt 10 title '$H_{z}$', \
     '../Fields_H_PW.dat_save' using 1:(1000*$2)  with lines lw 2 dt 2 lt -1 title '', \
     '../Fields_H_PW.dat_save' using 1:(1000*$3)  with lines lw 2 dt 2 lt -1 title '', \
     '../Fields_H_PW.dat_save' using 1:(1000*$4)  with lines lw 2 dt 2 lt -1 title ''