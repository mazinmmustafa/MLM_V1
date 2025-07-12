reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure_E.tex'

set key top right

set xrange [-2.0:+1.0]
set yrange [0.0:+1.8]

set xtics 0.5
set ytics 0.4

# set mxtics 5
# set mytics 5

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

set arrow from 0.0, 0.0 to 0.0, 1.8 nohead ls 1 dt 1 lc 'black' lw 1
set arrow from -0.2, 0.0 to -0.2, 1.8 nohead ls 1 dt 1 lc 'black' lw 1
set arrow from -0.5, 0.0 to -0.5, 1.8 nohead ls 1 dt 1 lc 'black' lw 1
set arrow from -1.0, 0.0 to -1.0, 1.8 nohead ls 1 dt 1 lc 'black' lw 1
set arrow from -1.3, 0.0 to -1.3, 1.8 nohead ls 1 dt 1 lc 'black' lw 1
set arrow from -1.5, 0.0 to -1.5, 1.8 nohead ls 1 dt 1 lc 'black' lw 1

set xlabel '$z$ [m]'
set ylabel '$|E|$ [V/m]'
set title ''

plot '../data_PW.txt' using 1:2 with lines lw 4 dt 1 lt 6 title '$E_{x}$', \
     '../data_PW.txt' using 1:3 with lines lw 4 dt 1 lt 7 title '$E_{y}$', \
     '../data_PW.txt' using 1:4 with lines lw 4 dt 1 lt 10 title '$E_{z}$', \
     '../Fields_E_PW.dat_save' using 1:2 with lines lw 2 dt 2 lt -1 title '', \
     '../Fields_E_PW.dat_save' using 1:3 with lines lw 2 dt 2 lt -1 title '', \
     '../Fields_E_PW.dat_save' using 1:4 with lines lw 2 dt 2 lt -1 title ''