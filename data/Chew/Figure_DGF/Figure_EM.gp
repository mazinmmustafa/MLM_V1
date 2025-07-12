reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure_EM.tex'

set key top right

set xrange [-3.0:+3.0]
set yrange [0.0:1.5]

set xtics 1.0
set ytics 0.5

# set mxtics 5
# set mytics 5

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

# set arrow from 0.0, 0.0 to 2.0, 0.0 nohead linestyle 1 lc 'black' lw 2

set xlabel '$x$ [m]'
set ylabel '$|E|$ [V/m]'
set title ''

plot '../data_EM.txt' using 1:2 with lines lw 4 dt 1 lt 6 title '$E_{x}$', \
     '../data_EM.txt' using 1:3 with lines lw 4 dt 1 lt 7 title '$E_{y}$', \
     '../data_EM.txt' using 1:4 with lines lw 4 dt 1 lt 10 title '$E_{z}$', \
     '../Fields_M.dat_save' using 1:2 with points pt 6 ps 0.5 lt 16 title '', \
     '../Fields_M.dat_save' using 1:3 with points pt 6 ps 0.5 lt 16 title '', \
     '../Fields_M.dat_save' using 1:4 with points pt 6 ps 0.5 lt 16 title ''