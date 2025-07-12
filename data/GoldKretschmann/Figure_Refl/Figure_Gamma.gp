reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure_Gamma.tex'

set key bottom right

set xrange [+0.0:+180.0]
set yrange [0.0:1.0]

set xtics 30.0
set ytics 0.2

# set mxtics 5
# set mytics 5

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

# set arrow from 0.0, 0.0 to 2.0, 0.0 nohead linestyle 1 lc 'black' lw 2

set xlabel '$\theta_{i}$ [deg]'
set ylabel '$|\Gamma|^{2}$'
set title ''

plot '../data_refl.txt' using 1:2 with lines lw 4 dt 1 lt 6 title '$e$', \
     '../data_refl.txt' using 1:3 with lines lw 4 dt 1 lt 7 title '$h$'