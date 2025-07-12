reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure_HM.tex'

set key top right

set xrange [-3.0:+3.0]
set yrange [0.0:5.0]

set xtics 1.0
set ytics 1.0

# set mxtics 5
# set mytics 5

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

# set arrow from 0.0, 0.0 to 2.0, 0.0 nohead linestyle 1 lc 'black' lw 2

set xlabel '$x$ [m]'
set ylabel '$|H|$ [mA/m]'
set title ''

plot '../data_HM.txt' using 1:(1000*$2) with lines lw 4 dt 1 lt 6 title '$H_{x}$', \
     '../data_HM.txt' using 1:(1000*$3) with lines lw 4 dt 1 lt 7 title '$H_{y}$', \
     '../data_HM.txt' using 1:(1000*$4) with lines lw 4 dt 1 lt 10 title '$H_{z}$', \
     '../Fields_M.dat_save' using 1:(1000*$5) with points pt 6 ps 0.5 lt 16 title '', \
     '../Fields_M.dat_save' using 1:(1000*$6) with points pt 6 ps 0.5 lt 16 title '', \
     '../Fields_M.dat_save' using 1:(1000*$7) with points pt 6 ps 0.5 lt 16 title ''