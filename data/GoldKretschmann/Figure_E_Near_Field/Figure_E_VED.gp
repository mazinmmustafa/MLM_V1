reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 16 cm, 10.0 cm
# size 12 cm, 7.5 cm # Default

set output 'Figure_E_VED.tex'

# set palette rgb 7,5,15; set title "traditional pm3d\n(black-blue-red-yellow)"
# set palette rgb 3,11,6; set title "green-red-violet"
# set palette rgb 23,28,3; set title "ocean (green-blue-white)\ntry also other permutations"
# set palette rgb 21,22,23; set title "hot (black-red-yellow-white)"
# set palette rgb 30,31,32; set title "color printable on gray\n(black-blue-violet-yellow-white)"
# set palette rgb 33,13,10; set title "rainbow (blue-green-yellow-red)"
# set palette rgb 34,35,36; set title "AFM hot (black-red-yellow-white)"
set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000') ; set title "MATLAB (jet)"

set size ratio 1

set xrange [-4:+4]
set yrange [-4:+4]

set cbrange [+200:+320]

set xtics 1
set ytics 1
set cbtics 20
set tics front

set arrow front from -4, 0 to +4, 0 nohead linestyle 1 lc 'black' lw 2
set arrow front from -4, -0.05to +4, -0.05 nohead linestyle 1 lc 'black' lw 2

set xlabel '$x$ [$\mu$m]'
set ylabel '$z$ [$\mu$m]'
set title '$20\log_{10}|\mathcal{R}e[E]|$ [dB]'

plot '../data_VED.txt' with image notitle