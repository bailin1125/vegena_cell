set terminal X11 persist
set title "Bandstructure" 
set xlabel "KPOINTS"   
set ylabel "Energy (eV)" 
set xrange [0.00:0.00]
set yrange [-0.9471080E+02:0.3194881E+02]
set grid
plot "bandstructure_1.txt" using 1:2 with  linespoints  pointtype 7 pointsize 0.4 linecolor rgb "blue" title "UP", "bandstructure_2.txt" using 1:2 with linespoints pointtype 2 pointsize 0.5 linecolor rgb "red" title "DOWN"
