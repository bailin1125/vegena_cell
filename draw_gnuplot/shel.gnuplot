set terminal postscript eps color solid  lw 2  "Times_New_Roman" 20
set title "test" 
set xlabel "test.x"   
set ylabel "test.y" 
set xrange [0:9]
set yrange [0:36]
set xtics 0,1,9
set output "test_shell.eps"
plot "metal_xuhao" with linespoints linecolor 3 linewidth 2 pointtype 7 pointsize 2
set output

set title "Be-O" 
set size ratio 0.8
unset key
set xlabel "test.x"   
set ylabel "test.y" 
set xrange [0:9]
set yrange [-0.5:60]
set xtics 0,1,7
set grid ytics
set grid xtics
set output "test_shell1.eps"
plot "Be-O.txt" using 1:2 with points linecolor -1 pointtype 7 pointsize 1,"" using 1:2 smooth csplines linecolor 7 linewidth 3  
set output
set term wxt