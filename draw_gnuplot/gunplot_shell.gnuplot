set terminal pdfcairo lw 2 font "Times_New_Roman,12"
set title "version1" 
set xlabel "real distance"   
set ylabel "metallic radius" 
set xrange [2:5]
set yrange [2:5]
set xtics 2,0.5,5
set output "all_version.pdf"
plot "all_newversion1"  using 1:2 with points lc "blue" pt 30 ps 0.05 title "data" , (x) with lines lw 2
set output
set term wxt