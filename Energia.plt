set term pngcairo size 900,700 enhanced font 'Verdana,10'
set output 'energia.png'
set multiplot layout 2,1
set title "H (no conservada)"
set ylabel "H"
plot "Energia.dat" u 1:2 with lines lw 1.5 lc rgb 'blue' notitle
set title "H' = H - {/Symbol w}p_{/Symbol f}  (constante del movimiento)"
set ylabel "H'"

stats "Energia.dat" u 1:3 nooutput


set yrange [STATS_min_y - (STATS_max_y * 0.1):STATS_max_y * 1.1]
plot "Energia.dat" u 1:3 with lines lw 1.5 lc rgb 'red' notitle
unset multiplot