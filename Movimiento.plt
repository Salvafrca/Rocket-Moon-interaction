set term gif animate delay 5
set output 'nave_animacion.gif'
set size square

# Fijamos los rangos para que el marco de la animación no salte
set xrange [-1.2:1.2]
set yrange [-1.2:1.2]

# Bucle para crear los fotogramas (ajusta el 300 según el número de líneas que quieras animar)
do for [a=1:3000:10] {
    plot "Posiciones.dat" every ::0::a u 1:2 with lines title 'Trayectoria', \
         "Posiciones.dat" every ::a::a u 1:2 with points pointtype 7 pointsize 1.5 title 'Nave',\
         "Luna.dat" every ::0::a u 1:2 with points pointtype 7 pointsize 1.5 title 'Luna'
}