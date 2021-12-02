set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Aproximacion por Cuadrados Minimos"
set xlabel "x"
set ylabel "P(x)"
set grid

plot "Variacion.dat" using 1:2 title 'Explicito' with lines,
