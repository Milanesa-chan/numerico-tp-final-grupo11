

set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Aproximacion por Cuadrados Minimos"
set xlabel "x"
set ylabel "P(x)"
set grid
plot "polinomio.dat" using 1:2 title 'P(x)=ao + a1 x + ... + an x^n' with lines, "puntos.dat" lt 3 title ''
