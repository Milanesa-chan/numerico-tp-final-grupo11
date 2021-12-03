# Ejemplo de script Gnuplot, para graficar los datos 
# que se encuentran en un archivo llamado "prueba.dat"


set   autoscale                        # escala los ejes automaticamente
unset log                              # quita la escala logaritmica (si la hubiera)
unset label                            # quita los titulos anteriores
set xtic auto                          # establece automaticamente las divisiones del eje x
set ytic auto                          # establece automaticamente las divisiones del eje y
set title "Euler Modificado"
set xlabel "EJE x"
set ylabel "EJE y"
set grid

plot "interp.dat" using 1:2 title 'puntos dato' with points,\
     "poli.dat" using 1:2 title 'polinomio interpolante' with lines



