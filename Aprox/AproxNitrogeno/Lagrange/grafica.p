 set xrange[   0.0000000000000000      :   1.9999999552965164E-002  ]
 unset log                              # quita la escala logaritmica (si la hubiera)
 unset label                            # quita los titulos anteriores
 set xtic auto                          # establece automaticamente las divisiones del eje x
 set ytic auto                          # establece automaticamente las divisiones del eje y
 set grid
 set title " Grafico de f(x)=0 para hallar raiz "
 set xlabel "x"
 set ylabel "y"
plot -73517126.479588*(x**4) + 6342523.0872* (x**3) -224327.4899 * (x**2) + 3792.8154*x + 73.89724 with lines
