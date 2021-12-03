!########### MODIFICACION DEL PROGRAMA ###########
!	-Tolerancia
!	-Maxima cantidad de iteraciones
!	-Metodo (primero graficar para ver los intervalos)
!	-Extremos a,b
!	-En subrutina "grafica" seleccionar
!		graficar y=x 
!		poner la funcion a graficar
!	-En la funcion "f(x)" poner la funcion a analizar
!	-En la funcion "fprima(x)" poner la derivada sacada a mano
!	-En la funcion "maxderivada" ajustar el paso
 
PROGRAM raices
IMPLICIT NONE
!#######################################################################
!# Este programa esta diseñado para encontrar raices de funciones.     #
!# Recordar colocar el metodo deseado para la resolucion en el         #
!# parametro "metodo".                                                 #
!#######################################################################

!####DECLARACION DE VARIABLES###########################################
REAL(8) :: a,b
REAL(8),PARAMETER :: tolerancia=0.0000001
REAL(8),PARAMETER :: e=2.718281828

!# En el paremetro tolerancia se coloca el maximo error.               #
!# Para el metodo de Biseccion, si no se especifica, usar 0.001        #
!#######################################################################
INTEGER,PARAMETER :: maxiter=100
!#######################################################################
!# En el paremetro maxiter se coloca la maxima cantidad iteraciones a  #
!# realizar.                                                           #
!#######################################################################
INTEGER,PARAMETER :: metodo=2
!#######################################################################
!# En el parametro metodo se selecciona lo que se desea realizar:      #
!#  1: Graficar la funcion.                                            #
!#  2: Metodo Biseccion.                                               #
!#  3: Metodo Punto fijo sistematico.                                  #
!#  4: Metodo Newton-Raphson.                                          #
!#######################################################################

!####SECCION EJECUTABLE#################################################
OPEN(UNIT=1,FILE='resolucion.dat',STATUS='REPLACE')
!Aca se colocan los extremos.
a=0.005
b=0.02
!·······································································
CALL verificaextremos(a,b)

SELECT CASE(metodo) 
    CASE(1)
        CALL grafica(a,b)
    CASE(2)
        CALL biseccion(a,b)
        WRITE(*,*)'La resolucion se encuentra en archivo resolucion.dat'
    CASE(3)
        CALL puntofijo(a,b)
        WRITE(*,*)'La resolucion se encuentra en archivo resolucion.dat'
    CASE(4)
        CALL newtonraphson(a)
        WRITE(*,*)'La resolucion se encuentra en archivo resolucion.dat'
END SELECT

CLOSE(UNIT=1,STATUS='KEEP')

CONTAINS
!####FUNCIONES Y SUBRUTINAS#############################################
SUBROUTINE verificaextremos(a,b)
!Esta subrutina verifica si los extremos se ingresaron correctamente.
REAL(8) :: a,b

IF (a>b) THEN
    WRITE (*,*)'a>b! Se ingresaron incorrectamente los extremos.'
END IF

END SUBROUTINE verificaextremos
!#######################################################################
SUBROUTINE grafica(a,b)
!Esta subrutina grafica una funcion que se coloca mas abajo.
REAL(8) :: a,b
OPEN(UNIT=2,FILE='grafica.p',STATUS='REPLACE')

WRITE(2,*)'set xrange[',a,':',b,' ]'
WRITE(2,*)'unset log                              # quita la escala logaritmica (si la hubiera)'
WRITE(2,*)'unset label                            # quita los titulos anteriores'
WRITE(2,*)'set xtic auto                          # establece automaticamente las divisiones del eje x'
WRITE(2,*)'set ytic auto                          # establece automaticamente las divisiones del eje y'
WRITE(2,*)'set grid'
WRITE(2,*)'set title " Grafico de f(x)=0 para hallar raiz "'
WRITE(2,*)'set xlabel "x"'
WRITE(2,*)'set ylabel "y"'
WRITE(2,'(a)',ADVANCE='NO')'plot '
!ACÁ AGREGUÉ PARA QUE GRAFIQUE Y=X, LO CANCELO EN CASO  DE NO SER NECESARIO
!WRITE(2,'(a)',ADVANCE='NO')'x'
!WRITE(2,'(a)',ADVANCE='NO')' with lines,\'
!Write(2,'(a)')
!.......................................................................
!En el siguiente renglon se escribe la funcion a graficar.
WRITE(2,'(a)',ADVANCE='NO')'-73517126.479588*(x**4) + 6342523.0872* (x**3) -224327.4899 * (x**2) + 3792.8154*x + 73.89724'

!·······································································
WRITE(2,'(a)',ADVANCE='NO')' with lines'
call system ("gnuplot -persist grafica.p")

CLOSE(UNIT=2,STATUS='KEEP')
END SUBROUTINE grafica
!#######################################################################
FUNCTION f(x)
!En esta funcion se coloca la funcion a analizar.
REAL(8) x,f

    f= -73517126.479588*(x**4) + 6342523.0872* (x**3) -224327.4899 * (x**2) + 3792.8154*x + 73.89724 - 95.
END FUNCTION f
!#######################################################################
FUNCTION fprima(x)
!En esta funcion se coloca la derivada de la funcion a analizar.
!Se tiene que derivar a mano.
REAL(8) x,fprima

    fprima= cos(x)-0.25*(1/sqrt(x))

END FUNCTION fprima
!#######################################################################
FUNCTION maxderivada(a,b)
!Esta funcion se encarga de encontrar el valor de la maxima derivada
!en el intervalo a analizar.
REAL(8) :: q,a,b,paso,maxderivada,t

q=a
!Este es el paso en el cual se recorre el intervalo, se debe cambiar
!dependiendo de la dimension del mismo.
paso=0.001
!·······································································
maxderivada=0.
DO WHILE (q<=b)
    t=fprima(q)
    IF (abs(t)>abs(maxderivada)) THEN
        maxderivada=t
    END IF
    q=q+paso
END DO
END FUNCTION maxderivada
!#######################################################################
SUBROUTINE biseccion(a,b)
!Esta subrutina utiliza el metodo de Biseccion para encontrar una raiz.
INTEGER :: i,n
REAL(8) :: a,b,a2,b2,m,errorenx,erroreny,raiz, xant, errorpaso

WRITE(1,'(A/)')' RESOLUCION MEDIANTE EL METODO DE BISECCION '
WRITE(1,*)
!n es la cantidad de iteraciones necesarias para asegurar el error deseado.
n=FLOOR((LOG(ABS(b-a)/tolerancia))/(LOG(2.0))+0.5)
!·······································································
WRITE(*,'(A,I4,A/)')'Se realizaran aproximadamente ',n,' iteraciones.'
a2=a
b2=b
erroreny=99999.0
i=0
m=(a2+b2)/2.0
errorpaso=1000.0
DO WHILE (tolerancia<ABS(erroreny))
!   DO WHILE (tolerancia < ABS(errorpaso))
    i=i+1
    xant=m
    errorenx=(abs(b2-a2)/2.0)
    erroreny=f((a2+b2)/2.0)
    WRITE(1,*)'------------------------------------------------------------------'
    WRITE(1,*)'Iteración: ',i
    WRITE(1,'(A,F20.7,A,F20.7,A)')'Se analiza intervalo: [ ',a2,' , ',b2,' ]'
    WRITE(1,'(A,F20.7)')'El error en x es: ',errorenx
    WRITE(1,'(A,F20.7)')'El error en y es: ',erroreny
    WRITE(1,*)
    IF (f(a2)*f(m)<0) THEN
        b2=m
        WRITE(1,*)'B pasara a ser: ',b2
        WRITE(1,*)'------------------------------------------------------------------'
    ELSE
        a2=m
        WRITE(1,*)'A pasara a ser: ',a2
        WRITE(1,*)'------------------------------------------------------------------'
    END IF
    m=(a2+b2)/2.0
    errorpaso= xant-m
END DO
IF ((ABS(erroreny)<tolerancia)) THEN
    errorenx=(abs(b2-a2)/2.0)
    erroreny=f(m)
    raiz=m
    WRITE(1,*) '------------------------------------------------------------------'
    WRITE(1,'(A,F20.7)')'La raiz encontrada es: ',raiz
    WRITE(1,'(A,F20.7)')'El error en x es: ',errorenx
    WRITE(1,'(A,F20.7)')'El error en y es: ',erroreny
ELSE
    WRITE(1,'(A/)')'No hay raices reales en el intervalo o se trata de una SINGULARIDAD'
END IF
WRITE(*,'(A,F20.7/)')'La raiz encontrada es: ',raiz

END SUBROUTINE biseccion
!#######################################################################
SUBROUTINE puntofijo(a,b)
!Esta subrutina utiliza el metodo de Punto Fijo Sistematico para
!encontrar una raiz.
REAL(8) :: a,b,x,lamda,error, xant, errorenx
INTEGER :: iter

WRITE(1,'(A/)')' RESOLUCION MEDIANTE PUNTO FIJO SISTEMATICO '
!Esto se hace para que la subrutina entre en el ciclo, NO BORRAR!
error=(tolerancia*2.0)
errorenx=(tolerancia*2.0)
!·······································································
iter=1
x=a
lamda=(1.0/maxderivada(a,b))
DO WHILE (error>tolerancia) 
!  DO WHILE (errorenx>tolerancia) 
    WRITE(1,*)'---------------------------------------'
    xant=x
    x=x-(lamda*f(x))
    WRITE(1,'(A,I4)')'Iteracion: ',iter
    WRITE(1,'(A,F20.7)')'Valor de x: ',x
    error=abs(f(x))
    WRITE(1,'(A,F20.7)')'Error o distancia a y: ',error
    iter=iter+1
    errorenx= abs(xant-x)
    WRITE(1,'(A,F20.7)')'Error en x: ',errorenx
END DO
WRITE(1,*)'---------------------------------------'
WRITE(1,'(A,F20.7)')'La raiz encontrada es: ',x
WRITE(*,'(A,F20.7)')'La raiz encontrada es: ',x

END SUBROUTINE puntofijo
!#######################################################################
SUBROUTINE newtonraphson (a)
!Esta subrutina utiliza el metodo de Newton-Raphson para encontrar una raiz.
REAL(8) :: a,x,error, errorenx, xant
INTEGER :: iter

WRITE(1,'(A/)') ' RESOLUCION MEDIANTE EL METODO DE NEWTONRAPSHON '
!Esto se hace para que la subrutina entre en el ciclo, NO BORRAR!
error=(tolerancia*2.0)
errorenx=(tolerancia*2.0)
!·······································································
iter=1
WRITE(*,*)'Ingresar el valor de "x0".'
WRITE(*,*)'En caso de no tenerlo, puede colocar el extremo "a".'
WRITE(*,*)'Recordar que el extremos "a" vale: ',a
READ(*,*)x
DO WHILE (error>=tolerancia)  
!   DO WHILE (errorenx> tolerancia)
    WRITE(1,*)'---------------------------------------'
    WRITE(1,*)
    xant=x
    x=x-(f(x)/fprima(x))
    WRITE(1,'(A,I4)')'Iteracion: ',iter
    WRITE(1,'(A,F20.7)')'Valor de x: ',x
    error=abs(f(x))
    WRITE(1,'(A,F20.7)')'Error o distancia a y: ',error
    iter=iter+1
    errorenx= abs(x-xant)
    WRITE(1,'(A,F20.7)')'Error en x : ',errorenx
END DO
WRITE(1,*)'---------------------------------------'
WRITE(1,*)
WRITE(1,'(A,F20.7/)')'La raiz encontrada es: ',x
WRITE(*,'(A,F20.7/)')'La raiz encontrada es: ',x

END SUBROUTINE newtonraphson
!#######################################################################
END PROGRAM

