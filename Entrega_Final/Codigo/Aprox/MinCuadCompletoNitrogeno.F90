PROGRAM MinCuadCompletoNitrogeno
	IMPLICIT NONE
	INTEGER g,k
	REAL(8), ALLOCATABLE :: x(:),y(:),m(:,:),a(:),e(:)
	REAL(8) var,rms,dur,tol,bis_a,bis_b
	
	!=== VARIABLES CONFIGURABLES ===
	g = 3		!GRADO DEL POLINOMIO
	k = 5		!CANTIDAD DE PUNTOS
	dur = 95.	!DUREZA BUSCADA
	tol = 0.000001	!TOLERANCIA PARA LA BISECCION
	bis_a = 0.0		!EXTREMOS PARA LA BISECCION
	bis_b = 0.02
	!===============================
	ALLOCATE(x(k),y(k),m(g+1,g+2),a(g+1),e(k))
	
	CALL ingresar_datos(k,x,y)
	CALL armar_sistema(k,g,x,y,m)
	CALL gauss(m,a,g+1)
	CALL errores(g,k,x,y,a,e)
	
	var = sum(e**2)/(k-g-1)
	rms=sqrt(sum(e**2)/k)
	
	CALL mostrar_aprox(g,a,var,rms)
	CALL biseccion(bis_a,bis_b,a,g,dur,tol)
CONTAINS

SUBROUTINE ingresar_datos(k,x,y)
	INTEGER i,k
	REAL(8) x(k),y(k)
	!INGRESA LOS PUNTOS DEL ARCHIVO "puntos_<gas>.txt"
	
	OPEN(UNIT=2,FILE='puntos_nitrogeno.txt',STATUS='old')
	
	DO i=1, k
		READ(2,*) x(i),y(i)
	END DO	
	
	CLOSE(2)
END SUBROUTINE

SUBROUTINE armar_sistema(k,g,x,y,m)
	INTEGER k,g,i,j
	REAL(8) x(k),y(k),m(g+1,g+2)
	!ARMA EL SISTEMA EN "Referencia_1.png"
	
	m(1,1)=k
	DO j=2,g+1
		m(1,j)=sum(x**(j-1))
	END DO
	
	DO i=2, g+1
		DO j=1, g+1
			m(i,j)=sum(x**(i+j-2))
		END DO
	END DO
	
	DO i=1, g+1
		m(i,g+2)=sum(y*x**(i-1))
	END DO
END SUBROUTINE

SUBROUTINE gauss(m,x,g)
	INTEGER i,k,g
	REAL(8) m(g,g+1),x(g)
	!APLICA TRIANGULACION DE GAUSS
	
	DO i=1,g
		IF(m(i,i).ne.0) THEN
			m(i,:)=m(i,:)/m(i,i)
			IF(i.ne.g) THEN
				DO k=i+1,g
					m(k,:)=m(k,:)-m(k,i)*m(i,:)
				END DO
			END IF
		END IF
	END DO
	
	CALL susreg(m,x,g)
END SUBROUTINE

SUBROUTINE susreg(m,x,g)
	INTEGER i,k,g
	REAL(8) m(g,g+1),x(g),s
	!APLICA SUSTITUCION REGRESIVA
	
	x(g)=m(g,g+1)/m(g,g)
	i=g-1
	DO WHILE (i>=1)
		s=0
		DO k=i+1,g
			s=m(i,k)*x(k)+s
		END DO
		x(i)=(m(i,g+1)-s)/m(i,i)
		i=i-1
	END DO
END SUBROUTINE

SUBROUTINE errores(g,k,x,y,a,e)
	INTEGER g,k,i,j
	REAL(8) x(k),y(k),a(g+1),e(k),p
	!CALCULA LOS ERRORES Y LOS DEVUELVE EN "e"
	
	DO j=1,k
		p=0
		DO i=1,g+1
			p=a(i)*(x(j)**(i-1))+p
		END DO
		e(j)=y(j)-p
	END DO
END SUBROUTINE

SUBROUTINE mostrar_aprox(g,a,var,rms)
	INTEGER g,i
	REAL(8) a(g+1), var, rms
	!MUESTRA EN PANTALLA LA APROXIMACION OBTENIDA
	
	WRITE(*,*) "Polinomio obtenido:"
	DO i=1,g+1
		WRITE(*,'(F15.5,A3,I2)',ADVANCE='NO') a(i),'*x^',i-1
		IF(i.ne.(g+1)) THEN
			WRITE(*,'(A3)',ADVANCE='NO') ' + '
		END IF
	END DO
	
	WRITE(*,*)
	WRITE(*,*) "Varianza: ",var," | RMS:",rms
END SUBROUTINE

SUBROUTINE biseccion(bis_a,bis_b,a,g,dur,tol)
	INTEGER g,n,i
	REAL(8) bis_a,bis_b,a(g+1),dur,tol,m
	REAL(8) aa,bb
	!CALCULA LA CONCENTRACION PARA LA DUREZA DADA POR BISECCION
	
	n = FLOOR((LOG(ABS(bis_b - bis_a)/tol))/LOG(2.) + 0.5)
	
	aa = bis_a
	bb = bis_b
	DO i=1,n
		m=(aa+bb)/2.0
		IF(f_raiz(aa,a,g,dur)*f_raiz(m,a,g,dur) < 0) THEN
			bb = m
		ELSE
			aa = m
		ENDIF
	END DO
	
	m=(aa+bb)/2.0
	WRITE(*,*) "Se realizaron ",n," iteraciones de biseccion"
	WRITE(*,*) "Se encontro un valor de dureza de ",dur," a una concentracion de ",m
END SUBROUTINE

FUNCTION f_raiz(x,a,g,dur)
	INTEGER g,i
	REAL(8) f_raiz,x,dur
	REAL(8) a(g+1)
	!OBTIENE EL VALOR DEL POLINOMIO CON RAIZ EN LA DUREZA BUSCADA
	
	f_raiz = 0.
	DO i=1,g+1
		f_raiz = a(i)*x**(i-1) + f_raiz
	END DO
	f_raiz = f_raiz - dur
END FUNCTION

END PROGRAM MinCuadCompletoNitrogeno
