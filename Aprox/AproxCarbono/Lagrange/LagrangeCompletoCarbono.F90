PROGRAM MinCuadCompletoCarbono
	IMPLICIT NONE
	INTEGER g,k
	REAL(8), ALLOCATABLE :: x(:),y(:),m(:,:),a(:)
	REAL(8) dur,tol,bis_a,bis_b
	
	!=== VARIABLES CONFIGURABLES ===
	k = 7		!CANTIDAD DE PUNTOS
	dur = 95.	!DUREZA BUSCADA
	tol = 0.0001	!TOLERANCIA PARA LA BISECCION
	bis_a = 0.		!EXTREMOS PARA LA BISECCION
	bis_b = 1.
	!===============================
	
	g = k-1		!GRADO DEL POLINOMIO
	ALLOCATE(x(k),y(k),m(g+1,g+2),a(g+1))
	
	CALL ingresar_datos(k,x,y)
	CALL lagrange(x,y,k,g,a)
	
	CALL mostrar_aprox(g,a)
	CALL biseccion(bis_a,bis_b,a,g,dur,tol)
CONTAINS

SUBROUTINE ingresar_datos(k,x,y)
	INTEGER i,k
	REAL(8) x(k),y(k)
	!INGRESA LOS PUNTOS DEL ARCHIVO "puntos_<gas>.txt"
	
	OPEN(UNIT=2,FILE='puntos_carbono.txt',STATUS='old')
	
	DO i=1, k
		READ(2,*) x(i),y(i)
	END DO	
	
	CLOSE(2)
END SUBROUTINE

SUBROUTINE lagrange(x, y, n, g, a)
	INTEGER g, n, c
	REAL(8) x(n),y(n),a(n),xk(n-1),pn(n),denom,mult,comb
	INTEGER k,i,cant_unos,vec(n-1)
	
	a = 0
	DO k=1,n
		xk = 0
		DO i=1,n
			IF(i<k) THEN
				xk(i)=x(i)
			ELSE IF (i>k) THEN
				xk(i-1)=x(i)
			END IF
		END DO
		
		denom = 1.0
		DO i=1,n-1
			denom = denom*(x(k)-xk(i))
		END DO
		
		mult = y(k)/denom
		
		vec = 0
		pn = 0.
		DO i=1,(2**g)
			cant_unos = sum(vec)
			comb = 1.
			DO c=1,g
				IF (vec(c).ne.0) THEN
					comb = comb*vec(c)*(-xk(c))
				ELSE
					comb = comb*1.0
				END IF
			END DO
			pn(n-cant_unos) = pn(n-cant_unos) + comb
			CALL unos(vec,g)
		END DO
		pn = pn*mult
		a = a + pn
	END DO
END SUBROUTINE

SUBROUTINE unos(vec,n)
	INTEGER n, i
	INTEGER vec(n)
	
	i=1
	DO WHILE (vec(i)==1 .AND. i<=n)
		vec(i)=0
		i = i+1
	END DO
	IF (i>n) THEN
		vec = 0
	ELSE
		vec(i) = 1
	END IF
END SUBROUTINE

SUBROUTINE mostrar_aprox(g,a)
	INTEGER g,i
	REAL(8) a(g+1)
	!MUESTRA EN PANTALLA LA APROXIMACION OBTENIDA
	
	WRITE(*,*) "Polinomio obtenido:"
	DO i=1,g+1
		!WRITE(*,'(F10.5,A3,I2)',ADVANCE='NO') a(i),'*x^',i-1
		WRITE(*,*) a(i),'*x^',i-1
		IF(i.ne.(g+1)) THEN
			WRITE(*,'(A3)',ADVANCE='NO') ' + '
		END IF
	END DO
	
	WRITE(*,*)
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

END PROGRAM MinCuadCompletoCarbono
