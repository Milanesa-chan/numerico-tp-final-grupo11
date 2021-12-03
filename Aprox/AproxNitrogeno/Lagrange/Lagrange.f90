!############### MODIFICACIONES DEL PROGRAMA #################
!	-Cambiar m (cantidad de puntos
!	-Ingresar x e y

PROGRAM nombreProg

  IMPLICIT NONE

! DECLARACION DE VARIABLES
integer yyy, i,t,k,g,OP
real(8), allocatable, dimension(:) :: ex, y, xr, coeficientes
real(8) Fm,xs, xx, yy, Error,xf,h
real(8), allocatable, dimension(:,:) :: datos
integer , parameter :: m=5 !cantidad de puntos
integer , parameter :: n=m-1 ! grado del polinomio


allocate (ex(m),y(m),coeficientes(0:m-1))
allocate (datos(0:m-1,2))


ex(1)=0.005 ; y(1)=88. 
ex(2)=0.01 ; y(2)=95.
ex(3)=0.015 ; y(3)=98.
ex(4)=0.02 ; y(4)=99.
ex(5)=0.024 ; y(5)=99.
do k=1,m
  datos(k-1,1)=ex(k)
  datos(k-1,2)=y(k)
enddo

call system ("clear") 

OP=1
DO WHILE (OP==1)
write (*,*) ' NECESARIAMENTE HAY QUE PASAR POR LA OPCION 1'
write (*,*)
write (*,*) ' Ingrese 0 si quiere obtener los coeficientes del polinomio'
write (*,*) ' Ingrese 1 si quiere obtener el polinomio evaluado en un X determinado. Tambien te da los Lk evaluados en el punto'
write (*,*) ' Ingrese 2 si quiere plotear'
write (*,*) ' Ingrese 3 si quiere estimar el error de interpolacion'

read (*,*) yyy

select case (yyy)

case(0)
    OPEN(33,file='polinomio.dat', status='replace')
    call lagrange(datos, coeficientes)
    print*,'El polinomio obtenido es:'
    write(*,'(A2)',advance='NO')'Y='
    write(33,'(A2)',advance='NO')'Y='
    do g=m-1,1,-1
           write(*,'(f20.6,A4,I2,A1)',advance='no')coeficientes(g),'* x^',g,'+'
           write(33,'(f20.6,A4,I2,A1)',advance='no')coeficientes(g),'* x^',g,'+'
    enddo
    write(*,*) coeficientes(0)
    write(33,*) coeficientes(0)
    close(33, status='keep')
case(1)
!LOS COEFICIENTES DEL POLINOMIO DEBEN SACARSE A MANO USANDO:
!L(k)=Multiplicatoria desde i=0,n de (x-xi)/(xk-xi)
!luego los L(k) se multiplican por el yk dato cada uno
!y luego se suman y resolviendo llego al polinomio de lagrange

!Es la iteracion K, en el punto Xk,Yk se hace la sumatoria con todos los
!xi de todas las iters

call system ("clear")

write (*,*) ' Ingrese el punto donde se evaluara el polinomio'
read (*,*) xx

!Graba tabla de valores
open (2, file='interp.dat', status='replace')
do i=1,m
 write (2,'(2F20.6)') ex(i), y(i)
end do
close (2, status='keep')

read(*,*)
call system ("clear")

call PoliLagrange(m,n,xx,yy,ex,y)
write(*,'(A28,F20.6, A4,F20.6)') 'El valor del polinomio en X=',xx,' es ',yy
write(*,'(A1,I1,A1,F20.6,A2,F20.6)') 'P',n,'(',xx,')=',yy
!*************************************************************************
case (2)
open (23, file='poli.dat', status='replace')
print*,'Ingrese un paso para graficar'
read*,h
xf=ex(1)
do while (xf.le.ex(m))
 call PoliLagrange(m,n,xf,yy,ex,y)
 write(23,'(2F20.6)')xf,yy
 xf=xf+h
end do
close (23, status='keep')
call system ("gnuplot -persist 'script.p'")

!*************************************************************************
case(3) ! ERROR DE INTERPOLACION 
 
  call system ("clear")
 write (*,*) ' Ingrese el numero de puntos a considerar'
 read (*,*) t
 allocate(xr(0:t-1))
  call system ("clear")
 write (*,'(A22,I2,A14)') ' Ingrese la derivada F',(t),' de la funcion'
 WRITE(*,*) 'LA DERIVADA QUE ME PIDE ES UN VALOR (NO UNA FUNCIÓN) YA QUE ESPECIALIZO EN EL PTO MEDIO'
 read (*,*) Fm
  call system ("clear")
 write (*,*) ' Ingrese los datos en X'
 do i=0,t-1
 write (*,'(A17,I2)') ' Ingrese el x', i
 read (*,*) xr(i)
 end do
   call system ("clear")
 write (*,*) ' Ingrese el x donde quiere estimar el error (NO TIENE POR QUÉ SER EL VALOR MEDIO DEL INTERVALO)'
 read (*,*) xs
   call system ("clear")
   Error=E(xs,t,xr,Fm)
Print*,'El error es ', Error    !xr son los x(i)
  read (*,*)
end select
END DO

contains

SUBROUTINE lagrange(datos, pol)
! Devuelve el polinomio de lagrange de grado (n-1) 
! como resultado de los productos de todos los (x - xi)
REAL(8), ALLOCATABLE :: polAux(:), xValues(:), polLk(:)
REAL(8) pol(0:), datos(0:, :), denom
INTEGER i, n, g, k

n = SIZE(datos, DIM=1)-1
ALLOCATE(polLk(0:n), polAux(0:n), xValues(0:n-1))

pol = 0
DO k=0, n

! Elimina el elemento k de la lista de valores de x
  xValues(:k-1) = datos(:k-1, 1)
  xValues(k:) = datos(k+1:, 1)

! Arma el primer monomio para comenzar los productos
  polAux(1) = 1.0
  polAux(0) = -xValues(0)
  denom = datos(k, 1) -xValues(0)

  DO g=1, n-1
    polLk(g+1) = 1.0
    DO i = g, 1, -1
      polLk(i) = polAux(i-1) - polAux(i)*xValues(g) !para hallar el coeficiente correspondiente al grado i, sumamos el de i-1 con el producto del i por el termino indep(del q agregamos(x-xk))
    END DO
    polLk(0) = -polAux(0)*xValues(g)
    polAux = polLk
    denom = denom*(datos(k, 1) -xValues(g))
  ENDDO
  polLk = polLk*datos(k, 2)/denom

!Suma cada polinomio polLk
  pol = pol + polLk
ENDDO

END SUBROUTINE

function E(x,t,xr,Fm1)
real(8) E, x, v, Fm1
real(8) xr(0:t-1)
integer i,t,fact,k

fact=1
do i=1,t
fact=fact*i
end do


v=1.0
 Do k=0,t-1
  v = v*(x-(xr(k)))
 end do 
print*,'Productoria ',v
print*,'Derivada ', Fm1
print*,'factorial ',fact

E= (v*Fm1)/fact

end function E

Subroutine PoliLagrange(iv,n,xx,yy,X,Y)  
  integer i,iv,k,n
  real(8)  X(iv), Y(iv)
  real(8)  xx,yy, U, S
  ! Se fija si el punto esta dentro del intervalo de interpolacion
  if ((xx<X(1)).or.(xx>X(iv))) then
   print *,' El punto elegido esta fuera del intervalo de interpolacion'
  else
  
  !interpolacion  
  yy=0
  do k = 1, n+1
  S=1  
  U=1
    do i = 1, n+1
      if (i/=k) then
      S=S*(X(k)-X(i))
      U=U*(xx-X(i))
      end if
    end do
    
    yy=yy+(U/S)*Y(k)
    !print*, 'El coeficiente L ',(k-1) ,'es ', (U/S) 
  end do
  !print*,'OJO!No son los coeficientes del polinomio!'
  end if
end subroutine


END PROGRAM


