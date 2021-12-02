!############### MODIFICACIONES DEL PROGRAMA #################
!	-Los puntos aca se ingresan por archivo.
!	-Grado de polinomio
!	-Cantidad de puntos
!	-Dx paso para graficar linea 158 (ACORDE A LA MAGNITUD DE LOS PUNTOS)
!	-Activar funcion de evaluar la funcion en un punto (en caso de ser necesario)
!	-Indicar el punto a evaluar en la subrutina linea 196

PROGRAM Cuad_Min
IMPLICIT NONE
integer n,k,i
real(8), allocatable :: x(:),y(:),m(:,:),a(:),e(:)
real(8) var,rms

!Sale un polinomio aproximante:  a0 + a1 x + a2 x^2 + ... + an x^n
n=3                               !GRADO POLINOMIO APROXIMANTE
k=7                                !CANTIDAD DE PUNTOS (enter al final del archivo)

allocate(x(k),y(k),m(n+1,n+2),a(n+1),e(k))  ! m es la matriz ampliada del sistema lineal

call ingreso_datos(k,x,y)
call sistema_lineal(k,n,x,y,m)
call gauss(m,a,n+1)
call grafica(n,k,x,a)
call errores(n,k,x,y,a,e)

write(*,*) 'Polinomio aproximante:'
write(*,*)
write(*,*) 'P(x)=a0 + a1 x + a2 x^2 + ... + an x^n'
write(*,*)
do i=1, n+1
   write(*,'(A2,I1,A2)',advance='no') 'a(',i-1,')='
   write(*,*) a(i)
end do

var=sum(e**2)/(k-n-1)
rms=sqrt(sum(e**2)/k)

write(*,*)
write(*,*) 'La varianza es', var
write(*,*) 
write(*,*) 'El error cuadratico medio (RMS) es', rms
!SI QUIERO LA FUNCION EN UN PUNTO, ACTIVO ACÁ:
!FUNCION PARA CALCULARLA FUNCION EN UN PUNTO
!WRITE(*,*)'la funcion en el punto ingresado es ' 
!write(*,'(F20.6)') calculopunto(a)
CONTAINS

subroutine ingreso_datos(k,x,y)

integer i,k
real(8) x(k),y(k)

open(unit=2,file='puntos_carbono.txt',status='old')

write(*,*) 'Ingresando datos desde el archivo puntos.dat'
write(*,*)

!PARA INGRESAR PUNTOS EN CASO DE QUE NO ANDE EL ARCHIVO
!MUTEAR LINEA 55 Y 77 (parece q estoy hablando de bondis)
!x(1)=0
!x(2)=.25
!x(3)=.5
!x(4)=.75
!x(5)=1.

!y(1)=1.
!y(2)=1.284
!y(3)=1.6487
!y(4)=2.117
!y(5)=2.7183
do i=1, k
   read(2,*) x(i),y(i)
end do

close(unit=2)

end subroutine
!######################################################
subroutine sistema_lineal(k,n,x,y,m)  !forma el sistema lineal a resolver

integer k,n,i,j
real(8) x(k),y(k),m(n+1,n+2)

!################# FORMA LA PRIMER FILA ##############

m(1,1)=k
do j=2, n+1
   m(1,j)=sum(x**(j-1))
end do

!######################################################

do i=2, n+1
   do j=1, n+1   
      m(i,j)=sum(x**(i+j-2))
   end do
end do

!######################################################

do i=1, n+1
   m(i,n+2)=sum(y*x**(i-1))
end do

end subroutine
!######################################################
subroutine gauss(m,x,n)   ! algoritmo Gauss

integer i,k,n
real(8) m(n,n+1), x(n)

  do i=1,n
    if (m(i,i).ne.0) then
       m(i,:)=m(i,:)/m(i,i)
      if (i.ne.n) then
       do K=i+1,n
         m(k,:)=m(k,:)-m(k,i)*m(i,:)
       end do
      end if
     end if
  end do


call regreg(m,x,n)

end subroutine
!######################################################
subroutine regreg(m,x,n)  !sustitucion por regresion regresiva

integer i,k,n
real(8) m(n,n+1), x(n), s

x(n)=m(n,n+1)/m(n,n)
i=n-1
do while (i>=1)
 s=0
 do k=i+1,n
  s=m(i,k)*x(k)+s
 end do
 x(i)=(m(i,n+1)-s)/m(i,i) 
 i=i-1
end do  

end subroutine
!######################################################
subroutine grafica(n,k,x,a)

integer n,i,k
real(8) x(k),a(n+1),xx,dx,p


open(unit=3,file='polinomio.dat',status='unknown')

xx=minval(x)
!ACÁ SE SETEA EL PASO PARA GRAFICAR
dx=0.1                                 !Paso para graficar

do while (xx<=maxval(x))
   p=0
   do  i=1, n+1
       p=a(i)*(xx**(i-1)) + p
   end do
   write(3,'(2F20.10)',advance='no') xx,p
   write(3,*)
   xx=xx + dx
end do



close(unit=3)


call system ("gnuplot -persist 'grafica_pol.p' ")

end subroutine
!######################################################
subroutine errores(n,k,x,y,a,e)

integer n,k,i,j
real(8) x(k),y(k),a(n+1),e(k),p

do j=1, k
   p=0
   do  i=1, n+1
      p=a(i)*(x(j)**(i-1)) + p
   end do
   e(j)=y(j)-p
end do

end subroutine
!########## ESTA FUNCION SOLO LA LLAMA SI ACTIVE PARA QUE ME DEVUELVA EL VALOR EN UN PTO.#############
FUNCTION calculopunto(a)
REAL(8) calculopunto, suma
REAL(8), PARAMETER:: punto=0.8         !ingreso el punto donde deseo calcular la funcion
REAL(8) a(n+1)
INTEGER i
suma=a(1)
DO i=2, n+1 
    suma=suma+ a(i)*(punto**(i-1))
END DO
calculopunto=suma
END FUNCTION 

END PROGRAM

