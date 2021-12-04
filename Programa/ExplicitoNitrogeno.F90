Program Explicito

IMPLICIT NONE
Real,Parameter:: Pi=3.141592654,Dab=1E-6
Real,Parameter:: largo=0.3
Real,Parameter :: conc_buscada = 0.0100546, punto_buscado = 0.005
Integer, Parameter :: max_iter = 5000
Integer :: nodos
Real :: z, dt, dx, tiempo
Real, ALLOCATABLE :: T(:)

!z=0.01
dt=8.
!dx=(dt*Dab)/z
dx=0.0025
!dt=(dx*z)/Dab
nodos = int(largo/dx + 1)
ALLOCATE(T(nodos))
Write(*,*) dx, nodos, z, dt, Dab

Open(1,File='Variacion.dat')
T(1)=0.022
T(2:nodos-1)=0.
T(nodos)=0.

Call Parabolicasexplicito(T, dt, dx, nodos, tiempo)
Close(1)
!call system ("gnuplot -persist grafica_pol.p")

Contains
Subroutine Parabolicasexplicito(T, dt, dx, nodos, tiempo)
Real T(nodos),Tant(nodos),Tf,dt,dx,r,tiempo
Integer i, nodos, i_buscado

i_buscado = int(punto_buscado/dx + 1.5)
Write(*,*) i_buscado
Read(*,*)
tf=0
do while (tf.lt.max_iter) !Tiempo final             !ESCRIBO EL VECTOR EN CADA TIEMPO
	Write(*,'(F10.2,F10.5)',Advance='yes') tf, t(i_buscado)
	if (t(i_buscado).gt.conc_buscada) then
		tiempo = tf
		!Write(*,*) tiempo
		!Write(*,'(F10.2,F10.5)',Advance='yes') tf, t(i_buscado)
		return
	endif
	Write(1,*)
	 tant=t 
	 do i=2,nodos-1
		r = dx*(i-1)
		t(i)=dt*Dab*((1./r)*((tant(i+1)-tant(i-1))/(2.*dx)) + ((tant(i+1)-2.*tant(i)+tant(i-1))/(dx**2.))) + tant(i)
		!t(i)=(z/r)*((tant(i+1)-tant(i-1))/(2.)) + z*((tant(i+1)-2.*tant(i)+tant(i-1))/(dx)) + tant(i)
	 end do
	 tf=tf+Dt
end do

Endsubroutine
EndProgram
