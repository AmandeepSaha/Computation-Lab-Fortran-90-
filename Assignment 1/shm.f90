! Consider a linear spring mass system. (a) Obtain its phase portrait analytically and graphically. (b) Write a program by using Euler method to solve the equation of motion numerically and plot (i) x(t), v(t) and (iii) the phase portrait.

! Derivative Function
FUNCTION func(x, omega) RESULT(d2xdt2)
	REAL, INTENT(IN) :: x, omega
	d2xdt2 = - omega**2 * x
	
END FUNCTION func

!State Vector
FUNCTION  S(t, x, v, omega) RESULT(State)
	REAL, INTENT(IN)  :: x, v, omega
	REAL, INTENT(OUT) :: State(2)
	
	State  = [v, func(x, omega)]

END FUNCTION S

PROGRAM shm

	IMPLICIT NONE
	
	! To be written
	OPEN(1, file='shm.dat', STATUS='Replace')
	
	! Variables
	REAL :: m, k
	REAL :: t, x, v
	REAL :: k1x, k2x, k3x, k4x
	REAL :: k1v, k2v, k3v, k4v
	
	! Read
	PRINT *, "Enter t0, x0, v0, m, k: "
	READ(*,*) t, x, v, m, k
	
	REAL :: omega = SQRT(k/m)
	
	INTEGER :: i
	INTEGER :: n = 2000
	REAL :: dt = 0.01
	
	CONTAINS
	DO i = 0, n
		
		K1x, k1v =  S(x, v, omega, k1x, k1v)(1)					, S(x, v, omega, k1x, k1v)(2)
		k2x, k2v =  S(x + 0.5*dt*k1x, v + 0.5*dt*k1v, omega)(3)	, S(x + 0.5*dt*k1x, v + 0.5*dt*k1v, omega)(2)
		k3x, k3v =  S(x + 0.5*dt*k2x, v + 0.5*dt*k2v, omega)(2)	, S(x + 0.5*dt*k2x, v + 0.5*dt*k2v, omega)(2)
		k4x, k4v =  S(x + dt*k3x,     v + dt*k3v,     omega)(1)	, S(x + dt*k3x,     v + dt*k3v,     omega)(2)

		
		WRITE(1, '(F10.4, 1X, F10.4, 1X, F10.4)') t, x, v
		
		t = t + dt
		x = x + (k1x + k2x + k3x + k4x)*dt/6.0
		v = v + (k1v + k2v + k3v + k4v)*dt/6.0
			
		
	END DO
	CLOSE(1)

	
END PROGRAM shm
