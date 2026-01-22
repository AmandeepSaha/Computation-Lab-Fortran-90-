PROGRAM shm
  IMPLICIT NONE

  ! Variables
  REAL :: m, k, b
  REAL :: t, x, v
  REAL :: k1x, k2x, k3x, k4x
  REAL :: k1v, k2v, k3v, k4v
  REAL :: omega
  REAL :: dt
  INTEGER :: i, n
  REAL :: svec(2)

  OPEN(1, file='shm.dat', STATUS='REPLACE')

  PRINT *, "Enter t0, x0, v0, m, k, b: "
  READ(*,*) t, x, v, m, k, b

  omega = SQRT(k/m)
  dt = 0.01
  n  = 2000

  DO i = 0, n

     svec = S(x, v, b, omega)
     k1x = svec(1)
     k1v = svec(2)

     svec = S(x + 0.5*dt*k1x, v + 0.5*dt*k1v, b, omega)
     k2x = svec(1)
     k2v = svec(2)

     svec = S(x + 0.5*dt*k2x, v + 0.5*dt*k2v, b, omega)
     k3x = svec(1)
     k3v = svec(2)

     svec = S(x + dt*k3x, v + dt*k3v, b, omega)
     k4x = svec(1)
     k4v = svec(2)

     WRITE(1, '(F10.4, 1X, F10.4, 1X, F10.4)') t, x, v

     t = t + dt
     x = x + (k1x + 2.0*k2x + 2.0*k3x + k4x)*dt/6.0
     v = v + (k1v + 2.0*k2v + 2.0*k3v + k4v)*dt/6.0

  END DO

  CLOSE(1)

CONTAINS

  FUNCTION func(x, xdot, b, omega) RESULT(d2xdt2)
    REAL, INTENT(IN) :: x, xdot, b, omega
    REAL :: d2xdt2

    d2xdt2 = -b*xdot -omega**2 * x
  END FUNCTION func

  FUNCTION S(x, v, b, omega) RESULT(State)
    REAL, INTENT(IN) :: x, v, b, omega
    REAL :: State(2)

    State(1) = v
    State(2) = func(x, v, b, omega)
  END FUNCTION S

END PROGRAM shm
