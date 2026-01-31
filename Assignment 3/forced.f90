PROGRAM forced
  IMPLICIT NONE

  ! Variables
  REAL :: m, k, b
  REAL :: t, x, v
  REAL :: a
  REAL :: k1x, k2x, k3x, k4x
  REAL :: k1v, k2v, k3v, k4v
  REAL :: omega
  REAL :: dt
  INTEGER :: i, n
  REAL :: svec(2)

  OPEN(1, file='datas.dat', STATUS='REPLACE')

  PRINT *, "Enter t0, x0, v0, m, k, b, a: "
  READ(*,*) t, x, v, m, k, b, a

  omega = SQRT(k/m)
  dt = 0.01
  n  = 2000

  DO i = 0, n

     svec = S(x, v, b, omega, a)
     k1x = svec(1)
     k1v = svec(2)

     svec = S(x + 0.5*dt*k1x, v + 0.5*dt*k1v, b, omega, a)
     k2x = svec(1)
     k2v = svec(2)

     svec = S(x + 0.5*dt*k2x, v + 0.5*dt*k2v, b, omega, a)
     k3x = svec(1)
     k3v = svec(2)

     svec = S(x + dt*k3x, v + dt*k3v, b, omega, a)
     k4x = svec(1)
     k4v = svec(2)

     WRITE(1, '(F10.4, 1X, F10.4, 1X, F10.4)') t, x, v

     t = t + dt
     x = x + (k1x + 2.0*k2x + 2.0*k3x + k4x)*dt/6.0
     v = v + (k1v + 2.0*k2v + 2.0*k3v + k4v)*dt/6.0

  END DO

! ---------- gnuplot ------------

CALL system( &
"gnuplot -e ""set multiplot layout 3,1; set grid; " // &
"plot 'datas.dat' u 1:2 w l t 'x(t)'; " // &
"plot 'datas.dat' u 1:3 w l t 'v(t)'; " // &
"plot 'datas.dat' u 2:3 w l t 'phase'; pause -1""")
  CLOSE(1)

CONTAINS

  FUNCTION func(x, xdot, b, omega, a) RESULT(d2xdt2)
    REAL, INTENT(IN) :: x, xdot, b, omega, a
    REAL :: d2xdt2

    d2xdt2 = -b*xdot -omega**2 * x + a*SIN(x)
  END FUNCTION func

  FUNCTION S(x, v, b, omega, a) RESULT(State)
    REAL, INTENT(IN) :: x, v, b, omega, a
    REAL :: State(2)

    State(1) = v
    State(2) = func(x, v, b, omega, a)
  END FUNCTION S

END PROGRAM forced
