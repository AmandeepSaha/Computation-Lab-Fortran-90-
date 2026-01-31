! Write the Lagrangian for a one-dimensional quartic oscillator and obtain its equation of motion. Use the plot of the potential function to obtain a schematic phase portrait. Write a program to obtain the phase portrait numerically and compare it with the schematic diagram you obtained.
! You may consult the book ' Nonlinear dynamics and chaos by Steven H. Strogatz for drawing the phase portrait.

PROGRAM cubic
  IMPLICIT NONE

  ! Variables
  REAL :: m, k
  REAL :: lambda
  REAL :: t, x, v
  REAL :: k1x, k2x, k3x, k4x
  REAL :: k1v, k2v, k3v, k4v
  REAL :: omega
  REAL :: dt
  INTEGER :: i, n
  REAL :: svec(2)
  
  
  OPEN(1, file='cubic_oscillator.dat', STATUS='REPLACE')

  PRINT *, "Enter t0, x0, v0, m, k, lambda: "
  READ(*,*) t, x, v, m, k, lambda

  omega = SQRT(k/m)
  dt = 0.01
  n  = 2000

  DO i = 0, n

     svec = S(x, v, k, m, lambda)
     k1x = svec(1)
     k1v = svec(2)

     svec = S(x + 0.5*dt*k1x, v + 0.5*dt*k1v, k, m, lambda)
     k2x = svec(1)
     k2v = svec(2)

     svec = S(x + 0.5*dt*k2x, v + 0.5*dt*k2v, k, m, lambda)
     k3x = svec(1)
     k3v = svec(2)

     svec = S(x + dt*k3x, v + dt*k3v, k, m, lambda)
     k4x = svec(1)
     k4v = svec(2)

     WRITE(1, '(F10.4, 1X, F10.4, 1X, F10.4)') t, x, v

     t = t + dt
     x = x + (k1x + 2.0*k2x + 2.0*k3x + k4x)*dt/6.0
     v = v + (k1v + 2.0*k2v + 2.0*k3v + k4v)*dt/6.0

  END DO

  CLOSE(1)

CALL SYSTEM("gnuplot -e ""set multiplot layout 3,1; " // &
    "set grid; " // &
     "plot 'cubic_oscillator.dat' u 1:2 w l t 'x(t)'; " // &
     "plot 'cubic_oscillator.dat' u 1:3 w l t 'v(t)'; " // &
     "plot 'cubic_oscillator.dat' u 2:3 w l t 'v(t)'; " // &
     "pause -1""")


CONTAINS

  FUNCTION func(x, v, k, m, lambda) RESULT(d2xdt2)
    REAL, INTENT(IN) :: x, v, k, m, lambda
    REAL :: d2xdt2

    d2xdt2 = (-k * x - 4.0 * lambda * x**3)/m
  END FUNCTION func

  FUNCTION S(x, v, k, m, lambda) RESULT(State)
    REAL, INTENT(IN) :: x, v, k, m, lambda
    REAL :: State(2)

    State(1) = v
    State(2) = func(x, v, k, m, lambda)
  END FUNCTION S
  

END PROGRAM cubic
