! nclude a cubic term in the Lagrangian for a simple harmonic oscillator to make it an anharmonic oscillator. Write the equation of motion and solve it using the method of successive approximation to show the shift in its equilibrium position. Consider the program you wrote for the SHM problem and modify it suitably to get x(t), v(t) and the phase portrait. Also, verify the expression for shifting in the equilibrium position as a function of the amplitude.
! Note: Consult  Course of Theoretical Physics-(Vol-1: Mechanics) by Landau and Lifshitz.

PROGRAM cubic
  IMPLICIT NONE

  ! Variables
  REAL :: m, k
  REAL :: a, b
  REAL :: t, x, v
  REAL :: k1x, k2x, k3x, k4x
  REAL :: k1v, k2v, k3v, k4v
  REAL :: omega
  REAL :: dt
  INTEGER :: i, n
  REAL :: svec(2)
  
  b = 0
  
  OPEN(1, file='cubic_oscillator.dat', STATUS='REPLACE')
  PRINT *, "We added cubic terms of all independent variables of the lagrangian, although we are solving it for 1D system."
  PRINT *, "Enter t0, x0, v0, m, k, a: "
  READ(*,*) t, x, v, m, k, a

  omega = SQRT(k/m)
  dt = 0.01
  n  = 2000

  DO i = 0, n

     svec = S(x, v, k, m, a, b)
     k1x = svec(1)
     k1v = svec(2)

     svec = S(x + 0.5*dt*k1x, v + 0.5*dt*k1v, k, m, a, b)
     k2x = svec(1)
     k2v = svec(2)

     svec = S(x + 0.5*dt*k2x, v + 0.5*dt*k2v, k, m, a, b)
     k3x = svec(1)
     k3v = svec(2)

     svec = S(x + dt*k3x, v + dt*k3v, k, m, a, b)
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

  FUNCTION func(x, v, k, m, a, b) RESULT(d2xdt2)
    REAL, INTENT(IN) :: x, v, k, m, a, b
    REAL :: d2xdt2

    d2xdt2 = (k * x - 3.0 * a * x**2)/(m - 6.0 * b * v)
  END FUNCTION func

  FUNCTION S(x, v, k, m, a, b) RESULT(State)
    REAL, INTENT(IN) :: x, v, k, m, a, b
    REAL :: State(2)

    State(1) = v
    State(2) = func(x, v, k, m, a, b)
  END FUNCTION S

END PROGRAM cubic
