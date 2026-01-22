PROGRAM analytic_shm
  IMPLICIT NONE

  REAL :: t, x, v
  REAL :: dt, tmax
  INTEGER :: i, n

  OPEN(1, file='analytic_shm.dat', STATUS='REPLACE')

  dt   = 0.01
  tmax = 10.0
  n    = INT(tmax/dt)

  DO i = 0, n
     t = i * dt

     x = 0.447 * SIN(2.23 * t)
     v = 0.99  * COS(2.23 * t)

     WRITE(1,'(F10.4,1X,F10.4,1X,F10.4)') t, x, v
  END DO

  CLOSE(1)

  PRINT *, "Data written to analytic_shm.dat"

END PROGRAM analytic_shm
