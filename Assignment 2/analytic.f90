PROGRAM DampedOscillationAnalytic
IMPLICIT NONE

REAL    :: t, x, v
REAL    :: m, b_, b, k, omega0, omegad
REAL    :: x0, v0
REAL    :: r1, r2, C1, C2
REAL    :: dt, t_max
INTEGER :: i, n

OPEN(1, FILE="underdamped.dat",  STATUS="REPLACE")
OPEN(2, FILE="overdamped.dat",   STATUS="REPLACE")
OPEN(3, FILE="criticaldamped.dat", STATUS="REPLACE")

PRINT *, "Enter m, b_, k, x0, v0 :"
READ(*,*) m, b_, k, x0, v0

b       = b_ / (2.0*m)
omega0 = SQRT(k/m)

dt    = 0.01
t_max = 20.0
n     = INT(t_max/dt)

IF (b < omega0) THEN
    ! -------- UNDERDAMPED --------
    omegad = SQRT(omega0*omega0 - b*b)

    DO i = 1, n
        t = i * dt

        x = EXP(-b*t) * ( x0*COS(omegad*t) + &
            (v0 + b*x0)/omegad * SIN(omegad*t) )

        v = EXP(-b*t) * ( -b*x0*COS(omegad*t) - &
            (b*(v0 + b*x0)/omegad)*SIN(omegad*t) - &
            x0*omegad*SIN(omegad*t) + &
            (v0 + b*x0)*COS(omegad*t) )

        WRITE(1,'(F10.4,1X,F10.4,1X,F10.4)') t, x, v
    END DO

ELSE IF (ABS(b - omega0) < 1.0E-6) THEN
    ! -------- CRITICALLY DAMPED --------
    DO i = 1, n
        t = i * dt

        x = EXP(-b*t) * ( x0 + (v0 + b*x0)*t )
        v = EXP(-b*t) * ( v0 - b*(v0 + b*x0)*t )

        WRITE(3,'(F10.4,1X,F10.4,1X,F10.4)') t, x, v
    END DO

ELSE
    ! -------- OVERDAMPED --------
    r1 = -b + SQRT(b*b - omega0*omega0)
    r2 = -b - SQRT(b*b - omega0*omega0)

    C1 = (v0 - r2*x0) / (r1 - r2)
    C2 = (r1*x0 - v0) / (r1 - r2)

    DO i = 1, n
        t = i * dt

        x = C1*EXP(r1*t) + C2*EXP(r2*t)
        v = C1*r1*EXP(r1*t) + C2*r2*EXP(r2*t)

        WRITE(2,'(F10.4,1X,F10.4,1X,F10.4)') t, x, v
    END DO

END IF

CALL SYSTEM("gnuplot -e ""set multiplot layout 3,1; " // &
            "set grid; " // &
            "plot 'underdamped.dat' u 1:2 w l t 'x(t)'; " // &
            "plot 'underdamped.dat' u 1:3 w l t 'v(t)'; " // &
            "plot 'underdamped.dat' u 2:3 w l t 'x–v'; " // &
            "pause -1""")


END PROGRAM DampedOscillationAnalytic
