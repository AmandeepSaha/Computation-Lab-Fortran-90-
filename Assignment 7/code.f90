program coupled_differential_eqn
implicit none

! ---------- declerations ------------

real :: m1, m2, k
real :: t0, x1, x2, v1, v2
real :: k1x1, k2x1, k3x1, k4x1
real :: k1x2, k2x2, k3x2, k4x2
real :: k1v1, k2v1, k3v1, k4v1
real :: k1v2, k2v2, k3v2, k4v2
real :: dt
real :: S(4)
integer :: i, n

open(1, file='datas.dat', status='replace')

print *, "Enter values of t0, m1, m2, k, x1, x2, v1, v2: "
read(*,*) t0, m1, m2, k, x1, x2, v1, v2
n = 1000000
dt = 1.0e-4
! ---------- executable ------------

do i = 1, n
	
	S = statevector(t, m1, m2, k, x1, x2, v1, v2)
	k1x1= S(1)
	k1x2 = S(2)
	k1v1 = S(3)
	k1v2 = S(4)
	
	S = statevector(t+ 0.5 * dt, m1, m2, k, x1 + 0.5 * dt * k1x1 , x2 + 0.5 * dt * k1x2, v1 + 0.5 * dt * k1v1, v2 + 0.5 * dt * k1v2)
	
	k2x1= S(1)
	k2x2 = S(2)
	k2v1 = S(3)
	k2v2 = S(4)
	
	S = statevector(t+ 0.5 * dt, m1, m2, k, x1 + 0.5 * dt * k2x1 , x2 + 0.5 * dt * k2x2, v1 + 0.5 * dt * k2v1, v2 + 0.5 * dt * k2v2)

	k3x1= S(1)
	k3x2 = S(2)
	k3v1 = S(3)
	k3v2 = S(4)
	
	S = statevector(t+  dt, m1, m2, k, x1 +  dt * k2x1 , x2 +  dt * k2x2, v1 +  dt * k2v1, v2 +  dt * k2v2)
	
	k4x1= S(1)
	k4x2 = S(2)
	k4v1 = S(3)
	k4v2 = S(4)
	
	write(1,'(5F12.6)') t, x1, v1, x2, v2

	t = t + dt
	x1 = x1 + (k1x1 + 2.0*k2x1 + 2.0*k3x1 + k4x1) * dt/6
	x2 = x2 + (k1x2 + 2.0*k2x2 + 2.0*k3x2 + k4x2) * dt/6
	v1 = v1 + (k1v1 + 2.0*k2v1 + 2.0*k3v1 + k4v1) * dt/6
	v2 = v2 + (k1v2 + 2.0*k2v2 + 2.0*k3v2 + k4v2) * dt/6


end do

CALL SYSTEM("gnuplot -e ""set multiplot layout 3,2; " // &
    "set grid; " // &
     "plot 'datas.dat' u 1:2 w l t 'x1(t)'; " // &
     "plot 'datas.dat' u 1:3 w l t 'v1(t)'; " // &
     "plot 'datas.dat' u 2:3 w l t 'v1(t)'; " // &
     "plot 'datas.dat' u 1:4 w l t 'x2(t)'; " // &
     "plot 'datas.dat' u 1:5 w l t 'v2(t)'; " // &
     "plot 'datas.dat' u 4:5 w l t 'v2(t)'; " // &
     "pause -1""")

contains
! ------------- d2fdt2 carries equation of motion in second order differential form, of two differentent mass systems ----------------
function d2fdt2(m1, m2, k, x1, x2) result diffs
	real intent(in) ::  m1, k, x1, x2
	real :: diffs(2)
	
	diffs(1) = -2*(k/m)*x1 + k*x2
	diffs(2) = -2*(k/m)*x2 + k*x1

function statevector(t, m1, m2, k, x1, x2, v1, v2) result S
	real intent(in) :: t, m1, m2, k, x1, x2, v1, v2
	real :: S(4)
	
	! ------instance ------
	acc_states = d2fdt2(m1, m2, k, x1, x2)
	
	
	S(1) = v1
	S(2) = v2
	S(3) = acc_states(1)
	S(4) = acc_states(2)
	
end program coupled_differential_eqn
