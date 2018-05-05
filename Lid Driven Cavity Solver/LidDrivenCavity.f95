! Upwind Difference Scheme with Gauss Seidel!
program LidDrivenCavity

implicit none
integer*8, parameter       			  :: ILX = 1, IUX = 129, ILY = 1, IUY = 129 
integer*8							  :: ICX, ICY 
real*8, dimension(:,:), ALLOCATABLE   :: UOLD, UNEW, VOLD, VNEW, POLD, PNEW 
real*8, dimension(:,:), ALLOCATABLE   :: UOLDP, VOLDP, POLDP, RESAVGM
integer*8                  			  :: I, J, ITER
real*8, parameter          			  :: pi = 4*ATAN(1.0_8)  
real*8, parameter        			  :: LX = 1.0, LY= 1.0, ULID = 1.0
real*8                  			  :: DX, DY, TIMESTOP, NDT
real*8, parameter   				  :: RHO= 1.0, MU= 1.0, RE= 1000.0
real*8, parameter					  :: DT= 0.0076	
real*8								  :: RESAVG, X, Y
real*8                                :: start, finish

RESAVG = 1.0
ITER = 1

DX = LX/(IUX-ILX)
DY = LY/(IUY-ILY)

ICX = (IUX+ILX)/2
ICY = (IUY+ILY)/2

ALLOCATE(UOLD(IUX, IUY))
ALLOCATE(UNEW(IUX, IUY))
ALLOCATE(VOLD(IUX, IUY))
ALLOCATE(VNEW(IUX, IUY))
ALLOCATE(POLD(IUX, IUY))
ALLOCATE(PNEW(IUX, IUY))
ALLOCATE(UOLDP(IUX, IUY))
ALLOCATE(VOLDP(IUX, IUY))
ALLOCATE(POLDP(IUX, IUY))
ALLOCATE(RESAVGM(IUX, IUY))

RESAVGM(:,:) = 0.0

call cpu_time(start)
call initAssign(UOLD, UNEW, VOLD, VNEW, POLD, PNEW, ILX, IUX, ILY, &    
& IUY, UOLDP, VOLDP, POLDP)

call assignBoundary(UOLD, UOLDP, UNEW, VOLD, VNEW, POLD, PNEW, ILX, IUX, ILY, IUY, ULID, UOLDP, VOLDP, POLD)

open(unit = 1, file = 'residuals.dat')
write(1,*) 'variables="Iteration""Non-dimensional Time Units""Residual"'
write(1,*) 'Zone T = "Residual"'
do while ((RESAVG .GE. 1e-8) .AND. (ITER .LE. 25000))

call xmomentum(UOLD, UNEW, UOLDP, VOLD, VNEW, VOLDP, POLD, DX, DY, DT, RHO, MU, RE, ILX, IUX, ILY, IUY, ITER)

!UOLDP = UOLD
!UOLD = UNEW

call ymomentum(UOLD, UNEW, UOLDP, VOLD, VNEW, VOLDP, POLD, DX, DY, DT, RHO, MU, RE, ILX, IUX, ILY, IUY, ITER)

! VOLD = VNEW??

call pressureequation(UOLD, UNEW, VOLD, VNEW, POLD, PNEW, DX, DY, DT, RHO, ILX, IUX, ILY, IUY)

call newvelocity(UOLD, UNEW, VOLD, VNEW, POLD, PNEW, DX, DY, RHO, MU, ILX, IUX, ILY, IUY, DT)

! call adamsBashforth(UOLD, UNEW, VOLD, VNEW, POLD, PNEW, UOLDP, VOLDP, POLDP, ILX, IUX, ILY, IUY, DT)

RESAVG = 0.0

call residual(UOLD, UNEW, VOLD, VNEW, ILX, IUX, ILY, IUY, RESAVG, RESAVGM)

NDT = ITER*DT
write(1,*) ITER ,' ', NDT, '  ', RESAVG
UOLDP = UOLD
VOLDP = VOLD
POLDP = POLD
UOLD = UNEW
VOLD = VNEW
POLD = PNEW
ITER = ITER +1

print *, ITER, ' ', RESAVG
end do
call cpu_time(finish)
print '("Time = ",f16.6," seconds.")',finish-start

open(unit = 2, file = 'solutionContour.dat')
write(2,*) 'TITLE="Error Field Data"'
write(2,*) 'variables="X(m)""Y(m)""Pressure""U- velocity""V- velocity""Residuals"'
write(2,*) 'Zone T = "n=1"'
write(2,*) 'I =', IUX, 'J =', IUY
write(2,*) 'DATAPACKING = POINT'
do J = ILY, IUY  
do I = ILX, IUX
X = DX*(I-1)
Y = DY*(J-1)
write(2,*) X, ' ', Y, ' ', PNEW(I, J), ' ', UNEW(I,J), ' ', VNEW(I,J), ' ', RESAVGM(I,J)
end do
end do

open(unit = 3, file = 'centerVelocity.dat')
write(3,*) 'variables="X(m)""U- velocity""V- velocity"'
write(3,*) 'Zone T = "CenterVelocity', IUX, '"'
do I=ILX, IUX
X = DX*(I-1)
write(3,*) X, ' ', UNEW(ICX,I), ' ', VNEW(I,ICY)
end do


end program LidDrivenCavity
