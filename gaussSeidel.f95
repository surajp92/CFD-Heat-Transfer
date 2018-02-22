program gaussSeidel 

USE assignArray
USE initialResidual
USE pointGaussSeidel

implicit none
real*8                             :: start, finish
real*8, parameter                  :: DX = 0.0125, DY = 0.0125
integer*8, parameter               :: ILX = 1,ILY = 1,ILAX = 1, ILAY = 1 
integer*8                          :: N, NCOUNT, IUX, IUY , IUAX , IUAY
integer*8                          :: I, J, ITER = 0, ntemp 
real*8                             :: RES = 0, REST, RESAVG, EPSIAVG = 0
real*8, dimension (:), allocatable :: PHIA, PHIN, PHIO, X, Y, ERROR, S
real*8                             :: AP, AE, AW, AN, AS


IUX = IDNINT(1/DX) + 1.0
IUY = IDNINT(1.0/DY) + 1.0
IUAX = IDNINT(1.0/DX) + 1.0
IUAY = IDNINT(1.0/DX) + 1.0
N = (IUY-1.0)*IUX + IUX
! print *, N

allocate(S(1:N))
allocate(PHIA(1:N))
allocate(PHIN(1:N))
allocate(PHIO(1:N))
allocate(X(1:N))
allocate(Y(1:N))
allocate(Error(1:N))

call cpu_time(start)
! Array Element Assignment
ITER = 0
RESAVG = 1
call arrayAssignmnt(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, PHIN, & 
                    PHIO, PHIA, X, Y, DX, DY, N)

! Initial Residual
call inRES(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, &
           PHIO, RES, REST, RESAVG, N)
           
! Successive Order Relaxation
call pointGS(ILX, IUX, ILY, IUY, ITER, AP, AE, AW, AN, AS, S, &
                   PHIO, PHIN, RES, REST, RESAVG, N, EPSIAVG)

! print *, ITER
             
call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start

! Writing th final solution
open(unit = 4, file = 'solution.dat')
write(4,*), 'Node', ' ', 'Analytical', ' ','Numerical' 
do I = 1, N
ERROR(I) = ABS(100 * (PHIA(I) - PHIN(I))/PHIA(I))
write(4,*), I, ' ', PHIA(I), ' ', PHIN(I), ' ', ERROR(I)
end do

! Writing the Solution at each node location (Contour Plot)
NCOUNT = N
open(unit = 5, file = 'solutionContour.dat')
write(5,*) 'TITLE="Error Field Data"'
write(5,*) 'variables="X(m)""Y(m)""Analytical""Numerical""Error"'
write(5,*) 'Zone T = "n=1"'
write(5,*) 'I =', IUX, 'J =', IUY
write(5,*) 'DATAPACKING = POINT'
do while (NCOUNT .GE. 1)
do I = IUX, ILX, -1
write(5,*) X(NCOUNT), ' ', Y(NCOUNT), ' ', PHIA(I), ' ', PHIN(I), ' ', ERROR(NCOUNT)
NCOUNT =NCOUNT - 1
end do
end do
end program gaussSeidel

!!!!!! Module that assigns the coefficients of System Matrix!!!!!!

module assignArray

implicit none
contains
subroutine arrayAssignmnt(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, &
                          PHIN, PHIO, PHIA, X, Y, DX, DY, N)

integer*8                         :: ILX, IUX, ILY, IUY, I, J, ntemp, N
real*8                            :: DX, DY
real*8, dimension(:), allocatable :: PHIN, PHIO, PHIA, X, Y, S
real*8                            :: AP, AE, AW, AN, AS

AP = 2.0*(1/(DX**2) + 1/(DY**2))
AE = 1.0/(DX**2)
AW = 1.0/(DX**2)
AN = 1.0/(DY**2)
AS = 1.0/(DY**2)

do J = ILY, IUY
do I = ILX, IUX
    ntemp = (J-1)*IUX + I
    X(ntemp) = DX*(I-1)
    Y(ntemp) = DY*(J-1)
    PHIA(ntemp) = 500*EXP(-50*(((1-X(ntemp))**2)+(Y(ntemp))**2)) + 100*(1-Y(ntemp))*X(ntemp)
    if (J == 1 .or. J == IUY .or. I == 1 .or. I == IUX) then
        PHIO(ntemp) =  PHIA(ntemp)
        S(ntemp) = 0            
    else 
        PHIO(ntemp) =  0
        S(ntemp) = -50000*(100*(((1-X(ntemp))**2+(Y(ntemp))**2))-2)*EXP(-50*(((1-X(ntemp))**2)+(Y(ntemp))**2))
    end if    
end do
end do

AP = 2*(1/(DX**2) + 1/(DY**2))
AE = 1/(DX**2)
AW = 1/(DX**2)
AN = 1/(DY**2)
AS = 1/(DY**2)
PHIN = PHIO

end subroutine arrayAssignmnt
end module assignArray

!!!!!! Module that calculates the residual at the start of iteration !!!!!!
module initialResidual
implicit none
contains 
subroutine inRES(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, &
                 PHIO, RES, REST, RESAVG, N)
                           
integer*8                          :: I, J, N, ILX, IUX, ILY, IUY, ntemp_P
real*8                             :: REST, RES, RESAVG 
real*8, dimension(:), allocatable  ::  PHIO, S
real*8                             :: AP, AE, AW, AN, AS

do J = ILY+1, IUY-1
do I = ILX+1, IUX-1
    ntemp_P = (J-1)*IUX + I
    REST = (-1*S(ntemp_p))-(AE*PHIO(ntemp_P+1)+AW*PHIO(ntemp_P-1) & 
                               +AN*PHIO(ntemp_P+IUX)+AS*PHIO(ntemp_P-IUX) &
                               - AP*PHIO(ntemp_P))
    RES = RES + REST**2
end do
end do
RESAVG = (RES/N)**0.5
end subroutine inRES
end module initialResidual

!!!!!! Solution to Linear System using Point Gauss Seidel Method !!!!!!
module pointGaussSeidel
implicit none
contains
subroutine pointGS(ILX, IUX, ILY, IUY, ITER, AP, AE, AW, AN, AS, S, &
                   PHIO, PHIN, RES, REST, RESAVG, N,EPSIAVG)
integer*8               :: I, J, ITER, N, ILX, IUX, ILY, IUY, ntemp_P
real*8                  :: REST, RES, RESAVG, EPSIAVGT = 0, EPSIAVG
real*8                  :: EPSIAVGK, EGV 
real*8, dimension(:), allocatable :: PHIO, PHIN, S
real*8                            :: AP, AE, AW, AN, AS

open(unit = 2, file = 'residuals.dat')
EPSIAVG = 1
! Point Gauss Seidel
do while (RESAVG .GE. 1e-8)
PHIO = PHIN
EPSIAVGT = 0
ITER = ITER+1
RES = 0 
do J = ILY+1, IUY-1
do I = ILX+1, IUX-1
    ntemp_P = (J-1)*IUX + I
    PHIN(ntemp_P) = (1/AP)*(AE*PHIO(ntemp_P+1) + AW*PHIN(ntemp_P-1) & 
                    + AN*PHIO(ntemp_P+IUX) + AS*PHIN(ntemp_P-IUX) + S(ntemp_P))
end do
end do
do J = ILY+1, IUY-1
do I = ILX+1, IUX-1
    ntemp_P = (J-1)*IUX + I
    REST = (-1*S(ntemp_p))-(AE*PHIN(ntemp_P+1)+AW*PHIN(ntemp_P-1) & 
                               +AN*PHIN(ntemp_P+IUX)+AS*PHIN(ntemp_P-IUX) &
                               - AP*PHIN(ntemp_P))
    RES = RES + REST**2
end do
end do
RESAVG = (RES/(N-2.0*IUX-2.0*IUY+4.0))**0.5
do I = 1,N
EPSIAVGT = EPSIAVGT+(PHIN(I) - PHIO(I))**2
end do
EPSIAVGK = EPSIAVG
EPSIAVG  = (EPSIAVGT/(N-2.0*IUX-2.0*IUY+4.0))**0.5
EGV = EPSIAVG/EPSIAVGK
write(2,*) ITER , '  ', RESAVG, ' ', EPSIAVG, ' ',EGV
end do 
end subroutine pointGS

end module pointGaussSeidel
