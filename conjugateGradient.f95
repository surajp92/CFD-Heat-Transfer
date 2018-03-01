program conjugateGradient 
USE assignArray
USE initialResidual
USE vecInnProd
USE matVecProd


implicit none
real*8                             :: start, finish
real*8, parameter                  :: DX = 0.0125, DY = 0.0125
integer*8, parameter               :: ILX = 1,ILY = 1,ILAX = 1, ILAY = 1 ! , IUX = 21 , IUY=21 , IUAX=21 , IUAY=21
integer*8                          :: N, NCOUNT, IUX, IUY , IUAX , IUAY, ntemp_P
integer*8                          :: I, J, ITER = 0, ntemp 
real*8                             :: REST, RESAVG, EPSIAVG = 0.0, RHOOLD, RHONEW, PTQ, EPSIAVGP=1.0
real*8, dimension (:), allocatable :: PHIA, PHIN, PHIO, X, Y, ERROR, S, RES, Z, P, Q
real*8                             :: AP, AE, AW, AN, AS, BETA, ALFA, EGV 


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
allocate(ERROR(1:N))
allocate(RES(1:N))
allocate(Z(1:N))
allocate(P(1:N))
allocate(Q(1:N))

do I = 1,N 
    RES(I)=0.0
    Q(I)=0.0
    Z(I)=0.0
end do

open(unit = 2, file = 'residuals.dat')
write(4,*), 'Ieration', ' ', 'AVG Resiual', ' ','AVG Error' 
EPSIAVGP = 1.0
call cpu_time(start)
! Array Element Assignment
ITER = 0
RESAVG = 0
call arrayAssignmnt(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, PHIN, & 
                    PHIO, PHIA, X, Y, DX, DY, N)

! Initial Residual
call inRES(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, &
           PHIO, RES, RESAVG, N)

do while (RESAVG .GE. 1e-8 )
! solce Mz = r, M = preconditioner Matrix, M= I for present code
    PHIO = PHIN
    Z = RES
    call vip(RHONEW, RES, Z, N)
    ITER = ITER +1
    if (ITER==1) then
            P = Z
! old = (i-2) iteration and new = (i-1) iteration
        RHOOLD = RHONEW
    else
        BETA = RHONEW/RHOOLD
        P = Z + BETA*P
        RHOOLD = RHONEW
    end if
! Calculate q = Ap by Matrix Vector Multiplication
    call mvp(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, P, Q)
! Calculate transpose(P)*Q
    call vip(PTQ, P, Q, N)
! Calculate alfa = rhoold/ (transpose(P)*Q)
    alfa = RHONEW/PTQ
! Calculate new solution
    PHIN = PHIO + alfa*P
! Calculate new residual
    RES = RES - alfa*Q
    
    RESAVG = 0.0 
    EPSIAVG = 0.0
    ERROR = PHIA - PHIN
    do I = 1, N 
        RESAVG = RESAVG + (RES(I))**2
        EPSIAVG = EPSIAVG + (ERROR(I))**2
    end do
    RESAVG = ((RESAVG/(N-2.0*IUX-2.0*IUY+4.0)))**0.5
    EPSIAVG = ((EPSIAVG/(N-2.0*IUX-2.0*IUY+4.0)))**0.5
    EGV = EPSIAVG/EPSIAVGP
    EPSIAVGP = EPSIAVG
    write(2,*) ITER , '  ', RESAVG, ' ', EPSIAVG, ' ',EGV
    print *, ITER, ' ', RESAVG, ' ', EPSIAVG, ' ', EGV
end do
        
call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start

! Writing th final solution
open(unit = 4, file = 'solution.dat')
write(4,*), 'Node', ' ', 'Analytical', ' ','Numerical' 
do I = 1, N
ERROR(I) = ABS(100 * (PHIA(I) - PHIN(I))/PHIA(I))
write(4,*), I, ' ', PHIA(I), ' ', PHIN(I), ' ', ERROR(I)
end do

!! Writing the Solution at each node location (Contour Plot)
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
end program conjugateGradient

!########## Module to assign all the Coefficients ##########!

module assignArray
implicit none
contains
subroutine arrayAssignmnt(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, &
                          PHIN, PHIO, PHIA, X, Y, DX, DY, N)

integer*8       :: ILX, IUX, ILY, IUY, I, J, ntemp, N
real*8          :: DX, DY
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
        S(ntemp) = 0.0            
    else 
        PHIO(ntemp) =  0.0
        S(ntemp) = -50000*(100*(((1-X(ntemp))**2+(Y(ntemp))**2))-2)*EXP(-50*(((1-X(ntemp))**2)+(Y(ntemp))**2))
    end if    
end do
end do
PHIN = PHIO
end subroutine arrayAssignmnt
end module assignArray

!########## Module to Calculate Initial Residual ##########!

module assignArray
implicit none
contains
subroutine arrayAssignmnt(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, &
                          PHIN, PHIO, PHIA, X, Y, DX, DY, N)

integer*8       :: ILX, IUX, ILY, IUY, I, J, ntemp, N
real*8          :: DX, DY
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
        S(ntemp) = 0.0            
    else 
        PHIO(ntemp) =  0.0
        S(ntemp) = -50000*(100*(((1-X(ntemp))**2+(Y(ntemp))**2))-2)*EXP(-50*(((1-X(ntemp))**2)+(Y(ntemp))**2))
    end if    
end do
end do
PHIN = PHIO
end subroutine arrayAssignmnt
end module assignArray

!########## Module to calculate inner product between two vectors ##########!
 
module vecInnProd
implicit none
contains 
subroutine vip(ATB, A, B, N)

real*8, dimension(:), allocatable :: A, B
real*8                            :: ATB
integer*8                         :: N, I

ATB = 0.0

do I = 1,N
        ATB = ATB + A(I)*B(I)
end do

end subroutine vip
end module vecInnProd

!########## Module to calculate Matrix- Vector Product ##########!

module matVecProd
implicit none
contains 
subroutine mvp(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, P, Q)

integer*8                          :: I, J, ILX, IUX, ILY, IUY, ntemp_P 
real*8, dimension(:), allocatable  :: P, Q
real*8                             :: AP, AE, AW, AN, AS

do J = ILY+1, IUY-1
do I = ILX+1, IUX-1
    ntemp_P = (J-1)*IUX + I
    Q(ntemp_P) = AP*P(ntemp_P)-AN*P(ntemp_P+IUX)-AS*P(ntemp_P-IUX)-AE*P(ntemp_P+1)-AW*P(ntemp_P-1)
end do
end do
end subroutine mvp
end module matVecProd
