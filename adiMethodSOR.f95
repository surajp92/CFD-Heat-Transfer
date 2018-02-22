program gaussSeidel 
USE assignArray
USE initialResidual
USE lineADI

implicit none
real*8                             :: start, finish
real*8, parameter                  :: DX = 1.0/80.0, DY = 1.0/80.0
integer*8, parameter               :: ILX = 1,ILY = 1,ILAX = 1, ILAY = 1
integer*8                          :: IUX , IUY , IUAX , IUAY , N, NCOUNT
integer*8                          :: I, J, ITER = 0, ntemp 
real*8                             :: RES = 0, REST, RESAVG, EPSIAVG = 0
real*8, dimension (:), allocatable :: PHIA, PHIN, PHIO, X, Y, ERROR, S
real*8                             :: AP, AE, AW, AN, AS, W, EGV


IUX = IDNINT(1/DX) + 1.0
IUY = IDNINT(1.0/DY) + 1.0
IUAX = IDNINT(1.0/DX) + 1.0
IUAY = IDNINT(1.0/DX) + 1.0
N = (IUY-1.0)*IUX + IUX

allocate(S(1:N))
allocate(PHIA(1:N))
allocate(PHIN(1:N))
allocate(PHIO(1:N))
allocate(X(1:N))
allocate(Y(1:N))
allocate(Error(1:N))

W = 1.0
open(unit = 3, file = 'iterationsSORADI.dat')
call cpu_time(start)
! Array Element Assignment
do while (W .LE. 2.0)
ITER = 0
RESAVG = 1
RES = 0
call arrayAssignmnt(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, PHIN, & 
                    PHIO, PHIA, X, Y, DX, DY, N)

! Initial Residual
call inRES(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, &
           PHIO, RES, REST, RESAVG, N)
           
! Successive Order Relaxation
call ADI(ILX, IUX, ILY, IUY, ITER, AP, AE, AW, AN, AS, S, &
                   PHIO, PHIN, RES, REST, RESAVG, N, EPSIAVG, W, EGV)

write(3,*), W , '  ', ITER, ' ', RESAVG, ' ', EGV
print *, W , '  ', ITER !, ' ', RESAVG
W = W+0.01
end do
             
call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start
end program gaussSeidel

!!!!!! Module to Assign coefficients of System matrix !!!!!!
module assignArray
implicit none
contains
subroutine arrayAssignmnt(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, &
                          PHIN, PHIO, PHIA, X, Y, DX, DY, N)

integer*8       :: ILX, IUX, ILY, IUY, I, J, ntemp, N
real*8          :: DX, DY
real*8, dimension(:), allocatable :: PHIN, PHIO, PHIA, X, Y, S
real*8                            :: AP, AE, AW, AN, AS

AP = 2*(1/(DX**2) + 1/(DY**2))
AE = 1/(DX**2)
AW = 1/(DX**2)
AN = 1/(DY**2)
AS = 1/(DY**2)

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

!!!!!! Module to calculate the initial residual at the start of iterations !!!!!!!
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

!!!!!! Module to solve linear system of equation using Alternate Direct Implicit method and Successive 
!!!!!! overrelaxation factor ranging from 1 to 2 !!!!!!
module lineADI
implicit none
contains
subroutine ADI(ILX, IUX, ILY, IUY, ITER, AP, AE, AW, AN, AS, S, &
                   PHIO, PHIN, RES, REST, RESAVG, N,EPSIAVG, W, EGV)
integer*8               :: I, J, ITER, N, ILX, IUX, ILY, IUY, ntemp_P
real*8                  :: REST, RES, RESAVG, EPSIAVGT = 0, EPSIAVG, W, EGV 
real*8                  :: EPSIAVGK 
real*8, dimension(:), allocatable :: PHIO, PHIN, S
real*8                            :: AP, AE, AW, AN, AS
real*8, dimension(ILX:IUX)    :: BB, DD, AA, CC
integer*8                         :: IL, IU
EPSIAVG = 1.0
! Alternate Direct Implicit Method
do while (RESAVG .GE. 1e-8 .AND. ITER .LE. 20000)
PHIO = PHIN
EPSIAVGT = 0
ITER = ITER+1
RES = 0 
! Sweep all J 
IL = ILX
IU = IUX
do J = ILY+1, IUY-1
do I = IL, IU
    ntemp_P = (J-1)*IUX + I
    if (I == IL) then
        BB(I) = 0.0
        DD(I) = 1.0
        AA(I) = 0.0
        CC(I) = PHIO(ntemp_P)
    else if (I == IU) then
        BB(I) = 0.0
        DD(I) = 1.0
        AA(I) = 0.0
        CC(I) = PHIO(ntemp_P)
    else
        BB(I) = -W*AW
        DD(I) = AP
        AA(I) = -W*AE
        CC(I) = AP*(1-W)*PHIO(ntemp_P) + W*AN*PHIO(ntemp_P+IUX)+W*AS*PHIN(ntemp_P-IUX) + W*S(ntemp_P)
    end if
end do
call tdmaSolver(IL, IU, BB, DD, AA, CC)

do I = IL, IU
    ntemp_P = (J-1)*IUX + I
    PHIN(ntemp_P) = CC(I)
end do
end do
PHIO = PHIN
! Sweep all I
IL = ILY
IU = IUY
do I = ILX+1, IUX-1
do J = IL, IU
    ntemp_P = (J-1)*IUX + I
    if (J == IL) then
        BB(J) = 0.0
        DD(J) = 1.0
        AA(J) = 0.0
        CC(J) = PHIO(ntemp_P)
    else if (J == IU) then
        BB(J) = 0.0
        DD(J) = 1.0
        AA(J) = 0.0
        CC(J) = PHIO(ntemp_P)
    else        
        BB(J) = -W*AS
        DD(J) = AP
        AA(J) = -W*AN
        CC(J) = AP*(1-W)*PHIO(ntemp_P) + W*AE*PHIO(ntemp_P+1)+W*AW*PHIN(ntemp_P-1) + W*S(ntemp_P)
    end if
end do
call tdmaSolver(IL, IU, BB, DD, AA, CC)
do J = IL, IU
    ntemp_P = (J-1)*IUX + I
    PHIN(ntemp_P) = CC(J)
end do
end do

! Calculation of Residual
do J = ILY+1, IUY-1
do I = ILX+1, IUX-1
    ntemp_P = (J-1)*IUX + I
    REST = (-1*S(ntemp_P))-(AE*PHIN(ntemp_P+1)+AW*PHIN(ntemp_P-1) & 
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
end do
end subroutine ADI
end module lineADI

!!!!!! TDMA solver using Thomas Algorithm !!!!!!
subroutine tdmaSolver(IL, IU, BB, DD, AA, CC)

implicit none
integer*8 :: IL, IU, I, J, LP
real*8 :: R
real*8, dimension(IL:IU) :: BB, DD, AA, CC
LP = IL + 1
! formation of upper triangular matrix
do I = LP, IU
    R = BB(I)/DD(I-1)
    DD(I)  = DD(I)-R*AA(I-1)
    CC(I) = CC(I)-R*CC(I-1)
end do
! Back substitution
CC(IU) = CC(IU)/ DD(IU)
do I = LP, IU
    J = IU-I+LP
    CC(J) = (CC(J)-AA(J)*CC(J+1))/DD(J)
end do

end subroutine tdmaSolver
