program multiGrid 
USE assignArray
USE calculateResidual
USE pointGaussSeidel
USE restrictionOper
USE prolongationOper

implicit none
real*8                             :: start, finish
real*8, parameter                  :: DX = 0.0125, DY = 0.0125
integer*8, parameter               :: ILXF = 1,ILYF = 1,ILAX = 1, ILAY = 1 
integer*8, parameter               :: ILXC = 1,ILYC = 1 
integer*8                          :: NF, NCOUNT, IUXF, IUYF , IUAX , IUAY, ntemp_P
integer*8                          :: NC, IUXC, IUYC
integer*8                          :: I, J, ITERFINE = 0, ITERCOARSE = 0,  ntemp, counter =0
real*8                             :: REST, RESAVGFINE, RESAVGCOARSE, EPSIAVG = 0.0, EPSIAVGP=1.0
real*8, dimension (:), allocatable :: PHIA, PHINFINE, PHIOFINE, X, Y, ERROR, S, RESFINE, RESCOARSE
real*8                             :: AP, AE, AW, AN, AS, BETA, ALFA, EGV, WU 
real*8                             :: APC, AEC, AWC, ANC, ASC 
real*8, dimension(:), allocatable  :: DPHINCOARSE, DPHIOCOARSE, DPHIFINE 

IUXF = IDNINT(1.0/DX) + 1.0
IUYF = IDNINT(1.0/DY) + 1.0
IUAX = IDNINT(1.0/DX) + 1.0
IUAY = IDNINT(1.0/DX) + 1.0
NF = (IUYF-1.0)*IUXF + IUXF

IUXC = IDNINT(0.5/DX) + 1.0
IUYC = IDNINT(0.5/DY) + 1.0
NC = (IUYC-1.0)*IUXC + IUXC

allocate(S(1:NF))
allocate(PHIA(1:NF))
allocate(PHINFINE(1:NF))
allocate(PHIOFINE(1:NF))
allocate(X(1:NF))
allocate(Y(1:NF))
allocate(ERROR(1:NF))
allocate(RESFINE(1:NF))
allocate(DPHIFINE(1:NF))

allocate(RESCOARSE(1:NC))
allocate(DPHINCOARSE(1:NC))
allocate(DPHIOCOARSE(1:NC))

    RESFINE(:)=0.0
    DPHIFINE(:)=0.0
    DPHINCOARSE(:)=0.0
    DPHIOCOARSE(:)=0.0
    RESCOARSE(:)=0.0
    PHIOFINE(:)=0.0
    PHINFINE(:)=0.0
    PHIA(:)=0.0

open(unit = 2, file = 'residuals.dat')
write(4,*), 'Ieration', ' ', 'AVG Resiual', ' ','AVG Error' 
EPSIAVGP = 1.0


! Array Element Assignment

ITERFINE = 0
RESAVGFINE = 0
ITERCOARSE = 0
call arrayAssignmnt(ILXF, IUXF, ILYF, IUYF, AP, AE, AW, AN, AS, S, PHINFINE, & 
                    PHIOFINE, PHIA, X, Y, DX, DY, NF)

APC = 0.25*AP
AEC = 0.25*AE
AWC = 0.25*AW
ANC = 0.25*AN
ASC = 0.25*AS
! Initial Residual
call RESIDUAL(ILXF, IUXF, ILYF, IUYF, AP, AE, AW, AN, AS, S, &
           PHIOFINE, RESAVGFINE, NF)
           

call cpu_time(start)
!######################################################################!
do while (RESAVGFINE .GE. 1e-8 .AND. counter .LE. 50000)
! solce Mz = r, M = preconditioner Matrix, M= I for present code
    PHIOFINE = PHINFINE
    counter = counter +1
    
! Call point gauss seidel method to solve A*x = b
! pointGS subroutine will change only PHINFINE
    call pointGSPHI(ILXF, IUXF, ILYF, IUYF, ITERFINE, AP, AE, AW, AN, AS, S, &
                   PHIOFINE, PHINFINE, RESFINE, NF, RESAVGFINE)

! Call restriction subroutine to transfer residual to coarse mesh
    call restriction(RESFINE, RESCOARSE, ILXC, ILYC, IUXC, IUYC, IUXF)
    
! Cal Residual function to calculate the L2 norm of residual R(k-1)    
     call RESIDUAL(ILXC, IUXC, ILYC, IUYC, APC, AEC, AWC, ANC, ASC, RESCOARSE, &
           DPHIOCOARSE, RESAVGCOARSE, NC)

! Call point gauss seidel method to solv A*dx = -res
    call pointGSDPHI(ILXC, IUXC, ILYC, IUYC, ITERCOARSE, APC, AEC, AWC, ANC, ASC, RESCOARSE, &
                   DPHIOCOARSE, DPHINCOARSE, NC, RESAVGCOARSE)


! Call prolongation subroutine to transfer solution error to fine mesh
    call prolongation(DPHINCOARSE, DPHIFINE, ILXF, ILYF, IUXF, IUYF, IUXC)
      
! Add solution error to initial guess 
	PHINFINE = PHINFINE + DPHIFINE
	DPHIOCOARSE(:) = 0.0

! Cal Residual function to calculate the L2 norm of residual R(k-1    
    RESAVGFINE = 0.0 
    call RESIDUAL(ILXF, IUXF, ILYF, IUYF, AP, AE, AW, AN, AS, S, &
           PHINFINE, RESAVGFINE, NF)
    EPSIAVG = 0.0
    ERROR = PHINFINE - PHIOFINE
    do I = 1, NF 
        EPSIAVG = EPSIAVG + (ERROR(I))**2
    end do
    EPSIAVG = ((EPSIAVG/(NF-2.0*IUXF-2.0*IUYF+4.0)))**0.5
    EGV = EPSIAVG/EPSIAVGP
    EPSIAVGP = EPSIAVG
    write(2,*) counter ,' ', ITERFINE, ' ', ITERCOARSE, '  ', RESAVGFINE, ' ', EPSIAVG, ' ',EGV    
end do
!########################################################################!        
call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start
WU = ITERFINE + 0.25*ITERCOARSE
! Writing th final solution
open(unit = 4, file = 'solution.dat')
write(4,*), 'Node', ' ', 'Analytical', ' ','Numerical' 
do I = 1, NF
ERROR(I) = ABS(100 * (PHIA(I) - PHINFINE(I))/PHIA(I))
write(4,*), I, ' ', PHIA(I), ' ', PHINFINE(I), ' ', ERROR(I)
end do
!! Writing the Solution at each node location (Contour Plot)
NCOUNT = NF
open(unit = 5, file = 'solutionContour.dat')
write(5,*) 'TITLE="Error Field Data"'
write(5,*) 'variables="X(m)""Y(m)""Analytical""Numerical""Error"'
write(5,*) 'Zone T = "n=1"'
write(5,*) 'I =', IUXF, 'J =', IUYF
write(5,*) 'DATAPACKING = POINT'
do while (NCOUNT .GE. 1)
do I = IUXF, ILXF, -1
write(5,*) X(NCOUNT), ' ', Y(NCOUNT), ' ', PHIA(I), ' ', PHINFINE(I), ' ', ERROR(NCOUNT)
NCOUNT =NCOUNT - 1
end do
end do
print*, WU, ' ', "Finish"
end program multiGrid

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module assignArray
implicit none
contains
subroutine arrayAssignmnt(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, &
                          PHIN, PHIO, PHIA, X, Y, DX, DY, N)

integer*8      					  :: ILX, IUX, ILY, IUY, I, J, ntemp, N
real*8         					  :: DX, DY
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
        S(ntemp) = 50000*(100*(((1-X(ntemp))**2+(Y(ntemp))**2))-2)*EXP(-50*(((1-X(ntemp))**2)+(Y(ntemp))**2))
    end if    
end do
end do
PHIN = PHIO
end subroutine arrayAssignmnt
end module assignArray
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module calculateResidual
implicit none
contains 
subroutine RESIDUAL(ILX, IUX, ILY, IUY, AP, AE, AW, AN, AS, S, &
                 PHIO, RESAVG, N)
                           
integer*8                          :: I, J, N, ILX, IUX, ILY, IUY, ntemp_P
real*8                             :: RESAVG 
real*8, dimension(:), allocatable  :: PHIO, S
real*8                             :: AP, AE, AW, AN, AS, REST

RESAVG = 0.0
do J = ILY+1, IUY-1
do I = ILX+1, IUX-1
    ntemp_P = (J-1)*IUX + I
    REST = (1.0*S(ntemp_P))-(AN*PHIO(ntemp_P+IUX)+AS*PHIO(ntemp_P-IUX) &
							  + AE*PHIO(ntemp_P+1)+AW*PHIO(ntemp_P-1) &
							  - AP*PHIO(ntemp_P))
    RESAVG = RESAVG + (REST)**2
end do
end do
RESAVG = ((RESAVG/(N-2.0*IUX-2.0*IUY+4.0)))**0.5
! print *, RESAVG
end subroutine RESIDUAL
end module calculateResidual
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module pointGaussSeidel
implicit none
contains

subroutine pointGSPHI(ILX, IUX, ILY, IUY, ITER, AP, AE, AW, AN, AS, S, &
                   PHIO, PHIN, RESFINE, N, RESAVGFINE)
integer*8                         :: I, J, ITER, N, ILX, IUX, ILY, IUY, ntemp_P
integer*8                         :: ITERLOOP
real*8                            :: RESNEW, RESPRE, RESRATIO, RES, REST, RESAVGFINE
real*8, dimension(:), allocatable :: PHIO, PHIN, S, RESFINE
real*8                            :: AP, AE, AW, AN, AS

RESNEW = 0.0
RESPRE = RESAVGFINE
RESRATIO = 0.0
ITERLOOP = 0
! Point Gauss Seidel
do while (RESRATIO.LE. 0.5 )!.AND. ITERLOOP .LE. 5)
! print *, "Point GS PHI", ILX, IUX, ILY, IUY
ITER = ITER+1
! ITERLOOP = ITERLOOP+1
RES = 0.0 
do J = ILY+1, IUY-1
do I = ILX+1, IUX-1
    ntemp_P = (J-1)*IUX + I
    PHIN(ntemp_P) = (1.0/AP)*(AE*PHIO(ntemp_P+1) + AW*PHIN(ntemp_P-1) & 
                    + AN*PHIO(ntemp_P+IUX) + AS*PHIN(ntemp_P-IUX) - S(ntemp_P))
end do
end do

do J = ILY+1, IUY-1
do I = ILX+1, IUX-1
    ntemp_P = (J-1)*IUX + I
    REST = (1.0*S(ntemp_p))-(AE*PHIN(ntemp_P+1)+AW*PHIN(ntemp_P-1) & 
                               +AN*PHIN(ntemp_P+IUX)+AS*PHIN(ntemp_P-IUX) &
                               - AP*PHIN(ntemp_P))
    RESFINE(ntemp_P) = REST
    RES = RES + REST**2
end do
end do
	RESNEW = (RES/(N-2.0*IUX-2.0*IUY+4.0))**0.5   
    RESRATIO = RESNEW/RESPRE
    RESPRE = RESNEW
    RESAVGFINE = RESNEW
    PHIO = PHIN
end do
end subroutine pointGSPHI
    
subroutine pointGSDPHI(ILX, IUX, ILY, IUY, ITER, AP, AE, AW, AN, AS, S, &
                   PHIO, PHIN, N, RESAVGCOARSE)
integer*8                         :: I, J, ITER, N, ILX, IUX, ILY, IUY, ntemp_P
integer*8                         :: ITERLOOP
real*8                            :: RESNEW, RESPRE, RESRATIO, RES, REST, RESAVGCOARSE
real*8, dimension(:), allocatable :: PHIO, PHIN, S
real*8                            :: AP, AE, AW, AN, AS

RESNEW = 0.0
RESPRE = RESAVGCOARSE
RESRATIO = 0.0
ITERLOOP = 0
! Point Gauss Seidel
do while (RESRATIO.LE. 0.5 )!.AND. ITERLOOP .LE. 5)
! print *, "Point GS PHI", ILX, IUX, ILY, IUY
ITER = ITER+1
! ITERLOOP = ITERLOOP+1
RES = 0.0 

do J = ILY+1, IUY-1    
do I = ILX+1, IUX-1
    ntemp_P = (J-1)*IUX + I
    PHIN(ntemp_P) = (1.0/AP)*(AE*PHIO(ntemp_P+1) + AW*PHIN(ntemp_P-1) & 
                    + AN*PHIO(ntemp_P+IUX) + AS*PHIN(ntemp_P-IUX) - S(ntemp_P))
end do
end do
RES = 0.0
do J = ILY+1, IUY-1
do I = ILX+1, IUX-1
    ntemp_P = (J-1)*IUX + I
    REST = (1.0*S(ntemp_p))-(AE*PHIN(ntemp_P+1)+AW*PHIN(ntemp_P-1) & 
                               +AN*PHIN(ntemp_P+IUX)+AS*PHIN(ntemp_P-IUX) &
                               - AP*PHIN(ntemp_P))
    RES = RES + REST**2
end do
end do
	RESNEW = (RES/(N-2.0*IUX-2.0*IUY+4.0))**0.5
    RESRATIO = RESNEW/RESPRE
    RESPRE = RESNEW
    RESAVGCOARSE = RESNEW
! print *, ITERLOOP, RESRATIO, RESNEW
PHIO = PHIN
end do 
end subroutine pointGSDPHI

end module pointGaussSeidel
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module restrictionOper
implicit none
contains
subroutine restriction(RESFINE, RESCOARSE, ILXC, ILYC, IUXC, IUYC, IUXF)
real*8, dimension(:), allocatable   :: RESFINE, RESCOARSE
integer*8                           :: ILXC, ILYC, IUXC, IUYC, IUXF
integer*8                           :: I, J, ITEMP, JTEMP, ntemp_PF, ntemp_PC


do J = ILYC+1, IUYC-1
do I = ILXC+1, IUXC-1
    ITEMP = 2*I-1
    JTEMP = 2*J-1
    ntemp_PF = (JTEMP-1)*IUXF+ITEMP
    ntemp_PC = (J-1)*IUXC+I
    RESCOARSE(ntemp_PC) = (RESFINE(ntemp_PF+1)+RESFINE(ntemp_PF-1)+RESFINE(ntemp_PF+IUXF) &
                           +RESFINE(ntemp_PF-IUXF)+4.0*RESFINE(ntemp_PF))/8.0
end do
end do
! print *, RESFINE, RESCOARSE

end subroutine restriction
end module restrictionOper
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module prolongationOper
implicit none
contains
subroutine prolongation(DPHINCOARSE, DPHIFINE, ILXF, ILYF, IUXF, IUYF, IUXC)
real*8, dimension(:), allocatable   :: DPHINCOARSE, DPHIFINE
integer*8                           :: ILXF, ILYF, IUXF, IUYF, IUXC
integer*8                           :: I, J, ITEMP, JTEMP, ntemp_PF, ntemp_PC

do J = ILYF+1, IUYF-1
do I = ILXF+1, IUXF-1

    if (MOD(I,2)==0 .AND. MOD(J,2)==0) then
        ITEMP = (I+2)/2         ! NE point's column location
        JTEMP = (J+2)/2         ! NE point's row location
        ntemp_PC = (JTEMP-1)*IUXC + ITEMP    ! NE point on Coarse Grid level
        ntemp_PF = (J-1)*IUXF + I    ! point under consideration on fine grid level
        ! fine = coarse(NE+NW+SE+SW)
        DPHIFINE(ntemp_PF) = (DPHINCOARSE(ntemp_PC)+DPHINCOARSE(ntemp_PC-1) &
                              +DPHINCOARSE(ntemp_PC-IUXC)+DPHINCOARSE(ntemp_PC-IUXC-1))/4.0
    
    else if (MOD(I,2) .NE. 0 .AND. MOD(J,2) .NE. 0) then
        JTEMP = (J+1)/2
        ITEMP = (I+1)/2
        ntemp_PC = (JTEMP-1)*IUXC + ITEMP    !  point on Coarse Grid level
        ntemp_PF = (J-1)*IUXF + I    ! point under consideration on fine grid level
        ! fine = coarse
        DPHIFINE(ntemp_PF) = DPHINCOARSE(ntemp_PC)
        
    else if (MOD(I,2) == 0 .AND. MOD(J,2) .NE. 0) then
        ITEMP = (I+2)/2                 ! East point's column location
        JTEMP = (J+1)/2                 ! East point's rows location
        ntemp_PC = (JTEMP-1)*IUXC + ITEMP    ! E point on Coarse Grid level
        ntemp_PF = (J-1)*IUXF + I    ! point under consideration on fine grid level
        ! fine = coarse(E+W)
        DPHIFINE(ntemp_PF) = (DPHINCOARSE(ntemp_PC) + DPHINCOARSE(ntemp_PC-1))/2.0
        
    else if (MOD(I,2) .NE. 0 .AND. MOD(J,2) == 0) then
        ITEMP = (I+1)/2                 ! North point's column location
        JTEMP = (J+2)/2                 ! North point's row location
        ntemp_PC = (JTEMP-1)*IUXC + ITEMP    ! N point on Coarse Grid level
        ntemp_PF = (J-1)*IUXF + I    ! point under consideration on fine grid level
        ! fine = coarse(N+S)
        DPHIFINE(ntemp_PF) = (DPHINCOARSE(ntemp_PC) + DPHINCOARSE(ntemp_PC-IUXC))/2.0        
    end if
    
end do
end do

end subroutine prolongation
end module prolongationOper
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

