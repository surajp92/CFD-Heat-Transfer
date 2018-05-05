subroutine pressureequation(UOLD, UNEW, VOLD, VNEW, POLD, PNEW, DX, DY, DT, RHO, ILX, IUX, ILY, IUY)

implicit none
integer*8							:: ILX, IUX, ILY, IUY, I, J
real*8, dimension(IUX,IUY)	        :: UOLD, UNEW, VOLD, VNEW, POLD, PNEW, S, POLDTEMP
real*8								:: DX, DY,RHO, DT
real*8, dimension(IUX, IUY)			:: AEM, AWM, ANM, ASM, APM 
real*8								:: RESAVGP, RESP

S(:,:) = 0.0
AEM(:,:) = 0.0
AWM(:,:) = 0.0
ANM(:,:) = 0.0
ASM(:,:) = 0.0
APM(:,:) = 0.0


RESAVGP = 1.0
RESP = 0.0

do J = 2, IUY-1
do I = 2, IUX-1

S(I,J) = ((RHO*UNEW(I+1,J)-RHO*UNEW(I,J))*DY + (RHO*VNEW(I,J+1)-RHO*VNEW(I,J))*DX)/DT
if (I == 2) then
	AEM(I,J) = DY/DX
	AWM(I,J) = 0.0
	ANM(I,J) = DX/DY
	ASM(I,J) = DX/DY
	APM(I,J) = (-2.0*DX/DY-1.0*DY/DX)
	if (J == 2) then
		ASM(I,J) = 0.0
		APM(I,J) = (-1.0*DX/DY-1.0*DY/DX)
	else if (J == IUY-1) then
		ANM(I,J) = 0.0
		APM(I,J) = (-1.0*DX/DY-1.0*DY/DX)
	end if
	
else if (I == IUX-1) then
	AEM(I,J) = 0.0
	AWM(I,J) = DY/DX
	ANM(I,J) = DX/Dy
	ASM(I,J) = DX/DY
	APM(I,J) = (-2.0*DX/DY-1.0*DY/DX)
	if (J == 2) then
		ASM(I,J) = 0.0
		APM(I,J) = (-1.0*DX/DY-1.0*DY/DX)
	else if (J == IUY-1) then
		ANM(I,J) = 0.0
		APM(I,J) = (-1.0*DX/DY-1.0*DY/DX)
	end if
	
else if (J == 2) then
	AEM(I,J) = DY/DX
	AWM(I,J) = DY/DX
	ANM(I,J) = DX/DY
	ASM(I,J) = 0.0
	APM(I,J) = (-1.0*DX/DY-2.0*DY/DX)
	if (I == 2) then
		AWM(I,J) = 0.0
		APM(I,J) = (-1.0*DX/DY-1.0*DY/DX)
	else if (I == IUX-1) then
		AEM(I,J) = 0.0
		APM(I,J) = (-1.0*DX/DY-1.0*DY/DX)
	end if
	
else if (J == IUY-1) then
	AEM(I,J) = DY/DX
	AWM(I,J) = DY/DX
	ANM(I,J) = 0.0
	ASM(I,J) = DX/DY
	APM(I,J) = (-1.0*DX/DY-2.0*DY/DX)
	if (I == 2) then
		AWM(I,J) = 0.0
		APM(I,J) = (-1.0*DX/DY-1.0*DY/DX)
	else if (I == IUX-1) then
		AEM(I,J) = 0.0
		APM(I,J) = (-1.0*DX/DY-1.0*DY/DX)
	end if	

else 	
	AEM(I,J) = DY/DX
	AWM(I,J) = DY/DX
	ANM(I,J) = DX/DY
	ASM(I,J) = DX/DY
	APM(I,J) = (-2.0*DX/DY-2.0*DY/DX)
end if

end do
end do

POLDTEMP = POLD
PNEW = POLD

! using Gauss Seidel method for linear solver
do while(RESAVGP .GE. 1E-5)

do J = 2, IUY-1
do I = 2, IUX-1

PNEW(I,J) = (-(AEM(I,J)*PNEW(I+1,J)+AWM(I,J)*PNEW(I-1,J)+ASM(I,J)*PNEW(I,J-1)+ANM(I,J)*PNEW(I,J+1))+S(I,J))/APM(I,J)

end do
end do

RESAVGP = 0.0
do J = 2, IUY-1
do I = 2, IUX-1

RESP = S(I,J)-(AEM(I,J)*PNEW(I+1,J)+AWM(I,J)*PNEW(I-1,J)+ASM(I,J)*PNEW(I,J-1)+ANM(I,J)*PNEW(I,J+1)+APM(I,J)*PNEW(I,J))
RESAVGP = RESAVGP+RESP**2
end do
end do
RESAVGP = (RESAVGP/((IUX-2)*(IUY-2)))**0.5
POLDTEMP = PNEW
!print *, RESAVGP
end do

do J = 2, IUY-1
do I = 2, IUX-1

! normalize with P(2,2) to limit the truncation error
PNEW(I,J) = PNEW(I,J)-POLDTEMP(2,2)

end do
end do

! update pressure at Boundaries
PNEW(1,:) =  PNEW(2,:)
PNEW(IUX,:) =  PNEW(IUX-1,:)
PNEW(:,1) =  PNEW(:,2)
PNEW(:,IUY) =  PNEW(:,IUY-1)


! linear solver to be included
! print *, "Pressure Correction Completed"
end subroutine pressureequation
