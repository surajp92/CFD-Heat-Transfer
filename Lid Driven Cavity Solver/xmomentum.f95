subroutine xmomentum(UOLD, UNEW, UOLDP, VOLD, VNEW, VOLDP, POLD, DX, DY, DT, RHO, MU, RE, ILX, IUX, ILY, IUY, ITER)

implicit none
integer*8							:: ILX, IUX, ILY, IUY, I, J, ITER
real*8, dimension(IUX,IUY)	        :: UOLD, UNEW, VOLD, VNEW, POLD, UOLDP, VOLDP
real*8								:: DX, DY, RHO, MU, RE, DT
real*8								:: RHOUE, RHOUW, RHOVN, RHOVS
real*8								:: RHOUEP, RHOUWP, RHOVNP, RHOVSP
real*8								:: AE, AW, AN, AS, AP, HN, HNP

do J = 2, IUY-1
do I = 3, IUX-1
! print *, I, " ", J
! using properties at n level
RHOUE = 0.5*(UOLD(I+1,J) + UOLD(I,J))
RHOUW = 0.5*(UOLD(I-1,J) + UOLD(I,J))
RHOVN = 0.5*(VOLD(I,J+1) + VOLD(I-1,J+1))
RHOVS = 0.5*(VOLD(I-1,J) + VOLD(I,J))

AE = (-0.5*RHOUE+1.0/(RE*DX))*DY
AW = (0.5*RHOUW+1.0/(RE*DX))*DY

AN = (-0.5*RHOVN+1.0/(RE*DY))*DX
AS = (0.5*RHOVS+1.0/(RE*DY))*DX

AP = (-0.5*RHOUE-1.0/(RE*DX))*DY + (0.5*RHOUW-1.0/(RE*DX))*DY + (-0.5*RHOVN-1.0/(RE*DY))*DX + (0.5*RHOVS-1.0/(RE*DY))*DX

HN = AE*UOLD(I+1,J) + AW*UOLD(I-1,J) + AN*UOLD(I,J+1) + AS*UOLD(I,J-1) + AP*UOLD(I, J)

! Using Properties at (n-1) level
RHOUEP = 0.5*(UOLDP(I+1,J) + UOLDP(I,J))
RHOUWP = 0.5*(UOLDP(I-1,J) + UOLDP(I,J))
RHOVNP = 0.5*(VOLDP(I,J+1) + VOLDP(I-1,J+1))
RHOVSP = 0.5*(VOLDP(I-1,J) + VOLDP(I,J))

!AE = (-0.5*RHOUEP+1.0/(RE*DX))*DY
!AW = (0.5*RHOUWP+1.0/(RE*DX))*DY

!AN = (-0.5*RHOVNP+1.0/(RE*DY))*DX
!AS = (0.5*RHOVSP+1.0/(RE*DY))*DX

!AP = (-0.5*RHOUEP-1.0/(RE*DX))*DY + (0.5*RHOUWP-1.0/(RE*DX))*DY + (-0.5*RHOVNP-1.0/(RE*DY))*DX + (0.5*RHOVSP-1.0/(RE*DY))*DX

HNP = AE*UOLDP(I+1,J) + AW*UOLDP(I-1,J) + AN*UOLDP(I,J+1) + AS*UOLDP(I,J-1) + AP*UOLDP(I, J)

! Calculate UNEW at (n+1/2) level
if (ITER == 1) then
UNEW(I,J) = UOLD(I,J) + DT*HN/(DX*DY)
else
UNEW(I,J) = UOLD(I,J) + DT*(3.0*HN/2.0 - HNP/2.0)/(DX*DY)
end if

end do
end do
! print *, "X- Momentum Completed"
end subroutine xmomentum
