subroutine assignBoundary(UOLD, UOLDP, UNEW, VOLD, VNEW, POLD, PNEW, ILX, IUX, ILY, IUY, ULID)

implicit none
integer*8									:: ILX, IUX, ILY, IUY
real*8									    :: ULID
real*8, dimension(IUX,IUY)			     	:: UOLD, UNEW, VOLD, VNEW, POLD, PNEW, UOLDP


UOLD(ILX+1,:)   		= 0.0
UOLD(IUX,:)     		= 0.0
UOLD(:,ILY)   		    = 0.0
UOLD(ILX:IUX,IUY)       = ULID
! assign BC to Unew also
UNEW = UOLD
UOLDP = UOLD

! print *, UOLD
end subroutine assignBoundary
