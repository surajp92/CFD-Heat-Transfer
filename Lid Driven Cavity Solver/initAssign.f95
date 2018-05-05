subroutine initAssign(UOLD, UNEW, VOLD, VNEW, POLD, PNEW, ILX, &        
& IUX, ILY, IUY, UOLDP, VOLDP, POLDP)

implicit none
integer*8							:: ILX, IUX, ILY, IUY
real*8, dimension(ILY,IUY)	        :: UOLD, UNEW, VOLD, VNEW, POLD, PNEW
real*8, dimension(IUX,IUY)			:: UOLDP, VOLDP, POLDP


UOLD(:,:) = 0.0
UNEW(:,:) = 0.0
VOLD(:,:) = 0.0
VNEW(:,:) = 0.0
POLD(:,:) = 0.0
PNEW(:,:) = 0.0
UOLDP(:,:) = 0.0
VOLDP(:,:) = 0.0
POLDP(:,:) = 0.0

end subroutine initAssign
