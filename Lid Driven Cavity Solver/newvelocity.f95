subroutine newvelocity(UOLD, UNEW, VOLD, VNEW, POLD, PNEW, DX, DY, RHO, MU, ILX, IUX, ILY, IUY, DT)

implicit none
integer*8							:: ILX, IUX, ILY, IUY, I, J
real*8, dimension(IUX,IUY)	        :: UOLD, UNEW, VOLD, VNEW, POLD, PNEW, S
real*8								:: DX, DY, RHO, MU, DT
real*8								:: AE, AW, AN, AS, AP


do J = 2, IUY-1
do I = 3, IUX-1

UNEW(I,J) = UNEW(I,J) - DT*(PNEW(I,J)-PNEW(I-1,J))/(DX)

end do
end do

do J = 3, IUY-1
do I = 2, IUX-1

VNEW(I,J) = VNEW(I,J) - DT*(PNEW(I,J)-PNEW(I,J-1))/(DY)

end do
end do

! add pressure correction to old pressure to get new pressure
!do J = 2, IUY-1
!do I = 2, IUX-1

!PNEW(I,J) = POLD(I,J) + PNEW(I,J)

!end do
!end do
! print *, "Velocity Updated"
end subroutine newvelocity
