subroutine setaux(mbc,mx,my,xlower,ylower,dxc,dyc,maux,aux,aux_copy_mask)
    ! ===========================================
    !  set auxiliary arrays for advection on a curvilinear grid
    !  on input, (xc(i),yc(j)) gives uniformly spaced computational grid.
    !  on output, 
    !    aux(1,i,j) is edge velocity at "left" boundary of grid point (i,j)
    !    aux(2,i,j) is edge velocity at "bottom" boundary of grid point (i,j)
    !    aux(3,i,j) = kappa  is ratio of cell area to dxc*dyc
    !       implicit double precision (a-h,o-z)

    implicit none
    integer, intent(in) :: mbc, mx, my, maux
    real(kind=8), intent(in) :: xlower, ylower, dxc, dyc
    real(kind=8), intent(inout) :: aux(3, 1-mbc:mx+mbc,1-mbc:my+mbc)
    integer(kind=1), intent(in) :: aux_copy_mask(1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j, ic
    real(kind=8) :: xccorn(5), yccorn(5), xpcorn(5), ypcorn(5)
    real(kind=8) :: area

    ! Stream function interface
    interface
        real(kind=8) function stream(x,y)
            implicit none
            real(kind=8), intent(in) :: x,y
        end function stream
    end interface

    do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            ! computational points (xc,yc) are mapped to physical
            ! coordinates (xp,yp) by mapc2p:

            ! lower left corner:
            xccorn(1) = xlower + (i-1)*dxc
            yccorn(1) = ylower + (j-1)*dyc
            call mapc2p(xccorn(1),yccorn(1),xpcorn(1),ypcorn(1))

            ! upper left corner:
            xccorn(2) = xccorn(1)
            yccorn(2) = yccorn(1) + dyc
            call mapc2p(xccorn(2),yccorn(2),xpcorn(2),ypcorn(2))

            ! upper right corner:
            xccorn(3) = xccorn(1) + dxc
            yccorn(3) = yccorn(1) + dyc
            call mapc2p(xccorn(3),yccorn(3),xpcorn(3),ypcorn(3))

            ! lower right corner:
            xccorn(4) = xccorn(1) + dxc
            yccorn(4) = yccorn(1)
            call mapc2p(xccorn(4),yccorn(4),xpcorn(4),ypcorn(4))
            
            ! compute edge velocities by differencing stream function:

	        aux(1,i,j) =  (stream(xpcorn(2),ypcorn(2)) - stream(xpcorn(1),ypcorn(1)))/ dyc

	        aux(2,i,j) = -(stream(xpcorn(4),ypcorn(4)) - stream(xpcorn(1),ypcorn(1)))/ dxc

            ! compute area of physical cell from four corners:
            xpcorn(5) = xpcorn(1)
            ypcorn(5) = ypcorn(1)
            area = 0.d0
	        do ic=1,4
	            area = area + 0.5d0 * (ypcorn(ic)+ypcorn(ic+1)) * (xpcorn(ic+1)-xpcorn(ic))
	        enddo
	    
	       aux(3,i,j) = area / (dxc*dyc)
       end do
   end do
end subroutine setaux
