subroutine setaux(mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,maux,aux,aux_copy_mask)
   ! set auxiliary arrays
   !
   ! advection
   !    aux(i,j,k,1) is u velocity on left face of cell
   !    aux(i,j,k,2) is v velocity on bottom face of cell
   !    aux(i,j,k,3) is w velocity on back face of cell

   implicit none

   integer, intent(in) :: mx, my, mz, mbc, maux
   real(kind=8), intent(in) :: xlower, ylower, zlower, dx, dy, dz
   real(kind=8), intent(in out) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
   integer(kind=1), intent(in) :: aux_copy_mask(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
   
   integer ::  i,j,k
   real(kind=8) :: pi, pi2, dx2, dy2, dz2
   real(kind=8) :: xc(1-mbc:mx+mbc)
   real(kind=8) :: yc(1-mbc:my+mbc)
   real(kind=8) :: zc(1-mbc:mz+mbc)

   common /piconstants/ pi, pi2

   pi = 4.d0*atan(1.d0)
   pi2 = 2.d0*pi

   do i = 1-mbc,mx+mbc
      xc(i) = xlower + (i-0.5d0)*dx
   enddo

   do j = 1-mbc,my+mbc
      yc(j) = ylower + (j-0.5d0)*dy
   enddo

   do k = 1-mbc,mz+mbc
      zc(k) = zlower + (k-0.5d0)*dz
   enddo

   dx2 = 0.5d0*dx
   dy2 = 0.5d0*dy
   dz2 = 0.5d0*dz

   do  k = 1-mbc,mz+mbc
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            aux(1,i,j,k) = compute_u(xc(i)-dx2, yc(j),zc(k))
            aux(2,i,j,k) = compute_v(xc(i),yc(j)-dy2, zc(k))
            aux(3,i,j,k) = compute_w(xc(i),yc(j), zc(k)-dz2)
         enddo
      enddo
   enddo

contains

   real(kind=8) function compute_u(x,y,z)
      implicit none
      real(kind=8), intent(in) :: x,y,z
      real(kind=8) :: pi, pi2
      common /piconstants/ pi, pi2
      compute_u = 2.d0 * dsin(pi*x)**2 * dsin(pi2*y) * dsin(pi2*z)
   end function compute_u

   real(kind=8) function compute_v(x,y,z)
      implicit none
      real(kind=8), intent(in) :: x,y,z
      real(kind=8) :: pi, pi2
      common /piconstants/ pi, pi2
      compute_v = -dsin(pi2*x) * dsin(pi*y)**2 * dsin(pi2*z)
   end function compute_v

   real(kind=8) function compute_w(x,y,z)
      implicit none
      real(kind=8), intent(in) :: x,y,z
      real(kind=8) :: pi, pi2
      common /piconstants/ pi, pi2
      compute_w = -dsin(pi2*x) * dsin(pi2*y) * dsin(pi*z)**2
   end function compute_w

end subroutine setaux

