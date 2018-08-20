module abl_module

  implicit none
  integer, parameter :: NO_ABL = 0
  integer, parameter :: TRIG = 1
  integer, parameter :: AC = 2 ! Appelo & Colonius
  integer, parameter :: PS = 3 ! Petersson & Sjogreen

  integer, parameter :: MAX_PARAM = 3

  real(kind=8), parameter :: PIH = 2.d0*datan(1.d0)


  integer :: abl_type = 0
  integer :: p_AC, q_AC
  real(kind=8) :: eps_AC, eps_PS
  real(kind=8) :: xdepth(2), ydepth(2), xpos(2), ypos(2)
  logical :: initialized = .false.

contains

  subroutine initialize()

    implicit none

    real(kind=8) :: xlower, xupper, ylower, yupper

    ! Read in ABL data
    ! open the unit with new routine from Clawpack 4.4 to skip over
    ! comment lines starting with #:
    call opendatafile(7, 'abl.data')

    ! Read in parameters:
    read(7,*) abl_type
    if (abl_type == TRIG) then
      print *, "ABL is trigonometric"
    else if (abl_type == AC) then
      read(7,*) eps_AC
      read(7,*) p_AC
      read(7,*) q_AC
      print *, "ABL is Appelo/Colonius with (eps,p,q)=", eps_AC, p_AC, q_AC
    else if (abl_type == PS) then
      read(7,*) eps_PS
      print *, "ABL is Petersson/Sjogreen with eps=", eps_PS
    end if
    read(7,*) xdepth(1), ydepth(1)
    read(7,*) xdepth(2), ydepth(2)
    close(7)

    ! Read in grid parameters
    ! open the unit with new routine from Clawpack 4.4 to skip over
    ! comment lines starting with #:
    call opendatafile(7, 'claw.data')

    read(7,*) xlower ! num_dim
    read(7,*) xlower, ylower ! lower
    read(7,*) xupper, yupper ! upper
    close(7)

    ! Compute ABL positions
    xpos(1) = xlower + xdepth(1)
    xpos(2) = xupper - xdepth(2)
    ypos(1) = ylower + ydepth(1)
    ypos(2) = yupper - ydepth(2)

    initialized = .true.

  end subroutine initialize

  subroutine set_scaling_factor(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)

    implicit none

    integer, intent(in) :: mbc, mx, my, maux
    real(kind=8), intent(in) :: xlower, ylower, dx, dy
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer :: i, j
    real(kind=8) :: xcenter, ycenter, lower(2)

    ! initialize module if needed
    if (.not. initialized) then
      call initialize()
    end if

    if (abl_type == 0) then
      return
    end if

    ! Loop over all cell edges
    do j=1-mbc,my+mbc

      ycenter = ylower + (j-0.5d0)*dy

      do i=1-mbc,mx+mbc

        xcenter = xlower + (i-0.5d0)*dx

        ! determine lower cell edge locations
        ! constrained to the computational domain without ghost cells
        lower(1) = max(xcenter-0.5d0*dx,xpos(1)-xdepth(1))
        lower(1) = min(lower(1),xpos(2)+xdepth(2)-dx)
        lower(2) = max(ycenter-0.5d0*dy,ypos(1)-ydepth(1))
        lower(2) = min(lower(2),ypos(2)+ydepth(2)-dy)

        ! Calculate inverse of Jacobian 1/(1 + g'(x_k-1/2)) & 1/(1 + g'(y_k-1/2))
        aux(maux-1,i,j) = 1.d0/(1.d0 + gp(lower(1),1))
        aux(maux,i,j) = 1.d0/(1.d0 + gp(lower(2),2))

      end do
    end do

  end subroutine set_scaling_factor

  function gp(z,direction)

    implicit none

    real(kind=8) :: gp, z
    integer :: direction

    real(kind=8) :: pos(2), depth(2), sig, val

    gp = 0.d0
    sig = 0.d0
    if (direction == 1) then
      pos = xpos
      depth = xdepth
    else if (direction == 2) then
      pos = ypos
      depth = ydepth
    end if

    if (z < pos(1)) then
      val = (pos(1)-z)/depth(1)
    else if (z > pos(2)) then
      val = (z-pos(2))/depth(2)
    else
      return
    end if

    if (abl_type == TRIG) then
      gp = dtan(PIH*val)**2
    else if (abl_type == AC) then
      sig = (1.0 - eps_AC)*(1.0 - (1.0 - val)**p_AC)**q_AC
      gp = sig/(1.0-sig)
    else if (abl_type == PS) then
      sig = (1.0 - eps_PS)* &
      val**6*(462.0 - 1980.0*val + 3465.0*val**2 - 3080.0*val**3 + 1386.0*val**4 - 252.0*val**5)
      gp = sig/(1.0-sig)
    end if

  end function gp

end module abl_module
