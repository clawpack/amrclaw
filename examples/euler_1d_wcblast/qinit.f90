subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.

    implicit none
    
    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    integer :: i
    real(kind=8) :: xcell, gamma1
 
    real(kind=8) :: gamma
    common /cparam/  gamma
 
    gamma1 = gamma - 1.d0

    do i=1,mx
        xcell = xlower + (i-0.5d0)*dx
        q(1,i) = 1.d0
        q(2,i) = 0.d0
        if (xcell .lt. 0.1d0) then
            q(3,i) = 1.d3/gamma1
          else if (xcell .lt. 0.9d0) then
            q(3,i) = 1.d-2/gamma1
          else
            q(3,i) = 1.d2/gamma1
          endif
    enddo

end subroutine qinit

