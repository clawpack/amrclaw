subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    ! Set initial conditions for q.
    ! Sample scalar equation with data that is piecewise constant with
    ! q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
    !     0.1  otherwise

    ! set concentration profile
    q(1, :, :) = 0.1d0

end subroutine qinit