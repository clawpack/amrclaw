! =========================================================
!> Check for NANs in solution q
subroutine check4nans(meqn,mbc,mx,my,q,t,grid_local_id,level)

    implicit none
    
    ! Input
    integer, intent(in) ::  meqn, mbc, mx, my, grid_local_id,level
    real(CLAW_REAL), intent(in) :: t, q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j, m

    ! print *,'Checking for NANs at grid_local_id = ',grid_local_id
    ! print *,'  mx = ',mx,'  my = ',my,'  meqn = ',meqn
      
    do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
            do m=1,meqn
                if (.not. (q(m,i,j) == q(m,i,j))) then
                    ! true if q(i,j,m) = NAN
                    print *, 'SOLUTION ERROR --- ABORTING CALCULATION'
                    print *, 'At grid_local_id = ',grid_local_id
                    print *, 'level = ',level
                    print *, '   mx,my,t:',mx,my,t
                    print *, '   m,i,j:',m,i,j
                    print *, '   q(m,i,j) = ',q(m,i,j)
                    stop
                endif
            enddo
        enddo
    enddo

    ! Uncomment the next line if desired when debugging:
    ! print *,'No NANs at grid_local_id = ',grid_local_id,' at t = ',t

end subroutine check4nans
