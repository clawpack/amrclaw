! =========================================================
! Check for NANs in solution q
subroutine check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,ichecknan)

    implicit none
    
    ! Input
    integer, intent(in) :: maxmx, maxmy, meqn, mbc, mx, my, ichecknan
    real(kind=8), intent(in) :: t, q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)

    ! Locals
    integer :: i, j, m

    ! print *,'Checking for NANs at ichecknan = ',ichecknan
    ! print *,'  maxmx = ',maxmx,'  maxmy = ',maxmy,'  meqn = ',meqn
      
    do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
            do m=1,meqn
                if (.not. (q(m,i,j) == q(m,i,j))) then
                    ! true if q(i,j,m) = NAN
	                print *, 'SOLUTION ERROR --- ABORTING CALCULATION'
	                print *, 'At ichecknan = ',ichecknan
                    print *, '   mx,my,t:',mx,my,t
                    print *, '   m,i,j:',m,i,j
                    print *, '   q(m,i,j) = ',q(m,i,j)
	                stop
                endif
            enddo
        enddo
    enddo

    ! Uncomment the next line if desired when debugging:
    ! print *,'No NANs at ichecknan = ',ichecknan,' at t = ',t

end subroutine check4nans
