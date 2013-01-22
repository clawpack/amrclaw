!  :::::::::::::: PREFILRECUR :::::::::::::::::::::::::::::::::::::::::::
!     For periodic boundary conditions more work needed to fill the
!     piece of the boundary. This routine was
!     called because the patch sticks out of the domain,
!     and has periodic bc.s preprocess the patch before calling
!     filpatch to shift the patch periodically back into the domain.
!
!     Inputs to this routine:
!     xl, xr, yb, yt = the location in physical space of
!     corners of a patch.
!     fill_indices = the location in index space of this patch.
!
!     Outputs from this routine:
!     The values around the border of the grid are inserted
!     directly into the enlarged valbig array for this piece.
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
recursive subroutine prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,nrowst,ncolst,fill_indices)

    use amr_module, only: iregsz, jregsz, nghost, xlower, ylower
    use amr_module, only: spheredom, hxposs, hyposs
    implicit none

    ! Input
    integer, intent(in) :: level, nvar, naux, mitot, mjtot, nrowst, ncolst, fill_indices(4)
    real(kind=8), intent(in) :: time

    ! Output
    real(kind=8), intent(in out) :: valbig(nvar,mitot,mjtot)
    real(kind=8), intent(in out) :: aux(naux,mitot,mjtot)
    
    ! Local storage
    integer :: i, j, ii, jj, ivar, nr, nc, ng, i1, i2, j1, j2, iputst, jputst
    integer :: jbump, iwrap1, iwrap2, jwrap1, tmp, locflip, rect(4)
    real(kind=8) :: xlwrap, ybwrap

    integer :: ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)
    real(kind=8) :: scratch(max(mitot,mjtot)*nghost*nvar)
    real(kind=8) :: scratchaux(max(mitot,mjtot)*nghost*naux)

!     # will divide patch into 9 possibilities (some empty): 
!       x sticks out left, x interior, x sticks out right
!       same for y. for example, the max. would be
!       i from (fill_indices(1),-1), (0,iregsz(level)-1), (iregsz(level),fill_indices(2))
        
!--        write(dbugunit,*)" entering prefilrecur with level ",level
!--        write(dbugunit,*)"     and patch indices fill_indices ",
!--     &             fill_indices

    ist(1)  = fill_indices(1)
    ist(2)  = 0
    ist(3)  = iregsz(level)
    iend(1) = -1
    iend(2) = iregsz(level)-1
    iend(3) = fill_indices(2)

    jst(1)  = fill_indices(3)
    jst(2)  = 0
    jst(3)  = jregsz(level)
    jend(1) = -1
    jend(2) = jregsz(level)-1
    jend(3) = fill_indices(4)

    ishift(1) = iregsz(level)
    ishift(2) = 0
    ishift(3) = -iregsz(level)
    jshift(1) = jregsz(level)
    jshift(2) = 0
    jshift(3) = -jregsz(level)

    do i = 1, 3
        i1 = max(fill_indices(1),  ist(i))
        i2 = min(fill_indices(2), iend(i))
        do j = 1, 3
            j1 = max(fill_indices(3),  jst(j))
            j2 = min(fill_indices(4), jend(j))

            ! part of patch in this region
            if ((i1 <= i2) .and. (j1 <= j2)) then 

                ! there is something to fill. j=2 case is interior, no special
                ! mapping needed even if spherical bc
                if (.not. spheredom .or. j == 2 ) then
                    iputst = (i1 - fill_indices(1)) + nrowst
                    jputst = (j1 - fill_indices(3)) + ncolst
                    rect = [i1+ishift(i),i2+ishift(i),j1+jshift(j),j2+jshift(j)]
                    call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot, &
                                  iputst,jputst,rect)
                else
                    nr = i2 - i1 + 1
                    nc = j2 - j1 + 1
                    ng = 0    ! no ghost cells in this little patch. fill everything.

                    jbump = 0
                    if (j1 < 0)   jbump = abs(j1)  ! bump up so new bottom is 0
                    if (j2 >= jregsz(level)) jbump = -(j2+1-jregsz(level)) ! bump down

                    ! next 2 lines would take care of periodicity in x
                    iwrap1 = i1 + ishift(i)
                    iwrap2 = i2 + ishift(i)
                    ! next 2 lines take care of mapped sphere bcs
                    iwrap1 = iregsz(level) - iwrap1 -1
                    iwrap2 = iregsz(level) - iwrap2 -1
                    ! swap so that smaller one is left index, etc since mapping reflects
                    tmp = iwrap1
                    iwrap1 = iwrap2
                    iwrap2 = tmp

                    jwrap1 = j1 + jbump
                    xlwrap = xlower + iwrap1*hxposs(level)
                    ybwrap = ylower + jwrap1*hyposs(level)
              
                    if (naux>0) then
                        call setaux(nr,nc,ng,nr,nc,xlwrap,ybwrap,hxposs(level),hyposs(level),naux,scratchaux)
                    endif
                    rect = [iwrap1,iwrap2,j1+jbump,j2+jbump]
                    call filrecur(level,nvar,scratch,scratchaux,naux,time,nr, &
                                  nc,1,1,rect)

                    ! copy back using weird mapping for spherical folding
                    do ii = i1, i2
                        do jj = j1, j2
                            ! write(dbugunit,'(" filling loc ",2i5," with ",2i5)') nrowst+ii-fill_indices(1),ncolst+jj-fill_indices(3),nr-(ii-i1),nc-jj+j1

                            do ivar = 1, nvar
                                valbig(ivar,nrowst+(ii-fill_indices(1)),ncolst+(jj-fill_indices(3))) = &
                                    scratch(iaddscratch(ivar,nr-(ii-i1),nc-(jj-j1)))
                            end do
                            ! write(dbugunit,'(" new val is ",4e15.7)')(valbig(ivar,nrowst+(ii-fill_indices(1)),ncolst+(jj-fill_indices(3))),ivar=1,nvar)
                        end do
                    end do 
                endif
            endif
        end do
    end do

contains

    integer pure function iadd(n,i,j)
        implicit none
        integer, intent(in) :: n, i, j
        iadd = locflip + n-1 + nvar*((j-1)*nr+i-1)
    end function iadd

    integer pure function iaddscratch(n,i,j)
        implicit none
        integer, intent(in) :: n, i, j
        iaddscratch = n + nvar*((j-1)*nr+i-1)  ! no subtract 1
    end function iaddscratch

end subroutine prefilrecur