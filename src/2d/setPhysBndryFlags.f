c ****************************************************************
!> If grid borders the physical domain then
!! set domain flags to 1 in buffer zone (which is outside the physical
!! domain). 
!! That way when shrink
!! by 1 to get proper nested domain, you wont lose the first
!! border cell of a grid
!!
!! flag array uses 0-based index space
!! 
!! if periodic, then have to look elsewhere to see if
!! last interior row/col that is flagged is ok.
!! this is done in the calling routine that
!! transfers flagged points to base grids
!!
!! \param iflags array that stores the flags for this grid
!! \param ilo global *i* index of left border of the grid
!! \param ihi global *i* index of right border of the grid
!! \param jlo global *j* index of lower border of the grid
!! \param jhi global *j* index of upper border of the grid
!! \param mbuff width of buffer zone
!! \param level AMR level of the grid being flagged
c ****************************************************************

c
c -----------------------------------------------------------
c
       subroutine setPhysBndryFlags(iflags,ilo,ihi,jlo,jhi,mbuff,level)

       use amr_module
       implicit double precision (a-h, o-z)

       integer*1  iflags(ilo-mbuff:ihi+mbuff, jlo-mbuff:jhi+mbuff)

       if (ilo-mbuff .lt. 0 .and. .not. xperdom) then  ! grid extends out of left side of domain
c       set left flagged points to be ok
          do j = jlo-mbuff, jhi+mbuff
            do i = ilo-mbuff, -1
c             iflags(i,j) = iflags(0,j)   ! extend using whatever 1st col inside is 
             iflags(i,j) = 1     ! set to 1 for bndry buffers that touch exterior domain
            end do
          end do
       endif

       if (ihi+mbuff .ge. iregsz(level) .and. .not. xperdom) then
c       set right flagged points to be ok
          do j = jlo-mbuff, jhi+mbuff
            do i = iregsz(level), ihi+mbuff
c             iflags(i,j) = iflags(iregsz(level)-1,j)
             iflags(i,j) = 1
            end do
          end do
       endif

 
       if (jlo-mbuff .lt. 0 .and. .not. yperdom) then
c       set bottom flagged points to be ok
         do j = jlo-mbuff, -1    
           do i = ilo-mbuff, ihi+mbuff
c             iflags(i,j) = iflags(i,0)
             iflags(i,j) = 1
            end do
          end do
       endif

       if (jhi+mbuff .ge. jregsz(level) .and. .not. yperdom) then
c       set top flagged points to be ok
         do j = jregsz(level), jhi+mbuff
           do i = ilo-mbuff, ihi+mbuff
c             iflags(i,j) = iflags(i,jregsz(level)-1) ! extend using  last flag in domain 
             iflags(i,j) = 1
            end do
          end do
       endif


       return
       end
