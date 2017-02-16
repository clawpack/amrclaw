c
c -----------------------------------------------------------
c
       subroutine setPhysBndryFlags(iflags,ilo,ihi,mbuff,level)

       use amr_module
       implicit double precision (a-h, o-z)

       integer*1  iflags(ilo-mbuff:ihi+mbuff)

c ****************************************************************
c  setPhysBndryFlags = if grid borders the physical domain then
c                 set domain flags to 1 in buffer zone. That way when shrink
c                 by 1 to get proper nested domain, you wont lose the first
c                 border cell of a grid
c
c                 flag array uses 0-based indexing
c              
c                 if periodic, then have to look elsewhere to see if
c                 last interior row/col that is flagged is ok.
c                 this is done in the calling routine that
c                 transfers flagged points to base grdis
c ****************************************************************

       if (ilo-mbuff .lt. 0 .and. .not. xperdom) then  ! grid extends out of left side of domain
c       set left flagged points to be ok
            do i = ilo-mbuff, -1
c             iflags(i) = iflags(0)   ! extend using whatever 1st col inside is
             iflags(i) = 1     ! set to 1 for bndry buffers that touch exterior domain
            end do
       endif

       if (ihi+mbuff .ge. iregsz(level) .and. .not. xperdom) then
c       set right flagged points to be ok
            do i = iregsz(level), ihi+mbuff
c             iflags(i) = iflags(iregsz(level)-1)
             iflags(i) = 1
            end do
       endif

       return
       end
