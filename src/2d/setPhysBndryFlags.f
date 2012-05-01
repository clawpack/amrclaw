c
c -----------------------------------------------------------
c
       subroutine setPhysBndryFlags(iflags,ilo,ihi,jlo,jhi,mbuff,level)

       use amr_module
       implicit double precision (a-h, o-z)

       integer*1  iflags(ilo-mbuff:ihi+mbuff, jlo-mbuff:jhi+mbuff)

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
c                 (NOT YET DONE)
c ****************************************************************

       if (ilo-mbuff .lt. 0 .and. .not. xperdom) then  ! grid extends out of left side of domain
c       set left flagged points to be ok
          do j = jlo-mbuff, jhi+mbuff
            do i = ilo-mbuff, -1
             iflags(i,j) = iflags(0,j)   ! extend using whatever 1st col inside is 
            end do
          end do
       endif

       if (ihi+mbuff .ge. iregsz(level) .and. .not. xperdom) then
c       set right flagged points to be ok
          do j = jlo-mbuff, jhi+mbuff
            do i = iregsz(level), ihi+mbuff
             iflags(i,j) = iflags(iregsz(level)-1,j)
            end do
          end do
       endif

 
       if (jlo-mbuff .lt. 0 .and. .not. yperdom) then
c       set bottom flagged points to be ok
         do j = jlo-mbuff, -1    
           do i = ilo-mbuff, ihi+mbuff
             iflags(i,j) = iflags(i,0)
            end do
          end do
       endif

       if (jhi+mbuff .ge. jregsz(level) .and. .not. yperdom) then
c       set top flagged points to be ok
         do j = jregsz(level), jhi+mbuff
           do i = ilo-mbuff, ihi+mbuff
             iflags(i,j) = iflags(i,jregsz(level)-1) ! extend using  last flag in domain 
            end do
          end do
       endif

       if (xperdom .or. yperdom) then
          write(*,*) "NOT SET YET FOR PERIODICITY IN SETPHYSBNDRYFLAGS"
       endif

       return
       end
