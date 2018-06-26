c
c -----------------------------------------------------------
c
       subroutine setPhysBndry(rectflags,ilo,ihi,jlo,jhi,mbuff,level)

       use amr_module
       implicit double precision (a-h, o-z)

       dimension rectflags(ilo-mbuff:ihi+mbuff, jlo-mbuff:jhi+mbuff)

c ****************************************************************
!> If grid borders the physical domain then
!! turn off any flagged points in buffer zone = those points
!! are not properly nested (and it doesnt matter).
!! But last row/col interior to grid if flagged is ok
!! 
!! if periodic, then have to look elsewhere to see if
!! last interior row/col that is flagged is ok.
!! (done in rest of colate2)
c ****************************************************************

       if (ilo .eq. 0 .and. .not. xperdom) then
c       set left flagged points to be ok
          do j = jlo-mbuff, jhi+mbuff
            do i = ilo-mbuff, ilo-1
             rectflags(i,j) = DONTFLAG
            end do
c           1st interior cell ok if on bndry. set back to pos if flagged
            rectflags(0,j) = abs(rectflags(0,j))
          end do
       endif

       if (ihi .eq. iregsz(level)-1 .and. .not. xperdom) then
c       set right flagged points to be ok
          do j = jlo-mbuff, jhi+mbuff
            do i = ihi+1, ihi+mbuff
             rectflags(i,j) = DONTFLAG
            end do
            rectflags(ihi,j) = abs(rectflags(ihi,j))
          end do
       endif

 
       if (jlo .eq. 0 .and. .not. yperdom) then
c       set bottom flagged points to be ok
          do i = ilo-mbuff, ihi+mbuff
            do j = jlo-mbuff, jlo-1    
             rectflags(i,j) = DONTFLAG
            end do
            rectflags(i,0) = abs(rectflags(i,0))
          end do
       endif

       if (jhi .eq. jregsz(level)-1 .and. .not. yperdom) then
c       set top flagged points to be ok
          do i = ilo-mbuff, ihi+mbuff
            do j = jhi+1, jhi+mbuff
             rectflags(i,j) = DONTFLAG
            end do
            rectflags(i,jhi) = abs(rectflags(i,jhi))
          end do
       endif

       return
       end
