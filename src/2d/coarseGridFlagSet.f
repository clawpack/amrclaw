c
!> Flag a whole subregion from (ixlo,ixhi) to (jxlo, jxhi) with integer
!! 1. 
!! The subregion is inside a grid described by (ilo_coarse, jlo_coarse)
!! and (ihi_coarse, jhi_coarse)
c -------------------------------------------------------------
c
       subroutine coarseGridFlagSet(iflags,ixlo,ixhi,jxlo,jxhi,
     .                              ilo_coarse,ihi_coarse,
     .                              jlo_coarse,jhi_coarse,mbuff)

                                     
       integer*1 iflags(ilo_coarse-mbuff:ihi_coarse+mbuff,
     .                  jlo_coarse-mbuff:jhi_coarse+mbuff)

c
c whole point of this routine is to index using the integer cords, not the
c  way its dimensioned and indexed in the calling routine setdomflags
c
       do 25 j = jxlo,jxhi         
       do 25 i = ixlo,ixhi         
          iflags(i,j) = 1
 25    continue

       return
       end
