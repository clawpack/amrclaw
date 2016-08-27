c
c -------------------------------------------------------------
c
       subroutine coarseGridFlagSet(iflags,ixlo,ixhi,
     .                              ilo_coarse,ihi_coarse,
     .                              mbuff)

                                     
       integer*1 iflags(ilo_coarse-mbuff:ihi_coarse+mbuff)

c
c whole point of this routine is to index using the integer cords, not the
c  way its dimensioned and indexed in the calling routine setdomflags
c
       do 25 i = ixlo,ixhi         
          iflags(i) = 1
 25    continue

       return
       end
