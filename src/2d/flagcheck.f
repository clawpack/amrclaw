!> Check if every cell in grid mptr is properly nested in base level
!! grids (base level in current refinement, usually represented by
!! **lbase**).
!! If not, turn off the flagging on that cell (mark as needing
!! no refinement).
c
c ----------------------------------------------------------------------------
c
      subroutine flagcheck(rectflags,ilo,ihi,jlo,jhi,mbuff,iflags,
     .                     imin,imax,jmin,jmax,mptr)

      use amr_module
      implicit double precision (a-h, o-z)

      dimension rectflags(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)
      integer*1 iflags(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff) 
c
c should really only check interior border cells here
c
      level = node(nestlevel,mptr)
c
c want to check all cells for flagging, including buffer
c but wrap if periodic when put on list
c if not periodic no need to check if outside domain
c

      do 10 j = jmin, jmax
      do 10 i = imin, imax
        if (rectflags(i,j) .le. DONTFLAG) go to 10
        if (iflags(i,j) .ne. 1) then  !point not nested. turn off
             rectflags(i,j) = DONTFLAG
             if (nprint) then
                write(outunit,100) i,j,mptr
 100            format("turning off point ",2i5," from grid ",i5)
             endif
        endif
 10   continue

      return
      end
