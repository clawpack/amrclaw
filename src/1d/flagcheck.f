c
c ----------------------------------------------------------------------------
c
      subroutine flagcheck(rectflags,ilo,ihi,mbuff,iflags,
     .                     imin,imax,mptr)

      use amr_module
      implicit double precision (a-h, o-z)

      dimension rectflags(ilo-mbuff:ihi+mbuff)
      integer*1 iflags(ilo-mbuff:ihi+mbuff)
c
c should really only check interior border cells here
c
      level = node(nestlevel,mptr)
c
c want to check all cells for flagging, including buffer
c but wrap if periodic when put on list
c if not periodic no need to check if outside domain
c

      do 10 i = imin, imax
        if (rectflags(i) .le. DONTFLAG) go to 10
        if (iflags(i) .ne. 1) then  !point not nested. turn off
             rectflags(i) = 0.
             if (nprint) then
                write(outunit,100) i,mptr
 100            format("turning off point ",2i5," from grid ",i5)
             endif
        endif
 10   continue

      return
      end
