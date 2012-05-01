c
c ----------------------------------------------------------------------------
c
      subroutine flagcheck(rectflags,ilo,ihi,jlo,jhi,mbuff,iflags,mptr)

      use amr_module
      implicit double precision (a-h, o-z)

      dimension rectflags(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)
      integer*1 iflags(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff) 
c
c should really only check interior border cells here
c
      level = node(nestlevel,mptr)
      imin = max(ilo-mbuff,0)
      jmin = max(jlo-mbuff,0)
      imax = min(ihi+mbuff,iregsz(level)-1)
      jmax = min(jhi+mbuff,jregsz(level))-1

c     do 10 j = jlo-mbuff, jhi+mbuff
c     do 10 i = ilo-mbuff, ihi+mbuff
      do 10 j = jmin, jmax
      do 10 i = imin, imax
        if (rectflags(i,j) .eq. goodpt) go to 10
        if (iflags(i,j) .ne. 1) then  !point not nested. turn off
             rectflags(i,j) = 0.
             write(outunit,100) i,j,mptr
 100         format("turning off point ",2i5," from grid ",i5)
        endif
 10   continue

      return
      end
