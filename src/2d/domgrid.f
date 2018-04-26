c
c ---------------------------------------------------------------------------------
c
      subroutine domgrid(lbase,lcheck)
c
      use amr_module
      implicit double precision (a-h,o-z)

      mbuff = max(nghost,ibuff+1)
c
!>  loop over base grids to get proper nesting domain for grids at level lcheck
!!  but only upsize to the lcheck grids dimensions
c

      mptr = lstart(lcheck)
 10   continue
      ilo    = node(ndilo,mptr)
      ihi    = node(ndihi,mptr)
      jlo    = node(ndjlo,mptr)
      jhi    = node(ndjhi,mptr)
      nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
      ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
c     up to mbuff cells on each side might be flagged and
c     buffered, so to allow shrinkage, need yet one additional
c     cell on each side to be set from base grids
      mibuff = nx + 2*mbuff
      mjbuff = ny + 2*mbuff
      ibytesPerDP = 8

c   bad names, for historical reasons. they are both smae size now
      locdomflags = igetsp( (mibuff*mjbuff)/ibytesPerDP+1)
      locdom2 = igetsp( (mibuff*mjbuff)/ibytesPerDP+1)


      node(domflags_base,mptr) = locdomflags
      node(domflags2,mptr) = locdom2
      call setdomflags(mptr,alloc(locdomflags),ilo,ihi,jlo,jhi,
     .                 mbuff,lbase,lcheck,mibuff,mjbuff)

      mptr = node(levelptr, mptr)
      if (mptr .ne. 0) go to 10

      return
      end
