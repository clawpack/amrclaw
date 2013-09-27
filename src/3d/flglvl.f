c
c -----------------------------------------------------------
c
      subroutine flglvl(nvar,naux,lcheck,nxypts,index,lbase,
     1                  i1flags,npts,t0,isize,jsize,ksize)
c
      use amr_module
      implicit double precision (a-h,o-z)

      integer*1 i1flags(isize+2,jsize+2,ksize+2)
      integer*1 dom1flags(isize+2,jsize+2,ksize+2)
c
c :::::::::::::::::::: FLGLVL :::::::::::::::::::::::::::::::::
c
c flglvl = controls the error estimation/flagging bad pts. for
c          an entire level of grids.  returns pointer into alloc
c          where the (x,y,z) coordinations of the flagged pts. are.
c input parameters:
c           lcheck = level to be flagged
c output parameters:
c           nxypts = no. of flagged pts. total
c           index  = starting index in alloc of the flagged pts.
c                    (which occupy numdim*nxypts locations).
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
      nxypts = 0
c
c   reserve space for entire domain worth of flagged points at
c   level lcheck. bits would be better, but integer will do
c   dom2 - holds domain flags
c   dom  - holds flagged pts.
c   dom3 - scratch
c
      isize = iregsz(lcheck)
      jsize = jregsz(lcheck)
      ksize = kregsz(lcheck)
c    ldom  = igetsp((isize+2)*(jsize+2)*(ksize+2))

c    
c   prepare domain in ldom2 (so can use ldom as scratch array before 
c   putting in the flags)
c
      idim = iregsz(lbase)
      jdim = jregsz(lbase)
      kdim = kregsz(lbase)
      call domprep(i1flags,lbase,idim,jdim,kdim)
      call domshrink(i1flags,dom1flags,idim,jdim,kdim)

      do 6 lev = lbase+1, lcheck
         call domup(i1flags,dom1flags,idim,jdim,kdim,
     1              intratx(lev-1)*idim,intraty(lev-1)*jdim,
     2              intratz(lev-1)*kdim,lev-1)
        idim = intratx(lev-1)*idim
        jdim = intraty(lev-1)*jdim
        kdim = intratz(lev-1)*kdim
         call domshrink(i1flags,dom1flags,idim,jdim,kdim)
 6    continue
c     # finish by transferring from iflags to iflags2
      call domcopy(i1flags,dom1flags,isize,jsize,ksize)
c
      numbad = 0
c     always call spest to set up stuff (initialize iflags, fill locbig)
      call spest(nvar,naux,lcheck,dom1flags,isize,jsize,ksize,
     1           tstart_thisrun)
      if (flag_richardson) call errest(nvar,naux,lcheck)

      call bufnst(nvar,naux,numbad,lcheck,dom1flags,isize,jsize,ksize)
      nxypts = nxypts + numbad
c
c  colate flagged pts into flagged points array
c
      if (nxypts .gt. 0) then
          index = igetsp(3*nxypts)
          call colate(alloc(index),nxypts,lcheck,
     1                dom1flags,i1flags,isize,jsize,ksize,npts)
      else 
          npts = nxypts
      endif
c
c      call reclam(ldom,  (isize+2)*(jsize+2)*(ksize+2))

      return
      end
