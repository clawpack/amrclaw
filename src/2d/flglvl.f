c
c -----------------------------------------------------------
c
      subroutine flglvl(nvar,naux,lcheck,nxypts,index,lbase,i1flags,
     .                  npts,t0,isize,jsize)
c        subroutine flglvl(nvar,naux,lcheck,nxypts,index,lbase,npts,t0)
c
      use amr_module
      implicit double precision (a-h,o-z)
      integer*1 i1flags(isize+2,jsize+2)
      integer*1 dom1flags(isize+2,jsize+2)


c
c :::::::::::::::::::: FLGLVL :::::::::::::::::::::::::::::::::
c
c flglvl = controls the error estimation/flagging bad pts. for
c          an entire level of grids.  returns pointer into alloc
c          where the (x,y) coordinations of the flagged pts. are.
c input parameters:
c           lcheck = level to be flagged
c output parameters:
c           nxypts = no. of flagged pts. total
c           index  = starting index in alloc of the flagged pts.
c                    (which occupy 2*nxypts locations).
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
c    
c   prepare domain in ldom2 (so can use ldom as scratch array before 
c   putting in the flags)
c
      idim = iregsz(lbase)
      jdim = jregsz(lbase)
      call domprep(i1flags,lbase,idim,jdim)

      call domshrink(i1flags,dom1flags,idim,jdim)

      do 6 lev = lbase+1, lcheck
         call domup(i1flags,dom1flags,idim,jdim,
     1              intratx(lev-1)*idim,intraty(lev-1)*jdim,lev-1)
         idim = intratx(lev-1)*idim
         jdim = intraty(lev-1)*jdim
         call domshrink(i1flags,dom1flags,idim,jdim)
 6    continue
c     # finish by transferring from iflags to iflags2
      call domcopy(i1flags,dom1flags,isize,jsize)
c
      numbad = 0
c     always call spest to set up stuff (initialize iflags, fill locbig)
c      call spest(nvar,naux,lcheck,dom1flags,isize,jsize,t0)
      call spest2(nvar,naux,lcheck,t0)
      if (tol .gt. 0.) call errest(nvar,naux,lcheck)

c      call bufnst(nvar,naux,numbad,lcheck,dom1flags,isize,jsize)
      call bufnst2(nvar,naux,numbad,lcheck)

      nxypts = nxypts + numbad
c
c  colate flagged pts into flagged points array
c  new version needs to check for proper nesting at this point
c  also needs to sort,  so can remove duplicates.
c
      if (nxypts .gt. 0) then  
c in new version, there are bad cells but nxypts isnt true count any longer
c since there are duplicates, and proper nesting not yet checked
           index = igetsp(2*nxypts)
c          call colate(alloc(index),nxypts,lcheck,
c     1                dom1flags,i1flags,isize,jsize,npts)
          call colate2(alloc(index),nxypts,lcheck,npts,lbase)
      else 
         npts = 0
      endif

      return
      end
