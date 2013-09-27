c
c --------------------------------------------------------------
c
      subroutine preicall(val,aux,nrow,ncol,nfil,nvar,naux,
     1                    ilo,ihi,jlo,jhi,klo,khi,level)
c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension val(nvar,nrow,ncol,nfil)
      dimension aux(naux,nrow,ncol,nfil)

      dimension ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)
      dimension kst(3), kend(3),                  kshift(3)

c
c  :::::::::::::: PREICALL :::::::::::::::::::::::::::::::::::::::::::
c     For periodic boundary conditions more work needed to initialize a
c     new grid that sticks out. This routine was
c     called because the patch sticks out of the domain,
c     and has periodic bc.s preprocess the patch before calling
c     icall to shift the patch periodically back into the domain.
c
c     Inputs to this routine:
c     ilo,ihi,jlo,jhi,klo,khi = the location in index space of this patch.
c
c     Outputs from this routine:
c     The values of the grid are inserted
c     directly into the enlarged val array for this piece.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c
c     # will divide patch into 9 possibilities (some empty): 
c       x sticks out left, x interior, x sticks out right
c       same for y. for example, the max. would be
c       i from (ilo,-1), (0,iregsz(level)-1), (iregsz(level),ihi)
        
        ist(1)    =  ilo
        ist(2)    =  0
        ist(3)    =  iregsz(level)
        iend(1)   = -1
        iend(2)   =  iregsz(level)-1
        iend(3)   =  ihi
        jst(1)    =  jlo
        jst(2)    =  0
        jst(3)    =  jregsz(level)
        jend(1)   = -1
        jend(2)   =  jregsz(level)-1
        jend(3)   =  jhi
        kst(1)    =  klo
        kst(2)    =  0
        kst(3)    =  kregsz(level)
        kend(1)   = -1
        kend(2)   =  kregsz(level)-1
        kend(3)   =  khi
        ishift(1) =  iregsz(level)
        ishift(2) =  0
        ishift(3) = -iregsz(level)
        jshift(1) =  jregsz(level)
        jshift(2) =  0
        jshift(3) = -jregsz(level)
        kshift(1) =  kregsz(level)
        kshift(2) =  0
        kshift(3) = -kregsz(level)

       do 30 i = 1, 3
          i1 = max(ilo,  ist(i))
          i2 = min(ihi, iend(i))
       do 20 j = 1, 3
          j1 = max(jlo,  jst(j))
          j2 = min(jhi, jend(j))
       do 10 k = 1, 3
          k1 = max(klo,  kst(k))
          k2 = min(khi, kend(k))


          if (((i1 .le. i2) .and. (j1 .le. j2)) .and. (k1 .le. k2)) then
            iputst = i1 - ilo + 1
            jputst = j1 - jlo + 1
            kputst = k1 - klo + 1
            call icall(val,aux,nrow,ncol,nfil,nvar,naux,
     1                    i1+ishift(i),i2+ishift(i),
     2                    j1+jshift(j),j2+jshift(j),
     3                    k1+kshift(k),k2+kshift(k),level,
     4                    iputst,jputst,kputst)
          endif

 10    continue
 20    continue
 30    continue
      
     


      return
      end
