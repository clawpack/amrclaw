c
c --------------------------------------------------------------
c
      subroutine preintcopy(val,mitot,nvar,ilo,ihi,
     1                      level,fliparray)
c
      use amr_module
      implicit double precision (a-h,o-z)

      dimension fliparray((mitot)*nghost*nvar)
      dimension val(nvar,mitot)
      dimension ist(3), iend(3), ishift(3)

c  NEW INDEXING ORDER SWITCHED
      iadd(ivar,i)  = locflip + ivar-1 + nvar*(i-1)
c
c  :::::::::::::: PREINTCOPY :::::::::::::::::::::::::::::::::::::::::::
c     For periodic boundary conditions more work needed to initialize a
c     new grid that sticks out. This routine was
c     called because the patch sticks out of the domain,
c     and has periodic bc.s preprocess the patch before calling
c     intcopy to shift the patch periodically back into the domain.
c
c     Inputs to this routine:
c     ilo,ihi = the location in index space of this patch.
c
c     Outputs from this routine:
c     The values of the grid are inserted
c     directly into the enlarged val array for this piece.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
        locflip = 1        

c
c     # will divide patch into 3 possibilities (some empty):
c       x sticks out left, x interior, x sticks out right
c       for example, the max. would be
c       i from (ilo,-1), (0,iregsz(level)-1), (iregsz(level),ihi)
        
        if (xperdom) then
           ist(1) = ilo
           ist(2) = 0
           ist(3) = iregsz(level)
           iend(1) = -1
           iend(2) = iregsz(level)-1
           iend(3) = ihi
           ishift(1) = iregsz(level)
           ishift(2) = 0
           ishift(3) = -iregsz(level)
        else    ! if not periodic, set vals to only have nonnull intersection for interior region
           ist(1)    = iregsz(level)
           ist(2)    = ilo
           ist(3)    = iregsz(level)
           iend(1)   = -1
           iend(2)   = ihi
           iend(3)   = -1
           ishift(1) = 0
           ishift(2) = 0
           ishift(3) = 0
        endif


       do 20 i = 1, 3
          i1 = max(ilo,  ist(i))
          i2 = min(ihi, iend(i))
          if (i1 .le. i2) then ! part of patch in this region
c
            iputst = (i1 - ilo) + 1
            call intcopy(val,mitot,nvar,
     2                    i1+ishift(i),i2+ishift(i),
     3                    level,iputst)

          endif

 20    continue
      
     
      return
      end
