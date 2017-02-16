c
c -----------------------------------------------------------------------------------
c
      subroutine setdomflags(mptr,igridflags,ilo,ihi,
     .                       mbuff,lbase,lcheck,mibuff)

      use amr_module

      integer*1 igridflags(ilo-mbuff:ihi+mbuff)
c icopy is dimensioned large enough, but will be used at several sizes
c and accessed using loi-mbuff:hi_i+mbuff, etc.      
      integer*1 icopy(mibuff)
      dimension ist(3), iend(3), ishift(3)
      dimension igridst(lcheck)
      dimension igridend(lcheck)

c
c set domain flags for this grid only, enlarged by buffer zone. check if any other base grids
c are in exterior or first interior border cell and mark ok.
c note that interior of base grids 1 away from edge are automatically ok for proper nesting
c  will shrink gridflags after setting to get proper nesting region
c
c  1. initialize this grids domain flags to 0, at lcheck
        do i = ilo-mbuff, ihi+mbuff
           igridflags(i) = 0
        end do
 
c
c
c    ... if lbase coarse than lcheck, set initial indices, before upscaling, for base transfer
c        so that dont have entire base grid upscaled
         igridst(lcheck)  = ilo
         igridend(lcheck) = ihi
         do lc = lcheck-1,lbase,-1  !NB: may be a 0 trip do loop, not old fortran
            ilo_coarse = floor(dfloat(igridst(lc+1))/intratx(lc))
            ihi_coarse = ceiling(dfloat(igridend(lc+1))/intratx(lc)) - 1
            if (ihi_coarse*intratx(lc) .lt. igridend(lc+1)) 
     .           ihi_coarse = ihi_coarse+1
            igridend(lc) = ihi_coarse
            igridst(lc) = ilo_coarse
         end do
         ilo_coarse = igridst(lbase)
         ihi_coarse = igridend(lbase)

      if (xperdom) then   ! set here once and for all, use coarsened "base" coords to set
         call setIndices(ist,iend,
     .              ilo_coarse-mbuff,ihi_coarse+mbuff,
     .              ishift,lbase)
       endif 

c  3.  loop over all other intersecting grids at base level staying fixed
c      set the buffer zone in igridflags to 1 if nested 
c      this is so when shrink by one you dont lose too much area.
        mbase = lstart(lbase)
 20     continue        
           iblo = node(ndilo,mbase)   ! if base grid coarser, need to scale up
           ibhi = node(ndihi,mbase)   ! if same grid will just mark interior cells as 1
c
c  3.5 if periodic bcs, then if grids buffer sticks out, will have to wrap the
c      coordinates and flag any intersecting base grids for wrapped  buffer.
c      do here instead of above since cant coarsen mbuff same way you can for regular grid
c      also grid itself (without enlarged mbuff zone) doesnt stick out
       if (xperdom) then
           do 25 i = 1, 3
c              i1 = max(ist(i)  + ishift(i)
c              i2 = iend(i) + ishift(i)
               i1 = max(ilo_coarse-mbuff,ist(i))
               i2 = min(ihi_coarse+mbuff,iend(i))

               if (.not. (i1 .le. i2)) go to 25 ! part of patch in this region
c
c              part of patch is in this region [i]
c              periodically wrap and fill if it intersects with grid mbase
c              note: this is done in two steps in hopes of greater clarity
               

c usual check would be ->  if (i1 .gt. i2) go to 25  ! no patch
c cant do that since have not yet included buffer zone - which is the part that would get wrapped

c             patch exist. does it intersect with mbase grid?
c             use wrapped coords of this grid to test if intersects with base grid
              ixlo = max(iblo,i1+ishift(i)) 
              ixhi = min(ibhi,i2+ishift(i))
c
              if (ixlo .gt. ixhi) go to 25  !this grid doesnt intersect
c
c             if wrapped region does intersect, be careful to set the INTERSECTED part of
c             the UNWRAPPED region  of original enlarged grid
              ixlo_unwrapped = ixlo - ishift(i)   
              ixhi_unwrapped = ixhi - ishift(i)
              call coarseGridFlagSet(igridflags,
     .                               ixlo_unwrapped,ixhi_unwrapped,
     .                               ilo_coarse,ihi_coarse,mbuff)

 25        continue

       else       
           ixlo = max(iblo,ilo_coarse-mbuff)
           ixhi = min(ibhi,ihi_coarse+mbuff)
c
c         does this patch intersect mbase grid?
           if (.not.(ixlo .le. ixhi)) go to 30  !this grid doesnt intersect
c
           call coarseGridFlagSet(igridflags,ixlo,ixhi,
     .                           ilo_coarse,ihi_coarse,
     .                           mbuff)
       endif

 30     mbase = node(levelptr,mbase)
        if (mbase .ne. 0) go to 20
c
c 3.5 set any part of grid buffer zone to 1 that is at physical boundary 
        call setPhysBndryFlags(igridflags,ilo_coarse,ihi_coarse,
     .                         mbuff,lbase)

c  4.  done setting flags. upscale to level needed then
c      shrink by 1 for actual nested region
c      shrink once - works if lcheck same as lbase
c      if going up 1 level each one needs to be nested, so still shrink first before upsizing
c
c  after loop above, dom flags in igridflags, copy to icopy
      call griddomcopy(icopy,igridflags,ilo_coarse,ihi_coarse,
     .                 mbuff)
c
c     shrink from icopy to dom2 flag array
      call griddomshrink(icopy,ilo_coarse,ihi_coarse,
     .                     mbuff,
     .                     alloc(node(domflags2,mptr)),lbase)

      do 40 lev = lbase+1, lcheck
c         ### for each level that upsize, calculate new coords starting from
c         ### actual fine grid and recoarsening down to needed level
c         ### cant take previous coarse coords and refine, since may be
c         ### too large. grid prob. not anchored at base grid corner.
         ilofine = igridst(lev)
         ihifine = igridend(lev)
c
c       flags in dom2, upsize to icopy array with finer dimensions
        call griddomup(alloc(node(domflags2,mptr)),icopy,
     .                 ilo_coarse,ihi_coarse,
     .                 mbuff,lev-1,
     .                 ilofine,ihifine)
c       flags in icopy, shrink one  back to dom2
        call griddomshrink(icopy,ilofine,ihifine,
     .                     mbuff,alloc(node(domflags2,mptr)),lev)
        ilo_coarse = ilofine
        ihi_coarse = ihifine
40    continue
c 
        return
        end

