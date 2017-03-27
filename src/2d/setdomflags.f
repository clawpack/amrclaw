c
!> set domain flags (not AMR flags) for grid **mptr** (only), enlarged by buffer zone. 
!!
!! \param mptr number (index) of the grid being processed
!! \param igridflags pointer to memory address storing domain flags (not
!! AMR flags)
!! \param ilo global *i* index of left-most cell of this grid
!! \param ihi global *i* index of right-most cell of this grid 
!! \param jlo global *j* index of lower-most cell of this grid 
!! \param jhi global *j* index of upper-most cell of this grid 
!! \param mbuff width of buffer zone
!! \param lbase base AMR level for current refinement, which stays
!!  fixed. Note that **lbase** is always less or equal to **lcheck**
!! \param lcheck AMR level of grid **mptr**
!! \param mibuff number of cells (including buffer zone) on grid
!! **mptr** in *i* direction 
!! \param mjbuff number of cells (including buffer zone) on grid
!! **mptr** in *j* direction 

c -----------------------------------------------------------------------------------
c
      subroutine setdomflags(mptr,igridflags,ilo,ihi,jlo,jhi,
     .                       mbuff,lbase,lcheck,mibuff,mjbuff)

      use amr_module

      integer*1 igridflags(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)
c icopy is dimensioned large enough, but will be used at several sizes
c and accessed using lo_i-mbuff:hi_i+mbuff, etc.      
c
      integer*1 icopy(mibuff,mjbuff) 
      dimension ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)
      dimension igridst(lcheck), igridend(lcheck)
      dimension jgridst(lcheck), jgridend(lcheck)
      logical borderx, bordery

c
c set domain flags for this grid only, enlarged by buffer zone. check if any other base grids
c are in exterior or first interior border cell and mark ok.
c note that interior of base grids 1 away from edge are automatically ok for proper nesting
c  will shrink gridflags after setting to get proper nesting region
c
c  1. initialize this grids domain flags to 0, at lcheck
        igridflags = 0
c
c    ... if lbase coarse than lcheck, set initial indices, before upscaling, for base transfer
c        so that dont have entire base grid upscaled
         igridst(lcheck)  = ilo
         igridend(lcheck) = ihi
         jgridst(lcheck)  = jlo
         jgridend(lcheck) = jhi
         do lc = lcheck-1,lbase,-1  !NB: may be a 0 trip do loop, not old fortran
            ilo_coarse = floor(dfloat(igridst(lc+1))/intratx(lc))
            jlo_coarse = floor(dfloat(jgridst(lc+1))/intraty(lc))
            ihi_coarse = ceiling(dfloat(igridend(lc+1))/intratx(lc)) - 1
            jhi_coarse = ceiling(dfloat(jgridend(lc+1))/intraty(lc)) - 1
            if (ihi_coarse*intratx(lc) .lt. igridend(lc+1)) 
     .           ihi_coarse = ihi_coarse+1
            if (jhi_coarse*intraty(lc) .lt. jgridend(lc+1)) 
     .           jhi_coarse = jhi_coarse+1
            igridend(lc) = ihi_coarse
            jgridend(lc) = jhi_coarse
            igridst(lc) = ilo_coarse
            jgridst(lc) = jlo_coarse
         end do
!       get out coarsened indices in case level lbase == lcheck (zero trip loop)
         ilo_coarse = igridst(lbase)
         ihi_coarse = igridend(lbase)
         jlo_coarse = jgridst(lbase)
         jhi_coarse = jgridend(lbase)
    ! If we project grid mptr to level "lbase" to get a grid on that
    ! level, ilo_coarse is global index of its left border and similarly
    ! for other three

c
c  3.  loop over all intersecting grids at base level staying fixed
c      to make the proper nesting dodmain.
c      set the buffer zone in igridflags to 1 if nested 
c      this is so when shrink by one you dont lose too much area.
c
        mbase = lstart(lbase)
 20     continue        
        
           iblo = node(ndilo,mbase)   ! if base grid coarser, need to scale up
           ibhi = node(ndihi,mbase)   ! if same grid will just mark interior cells as 1
           jblo = node(ndjlo,mbase)
           jbhi = node(ndjhi,mbase)
c
c  3.5 if periodic bcs, then if grids buffer sticks out, will have to wrap the
c      coordinates and flag any intersecting base grids for wrapped  buffer.
c      do here instead of above since cant coarsen mbuff same way you can for regular grid
c      also grid itself (without enlarged mbuff zone) doesnt stick out

        borderx = (ilo_coarse.le. 0 .or. ihi_coarse.ge.iregsz(lbase)-1)
        bordery = (jlo_coarse.le. 0 .or. jhi_coarse.ge.jregsz(lbase)-1)
        if ((xperdom .and. borderx) .or. (yperdom .and. bordery)) then   
           call setIndices(ist,iend,jst,jend,
     .                     ilo_coarse-mbuff,ihi_coarse+mbuff,
     .                     jlo_coarse-mbuff,jhi_coarse+mbuff,
     .                     ishift,jshift,lbase)
 
           do 25 i = 1, 3
               i1 = max(ilo_coarse-mbuff,ist(i))
               i2 = min(ihi_coarse+mbuff,iend(i))
           do 24 j = 1, 3
               j1 = max(jlo_coarse-mbuff,jst(j))
               j2 = min(jhi_coarse+mbuff, jend(j))

               if (.not. ((i1 .le. i2) .and. (j1 .le. j2))) go to 24 ! part of patch in this region
c
c              part of patch is in this region [i,j]
c              periodically wrap and fill if it intersects with grid mbase
c              note: this is done in two steps in hopes of greater clarity
               

c usual check would be ->  if ((i1 .gt. i2) .or. (j1 .gt. j2)) go to 24  ! no patch
c cant do that since have not yet included buffer zone - which is the part that would get wrapped

c             patch exist. does it intersect with mbase grid?
c             use wrapped coords of this grid to test if intersects with base grid
              ixlo = max(iblo,i1+ishift(i)) 
              ixhi = min(ibhi,i2+ishift(i))
              jxlo = max(jblo,j1+jshift(j))
              jxhi = min(jbhi,j2+jshift(j))
c
              if ((ixlo .gt. ixhi) .or. (jxlo .gt. jxhi)) go to 24  !this grid doesnt intersect
c
c             if wrapped region does intersect, be careful to set the INTERSECTED part of
c             the UNWRAPPED region  of original enlarged grid
              ixlo_unwrapped = ixlo - ishift(i)   
              ixhi_unwrapped = ixhi - ishift(i)
              jxlo_unwrapped = jxlo - jshift(j)
              jxhi_unwrapped = jxhi - jshift(j)
              call coarseGridFlagSet(igridflags,
     .                               ixlo_unwrapped,ixhi_unwrapped,
     .                               jxlo_unwrapped,jxhi_unwrapped,
     .                               ilo_coarse,ihi_coarse,
     .                               jlo_coarse,jhi_coarse,mbuff)

 24        continue
 25        continue

       else       
           ixlo = max(iblo,ilo_coarse-mbuff)
           ixhi = min(ibhi,ihi_coarse+mbuff)
           jxlo = max(jblo,jlo_coarse-mbuff)
           jxhi = min(jbhi,jhi_coarse+mbuff)
c
c         does this patch intersect mbase grid?
           if (.not.((ixlo .le. ixhi) .and. (jxlo .le. jxhi))) go to 30  !this grid doesnt intersect
c
c          use subroutine call since dimension of igridflags not same as above declaration
c          when on coarser grids
           call coarseGridFlagSet(igridflags,ixlo,ixhi,jxlo,jxhi,
     .                           ilo_coarse,ihi_coarse,
     .                           jlo_coarse,jhi_coarse,mbuff)
       endif

 30     mbase = node(levelptr,mbase)
        if (mbase .ne. 0) go to 20
c
c 3.5 set any part of grid buffer zone to 1 that is at physical boundary 
        call setPhysBndryFlags(igridflags,ilo_coarse,ihi_coarse,
     .                         jlo_coarse,jhi_coarse,mbuff,lbase)

c  4.  done setting flags on base level. next step is to transfer the
c      properly nested domain flags to lcheck - i.e. upscale to level needed 
c      first shrink by 1 for actual nested region.
c      always shrink once - so works if lcheck same as lbase
c      if going up 1 level each one needs to be nested, so still shrink first before upsizing
c
c      after loop above, dom flags in igridflags, copy to icopy (in subr for dimensioning reasons)
      call griddomcopy(icopy,igridflags,ilo_coarse,ihi_coarse,
     .                 jlo_coarse,jhi_coarse,mbuff)
c
c     shrink from icopy to dom2 flag array. This is where shrinking occurs if
c     lbase = lcheck, for proper nesting
      call griddomshrink(icopy,ilo_coarse,ihi_coarse,jlo_coarse,
     .                     jhi_coarse,mbuff,
     .                     alloc(node(domflags2,mptr)),lbase)

      do 40 lev = lbase+1, lcheck
c         ### for each level that upsize, calculate new coords starting from
c         ### actual fine grid and recoarsening down to needed level
c         ### cant take previous coarse coords and refine, since may be
c         ### too large. grid prob. not anchored at base grid corner.
         ilo_fine = igridst(lev)
         ihi_fine = igridend(lev)
         jlo_fine = jgridst(lev)
         jhi_fine = jgridend(lev)
c
c       flags in dom2, upsize to icopy array with finer dimensions
        call griddomup(alloc(node(domflags2,mptr)),icopy,
     .                 ilo_coarse,ihi_coarse,jlo_coarse,jhi_coarse,
     .                 mbuff,lev-1,
     .                 ilo_fine,ihi_fine,jlo_fine,jhi_fine)
c       flags in icopy, shrink one  back to dom2
        call griddomshrink(icopy,ilo_fine,ihi_fine,jlo_fine,jhi_fine,
     .                     mbuff,alloc(node(domflags2,mptr)),lev)
        ilo_coarse = ilo_fine
        ihi_coarse = ihi_fine
        jlo_coarse = jlo_fine
        jhi_coarse = jhi_fine
40    continue
c 
        return
        end

