c
c -----------------------------------------------------------------------------------
c
      subroutine setdomflags(mptr,igridflags,ilo,ihi,jlo,jhi,
     .                       mbuff,lbase,lcheck,mibuff,mjbuff)

      use amr_module

      integer*1 igridflags(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)
c icopy is dimensioned large enough, but will be used at several sizes
c and accessed using loi-mbuff:hi_i+mbuff, etc.      
      integer*1 icopy(mibuff,mjbuff) 
c
c set domain flags for this grid only, enlarged by buffer zone. check if any other base grids
c are in exterior or first interior border cell and mark ok.
c note that interior of base grids 1 away from edge are automatically ok for proper nesting
c  will shrink gridflags after setting to get proper nesting region
c
c  1. initialize this grids domain flags to 0, at lcheck
        do j = jlo-mbuff, jhi+mbuff
        do i = ilo-mbuff, ihi+mbuff
           igridflags(i,j) = 0
        end do
        end do
c
c    ... and figure ratios from coarser levels
      ratiox = 1.d0   ! needs to be real for proper ceiling/floor behavior below
      ratioy = 1.d0   ! so that coarser grid completely enclosed projected fine grid
      do lev = lbase+1,lcheck   ! difficulty is that if grid more than 1 level away
         ratiox = ratiox * intratx(lev-1)   ! may not hook up nicely at corners
         ratioy = ratioy * intraty(lev-1)
      end do
c
c    ... if lbase coarse than lcheck, set initial indices, before upscaling, for base transfer
c        so that dont have entire base grid upscaled
      ilo_coarse = floor(ilo/ratiox )
      ihi_coarse = ceiling((ihi+1)/ratiox) - 1
      jlo_coarse = floor(jlo/ratioy)
      jhi_coarse = ceiling((jhi+1)/ratioy) - 1
      

c  3.  loop over all other intersecting grids at base level staying fixed
c      set the buffer zone in igridflags to 1 if nested 
c      this is so when shrink by one you dont lose too much area.
        mbase = lstart(lbase)
 20     continue        
           iblo = node(ndilo,mbase)   ! if base grid coarser, need to scale up
           ibhi = node(ndihi,mbase)   ! if same grid will just mark interior cells as 1
           jblo = node(ndjlo,mbase)
           jbhi = node(ndjhi,mbase)
           ixlo = max(iblo,ilo_coarse-mbuff)
           ixhi = min(ibhi,ihi_coarse+mbuff)
           jxlo = max(jblo,jlo_coarse-mbuff)
           jxhi = min(jbhi,jhi_coarse+mbuff)
           if (.not.((ixlo .le. ixhi) .and. (jxlo .le. jxhi))) go to 30  !this grid doesnt intersect
c
           call coarseGridFlagSet(igridflags,ixlo,ixhi,jxlo,jxhi,
     .                           ilo_coarse,ihi_coarse,
     .                           jlo_coarse,jhi_coarse,mbuff)
 30     mbase = node(levelptr,mbase)
        if (mbase .ne. 0) go to 20
c
c 3.5 set any part of grid buffer zone to 1 that is at physical boundary 
        call setPhysBndryFlags(igridflags,ilo_coarse,ihi_coarse,
     .                         jlo_coarse,jhi_coarse,mbuff,lbase)

c  4.  done setting flags. upscale to level needed then
c      shrink by 1 for actual nested region
c      shrink once - works if lcheck same as lbase
c      if going up 1 level each one needs to be nested, so still shrink first before upsizing
c
c  after loop above, dom flags in igridflags, copy to icopy
      call griddomcopy(icopy,igridflags,ilo_coarse,ihi_coarse,
     .                 jlo_coarse,jhi_coarse,mbuff)
c
c     shrink from icopy to upsized flag array
      call griddomshrink(icopy,ilo_coarse,ihi_coarse,jlo_coarse,
     .                     jhi_coarse,mbuff,
     .                     alloc(node(domflags_upsized,mptr)),lbase)

      do 40 lev = lbase+1, lcheck
c          cant scale up in a simple way, since might  end up with too large grid
c          if all grids not exactly anchored at base grid corner.  instead
c          need to recalculate from the given (fine)  grid coords at each level
           rx = 1.d0
           ry = 1.d0
              do lc = lcheck,lev+1,-1  !NB: ,ay be a 0 trip do loop, not old fortran
                rx = rx * intratx(lc-1)
                ry = ry * intraty(lc-1)
              end do
           ilofine = floor(ilo/rx)
           jlofine = floor(jlo/ry)
           ihifine = ceiling((ihi+1)/rx) - 1
           jhifine = ceiling((jhi+1)/ry) - 1

c       flags in upsized, upsize to icopy array with finer dimensions
        call griddomup(alloc(node(domflags_upsized,mptr)),icopy,
     .                 ilo_coarse,ihi_coarse,jlo_coarse,jhi_coarse,
     .                 mbuff,lev-1,
     .                 ilofine,ihifine,jlofine,jhifine)
c       flags in icopy, shrink one  back to upsized
        call griddomshrink(icopy,ilofine,ihifine,jlofine,jhifine,
     .                     mbuff,alloc(node(domflags_upsized,mptr)),lev)
        ilo_coarse = ilofine
        ihi_coarse = ihifine
        jlo_coarse = jlofine
        jhi_coarse = jhifine
40    continue
c 
        return
        end

