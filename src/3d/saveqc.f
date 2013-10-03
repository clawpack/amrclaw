c
c  ================================================================
      subroutine saveqc(level,nvar,naux)
c  ================================================================
c
      use amr_module
      implicit double precision (a-h,o-z)

      logical sticksout
c
c ::::::::::::::::::::::::: SAVEQC :::::::::::::::::::::::::::::::::
c prepare new fine grids to save fluxes after each integration step
c for future conservative fix-up on coarse grids.
c save all boundary fluxes of fine grid (even if on a  phys. bndry.) -
c but only save space for every intrat of them.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      levc = level - 1
      hxc  = hxposs(levc)
      hyc  = hyposs(levc)
      hzc  = hzposs(levc)

      mkid = lstart(level)
 10   if (mkid .eq. 0) go to 99
       nx    = node(ndihi,mkid)-node(ndilo,mkid) + 1
       ny    = node(ndjhi,mkid)-node(ndjlo,mkid) + 1
       nz    = node(ndkhi,mkid)-node(ndklo,mkid) + 1
       ikeep = nx/intratx(level-1)
       jkeep = ny/intraty(level-1)
       kkeep = nz/intratz(level-1)
       lenbc = 2*(ikeep*jkeep+jkeep*kkeep+kkeep*ikeep)
       ist   = node(ffluxptr,mkid)
       time  = rnode(timemult,mkid)

c         make coarsened enlarged patch for conservative fixup
       ilo = node(ndilo,mkid)
       jlo = node(ndjlo,mkid)
       klo = node(ndklo,mkid)
       ihi = node(ndihi,mkid)
       jhi = node(ndjhi,mkid)
       khi = node(ndkhi,mkid)

       iclo = ilo/intratx(level-1) - 1
       jclo = jlo/intraty(level-1) - 1
       kclo = klo/intratz(level-1) - 1
       ichi = (ihi+1)/intratx(level-1)
       jchi = (jhi+1)/intraty(level-1)
       kchi = (khi+1)/intratz(level-1)

       nrow = ichi-iclo+1
       ncol = jchi-jclo+1
       nfil = kchi-kclo+1

       xl   = rnode(cornxlo,mkid) - hxc
       yf   = rnode(cornylo,mkid) - hyc
       zb   = rnode(cornzlo,mkid) - hzc
       xr   = rnode(cornxhi,mkid) + hxc
       yr   = rnode(cornyhi,mkid) + hyc
       zt   = rnode(cornzhi,mkid) + hzc

       loctmp = igetsp(nrow*ncol*nfil*(nvar+naux))
       loctx  = loctmp + nrow*ncol*nfil*nvar
       locaux = node(storeaux,mkid)

       if (iclo .lt. 0 .or. ichi .eq. iregsz(levc) .or.
     1     jclo .lt. 0 .or. jchi .eq. jregsz(levc) .or.
     2     kclo .lt. 0 .or. kchi .eq. kregsz(levc)) then
           sticksout = .true.
       else
           sticksout = .false.
       endif

       if (xperdom .and. yperdom .and. zperdom .and. sticksout) then
         call preicall(alloc(loctmp),alloc(loctx),nrow,ncol,nfil,
     .                    nvar,naux,
     .                 iclo,ichi,jclo,jchi,kclo,kchi,level-1)
       else
         call icall(alloc(loctmp),alloc(loctx),nrow,ncol,nfil,
     .                 nvar,naux,
     .                   iclo,ichi,jclo,jchi,kclo,kchi,level-1,1,1,1)
          endif
c
          call bc3amr(alloc(loctmp),alloc(loctx),nrow,ncol,nfil,
     &          nvar,naux,hxc,hyc,hzc,level,time,
     .                xl,xr,yf,yr,zb,zt,
     .                xlower,ylower,zlower,xupper,yupper,zupper,
     .                xperdom,yperdom,zperdom)

       call cstore(alloc(loctmp),nrow,ncol,nfil,nvar,
     .                alloc(ist+nvar*lenbc),lenbc,naux,alloc(loctx),
     .                alloc(ist+2*nvar*lenbc))
       call reclam(loctmp,nrow*ncol*nfil*(nvar+naux))

          mkid = node(levelptr,mkid)
          go to 10
 99    return
       end
