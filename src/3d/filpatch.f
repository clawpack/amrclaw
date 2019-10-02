c
c ---------------------------------------------------------------
c
        recursive subroutine filrecur(level,nvar,valbig,aux,naux,
     1                                time,mitot,mjtot,mktot,
     2                                nrowst,ncolst,nfilst,ilo,ihi,
     3                                jlo,jhi,klo,khi,patchOnly,msrc)

c :::::::::::::::::::::::::::: FILPATCH ::::::::::::::::::::::::::
c
c  fill the portion of valbig from rows  nrowst
c                             and  cols  ncolst
c                             and  fils  nfilst
c  the patch can also be described by the corners
c  (xlp,yfp,zbp) by (xrp,yrp,ztp).
c  vals are needed at time time, and level level,
c
c  first fill with  values obtainable from the level level
c  grids. if any left unfilled, then enlarge remaining rectangle of
c  unfilled values by 1 (for later linear interp), and recusively
c  obtain the remaining values from  coarser levels.
c
c  Adapted from 2D recursive routine, 10/22/2012.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      use amr_module
      implicit double precision (a-h,o-z)


      logical   set, sticksout, patchOnly
      dimension valbig(nvar,mitot,mjtot,mktot)
      dimension    aux(naux,mitot,mjtot,mktot)
      
      !for setaux timing
      integer :: clock_start, clock_finish, clock_rate
      real(kind=8) :: cpu_start, cpu_finish


c     Use stack-based scratch arrays instead of alloc for better
c     OpenMP-friendliness and dynamic memory resizing.  These scratch
c     arrays are 1D instead of 3D because when we pass them to subroutines
c     they treat it as having different dimensions than the max size
c     that we need to allocate here.  The +3 is to expand the coarse
c     grid to enclose the fine, including an offset of 1.

      dimension valcrse((ihi-ilo+3)*(jhi-jlo+3)*(khi-klo+3)*nvar)
      dimension auxcrse((ihi-ilo+3)*(jhi-jlo+3)*(khi-klo+3)*naux)

      dimension flaguse(ihi-ilo+1, jhi-jlo+1, khi-klo+1)

c     Some helper functions
      ivalc(ivar, i, j, k) = ivar + nvar*(i-1) + nvar*nrowc*(j-1)
     &                            + nvar*nrowc*ncolc*(k-1)

      sticksout(iplo, iphi, jplo, jphi, kplo, kphi) =
     &     (iplo < 0 .or. jplo < 0 .or. kplo < 0 .or.
     &      iphi >= iregsz(levc) .or. jphi >= jregsz(levc)
     &      .or. kphi >= kregsz(levc))
      ! *** NOTE *** levc needs to be defined later as level-1

c     Begin by filling values for grids at level "level".  If all values
c     can be filled this way, return.

      hxf     = hxposs(level)
      hyf     = hyposs(level)
      hzf     = hzposs(level)
      xlp     = xlower + dble(ilo  )*hxf
      xrp     = xlower + dble(ihi+1)*hxf
      ybp     = ylower + dble(jlo  )*hyf
      ytp     = ylower + dble(jhi+1)*hyf
      zfp     = zlower + dble(klo  )*hzf
      zrp     = zlower + dble(khi+1)*hzf
      nrowp   = ihi - ilo + 1
      ncolp   = jhi - jlo + 1
      nfilp   = khi - klo + 1

      call intfil
     & (valbig, mitot, mjtot, mktot, time, flaguse,
     &  nrowst, ncolst, nfilst, ilo, ihi, jlo, jhi, klo, khi,
     &  level, nvar, naux,msrc)

c Trimbd returns set = true if all of the entries are filled (=1.).
c set = false, otherwise. If set = true, then no other levels are
c are required to interpolate, and we return. If false, the minimum
c enclosing rectangle is returned in (il,jl,kl) x (ir,jr,kr)
c
c Note that the used array is filled entirely in intfil, i.e. the
c marking done there also takes into account the points filled by
c the boundary conditions. physbd will be called later (from bound), after
c all 4 boundary pieces filled.

      call trimbd(flaguse, set, il, ir, jl, jr, kl, kr,
     &            ilo, ihi, jlo, jhi, klo, khi)

      if (set) go to 90    ! All done except for BCs

c
c otherwise make recursive calls to coarser levels to fill remaining unset points
c
      if (level == 1) then
           write(outunit,*)" error in filrecur - level 1 not set"
           write(outunit,900) nrowst,ncolst,nfilst
           write(*,*)" error in filrecur - level 1 not set"
           write(*,*)" should not need more recursion "
           write(*,*)" to set patch boundaries"
           write(*,900) nrowst,ncolst,nfilst
900        format("start at row: ",i4," col ",i4," file ",i4)
           stop
      end if

c set = false. we will have to interpolate some values from coarser
c levels. We begin by initializing the level level arrays, so that we can use
c purely recursive formulation for interpolating.

      levc = level - 1
      hxc  = hxposs(levc)
      hyc  = hyposs(levc)
      hzc  = hzposs(levc)

c
c     coarsen
      lratiox = intratx(levc)
      lratioy = intraty(levc)
      lratioz = intratz(levc)
      iplo   = (il-lratiox+nghost*lratiox)/lratiox - nghost
      jplo   = (jl-lratioy+nghost*lratioy)/lratioy - nghost
      kplo   = (kl-lratioz+nghost*lratioz)/lratioz - nghost
      iphi   = (ir+lratiox)/lratiox
      jphi   = (jr+lratioy)/lratioy
      kphi   = (kr+lratioz)/lratioz

      xlc  =  xlower + iplo*hxc
      ybc  =  ylower + jplo*hyc
      zfc  =  zlower + kplo*hzc

      nrowc = iphi - iplo + 1
      ncolc = jphi - jplo + 1
      nfilc = kphi - kplo + 1
      ntot  = nrowc*ncolc*nfilc*(nvar+naux)

      if (nrowc > ihi-ilo+3 .or. ncolc > jhi-jlo+3
     &    .or. nfilc > khi-klo+3) then
         write(*,*)" did not make big enough work space in filrecur "
         write(*,*)" need coarse space with nrowc,ncolc,nfilc ",
     &             nrowc,ncolc,nfilc
         write(6,*)" made space for ilo,ihi,jlo,jhi,klo,khi ",
     &             ilo,ihi,jlo,jhi,klo,khi
         stop
      endif

      if (naux > 0) then
         mx = nrowc - 2*nghost
         my = ncolc - 2*nghost
         mz = nfilc - 2*nghost
         xl = xlc + nghost*hxc
         yb = ybc + nghost*hyc
         zf = zfc + nghost*hzc
         
         call system_clock(clock_start, clock_rate)
         call cpu_time(cpu_start)
         call setaux(nghost,mx,my,mz,xl,yb,zf,
     &               hxc,hyc,hzc,naux,auxcrse)
         call system_clock(clock_finish, clock_rate)
         call cpu_time(cpu_finish)
      endif

      if ((xperdom .or. yperdom .or. zperdom) .and.
     &     sticksout(iplo,iphi,jplo,jphi,kplo,kphi)) then
         call prefilrecur(levc,nvar,valcrse,auxcrse,
     1                    naux,time,nrowc,ncolc,nfilc,1,1,1,
     2                    iplo,iphi,jplo,jphi,kplo,kphi,
     3                    iplo,iphi,jplo,jphi,kplo,kphi,
     4                    .true.)
      else  !-1 indicates not a real grid in this next call to filpatch
         call filrecur(levc,nvar,valcrse,auxcrse,naux,
     1                 time,nrowc,ncolc,nfilc,1,1,1,
     2                 iplo,iphi,jplo,jphi,kplo,kphi,.true.,-1)
      endif

      ! convert to real for use in floor function below
      ratiox = float(lratiox)
      ratioy = float(lratioy)
      ratioz = float(lratioz)

      do 100 iff = 1,nrowp
         !ic = 2 + (iff-(il-ilo)-1)/lratiox
         ic =floor((iff+ilo-1)/ratiox) - iplo + 1
         !eta1 = (-0.5d0+dble(mod(iff-1,lratiox)))/dble(lratiox)
         xcent_coarse = xlc + (ic-.5d0)*hxc
         xcent_fine =  xlower + (iff-1+ilo + .5d0)*hxf
         eta1 = (xcent_fine-xcent_coarse)/hxc
         if (abs(eta1) .gt. .5) then
            write(*,*)" filpatch x indices wrong: eta1 = ",eta1
         endif


         do 100 jf  = 1,ncolp
          !jc = 2 + (jf -(jl-jlo)-1)/lratioy
          jc =floor((jf+jlo-1)/ratioy) - jplo + 1
          !eta2 = (-0.5d0+dble(mod(jf -1,lratioy)))/dble(lratioy)
          ycent_coarse = ybc + (jc-.5d0)*hyc
          ycent_fine =  ylower + (jf-1+jlo + .5d0)*hyf
          eta2 = (ycent_fine-ycent_coarse)/hyc
         if (abs(eta2) .gt. .5) then
            write(*,*)" filpatch y indices wrong: eta2 = ",eta2
         endif

          do 100 kf = 1,nfilp
           !kc = 2 + (kf - (kl-klo)-1)/lratioz
           kc =floor((kf+klo-1)/ratioz) - kplo + 1
           !eta3 = (-0.5d0+dble(mod(kf-1, lratioz)))/dble(lratioz)
           zcent_coarse = zfc + (kc-.5d0)*hzc
           zcent_fine =  zlower + (kf-1+klo + .5d0)*hzf
           eta3 = (zcent_fine-zcent_coarse)/hzc
           if (abs(eta3) .gt. .5) then
              write(*,*)" filpatch z indices wrong: eta3 = ",eta3
           endif

           flag = flaguse(iff,jf,kf)
           if (flag .eq. 0.0) then

            do 101 ivar = 1,nvar

               valp100 = valcrse(ivalc(ivar,ic+1,jc  ,kc   ))
               valm100 = valcrse(ivalc(ivar,ic-1,jc  ,kc   ))
               valc    = valcrse(ivalc(ivar,ic  ,jc  ,kc   ))
               valp010 = valcrse(ivalc(ivar,ic  ,jc+1,kc   ))
               valm010 = valcrse(ivalc(ivar,ic  ,jc-1,kc   ))
               valp001 = valcrse(ivalc(ivar,ic  ,jc  ,kc+1 ))
               valm001 = valcrse(ivalc(ivar,ic  ,jc  ,kc-1 ))

               ! Use monotonized centered limiter to reconstruct in all axes
               dupc = valp100 - valc
               dumc = valc    - valm100
               ducc = valp100 - valm100
               du   = dmin1(dabs(dupc),dabs(dumc))
               du   = dmin1(2.d0*du,.5d0*dabs(ducc))
               fu   = dmax1(0.d0,dsign(1.d0,dupc*dumc))
               uslope = du*sign(1.d0,ducc)*fu ! not really-should divide by h

               dvpc = valp010 - valc
               dvmc = valc    - valm010
               dvcc = valp010 - valm010
               dv   = dmin1(dabs(dvpc),dabs(dvmc))
               dv   = dmin1(2.d0*dv,.5d0*dabs(dvcc))
               fv   = dmax1(0.d0,dsign(1.d0,dvpc*dvmc))
               vslope = dv*sign(1.d0,dvcc)*fv  ! but faster to put with eta above

               dwpc = valp001 - valc
               dwmc = valc    - valm001
               dwcc = valp001 - valm001
               dw   = dmin1(dabs(dwpc),dabs(dwmc))
               dw   = dmin1(2.d0*dw,.5d0*dabs(dwcc))
               fw   = dmax1(0.d0,dsign(1.d0,dwpc*dwmc))
               wslope = dw*sign(1.d0,dwcc)*fw ! instead of recomputing

               valint = valc + eta1*uslope
     &                       + eta2*vslope
     &                       + eta3*wslope

               valbig(ivar,iff+nrowst-1,jf+ncolst-1,kf+nfilst-1)
     &                = valint

 101        continue

         endif

 100  continue

 90   continue
c     
c    ## only call if a small coarser recursive patch
c    ## otherwise whole grid bcs done from bound
      if (patchOnly) then
         call bc3amr(valbig,aux,mitot,mjtot,mktot,nvar,naux,
     1               hxf,hyf,hzf,level,time,
     2               xlp,xrp,ybp,ytp,zfp,zrp)
      endif

      return
      end
