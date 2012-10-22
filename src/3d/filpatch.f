c
c ---------------------------------------------------------------
c
      subroutine filpatch(level,nvar,valbig,aux,naux,
     1                    time,mitot,mjtot,mktot,
     2                    nrowst,ncolst,nfilst,ilo,ihi,jlo,jhi,klo,khi)

c :::::::::::::::::::::::::::: FILPATCH :::::::::::::::::::::::::;
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
c :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::;

      implicit double precision (a-h,o-z)

      include  "call.i"

      logical   set, sticksout
      dimension valbig(mitot,mjtot,mktot,nvar)
      dimension aux   (mitot,mjtot,mktot,naux)

      iadflag(i,j,k)    =  locuse  +    (i-1)
     &                             +    (j-1)*nrowp
     &                             +    (k-1)*nrowp*ncolp
      ivalc(i,j,k,ivar) =  loccrse +    (i-1)
     &                             +    (j-1)*nrowc
     &                             +    (k-1)*nrowc*ncolc
     &                             + (ivar-1)*nrowc*ncolc*nfilc
      sticksout(iplo,iphi,jplo,jphi,kplo,kphi)  =
     &          (iplo .lt. 0 .or. jplo .lt. 0 .or. kplo .lt. 0 .or.
     &           iphi .ge. iregsz(levc) .or. jphi .ge. jregsz(levc) .or.
     &           kphi .ge. kregsz(levc))
c
c We begin by filling values for grids at level level. If all values can be
c filled in this way, we return;

        hxf     = hxposs(level)
        hyf     = hyposs(level)
        hzf     = hzposs(level)
	xlp     = xlower + dble(ilo  )*hxf
	xrp     = xlower + dble(ihi+1)*hxf
	yfp     = ylower + dble(jlo  )*hyf
	yrp     = ylower + dble(jhi+1)*hyf
        zbp     = zlower + dble(klo  )*hzf
        ztp     = zlower + dble(khi+1)*hzf
	nrowp   = ihi - ilo + 1
	ncolp   = jhi - jlo + 1
        nfilp   = khi - klo + 1
        locuse  = igetsp(nrowp*ncolp*nfilp)

        call intfil
     &  (valbig,mitot,mjtot,mktot,time,locuse,nrowst,ncolst,nfilst,
     &   ilo,ihi,jlo,jhi,klo,khi,level,nvar,naux)
c
c Trimbd returns set = true if all of the entries are filled (=1.).
c set = false, otherwise. If set = true, then no other levels are
c are required to interpolate, and we return.
c
c Note that the used array is filled entirely in intfil, i.e. the
c marking done there also takes into account the points filled by
c the boundary conditions. physbd will be called later (from bound), after
c all 4 boundary pieces filled.

        call trimbd(alloc(locuse),set,isl,isr,jsf,jsr,ksb,kst,
     &                                ilo,ihi,jlo,jhi,klo,khi)
        if (set) then
           call reclam(locuse,nrowp*ncolp*nfilp)
           return
	else if (level .eq. 1) then
	   write(outunit,*)" error in filpatch - level 1 not set"
	   write(outunit,900) nrowst,ncolst,nfilst
	   write(*,*)" error in filpatch - level 1 not set"
	   write(*,900) nrowst,ncolst,nfilst
900        format("starting at row: ",i4," col ",i4," fil ",i4)
	   stop
        endif

c set = false. we will have to interpolate some values from coarser
c levels. We begin by initializing the level level arrays, so that we can use
c recursive formulation for interpolating.
c IS THIS TRUE ? - the fine grid patch remaining unfilled is always
c anchored to a coarse cell.

        levc = level - 1

c 2/27/02 : uncommented these three lines.
        hxc  = hxposs(levc)
        hyc  = hyposs(levc)
        hzc  = hzposs(levc)

c       isl  = il + ilo - 1
c       isr  = ir + ilo - 1
c       jsf  = jf + jlo - 1
c       jsr  = jr + jlo - 1
c       ksb  = kb + klo - 1
c       kst  = kt + klo - 1

c       coarsen
	lratiox = intratx(levc)
	lratioy = intraty(levc)
	lratioz = intratz(levc)
        iplo   = (isl-lratiox  + nghost*lratiox)/lratiox - nghost
        jplo   = (jsf-lratioy  + nghost*lratioy)/lratioy - nghost
        kplo   = (ksb-lratioz  + nghost*lratioz)/lratioz - nghost
	iphi   = (isr+lratiox  )/lratiox
	jphi   = (jsr+lratioy  )/lratioy
	kphi   = (kst+lratioz  )/lratioz

*
        xlc  =  xlower + dble(iplo  )*hxc
        yfc  =  ylower + dble(jplo  )*hyc
        zbc  =  zlower + dble(kplo  )*hzc
        xrc  =  xlower + dble(iphi+1)*hxc
        yrc  =  ylower + dble(jphi+1)*hyc
        ztc  =  zlower + dble(kphi+1)*hzc

	nrowc   =  iphi - iplo + 1
	ncolc   =  jphi - jplo + 1
        nfilc   =  kphi - kplo + 1
        ntot    = nrowc*ncolc*nfilc*(nvar+naux)
        loccrse = igetsp(ntot)
	locauxc = loccrse + nrowc*ncolc*nfilc*nvar
	if (naux.gt.0) then
           mx = nrowc
           my = ncolc
           mz = nfilc
           mbc = 0
           call setaux(mx,my,mz,mbc,mx,my,mz,xlc,yfc,zbc,
     &           hxc,hyc,hzc, naux,alloc(locauxc))
	endif


	if ((xperdom .or. yperdom .or. zperdom) .and.
     &       sticksout(iplo,iphi,jplo,jphi,kplo,kphi)) then
	     call prefil2(levc,nvar,alloc(loccrse),alloc(locauxc),naux,
     1                    time,nrowc,ncolc,nfilc,1,1,1,
     2                    iplo,iphi,jplo,jphi,kplo,kphi)
	else
          call filpatch2(levc,nvar,alloc(loccrse),alloc(locauxc),naux,
     1                   time,nrowc,ncolc,nfilc,1,1,1,
     2                   iplo,iphi,jplo,jphi,kplo,kphi)
	endif


c       interpolate back up

20      continue

        do 100 iff = 1,nrowp
          ic = 2 + (iff-(isl-ilo)-1)/lratiox
          eta1 = (-0.5d0+dble(mod(iff-1,lratiox)))/dble(lratiox)
        do 100 jf  = 1,ncolp
          jc = 2 + (jf -(jsf-jlo)-1)/lratioy
          eta2 = (-0.5d0+dble(mod(jf -1,lratioy)))/dble(lratioy)
        do 100 kf  = 1,nfilp
          kc = 2 + (kf -(ksb-klo)-1)/lratioz
          eta3 = (-0.5d0+dble(mod(kf -1,lratioz)))/dble(lratioz)

   	flag = alloc(iadflag(iff,jf,kf))
	if (flag .eq. 0.0) then

*         xif = xlp + (.5d0 + dble(iff-1))*hxf
*         yjf = yfp + (.5d0 + dble(jf -1))*hyf
*         zkf = zbp + (.5d0 + dble(kf -1))*hzf

*         ic=idint((xif-xlc+.5d0*hxc)/hxc)
*         jc=idint((yjf-yfc+.5d0*hyc)/hyc)
*         kc=idint((zkf-zbc+.5d0*hzc)/hzc)

*         xc = xlc + (dble(ic) - .5d0)*hxc
*         yc = yfc + (dble(jc) - .5d0)*hyc
*         zc = zbc + (dble(kc) - .5d0)*hzc

*         eta1 = (xif - xc)/hxc
*         eta2 = (yjf - yc)/hyc
*         eta3 = (zkf - zc)/hzc

          do 101 ivar = 1,nvar

          val111 = alloc(ivalc(ic  ,jc  ,kc  ,ivar))
          valp11 = alloc(ivalc(ic+1,jc  ,kc  ,ivar))
          valm11 = alloc(ivalc(ic-1,jc  ,kc  ,ivar))
          val1p1 = alloc(ivalc(ic  ,jc+1,kc  ,ivar))
          val1m1 = alloc(ivalc(ic  ,jc-1,kc  ,ivar))
          val11p = alloc(ivalc(ic  ,jc  ,kc+1,ivar))
          val11m = alloc(ivalc(ic  ,jc  ,kc-1,ivar))

          dupc =   valp11 - val111
          dumc =   val111 - valm11
          ducc =   valp11 - valm11
          du   = dmin1(dabs(dupc),dabs(dumc))
          du   = dmin1(2.d0*du, .5d0*dabs(ducc))
          fu   = dmax1(0.d0,dsign(1.d0, dupc*dumc))

          dvpc =   val1p1 - val111
          dvmc =   val111 - val1m1
          dvcc =   val1p1 - val1m1
          dv   = dmin1(dabs(dvpc),dabs(dvmc))
          dv   = dmin1(2.d0*dv, .5d0*dabs(dvcc))
          fv   = dmax1(0.d0,dsign(1.d0, dvpc*dvmc))
          dwpc =   val11p - val111
          dwmc =   val111 - val11m
          dwcc =   val11p - val11m
          dw   = dmin1(dabs(dwpc),dabs(dwmc))
          dw   = dmin1(2.d0*dw, .5d0*dabs(dwcc))
          fw   = dmax1(0.d0,dsign(1.d0, dwpc*dwmc))

          valint  = val111 + eta1*du*dsign(1.d0,ducc)*fu
     2                     + eta2*dv*dsign(1.d0,dvcc)*fv
     3                     + eta3*dw*dsign(1.d0,dwcc)*fw

          valbig(iff+nrowst-1,jf+ncolst-1,kf+nfilst-1,ivar) = valint

101       continue

        endif

100     continue

        call reclam(loccrse,ntot)

        call reclam(locuse,nrowp*ncolp*nfilp)

        return
        end
