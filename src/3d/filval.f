c
c ------------------------------------------------------------------
c
      subroutine filval(val,mitot,mjtot,mktot,hx,hy,hz,lev,time,
     1                  valc,auxc,mic,mjc,mkc,
     2                  xleft,xright,yfront,yrear,zbot,ztop,nvar,
     3                  mptr,ilo,ihi,jlo,jhi,klo,khi,aux,naux)

      use amr_module
      implicit double precision (a-h,o-z)


      dimension   val(nvar,mitot,mjtot,mktot), valc(nvar,mic,mjc,mkc)
      dimension   aux(nvar,mitot,mjtot,mktot), auxc(nvar,mic,mjc,mkc)
      dimension   dudx(max1d), dudy(max1d), dudz(max1d)
c
c :::::::::::::::::::::::::::::: FILVAL ::::::::::::::::::::::::::
c
c create and fill coarser (lev-1) patch with one extra coarse cell all
c around, plus the ghost cells . will interpolate from this patch to grid mptr
c without needing special boundary code.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      levc    = lev - 1
      lratiox = intratx(levc)
      lratioy = intraty(levc)
      lratioz = intratz(levc)
      hxcrse  = hx*lratiox
      hycrse  = hy*lratioy
      hzcrse  = hz*lratioz
      xl      = xleft  - hxcrse
      xr      = xright + hxcrse
      yf      = yfront - hycrse
      yr      = yrear  + hycrse
      zb      = zbot   - hzcrse
      zt      = ztop   + hzcrse
c
c     set integer indices for coarser patch enlarged by 1 cell
c     (can stick out of domain). proper nesting will insure this one
c     call is sufficient.
      iclo   = ilo/lratiox - 1
      jclo   = jlo/lratioy - 1
      kclo   = klo/lratioz - 1
      ichi   = (ihi+1)/lratiox - 1 + 1
      jchi   = (jhi+1)/lratioy - 1 + 1
      kchi   = (khi+1)/lratioz - 1 + 1
      ng     = 0

c    :::  mcapa  is the capacity function index

      if (mcapa .eq. 0) then
        if (xperdom .or. yperdom .or. zperdom) then
          call preintcopy(valc,mic,mjc,mkc,nvar,iclo,ichi,jclo,jchi,
     &                                          kclo,kchi,levc)
        else
          call intcopy(valc,mic,mjc,mkc,nvar,iclo,ichi,jclo,jchi,
     &                                        kclo,kchi,levc,1,1,1)
        endif
      else
        if (xperdom .or. yperdom .or. zperdom) then
          call preicall(valc,auxc,mic,mjc,mkc,nvar,naux,iclo,ichi,jclo,
     &                  jchi,kclo,kchi,levc)
        else
          call icall(valc,auxc,mic,mjc,mkc,nvar,naux,iclo,ichi,jclo,
     &               jchi,kclo,kchi,levc,1,1,1)
        endif
      endif
c      call physbd(valc,auxc,mic,mjc,mkc,nvar,naux,
c     1            hxcrse,hycrse,hzcrse,levc,time,
c     2            xl,xr,yf,yr,zb,zt,
c     3            xlower,ylower,zlower,xupper,yupper,zupper,
c     4            xperdom,yperdom,zperdom)

c 2/28/02 : new call to bc3amr
      call bc3amr(valc,auxc,mic,mjc,mkc,nvar,naux,
     1            hxcrse,hycrse,hzcrse,levc,time,
     2            xl,xr,yf,yr,zb,zt)


c
c     prepare slopes - use min-mod limiters
c
      do 30 k    = 2, mkc-1
      do 30 j    = 2, mjc-1
      do 30 ivar = 1, nvar
      do 10 i    = 2, mic-1

         slp = valc(ivar,i+1,j,k) - valc(ivar,i  ,j,k)
         slm = valc(ivar,i  ,j,k) - valc(ivar,i-1,j,k)
         slopex = dmin1(dabs(slp),dabs(slm))*
     .            dsign(1.0d0,valc(ivar,i+1,j,k) - valc(ivar,i-1,j,k))
c        # if theres a sign change, set slope to 0.
         if ( slm*slp .gt. 0.d0) then
           dudx(i) = slopex
         else
           dudx(i) = 0.d0
         endif

         slp = valc(ivar,i,j+1,k) - valc(ivar,i,j  ,k)
         slm = valc(ivar,i,j  ,k) - valc(ivar,i,j-1,k)
         slopey = dmin1(dabs(slp),dabs(slm))*
     .            dsign(1.0d0,valc(ivar,i,j+1,k) - valc(ivar,i,j-1,k))
         if ( slm*slp .gt. 0.d0) then
           dudy(i) = slopey
         else
           dudy(i) = 0.d0
         endif

         slp = valc(ivar,i,j,k+1) - valc(ivar,i,j,k  )
         slm = valc(ivar,i,j,k  ) - valc(ivar,i,j,k-1)
         slopez = dmin1(dabs(slp),dabs(slm))*
     .            dsign(1.0d0,valc(ivar,i,j,k+1) - valc(ivar,i,j,k-1))
         if ( slm*slp .gt. 0.d0) then
           dudz(i) = slopez
         else
           dudz(i) = 0.d0
         endif

 10   continue
c
c     interp. from coarse cells to fine grid
c
      do 20 kco = 1,lratioz
      kfine = (k-2)*lratioz + nghost + kco
      zoff  = (dble(kco) - .5d0)/lratioz - .5d0
         do 20 jco = 1,lratioy
         jfine = (j-2)*lratioy + nghost + jco
         yoff  = (dble(jco) - .5d0)/lratioy - .5d0
            do 20 ico = 1,lratiox
            xoff = (dble(ico) - .5d0)/lratiox - .5d0
               do 20 i = 2, mic-1
               ifine   = (i-2)*lratiox + nghost + ico
               val(ivar,ifine,jfine,kfine) = valc(ivar,i,j,k)
     1                                   + xoff*dudx(i) + yoff*dudy(i)
     2                                                  + zoff*dudz(i)
 20   continue
c
 30   continue

      if (mcapa .ne. 0) then
        call fixcapaq(val,aux,mitot,mjtot,mktot,valc,auxc,mic,mjc,mkc,
     &                nvar,naux,levc)
      endif
c
c  overwrite interpolated values with fine grid values, if available.
c
      call intcopy(val,mitot,mjtot,mktot,nvar,ilo-nghost,ihi+nghost,
     &             jlo-nghost,jhi+nghost,klo-nghost,khi+nghost,
     &             lev,1,1,1)

      return
      end
