c
c --------------------------------------------------------------
c
      subroutine bound(time,nvar,ng,valbig,mitot,mjtot,mktot,mptr,
     1                 aux,naux)

c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension valbig(nvar,mitot,mjtot,mktot)
      dimension aux   (naux,mitot,mjtot,mktot)
      logical   xsticksout, ysticksout
c
c  :::::::::::::: BOUND :::::::::::::::::::::::::::::::::::::::::::
c     This routine sets the boundary values for a given grid
c     at level level.
c     We are setting the values for a strip ng zones wide all
c     the way around the border, in 4 rectangular strips.
c
c     Outputs from this routine:
c     The values around the border of the grid are inserted
c     directly into the enlarged valbig array.
c
c     This routine calls the routine filpatch
c     which for any block of mesh points on a given level,
c     intersects that block with all grids on that level and with
c     the physical boundaries, copies the values into the
c     appropriate intersecting regions, and interpolates the remaining
c     cells from coarser grids as required.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      xlo    = rnode(cornxlo, mptr)
      xhi    = rnode(cornxhi, mptr)
      ylo    = rnode(cornylo, mptr)
      yhi    = rnode(cornyhi, mptr)
      zlo    = rnode(cornzlo, mptr)
      zhi    = rnode(cornzhi, mptr)
      ilo    = node(ndilo, mptr)
      ihi    = node(ndihi, mptr)
      jlo    = node(ndjlo, mptr)
      jhi    = node(ndjhi, mptr)
      klo    = node(ndklo, mptr)
      khi    = node(ndkhi, mptr)
      level  = node(nestlevel, mptr)
      hx     = hxposs(level)
      hy     = hyposs(level)
      hz     = hzposs(level)

c     left boundary

      xl = xlo - ng*hx
      xr = xlo
      yf = ylo - ng*hy
      yr = yhi + ng*hy
      zb = zlo
      zt = zhi
      ysticksout =  ((yf .lt. ylower) .or. (yr .gt. yupper))

      if ((xperdom .and. xl .lt. xlower) .or. 
     1    (yperdom .and. ysticksout)) then
	      call  prefilrecur(level,nvar,valbig,aux,naux,time,
     1                    mitot,mjtot,mktot,
     2                    1,1,ng+1,
     3                    ilo-ng,ilo-1,jlo-ng,jhi+ng,klo,khi)
      else
      	call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,mktot,
     1                1,1,ng+1,
     2                ilo-ng,ilo-1,jlo-ng,jhi+ng,klo,khi)
      endif


c     right boundary

      xl = xhi
      xr = xhi + ng*hx
      yf = ylo - ng*hy
      yr = yhi + ng*hy
      zb = zlo
      zt = zhi
      ysticksout =  ((yf .lt. ylower) .or. (yr .gt. yupper))

      if ((xperdom .and. xr .gt. xupper) .or. 
     1    (yperdom .and. ysticksout)) then
      	call  prefilrecur(level,nvar,valbig,aux,naux,time,
     1                    mitot,mjtot,mktot,
     2                    mitot-ng+1,1,ng+1,
     3                    ihi+1,ihi+ng,jlo-ng,jhi+ng,klo,khi)
      else
      	call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,mktot,
     1                mitot-ng+1,1,ng+1,
     2                ihi+1,ihi+ng,jlo-ng,jhi+ng,klo,khi)
      endif


c     front boundary
      xl = xlo
      xr = xhi
      yf = ylo - ng*hy
      yr = ylo
      zb = zlo
      zt = zhi
      if (yperdom .and. yf .lt. ylower) then
        call prefilrecur(level,nvar,valbig,aux,naux,time,
     1                   mitot,mjtot,mktot,
     2                   ng+1,1,ng+1,
     3                   ilo,ihi,jlo-ng,jlo-1,klo,khi)
      else
        call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,mktot,
     1                ng+1,1,ng+1,
     2                ilo,ihi,jlo-ng,jlo-1,klo,khi)
      endif

c     rear boundary
      xl = xlo
      xr = xhi
      yf = yhi
      yr = yhi + ng*hy
      zb = zlo
      zt = zhi
      if (yperdom .and. yr .gt. yupper) then
      	call prefilrecur(level,nvar,valbig,aux,naux,time,
     1                   mitot,mjtot,mktot,
     2                   ng+1,mjtot-ng+1,ng+1,
     3                   ilo,ihi,jhi+1,jhi+ng,klo,khi)
      else
        call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,mktot,
     1                ng+1,mjtot-ng+1,ng+1,
     2                ilo,ihi,jhi+1,jhi+ng,klo,khi)
      endif

c     bottom boundary
      xl = xlo - ng*hx
      xr = xhi + ng*hx
      yf = ylo - ng*hy
      yr = yhi + ng*hy
      zb = zlo - ng*hz
      zt = zlo
      xsticksout = (xl .lt. xlower) .or. (xr .gt. xupper)
      ysticksout = (yf .lt. ylower) .or. (yr .gt. yupper)

      if ((zperdom .and. zb .lt. zlower) .or. 
     1    (xperdom .and. xsticksout) .or.
     2    (yperdom .and. ysticksout)) then
        call prefilrecur(level,nvar,valbig,aux,naux,time,
     1                   mitot,mjtot,mktot,
     2                   1,1,1,
     3                   ilo-ng,ihi+ng,jlo-ng,jhi+ng,klo-ng,klo-1)
      else
        call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,mktot,
     1                1,1,1,
     2                ilo-ng,ihi+ng,jlo-ng,jhi+ng,klo-ng,klo-1)
      end if

c     top boundary
      xl = xlo - ng*hx
      xr = xhi + ng*hx
      yf = ylo - ng*hy
      yr = yhi + ng*hy
      zb = zhi
      zt = zhi + ng*hz
      xsticksout = (xl .lt. xlower) .or. (xr .gt. xupper)
      ysticksout = (yf .lt. ylower) .or. (yr .gt. yupper)

      if ((zperdom .and. zt .gt. zupper) .or.
     1    (xperdom .and. xsticksout) .or.
     2    (yperdom .and. ysticksout)) then
        call prefilrecur(level,nvar,valbig,aux,naux,time,
     1                   mitot,mjtot,mktot,
     2                   1,1,mktot-ng+1,
     3                   ilo-ng,ihi+ng,jlo-ng,jhi+ng,khi+1,khi+ng)
      else
        call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,mktot,
     1                1,1,mktot-ng+1,
     2                ilo-ng,ihi+ng,jlo-ng,jhi+ng,khi+1,khi+ng)
      end if

c
c
      return
      end
