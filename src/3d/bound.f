c
c --------------------------------------------------------------
c
      subroutine bound(time,nvar,ng,valbig,mitot,mjtot,mktot,mptr,
     1                 aux,naux)

c
      implicit double precision (a-h,o-z)

      include  "call.i"

      dimension valbig(mitot,mjtot,mktot,nvar)
      dimension aux   (mitot,mjtot,mktot,naux)
      logical periodic

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

      periodic = .false.
c     left boundary
      xl = xlo - ng*hx
      xr = xlo
      yf = ylo - ng*hy
      yr = yhi + ng*hy
      zb = zlo
      zt = zhi
      periodic = (periodic) .or. (xl .lt. xlower .and. xperdom)
c     right boundary
      xl = xhi
      xr = xhi + ng*hx
      yf = ylo - ng*hy
      yr = yhi + ng*hy
      zb = zlo
      zt = zhi
      periodic = (periodic) .or. (xr .gt. xupper .and. xperdom)
c     front boundary
      xl = xlo
      xr = xhi
      yf = ylo - ng*hy
      yr = ylo
      zb = zlo
      zt = zhi
      periodic = (periodic) .or. (yf .lt. ylower .and. yperdom)
c     rear boundary
      xl = xlo
      xr = xhi
      yf = yhi
      yr = yhi + ng*hy
      zb = zlo
      zt = zhi
      periodic = (periodic) .or. (yr .gt. yupper .and. yperdom)
c     bottom boundary
      xl = xlo - ng*hx
      xr = xhi + ng*hx
      yf = ylo - ng*hy
      yr = yhi + ng*hy
      zb = zlo - ng*hz
      zt = zlo
      periodic = (periodic) .or. (zb .lt. zlower .and. zperdom)
c     top boundary
      xl = xlo - ng*hx
      xr = xhi + ng*hx
      yf = ylo - ng*hy
      yr = yhi + ng*hy
      zb = zhi
      zt = zhi + ng*hz
      periodic = (periodic) .or. (zt .gt. zupper .and. zperdom)

c     left boundary

      xl = xlo - ng*hx
      xr = xlo
      yf = ylo - ng*hy
      yr = yhi + ng*hy
      zb = zlo
      zt = zhi

      if (periodic) then
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

      if (periodic) then
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
      if (periodic) then
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
      if (periodic) then
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
      if (periodic) then
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
      if (periodic) then
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
c external boundary conditions
c
c$$$      if (.not. ((xperdom .and. yperdom) .and. zperdom)) then
c$$$	xl = xlo    - ng*hx
c$$$	yf = ylo    - ng*hy
c$$$        zb = zlo    - ng*hz
c$$$	xr = xhi    + ng*hx
c$$$	yr = yhi    + ng*hy
c$$$        zt = zhi    + ng*hz
c$$$
c$$$	call bc3amr( valbig,aux,mitot,mjtot,mktot,nvar,naux,
c$$$     1        	     hx,hy,hz,level,time,xl,xr,yf,yr,zb,zt,
c$$$     3               xlower,ylower,zlower,xupper,yupper,zupper,
c$$$     4               xperdom,yperdom,zperdom)
c$$$      endif
c
      return
      end
