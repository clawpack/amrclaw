      subroutine src3(meqn,mbc,mx,my,mz,
     &      xlower, ylower, zlower, dx, dy, dz, q,maux,aux,t,dt)

      implicit double precision (a-h,o-z)
      dimension    q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      dimension  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

c
c      # dummy subroutine for use when equation has no source term.
c      # If method(5)=0 then this routine is never called, but its
c      # existence may be required by some compilers.
c
      return
      end
