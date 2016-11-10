c
c ----------------------------------------------------------
c
      subroutine copysol(valbig,val,nvar,mitot,nghost,
     1                   midub,ngbig)
c
      implicit double precision (a-h,o-z)

      dimension  valbig(nvar,midub), val(nvar,mitot)
c
c copy solution into grid with different number ghost cells
c
       do 10 i = nghost+1, mitot-nghost
       do 10 ivar = 1, nvar
          valbig(ivar,i-nghost+ngbig) = val(ivar,i)
 10    continue
c
       return
       end
