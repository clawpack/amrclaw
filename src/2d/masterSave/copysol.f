c
c ----------------------------------------------------------
c
      subroutine copysol(valbig,val,nvar,mitot,mjtot,nghost,
     1                   midub,mjdub,ngbig)
c
      implicit double precision (a-h,o-z)

      dimension  valbig(nvar,midub,mjdub), val(nvar,mitot,mjtot)
c
c copy solution into grid with different number ghsot cells
c
       do 10 j = nghost+1, mjtot-nghost
       do 10 i = nghost+1, mitot-nghost
       do 10 ivar = 1, nvar
          valbig(ivar,i-nghost+ngbig,j-nghost+ngbig) = val(ivar,i,j)
 10    continue
c
       return
       end
