c
c ----------------------------------------------------------
c
      subroutine copysol(valbig,val,nvar,mitot,mjtot,mktot,nghost,
     1                                   midub,mjdub,mkdub,ngbig)
c
      implicit double precision (a-h,o-z)

      dimension  valbig(nvar,midub,mjdub,mkdub)
      dimension  val   (nvar,mitot,mjtot,mktot)
c
c copy solution into grid with different number of ghost cells
c
       do 10 k = nghost+1, mktot-nghost
       do 10 j = nghost+1, mjtot-nghost
       do 10 i = nghost+1, mitot-nghost
       do 10 ivar = 1, nvar
          valbig(ivar,i-nghost+ngbig,j-nghost+ngbig,k-nghost+ngbig) =
     &                                                val(ivar,i,j,k)
 10    continue
c
       return
       end
