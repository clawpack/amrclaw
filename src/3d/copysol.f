c
c ----------------------------------------------------------
c
      subroutine copysol(valbig,val,nvar,mitot,mjtot,mktot,nghost,
     1                                   midub,mjdub,mkdub,ngbig)
c
      implicit double precision (a-h,o-z)

      dimension  valbig(midub,mjdub,mkdub,nvar)
      dimension  val   (mitot,mjtot,mktot,nvar)
c
c copy solution into grid with different number of ghost cells
c
       do 10 ivar = 1, nvar
       do 10 k = nghost+1, mktot-nghost
       do 10 j = nghost+1, mjtot-nghost
       do 10 i = nghost+1, mitot-nghost
	  valbig(i-nghost+ngbig,j-nghost+ngbig,k-nghost+ngbig,ivar)
     &  = val   (i             ,j             ,k             ,ivar)
 10    continue
c
       return
       end
