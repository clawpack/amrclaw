c
c -----------------------------------------------------------
c
      subroutine conck(level, nvar)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      iadd(i,j,k,ivar)   = loc    +    (i-1)
     &                            +    (j-1)*mitot
     &                            +    (k-1)*mitot*mjtot
     &                            + (ivar-1)*mitot*mjtot*mktot
      iaddaux(i,j,k)     = locaux +    (i-1)
     &                            +    (j-1)*mitot
     &                            +    (k-1)*mitot*mjtot
     &                            +(mcapa-1)*mitot*mjtot*mktot
c
c ******************************************************************
c conck - conservation check  for specified level
c         mostly a debugging tool
c         this assumes grids don't overlap
c ******************************************************************
c
c
c  grid loop for given level
c
      hx      = hxposs(level)
      hy      = hyposs(level)
      hz      = hzposs(level)
      dt      = possk(level)
      totmass = 0.d0

      mptr = lstart(level)
 20   if (mptr .eq. 0) go to 85
         loc    = node(store1,mptr)
	 locaux = node(storeaux,mptr)
	 nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
	 ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         nz     = node(ndkhi,mptr) - node(ndklo,mptr) + 1
	 mitot  = nx + 2*nghost
	 mjtot  = ny + 2*nghost
         mktot  = nz + 2*nghost
c
	 if (mcapa .eq. 0) then
           do 50 k  = nghost+1, mktot-nghost
           do 50 j  = nghost+1, mjtot-nghost
           do 50 i  = nghost+1, mitot-nghost
              totmass = totmass + alloc(iadd(i,j,k,1)) 
 50           continue
	  else
c          # with capa array:
           do 60 k  = nghost+1, mktot-nghost
           do 60 j  = nghost+1, mjtot-nghost
           do 60 i  = nghost+1, mitot-nghost
              totmass = totmass + alloc(iadd(i,j,k,1))
     &                          * alloc(iaddaux(i,j,k)) 
 60           continue
	  endif
c
       mptr = node(levelptr,mptr)
       go to 20
c
 85    totmass = totmass * hx * hy * hz
       write(outunit,*)" total mass = ", totmass
c
 99   return
      end
