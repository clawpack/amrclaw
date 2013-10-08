c
c -----------------------------------------------------------
c
      subroutine conck(level, nvar, naux, time, rest)
c
      use amr_module
      implicit double precision (a-h,o-z)

      logical  rest


!--      iadd(i,j,k,ivar)   = loc    +    (i-1)
!--     &                            +    (j-1)*mitot
!--     &                            +    (k-1)*mitot*mjtot
!--     &                            + (ivar-1)*mitot*mjtot*mktot
!--      iaddaux(i,j,k)     = locaux +    (i-1)
!--     &                            +    (j-1)*mitot
!--     &                            +    (k-1)*mitot*mjtot
!--     &                            +(mcapa-1)*mitot*mjtot*mktot
       iadd(ivar,i,j,k) = loc  + (ivar-1)
     &                         + (i-1)*nvar
     &                         + (j-1)*nvar*mitot
     &                         + (k-1)*nvar*mitot*mjtot
       iaddaux(i,j,k) = locaux + (mcapa-1)
     &                         + (i-1)*naux
     &                         + (j-1)*naux*mitot
     &                         + (k-1)*naux*mitot*mjtot
         
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
              totmass = totmass + alloc(iadd(1,i,j,k)) 
c             write(66,999) i,j,k,alloc(iadd(1,i,j,k))
c999          format(3i5,e30.20)
 50           continue
         else
c          # with capa array:
           do 60 k  = nghost+1, mktot-nghost
           do 60 j  = nghost+1, mjtot-nghost
           do 60 i  = nghost+1, mitot-nghost
              totmass = totmass + alloc(iadd(1,i,j,k))
     &                          * alloc(iaddaux(i,j,k)) 
 60           continue
         endif
c
       mptr = node(levelptr,mptr)
       go to 20
c
 85    totmass = totmass * hx * hy * hz
       if (time.eq.t0 .and. (level.eq.1) .and. .not. rest) then
           tmass0 = totmass
           write(6,*) 'Total mass at initial time: ',tmass0
           endif
       write(outunit,777) time, totmass, totmass-tmass0
 777   format('time t = ',e12.5,',  total mass = ',e22.15, '  diff = ',
     &         e11.4)
c
 99   return
      end
