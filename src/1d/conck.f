c
c -----------------------------------------------------------
c
      subroutine conck(level, nvar, naux, time, rest)
c
      use amr_module
      implicit double precision (a-h,o-z)

      logical  rest

c      iadd(i,j,ivar)  = loc + i - 1 + mitot*((ivar-1)*mjtot+j-1) OLD INDEXING
c      iaddaux(i,j) = locaux + i - 1 + mitot*(j-1) +
c     .                        mitot*mjtot*(mcapa-1)

c ## indexing into mcapa assumes cell volume is in mcapa location
      iadd(ivar,i)  = loc + ivar - 1 + nvar*(i-1)
      iaddaux(i) = locaux + mcapa - 1 + naux*(i-1)
c
c
c ******************************************************************
c conck - conservation check  for specified level
c         mostly a debugging tool
c         this assumes grids don't overlap
c
c ******************************************************************
c
c
c  grid loop for given level
c
      hx      = hxposs(level)
      dt      = possk(level)
      totmass = 0.d0

      mptr = lstart(level)
 20   if (mptr .eq. 0) go to 85
         loc    = node(store1,mptr)
         locaux = node(storeaux,mptr)
         nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
         mitot  = nx + 2*nghost
c
         if (mcapa .eq. 0) then
           do 50 i  = nghost+1, mitot-nghost
              totmass = totmass + alloc(iadd(1,i))
 50           continue
          else
c          # with capa array:
           do 60 i  = nghost+1, mitot-nghost
              totmass = totmass + alloc(iadd(1,i))*alloc(iaddaux(i))
 60           continue
          endif
c
       mptr = node(levelptr,mptr)
       go to 20
c
 85    totmass = totmass * hx
       if (time.eq. t0 .and. (level.eq.1) .and. .not. rest) then
           tmass0 = totmass
           write(6,*) 'Total mass at initial time: ',tmass0
           endif
       write(outunit,777) time, totmass, totmass-tmass0
 777   format('time t = ',e12.5,',  total mass = ',e22.15, '  diff = ',
     &         e11.4)
c
 99   return
      end
