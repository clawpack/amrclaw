c
c -----------------------------------------------------------
c
!> Conservation check for specified level.
!! This is mostly a debugging tool and assumes grids don't overlap
      subroutine conck(level, nvar, naux, time, rest)
c
      use amr_module
      implicit double precision (a-h,o-z)

      logical  rest

c      iadd(i,j,ivar)  = loc + i - 1 + mitot*((ivar-1)*mjtot+j-1) OLD INDEXING
c      iaddaux(i,j) = locaux + i - 1 + mitot*(j-1) +
c     .                        mitot*mjtot*(mcapa-1)

c ## indexing into mcapa assumes cell volume is in mcapa location
      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
      iaddaux(i,j) = locaux + mcapa - 1 + naux*(i-1) +
     .                                    naux*mitot*(j-1)
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
      hy      = hyposs(level)
      dt      = possk(level)
      totmass = 0.d0

      mptr = lstart(level)
 20   if (mptr .eq. 0) go to 85
         loc    = node(store1,mptr)
         locaux = node(storeaux,mptr)
         nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         mitot  = nx + 2*nghost
         mjtot  = ny + 2*nghost
c
         if (mcapa .eq. 0) then
           do j  = nghost+1, mjtot-nghost
           do i  = nghost+1, mitot-nghost
              totmass = totmass + alloc(iadd(1,i,j)) 
           end do
           end do
          else
c          # with capa array:
           do j  = nghost+1, mjtot-nghost
           do i  = nghost+1, mitot-nghost
              totmass = totmass + alloc(iadd(1,i,j))*alloc(iaddaux(i,j)) 
           end do
           end do
          endif
c
       mptr = node(levelptr,mptr)
       go to 20
c
 85    totmass = totmass * hx * hy
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
