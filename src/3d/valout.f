c
c -----------------------------------------------------
c
      subroutine valout (lst, lend, time, nvar, naux)
c
c     # Output the results for a general system of conservation laws
c     # in 3 dimensions
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw3.m
c
c     # set outaux = .true. to also output the aux arrays to fort.a<iframe>

      implicit double precision (a-h,o-z)
      character*10  matname1, matname2, matname3
      logical outaux

      include  "call.i"

      iadd(i,j,k,ivar)   = loc     +     (i-1)
     &                             +     (j-1)*mitot
     &                             +     (k-1)*mitot*mjtot
     &                             +  (ivar-1)*mitot*mjtot*mktot
      iaddaux(i,j,k,ivar) = locaux +     (i-1)
     &                             +     (j-1)*mitot
     &                             +     (k-1)*mitot*mjtot
     &                             +  (ivar-1)*mitot*mjtot*mktot

      outaux = .false.
c
c

c     ### NCAR graphics output

      if (ncarout) then

        call basic (time, lst, lend )
c
        write(pltunit1,100)  nvar
100     format(10h*VALS     ,i10)
c
        level = lst
10      if (level .gt. lend) go to 60
            mptr = lstart(level)
20          if (mptr .eq. 0) go to 50
                nx = node(ndihi,mptr)-node(ndilo,mptr) + 1
                ny = node(ndjhi,mptr)-node(ndjlo,mptr) + 1
                nz = node(ndkhi,mptr)-node(ndklo,mptr) + 1
                mitot = nx + 2*nghost
                mjtot = ny + 2*nghost
                mktot = nz + 2*nghost
                loc = node(store1,mptr)
                call outvar(alloc(loc),mitot,mjtot,mktot,
     &                nvar,mptr,nghost)
                mptr = node(levelptr,mptr)
            go to 20
50          continue
            level = level + 1
        go to 10
c
60    endif


c     ### MATLAB graphics output
c
c

      if (matlabout) then
c        ###  make the file names and open output files
         matname1 = 'fort.qxxxx'
         matname2 = 'fort.txxxx'
         matname3 = 'fort.axxxx'
         matunit1 = 50
         matunit2 = 60
         matunit3 = 70
         nstp     = matlabu
         do 55 ipos = 10, 7, -1
            idigit = mod(nstp,10)
            matname1(ipos:ipos) = char(ichar('0') + idigit)
            matname2(ipos:ipos) = char(ichar('0') + idigit)
            matname3(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue
         open(unit=matunit1,file=matname1,status='unknown',
     .       form='formatted')

         level = lst
         ngrids = 0
c65      if (level .gt. lfine) go to 90
 65      if (level .gt. lend) go to 90
            mptr = lstart(level)
 70         if (mptr .eq. 0) go to 80
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              nz      = node(ndkhi,mptr) - node(ndklo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              mktot   = nz + 2*nghost
              write(matunit1,1001) mptr, level, nx, ny, nz
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my',/,
     &       i5,'                 mz')


              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              zlow = rnode(cornzlo,mptr)
              write(matunit1,1002) xlow,ylow,zlow,hxposs(level),
     &              hyposs(level),hzposs(level)
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    ylow', /,
     &       e18.8,'    zlow', /,
     &       e18.8,'    dx', /,
     &       e18.8,'    dy', /,
     &       e18.8,'    dz',/)

      do 75 k = nghost+1, mktot-nghost
         do 76 j = nghost+1, mjtot-nghost
            do 77 i = nghost+1, mitot-nghost
               do ivar=1,nvar
                  if (dabs(alloc(iadd(i,j,k,ivar))) .lt. 1d-90) then
                     alloc(iadd(i,j,k,ivar)) = 0.d0
                  endif
               enddo
               write(matunit1,109)
     &               (alloc(iadd(i,j,k,ivar)), ivar=1,nvar)
  109          format(10e26.16)
   77       continue
            write(matunit1,*) ' '
   76    continue
         write(matunit1,*) ' '
   75 continue
      write(matunit1,*) ' '

      mptr = node(levelptr, mptr)
      go to 70
   80 level = level + 1
      go to 65

   90 continue

      if (outaux) then
c        # also output aux array to fort.aXXXX
         open(unit=matunit3,file=matname1,status='unknown',
     .       form='formatted')
         level = lst
         ngrids = 0
 165     if (level .gt. lfine) go to 190
            mptr = lstart(level)
 170        if (mptr .eq. 0) go to 180
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              nz      = node(ndkhi,mptr) - node(ndklo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              mktot   = nz + 2*nghost
              write(matunit3,1001) mptr, level, nx, ny, nz

              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              zlow = rnode(cornzlo,mptr)
              write(matunit3,1002) xlow,ylow,zlow,hxposs(level),
     &              hyposs(level),hzposs(level)

      do 175 k = nghost+1, mktot-nghost
         do 176 j = nghost+1, mjtot-nghost
            do 177 i = nghost+1, mitot-nghost
               do ivar=1,naux
                  if (dabs(alloc(iaddaux(i,j,k,ivar))) .lt. 1d-90) then
                     alloc(iaddaux(i,j,k,ivar)) = 0.d0
                  endif
               enddo
               write(matunit3,109)
     &               (alloc(iaddaux(i,j,k,ivar)), ivar=1,naux)
  177       continue
            write(matunit3,*) ' '
  176    continue
         write(matunit3,*) ' '
  175 continue
      write(matunit3,*) ' '

      mptr = node(levelptr, mptr)
      go to 170
  180 level = level + 1
      go to 165
      
  190 continue
      close(unit=matunit3)
      endif  !# end output aux array


      open(unit=matunit2,file=matname2,status='unknown',
     .       form='formatted')
      write(matunit2,1000) time,nvar,ngrids,naux
 1000 format(e18.8,'    time', /,
     &      i5,'                 meqn'/,
     &      i5,'                 ngrids'/,
     &      i5,'                 naux'/,/)
c

      write(6,601) matlabu,time
  601 format('AMRCLAW: Frame ',i4,
     &       ' matlab plot files done at time t = ', d12.5,/)

      matlabu = matlabu + 1

      close(unit=matunit1)
      close(unit=matunit2)
      endif

      return
      end
