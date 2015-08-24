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

      use amr_module
      implicit double precision (a-h,o-z)
      character*10  fname1, fname2, fname3, fname4, fname5

      logical outaux
      integer output_aux_num


      iadd(ivar,i,j,k)   = loc     +    (ivar-1)
     &                             +    (i-1)*nvar
     &                             +    (j-1)*nvar*mitot
     &                             +    (k-1)*nvar*mitot*mjtot
      iaddaux(iaux,i,j,k) = locaux +    (iaux-1)
     &                             +    (i-1)*naux
     &                             +    (j-1)*naux*mitot
     &                             +    (k-1)*naux*mitot*mjtot

      errmax = 0.d0

      output_aux_num = 0
      do i = 1, naux
        output_aux_num = output_aux_num + output_aux_components(i)
      end do

c     # Currently outputs all aux components if any are requested!
      outaux = ((output_aux_num > 0) .and. 
     .         ((.not. output_aux_onlyonce) .or. (time==t0)))

c     open(unit=77,file='fort.b',status='unknown',access='stream')


c     ### Python graphics output
c
c        ###  make the file names and open output files
      fname1 = 'fort.qxxxx'
      fname2 = 'fort.txxxx'
      fname3 = 'fort.axxxx'
      fname4 = 'fort.bxxxx'
      matunit1 = 50
      matunit2 = 60
      matunit3 = 70
      matunit4 = 71
      nstp     = matlabu
      do 55 ipos = 10, 7, -1
         idigit = mod(nstp,10)
         fname1(ipos:ipos) = char(ichar('0') + idigit)
         fname2(ipos:ipos) = char(ichar('0') + idigit)
         fname3(ipos:ipos) = char(ichar('0') + idigit)
         fname4(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
 55   continue

      open(unit=matunit1,file=fname1,status='unknown',
     .       form='formatted')

      if (output_format == 3) then
c         # binary output          
          open(unit=matunit4,file=fname4,status='unknown',
     &            access='stream')
          endif

      level = lst
      ngrids = 0
 65   if (level .gt. lend) go to 90
         mptr = lstart(level)
 70      if (mptr .eq. 0) go to 80
           ngrids  = ngrids + 1
           nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
           ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
           nz      = node(ndkhi,mptr) - node(ndklo,mptr) + 1
           loc     = node(store1, mptr)
           locaux  = node(storeaux,mptr)
           mitot   = nx + 2*nghost
           mjtot   = ny + 2*nghost
           mktot   = nz + 2*nghost
           call lookat(alloc(loc),nx,ny,nz,nghost)
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
     &           hyposs(level),hzposs(level)
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    ylow', /,
     &       e18.8,'    zlow', /,
     &       e18.8,'    dx', /,
     &       e18.8,'    dy', /,
     &       e18.8,'    dz',/)



         if (output_format == 1) then
            do 75 k = nghost+1, mktot-nghost
            z = zlow + (k-.5d0-nghost)*hzposs(level)
            do 76 j = nghost+1, mjtot-nghost
            y = ylow + (j-.5d0-nghost)*hyposs(level)
            do 77 i = nghost+1, mitot-nghost
            x = xlow + (i-.5d0-nghost)*hxposs(level)
            do ivar=1,nvar
               if (dabs(alloc(iadd(ivar,i,j,k))) < 1d-90) then
                  alloc(iadd(ivar,i,j,k)) = 0.d0
               endif
            enddo
            exact = qtrue(x,y,z,time)
            approx = alloc(iadd(1,i,j,k))
            err =  abs(approx-exact)  ! error plot instead
            errmax = max(err, errmax)
            write(matunit1,109)
!    &            (alloc(iadd(ivar,i,j,k)), ivar=1,nvar)
     &            err ! error plot instead
  109       format(50e26.16)
   77       continue
            write(matunit1,*) ' '
   76       continue
            write(matunit1,*) ' '
   75       continue
            write(matunit1,*) ' '
         endif

         if (output_format == 3) then
c            # binary output          
             i1 = iadd(1,1,1,1)
             i2 = iadd(nvar,mitot,mjtot,mktot)
c            # NOTE: we are writing out ghost cell data also, unlike ascii
             write(matunit4) alloc(i1:i2)
         endif


         mptr = node(levelptr, mptr)
         go to 70
   80 level = level + 1
      go to 65

   90 continue
      write(*,*) "max error found at time", time," is ",errmax

c       -------------------
c       # output aux arrays
c       -------------------

      if (outaux) then
c        #  output aux array to fort.aXXXX

         level = lst
         ngrids = 0
 165     if (level .gt. lend) go to 190
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




     	       if (output_format == 1) then
                open(unit=matunit3,file=fname3,status='unknown',
     .               form='formatted')
                write(matunit3,1001) mptr,level,nx,ny,nz
                xlow = rnode(cornxlo,mptr)
                ylow = rnode(cornylo,mptr)
                zlow = rnode(cornzlo,mptr)
                write(matunit3,1002) xlow,ylow,zlow,hxposs(level),
     &                               hyposs(level),hzposs(level)

                do 175 k = nghost+1, mktot-nghost
                do 176 j = nghost+1, mjtot-nghost
                do 177 i = nghost+1, mitot-nghost
                do ivar=1,naux
                   if (dabs(alloc(iaddaux(ivar,i,j,k))) .lt. 1d-90) then
                            alloc(iaddaux(ivar,i,j,k)) = 0.d0
                   endif
                enddo
                write(matunit3,109)
     &               (alloc(iaddaux(ivar,i,j,k)), ivar=1,naux)
 177            continue
                write(matunit3,*) ' '
 176            continue
                write(matunit3,*) ' '
 175            continue
                write(matunit3,*) ' '
             endif

             if (output_format == 3) then
c            # binary output          
                open(unit=matunit3,file=fname3,status='unknown',
     &               access='stream')
                i1 = iaddaux(1,1,1,1)
                i2 = iaddaux(naux,mitot,mjtot,mktot)
c               # NOTE: we are writing out ghost cell data also, unlike ascii
                write(matunit3) alloc(i1:i2)
             endif

            mptr = node(levelptr, mptr)
            go to 170
  180 level = level + 1
      go to 165
      
  190 continue

      close(unit=matunit3)
      endif  !# end outputting aux array


c     --------------
c     # fort.t file:
c     --------------

      open(unit=matunit2,file=fname2,status='unknown',
     .       form='formatted')

      ndim = 3
c     # NOTE: we need to print out nghost too in order to strip
c     #       ghost cells from q when reading in pyclaw.io.binary
      write(matunit2,1000) time,nvar,ngrids,naux,ndim,nghost
 1000 format(e18.8,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 naux'/,
     &       i5,'                 ndim'/,
     &       i5,'                 nghost'/,/)
c

      write(6,601) matlabu,time
  601 format('AMRCLAW: Frame ',i4,
     &       ' output files done at time t = ', d12.6,/)

      matlabu = matlabu + 1

      close(unit=matunit1)
      close(unit=matunit2)
      if (output_format == 3) then
          close(unit=matunit4)
          endif

      return
      end

c ---------------

      subroutine lookat(val,mx,my,mz,nghost)
      implicit double precision (a-h,o-z)
      
      dimension val(1-nghost:mx+nghost,1-nghost:my+nghost,
     .          1-nghost:mz+nghost)
      i=10
      J=10
      K=10

      !write(*,*) val(i,j,k)
      return
      end
