c 
c -----------------------------------------------------
c This is the updated valout subroutine for amrclaw.  It creates a netcdf file for each time step.
c -----------------------------------------------------

c
c -----------------------------------------------------
c     Routine to write netcdf files in the classic format
!        #jj-2011.03.29
!        # Each file written by the fortran code has 
!        # Dimensions:
!        #           timedimension : UNLIMITED
!        #           meqn          : The number of equations
!        #           dimx_<gridno> : X dimension for grid number <gridno>
!        #           dimy_<gridno> : Y dimension for grid number <gridno>
!        # Variables:
!        #           timedimension : Stores the time of the frame
!        #           ngrids        : Number of grids in this frame
!        #           naux          : Number of Auxilary Variables
!        #           ndim          : Number of Dimensions in the frame
!        #           grid_<gridno> : A grid of (dimx,dimy,meqn)
!        # Attributes:
!        # (grid_<no>) gridno      : The number of this grid <grid_no>
!        #           level         : The AMR level
!        #           dim_names     : a list of dimensions [dimx,dimy]
!        #           dim<x,y>.low  : The lowest dimension value 
!        #           dim<x,y>.d    : The distance between grid points 
c -----------------------------------------------------
      subroutine valout (lst, lend, time, nvar, naux)
c
      implicit double precision (a-h,o-z)
      character*10  matname1, matname2, matname3

c     # Output the results for a general system of conservation laws
c     # in 2 dimensions
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw2.m or Python tools
c
c     # set outaux = .true. to also output the aux arrays to fort.a<iframe>

      include  "call.i"

      include 'netcdf.inc'
      real(kind=8) time
      integer ncid,rcode
      integer timeid,tVarID,meqnID,ngridsVarID,nauxVarID,ndimVarID
      integer dimxid,dimyid,xlowid,ylowid,dxid,dyid
      integer gridid
      integer ntimes
      character*2 gridstr
      character*40 dim_names
      REAL, ALLOCATABLE  ::grid(:,:,:)
      real dx,dy,xlow,ylow      
      
      logical outaux

      iadd(i,j,ivar) = loc + i - 1 + mitot*((ivar-1)*mjtot+j-1)
      iaddaux(i,j,ivar) = locaux + i - 1 + mitot*((ivar-1)*mjtot+j-1)
c
      outaux = .false.

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
                mitot = nx + 2*nghost
                mjtot = ny + 2*nghost
                loc = node(store1,mptr)
                call outvar(alloc(loc),mitot,mjtot,nvar,mptr,nghost)
                mptr = node(levelptr,mptr)
            go to 20
50          continue
            level = level + 1
        go to 10
c
      endif
60    continue


c     ### netcdf/Python graphics output--jj-3/30/2011
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
         !!!!Define netcdf file
         rcode=NF_CREATE(matname1//'.nc',NF_NOCLOBBER,ncid)
         if(rcode.ne.NF_NOERR) print *,'ERROR OPENING NETCDF FILE'
         rcode=NF_DEF_DIM(ncid,'timedimension',NF_UNLIMITED,timeid)
         rcode=NF_DEF_VAR(ncid,'timedimension',NF_FLOAT,1,timeid,tVarID)
         rcode=NF_DEF_DIM(ncid,'meqn',nvar,meqnid)
         rcode=NF_DEF_VAR(ncid,'ngrids',NF_INT,0,0,ngridsVarID)
         rcode=NF_DEF_VAR(ncid,'naux',NF_INT,0,0,nauxVarID)
         rcode=NF_DEF_VAR(ncid,'ndim',NF_INT,0,0,ndimVarID)
         rcode=NF_ENDDEF(ncid)

         level = lst
         ngrids = 0
c65      if (level .gt. lfine) go to 90
 65      if (level .gt. lend) go to 90
            mptr = lstart(level)
 70         if (mptr .eq. 0) go to 80
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              rcode=NF_REDEF(ncid)
              if(rcode.ne.NF_NOERR) print *,'ERROR  REDEFINE MODE'
              write(gridstr,67) mptr
67            format(I2.2)              
              rcode=NF_DEF_DIM(ncid,'dimx_'//trim(gridstr),nx,dimxid)
              if(rcode.ne.NF_NOERR) print *,'ERROR  DEFINE DIMS'
              rcode=NF_DEF_DIM(ncid,'dimy_'//trim(gridstr),ny,dimyid)
              if(rcode.ne.NF_NOERR) print *,'ERROR  DEFINE DIMS'
              

              rcode=NF_DEF_Var(ncid,'grid_'//trim(gridstr),NF_FLOAT,4,
     &              (/dimxid,dimyid,meqnid,timeid/),gridid)
               if(rcode.ne.NF_NOERR) print *,'ERROR  DEFINE VAR'
               
              rcode=NF_PUT_ATT_INT(ncid,gridid,'gridno',NF_INT,1,
     &              mptr)
     
              rcode=NF_PUT_ATT_INT(ncid,gridid,'level',NF_INT,1,level)
              
              dim_names="['dimx','dimy']"
              rcode=NF_PUT_ATT_TEXT(ncid,gridid,'dim_names',
     &         LEN_TRIM(dim_names),TRIM(dim_names))
     
              rcode=NF_PUT_ATT_REAL(ncid,gridid,'dimx.lower',NF_DOUBLE,
     &              1,xlow)     
              rcode=NF_PUT_ATT_REAL(ncid,gridid,'dimy.lower',NF_DOUBLE,
     &              1,ylow)
     
              dx=hxposs(level)
              dy=hyposs(level)
              rcode=NF_PUT_ATT_REAL(ncid,gridid,'dimx.d',NF_FLOAT,1,
     &          dx)
              rcode=NF_PUT_ATT_REAL(ncid,gridid,'dimy.d',NF_FLOAT,1,
     &          dy) 
     
              rcode=NF_ENDDEF(ncid)
              if(rcode.ne.NF_NOERR) print *,'ERROR  REDEFINE MODE'
              allocate(grid(nx,ny,nvar))
              grid=0.d0            

         do j = nghost+1, mjtot-nghost
            do i = nghost+1, mitot-nghost
               do ivar=1,nvar
                  if (dabs(alloc(iadd(i,j,ivar))) .lt. 1d-90) then
                     alloc(iadd(i,j,ivar)) = 0.d0
                  endif
                  grid(i-nghost,j-nghost,ivar)=alloc(iadd(i,j,ivar))
               enddo
            
!              surface = alloc(iadd(i,j,1)) + alloc(iaddaux(i,j,1))
!              grid(i-nghost,j-nghost,4)=surface
            enddo
         enddo

         rcode=NF_PUT_VARA_REAL(ncid,gridid,(/1,1,1,1/), 
     &      (/nx,ny,nvar,1/),grid)
         !rcode=NF_SYNC(NCID)

         deallocate(grid)
            mptr = node(levelptr, mptr)
            go to 70
 80      level = level + 1
         go to 65

 90     continue

      rcode=NF_PUT_VAR_DOUBLE(ncid,tVarID,time)
      if(rcode.ne.NF_NOERR) print *,'ERROR  Write Time'
      rcode=NF_PUT_VAR_INT(ncid,ngridsVarID,int(ngrids))
      if(rcode.ne.NF_NOERR) print *,'ERROR  Write GridNo'
      rcode=NF_PUT_VAR_INT(ncid,nauxVarID,3)
      rcode=NF_PUT_VAR_INT(ncid,ndimVarID,2)
      rcode=NF_CLOSE(ncid)

 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')
 1003 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx')
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    ylow', /,
     &       e18.8,'    dx', /,
     &       e18.8,'    dy',/)
 1004 format(e18.8,'    xlow', /,
     &       e18.8,'    dx', /)
  109       format(4e26.16)
       
        if (outaux) then
c        # output aux array to fort.aXXXX
         open(unit=matunit3,file=matname3,status='unknown',
     .       form='formatted')
         level = lst
         ngrids = 0
 165     if (level .gt. lfine) go to 190
            mptr = lstart(level)
 170        if (mptr .eq. 0) go to 180
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              if (ny.gt.1) then
                  write(matunit3,1001) mptr, level, nx, ny
                else
c                 # output in 1d format if ny=1:
                  write(matunit3,1003) mptr, level, nx
                endif
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              if (ny.gt.1) then
                  write(matunit3,1002)
     &              xlow,ylow,hxposs(level),hyposs(level)
                else
                  write(matunit3,1004)
     &              xlow,hxposs(level)
                endif

         do j = nghost+1, mjtot-nghost
            do i = nghost+1, mitot-nghost
               do ivar=1,naux
                  if (dabs(alloc(iaddaux(i,j,ivar))) .lt. 1d-90) then
                     alloc(iaddaux(i,j,ivar)) = 0.d0
                  endif
               enddo
               write(matunit3,109) (alloc(iaddaux(i,j,ivar)), 
     &                              ivar=1,naux)
            enddo
            write(matunit3,*) ' '
         enddo

            mptr = node(levelptr, mptr)
            go to 170
 180     level = level + 1
         go to 165

 190    continue
        close(unit=matunit3)
        endif !# end outputting aux array

      open(unit=matunit2,file=matname2,status='unknown',
     .       form='formatted')
      if (ny.gt.1) then 
          ndim = 2
        else
c         # special case where 2d AMR is used for a 1d problem
c         # and we want to use 1d plotting routines
          ndim = 1
        endif

      write(matunit2,1000) time,nvar,ngrids,naux,ndim
 1000 format(e18.8,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 naux'/,
     &       i5,'                 ndim'/,/)
c

      write(6,601) matlabu,time
  601 format('AMRCLAW: Frame ',i4,
     &       ' output files done at time t = ', d12.6,/)

      matlabu = matlabu + 1

      close(unit=matunit2)
      endif

      return
      end
