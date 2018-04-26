c
c --------------------------------------------------------------
c
!> Output a subtree of the grids. All grids from the grid level of 
!! grid **mlev** to the finest level are output.
!! \param[in] mlev root grid for current output.
!! \param[in] outgrd If true, output value on grid
!! \param[in] nvar number of equations for the system
!! \param[in] naux number of auxiliary variables
      subroutine outtre(mlev,outgrd,nvar,naux)
c
      use amr_module
      implicit double precision (a-h,o-z)
      logical  outgrd

c
c ::::::::::::::::::::::: OUTTRE :::::::::::::::::::::::::::::::::::
c
c outtre - output subtree
c input parameters:
c    mlev   - ptr to subtree to output i.e., start at level(mlev)
c    outgrd - if true, output the values on the grid
c tree is output from 'level' to finest level.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      write (outunit,1)
1     format(1x,14hthe subtree is)
c
      level = node(nestlevel, mlev)
10    if (level .gt. lfine) go to 99
          mptr    = lstart(level)
 20       if (mptr .eq. 0) go to 30
              call outmsh(mptr,outgrd,nvar,naux)
              mptr = node(levelptr, mptr)
          go to 20
 30       continue
          write(outunit,2) numgrids(level), level,iregst(level),
     1                     jregst(level),iregend(level),jregend(level)
 2        format(/,i5," grids at level ",i5," go from ",2i9," to",2i9,/)
          level = level + 1
      go to 10
c
 99   return
      end
c
c --------------------------------------------------------------
c
!> Output all grids on the same level as grid **mlev**.
!! \param[in] mlev representative grid of the grid level being output
!! \param[in] outgrd If true, output value on grid
!! \param[in] nvar number of equations for the system
!! \param[in] naux number of auxiliary variables
      subroutine outlev(mlev,outgrd,nvar,naux)
c
      use amr_module
      implicit double precision (a-h,o-z)
      logical  outgrd

c
c ::::::::::::::::::::::: OUTTRE :::::::::::::::::::::::::::::::::::
c
c outtre - output subtree
c input parameters:
c    mlev   - ptr to subtree to output i.e., start at level(mlev)
c    outgrd - if true, output the values on the grid
c tree is output from 'level' to finest level.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      write (outunit,1)
1     format(1x,"the partially built new subtree is:")
c
          mptr = mlev
 20       if (mptr .eq. 0) go to 30
              call outmsh(mptr,outgrd,nvar,naux)
              mptr = node(levelptr, mptr)
          go to 20
 30       continue
c
 99   return
      end
