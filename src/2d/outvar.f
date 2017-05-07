c
c -------------------------------------------------------------
c
      subroutine outvar(rect,mitot,mjtot,nvar,mptr,ng)
c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension rect(nvar,mitot,mjtot)

c ::::::::::::::: OUTVAR ::::::::::::::::::::::::::::::::::
c
!>  dump soln for graphics 
c
c  only output max - 1 rows and cols, since with cell centered
c  variables there is one extra cell outside the grid.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      write(pltunit1,100) mptr
 100  format('*SOLN     ',i10,' is the grid - all variables')
c
      do 20 ivar = 1, nvar
         write(pltunit1,101) ((rect(ivar,i,j),i=ng+1,mitot-ng),
     .                                        j=ng+1,mjtot-ng)
 101     format(5e13.6)
 20   continue
c
      return
      end
