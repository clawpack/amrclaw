c
c =======================================================================
      subroutine outval(val,nvar,mitot,mptr,outgrd,naux,aux)
c =======================================================================
c
      use amr_module
      implicit double precision (a-h,o-z)

      dimension  val(nvar,mitot)
      dimension  aux(naux,mitot)
      logical    outgrd


c ::::::::::::::::::::::OUTVAL :::::::::::::::::::::::::::::::
c print solution and aux. variables to output. 
c if cell outside domain, don't print soln. value - nothing
c currently in ghost cells.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      if (.not. outgrd) go to 99
      level = node(nestlevel,mptr)
      hx    = hxposs(level)
      cornx = rnode(cornxlo,mptr) -  nghost*hx
c
      do 20 i=nghost+1,mitot-nghost

          x  = cornx + hx*(dble(i)-.5d0)
          write(outunit,107) x,i,(val(ivar,i),ivar=1,nvar)
 107      format(2hx=,f6.3,3h,i=,i3,' a=',e25.15)
c    *           5(e9.3,1x))
          if (naux.gt.0) write(outunit,108) (aux(iaux,i),iaux=1,naux)
 108      format(1x,'aux = ',7(e9.3,1x))

 20   continue

 99   return
      end
