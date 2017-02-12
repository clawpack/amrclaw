c
c =======================================================================
      subroutine outval(val,nvar,mitot,mjtot,mktot,mptr,outgrd,naux,aux)
c =======================================================================
c
      use amr_module
      implicit double precision (a-h,o-z)

      dimension  val(nvar,mitot,mjtot,mktot)
      dimension  aux(naux,mitot,mjtot,mktot)
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
      hy    = hyposs(level)
      hz    = hzposs(level)
      cornx = rnode(cornxlo,mptr) - nghost*hx
      corny = rnode(cornylo,mptr) - nghost*hy
      cornz = rnode(cornzlo,mptr) - nghost*hz
c
      do 30 k=nghost+1,mktot-nghost
      do 25 j=nghost+1,mjtot-nghost
      do 20 i=nghost+1,mitot-nghost

         x  = cornx + hx*(dble(i)-.5d0)
         y  = corny + hy*(dble(j)-.5d0)
         z  = cornz + hz*(dble(k)-.5d0)
         write(outunit,107) x,y,z,i,j,k,(val(ivar,i,j,k),ivar=1,nvar)
 107     format(2hx=,f7.3,3h y=,f7.3,3h z=,f7.3,
     *          4h, i=,i3,4h, j=,i3 ,4h, k=,i3 ,
     *          ' a= ',5(e11.4,1x))
         if (naux.gt.0) write(outunit,108) (aux(iaux,i,j,k),iaux=1,naux)
 108     format(1x,'aux = ',7(e10.3,1x))

 20   continue
 25   continue
 30   continue

 99   return
      end
