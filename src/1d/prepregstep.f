c
c-------------------------------------------------------------------------------------
c
       subroutine prepregstep(nvar,naux,lcheck,mptr,nx,mitot,
     .                        valbig,auxbig)

       use amr_module
       implicit double precision (a-h,o-z)

       dimension fp(nvar,mitot)
       dimension fm(nvar,mitot)

       hx   = hxposs(lcheck)
       dt     = possk(lcheck)
       time   = rnode(timemult,mptr)

       xlow   = rnode(cornxlo,mptr) - nghost*hx

c
       call stepgrid(valbig,fm,fp,
     1               mitot,nghost,
     2               dt,dtnew,hx,nvar,
     3               xlow,time,mptr,naux,auxbig)
c
c     update counts for error estimation step
       evol   = evol + nx

       return
       end
