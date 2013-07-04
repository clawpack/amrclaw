c
c-------------------------------------------------------------------------------------
c
       subroutine prepregstep(nvar,naux,lcheck,mptr,nx,ny,mitot,mjtot,
     .                        valbig,auxbig)

       use amr_module
       implicit double precision (a-h,o-z)

       dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
       dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)


       hx   = hxposs(lcheck)
       hy   = hyposs(lcheck)
       dt     = possk(lcheck)
       time   = rnode(timemult,mptr)

       xlow   = rnode(cornxlo,mptr) - nghost*hx
       ylow   = rnode(cornylo,mptr) - nghost*hy

c
       call stepgrid(valbig,fm,fp,gm,gp,
     1               mitot,mjtot,nghost,
     2               dt,dtnew,hx,hy,nvar,
     3               xlow,ylow,time,mptr,naux,auxbig)
c
c     update counts for error estimation step
       evol   = evol + nx * ny

       return
       end
