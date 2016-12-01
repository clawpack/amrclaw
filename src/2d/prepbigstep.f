c
c-------------------------------------------------------------------------------------
c
       subroutine prepbigstep(nvar,naux,lcheck,mptr,nx,ny,midub,mjdub,
     .                       valbgc,auxbgc,mi2tot,mj2tot)

       use amr_module
       implicit double precision (a-h,o-z)

       double precision valdub(nvar,midub,mjdub)
       double precision auxdub(naux,midub,mjdub)
       double precision valbgc(nvar,mi2tot,mj2tot)
       double precision auxbgc(naux,mi2tot,mj2tot)
       dimension fp(nvar,mi2tot,mj2tot),gp(nvar,mi2tot,mj2tot)
       dimension fm(nvar,mi2tot,mj2tot),gm(nvar,mi2tot,mj2tot)
       
          hx  = hxposs(lcheck)
          hy  = hyposs(lcheck)
          hx2 = 2.d0*hx
          hy2 = 2.d0*hy
          dt  = possk(lcheck)
          dt2 = 2. * dt
          time  = rnode(timemult,mptr)
          tpre  = time - dt

          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          ng2    = 2*nghost
          locold = node(store2,mptr)
          xlow   = rnode(cornxlo,mptr) - nghost*hx2
          ylow   = rnode(cornylo,mptr) - nghost*hy2
c

c         # transfer soln. into grid with twice the ghost cells
          call copysol(valdub,alloc(locold),nvar,mitot,mjtot,
     1              nghost,midub,mjdub,ng2)

c
          if (naux .gt. 0) then
              xl     = rnode(cornxlo, mptr)
              yb     = rnode(cornylo, mptr)
              mx = midub - 4*nghost
              my = mjdub - 4*nghost
              auxdub = NEEDS_TO_BE_SET  ! signal that needs a val
              
              call setaux(2*nghost,mx,my,xl,yb,hx,hy,
     &                    naux,auxdub)
              call auxcoarsen(auxdub,midub,mjdub,
     1                     auxbgc,mi2tot,mj2tot,naux,auxtype)
          endif

c         # fill it - use enlarged (before coarsening) aux arrays
          call bound(tpre,nvar,ng2,valdub,midub,mjdub,mptr,
     1               auxdub,naux)

c         coarsen by 2 in every direction
          call coarsen(valdub,midub,mjdub,
     1                 valbgc,mi2tot,mj2tot,nvar)

          call stepgrid(valbgc,fm,fp,gm,gp,
     1                mi2tot,mj2tot,nghost,
     2                dt2,dtnew2,hx2,hy2,nvar,
     3                xlow,ylow,tpre,mptr,naux,auxbgc)

c         update counts for error estimation work
          evol = evol + (nx/2)*(ny/2)

           return
           end
