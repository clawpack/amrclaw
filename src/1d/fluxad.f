c
c -------------------------------------------------------
c
      subroutine fluxad(xfluxm,xfluxp,
     1                  svdflx,mptr,mitot,
     2                   nvar,lenbc,lratiox,ng,dtf,dx)
c

      use amr_module
      implicit double precision (a-h,o-z)


c :::::::::::::::::::: FLUXAD ::::::::::::::::::::::::::::::::::
c  save fine grid fluxes  at the border of the grid, for fixing
c  up the adjacent coarse cells. at each edge of the grid, only
c  save the plus or minus fluxes, as necessary. For ex., on
c  left edge of fine grid, it is the minus xfluxes that modify the
c  coarse cell.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      dimension xfluxm(nvar,mitot)
      dimension xfluxp(nvar,mitot)
      dimension svdflx(nvar,lenbc)
 
c ::::: left side saved first
      lind = 0

      lind = lind + 1
      do 110 ivar = 1, nvar
            svdflx(ivar,lind) = svdflx(ivar,lind) +
     1                             xfluxm(ivar,ng+1)*dtf
c           write(dbugunit,900)lind,xfluxm(ivar,1),
c     .                           xfluxp(ivar,1)
 900        format(' lind ', i4,' m & p ',2e15.7,' svd ',e15.7)
110    continue
 
c ::::: right side
       lind = lind + 1
       do 310 ivar = 1, nvar
           svdflx(ivar,lind) = svdflx(ivar,lind) +
     1                    xfluxp(ivar,mitot-ng+1)*dtf
c          write(dbugunit,900)lind,xfluxm(ivar,mitot-ng+1
c                                 ),xfluxp(ivar,mitot-ng+1)
310    continue
 
      return
      end
