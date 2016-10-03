
c ------------------------------------------------------------------

        subroutine fixcapaq(val,aux,mitot,valc,auxc,mic,
     &                      nvar,naux,levc,setflags)

      use amr_module
      implicit double precision (a-h,o-z)

c
c :::::::::::::::::::::::  FIXCAPAQ ::::::::::::::::::::::::::::::
c  new fine grid solution q was linearly interpolated. but want
c  to conserve kappa*q, not q. calculate the discrepancy
c  in kappa*q using this q, and modify q to account for it.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      dimension   val(nvar,mitot), valc(nvar,mic)
      dimension   aux(naux,mitot), auxc(naux,mic)
      dimension   setflags(mitot)

      dcapamax = 0.d0
      lratiox  = intratx(levc)

      do 10 ic = 2, mic-1

       do 15 ivar = 1, nvar

       capaqfine = 0.d0

       do 20 ico = 1, lratiox
         ifine = (ic-2)*lratiox + nghost + ico
         capaqfine = capaqfine + aux(mcapa,ifine)*
     &                           val(ivar,ifine)
20     continue

       dcapaq = auxc(mcapa,ic)*valc(ivar,ic)-
     &          capaqfine/(lratiox)
       dcapamax = dmax1(dcapamax,dabs(dcapaq))
      
       do 30 ico = 1, lratiox
         ifine = (ic-2)*lratiox + nghost + ico

         if (setflags(ifine) .eq. NEEDS_TO_BE_SET) then
         ! was set by coarsegrid, need to check for adjustment
           val(ivar,ifine) = val(ivar,ifine) +
     &                             dcapaq/aux(mcapa,ifine)
         endif
30     continue

15     continue

10     continue

c      write(6,*)" max discrepancy ", dcapamax

       return
       end
