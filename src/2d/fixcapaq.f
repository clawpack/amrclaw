
c ------------------------------------------------------------------

        subroutine fixcapaq(val,aux,mitot,mjtot,valc,auxc,mic,mjc,
     &                      nvar,naux,levc,setflags)

      use amr_module
      implicit double precision (a-h,o-z)

c
c :::::::::::::::::::::::  FIXCAPAQ ::::::::::::::::::::::::::::::
c  new fine grid solution q was linearly interpolated. but want
c  to conserve kappa*q, not q. calculate the discrepancy
c  in kappa*q using this q, and modify q to account for it.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      dimension   val(nvar,mitot,mjtot), valc(nvar,mic,mjc)
      dimension   aux(naux,mitot,mjtot), auxc(naux,mic,mjc)
      dimension   setflags(mitot,mjtot)

      dcapamax = 0.d0
      lratiox  = intratx(levc)
      lratioy  = intraty(levc)

      do ic = 2, mic-1
      do jc = 2, mjc-1


       do 15 ivar = 1, nvar

       capaqfine = 0.d0

       do ico = 1, lratiox
       ifine = (ic-2)*lratiox + nghost + ico
       do jco = 1, lratioy
         jfine = (jc-2)*lratioy + nghost + jco
         capaqfine = capaqfine + aux(mcapa,ifine,jfine)*
     &                           val(ivar,ifine,jfine)
       end do
       end do

       dcapaq = auxc(mcapa,ic,jc)*valc(ivar,ic,jc)-
     &          capaqfine/(lratiox*lratioy)
       dcapamax = dmax1(dcapamax,dabs(dcapaq))
      
       do ico = 1, lratiox
       ifine = (ic-2)*lratiox + nghost + ico
       do jco = 1, lratioy
         jfine = (jc-2)*lratioy + nghost + jco

         if (setflags(ifine,jfine) .eq. NEEDS_TO_BE_SET) then 
         ! was set by coarsegrid, need to check for adjustment
           val(ivar,ifine,jfine) = val(ivar,ifine,jfine) +
     &                             dcapaq/aux(mcapa,ifine,jfine)
         endif
       end do
       end do

15     continue

       end do
       end do

c      write(6,*)" max discrepancy ", dcapamax

       return
       end
