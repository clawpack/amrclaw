c
        subroutine fixcapaq(val ,aux ,mitot,mjtot,mktot,
     &                      valc,auxc,mic  ,mjc  ,mkc  ,
     &                      nvar,naux,levc)

      implicit double precision (a-h,o-z)

      include "call.i"
c
c :::::::::::::::::::::::  FIXCAPAQ ::::::::::::::::::::::::::::::
c  new fine grid solution q was linearly interpolated. but want
c  to conserve kappa*q, not q. calculate the discrepancy
c  in kappa*q using this q, and modify q to account for it.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      dimension   val(mitot,mjtot,mktot,nvar), valc(mic,mjc,mkc,nvar)
      dimension   aux(mitot,mjtot,mktot,naux), auxc(mic,mjc,mkc,naux)

      dcapamax = 0.d0
      lratiox  = intratx(levc)
      lratioy  = intraty(levc)
      lratioz  = intratz(levc)

      do 10 kc = 2, mkc-1
      do 10 jc = 2, mjc-1
      do 10 ic = 2, mic-1


       do 15 ivar = 1, nvar

       capaqfine = 0.d0

       do 20 kco = 1, lratioz
         kfine = (kc-2)*lratioz + nghost + kco
         do 20 jco = 1, lratioy
           jfine = (jc-2)*lratioy + nghost + jco
           do 20 ico = 1, lratiox
             ifine = (ic-2)*lratiox + nghost + ico
	     capaqfine = capaqfine + aux(ifine,jfine,kfine,mcapa)
     &                             * val(ifine,jfine,kfine,ivar)
20     continue

       dcapaq = auxc(ic,jc,kc,mcapa)*valc(ic,jc,kc,ivar)-
     &          capaqfine/(lratiox*lratioy*lratioz)
       dcapamax = dmax1(dcapamax,dabs(dcapaq))
      
       do 30 kco = 1, lratioz
         kfine = (kc-2)*lratioz + nghost + kco
         do 30 jco = 1, lratioy
           jfine = (jc-2)*lratioy + nghost + jco
           do 30 ico = 1, lratiox
             ifine = (ic-2)*lratiox + nghost + ico
             val(ifine,jfine,kfine,ivar)  = val(ifine,jfine,kfine,ivar)
     &                             + dcapaq/aux(ifine,jfine,kfine,mcapa)
30     continue

15     continue

10     continue

c      write(6,*)" max discrepancy ", dcapamax

       return
       end
