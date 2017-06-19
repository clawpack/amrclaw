c
c --------------------------------------------------------------
c
      subroutine setflags(iflags,isize,jsize,
     1                    rctold,idim3,mitot,mjtot,mptr)
c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension rctold(idim3,mitot,mjtot)
      integer(kind=1) iflags(0:isize+1,0:jsize+1)

c :::::::::::::::::::::: SETFLAGS ::::::::::::::::::::::::::::::::::
!> transfer flagged arrays into 1 large array of entire domain
!! makes buffering, projecting, etc. easier without searching 
!! through all kinds of grids
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      ibeg = node(ndilo,mptr) - nghost
      jbeg = node(ndjlo,mptr) - nghost

      do 10 j = nghost+1, mjtot-nghost
      do 10 i = nghost+1, mitot-nghost
        iflags(ibeg+i,jbeg+j) = iflags(ibeg+i,jbeg+j) + rctold(1,i,j)
 10   continue
c
 99   return
      end
