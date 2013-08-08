c
c --------------------------------------------------------------
c
      subroutine setflags(iflags,isize,jsize,ksize,
     1                    rctold,idim3,mitot,mjtot,mktot,mptr)
c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension rctold(idim3,mitot,mjtot,mktot)
      integer*1 iflags(0:isize+1,0:jsize+1,0:ksize+1)

c :::::::::::::::::::::: SETFLAGS ::::::::::::::::::::::::::::::::::
c transfer flagged arrays into 1 large array of entire domain
c makes buffering, projecting, etc. easier without searching 
c through all kinds of grids
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      ibeg = node(ndilo,mptr) - nghost
      jbeg = node(ndjlo,mptr) - nghost
      kbeg = node(ndklo,mptr) - nghost

      do 10 k = nghost+1, mktot-nghost
      do 10 j = nghost+1, mjtot-nghost
      do 10 i = nghost+1, mitot-nghost
         iflags(ibeg+i,jbeg+j,kbeg+k) =  iflags(ibeg+i,jbeg+j,kbeg+k)
     .                                   + rctold(1,i,j,k)
 10   continue
c
 99   return
      end
