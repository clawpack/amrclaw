c
c ----------------------------------------------------
c
      subroutine domprep(domflags,lbase,ibase,jbase,kbase)

      use amr_module
      implicit double precision (a-h, o-z)


      integer*1 domflags(0:ibase+1,0:jbase+1,0:kbase+1)

c
c ::::::::::::::::::::::::::: PREPDOM :::::::::::::::::::::
c 
c  prepare 3 dimensional array of domain for proper nesting
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::


      do 10 k = 0, kbase+1
      do 10 j = 0, jbase+1
      do 10 i = 0, ibase+1
        domflags(i,j,k) = 0
 10   continue

      mptr = lstart(lbase)
 15   continue
      do 20 k = node(ndklo,mptr) + 1, node(ndkhi,mptr) + 1
      do 20 j = node(ndjlo,mptr) + 1, node(ndjhi,mptr) + 1
      do 20 i = node(ndilo,mptr) + 1, node(ndihi,mptr) + 1
        domflags(i,j,k) = 1
 20   continue
      mptr = node(levelptr, mptr)
      if (mptr .ne. 0) go to 15

c
c take care of periodic domains or if border of domain touches a 
c  physical boundary then set domain in ghost cell as well
c
      if (xperdom) then
         do 25 k = 0, kbase+1
         do 25 j = 0, jbase+1
          domflags(      0,j,k) = domflags(ibase,j,k)
          domflags(ibase+1,j,k) = domflags(    1,j,k)
 25      continue
       else
         do 65 k = 1, kbase
         do 65 j = 1, jbase
           if (domflags(    1,j,k) .eq. 1) domflags(      0,j,k) = 1
           if (domflags(ibase,j,k) .eq. 1) domflags(ibase+1,j,k) = 1
 65      continue
      endif
      if (yperdom) then
         do 35 k = 0, kbase+1
         do 35 i = 0, ibase+1
           domflags(i,      0,k) = domflags(i,jbase,k)
           domflags(i,jbase+1,k) = domflags(i,    1,k)
 35      continue
       else
         do 55 k = 1, kbase
         do 55 i = 1, ibase
            if (domflags(i,    1,k) .eq. 1) domflags(i,      0,k) = 1
            if (domflags(i,jbase,k) .eq. 1) domflags(i,jbase+1,k) = 1
 55      continue
      endif
      if (zperdom) then
         do 45 j = 0, jbase+1
         do 45 i = 0, ibase+1
           domflags(i,j,      0) = domflags(i,j,kbase)
           domflags(i,j,kbase+1) = domflags(i,j,    1)
 45      continue
      else
         do 75 j = 1, jbase
         do 75 i = 1, ibase
           if (domflags(i,j,    1) .eq. 1) domflags(i,j,      0) = 1
           if (domflags(i,j,kbase) .eq. 1) domflags(i,j,kbase+1) = 1
 75       continue
      endif
c
c     corner stuff
c
      do k = 1, kbase
        if (domflags(0,1,k)+domflags(1,0,k).eq.2)  domflags(0,0,k) = 1
        if (domflags(0,jbase,k)+domflags(1,jbase+1,k).eq.2)
     .                                       domflags(0,jbase+1,k) = 1
        if (domflags(ibase,0,k)+domflags(ibase+1,1,k).eq.2)
     .                                       domflags(ibase+1,0,k) = 1
        if (domflags(ibase,jbase+1,k)+domflags(ibase+1,jbase,k).eq.2)
     .                                 domflags(ibase+1,jbase+1,k) = 1
      end do
      do i = 1, ibase
        if (domflags(i,0,1)+domflags(i,1,0).eq.2)  domflags(i,0,0) = 1
        if (domflags(i,jbase,0)+domflags(i,jbase+1,1).eq.2)
     .                                       domflags(i,jbase+1,0) = 1
        if (domflags(i,0,kbase)+domflags(i,1,kbase+1).eq.2)
     .                                       domflags(i,0,kbase+1) = 1
        if (domflags(i,jbase+1,kbase)+domflags(i,jbase,kbase+1).eq.2)
     .                                 domflags(i,jbase+1,kbase+1) = 1
      end do
      do j = 1, jbase
        if (domflags(0,j,1)+domflags(1,j,0).eq.2)  domflags(0,j,0) = 1
        if (domflags(ibase,j,0)+domflags(ibase+1,j,1).eq.2)
     .                                       domflags(ibase+1,j,0) = 1
        if (domflags(0,j,kbase)+domflags(1,j,kbase+1).eq.2)
     .                                       domflags(0,j,kbase+1) = 1
        if (domflags(ibase+1,j,kbase)+domflags(ibase,j,kbase+1).eq.2)
     .                                 domflags(ibase+1,j,kbase+1) = 1
      end do
c
c     The neigbors of a corner are the three locations each having an
c     address which differs from that of the corner by 1 in just one of
c     the three indexes. If a corner index is 0, the corresponding neighbor
c     index may be 0 or 1. If a corner index is base+1, the corresponding
c     neighbor index may be base+1 or base (where base => [i,j,k]base).
c
      if (( domflags(      0,      0,      1)
     &     +domflags(      0,      1,      0)
     &     +domflags(      1,      0,      0)) .eq. 3)
     .      domflags(      0,      0,      0)   =   1
 
      if (( domflags(ibase+1,      0,      1)
     &     +domflags(ibase+1,      1,      0)
     &     +domflags(ibase  ,      0,      0)) .eq. 3)
     .      domflags(ibase+1,      0,      0)   =   1
 
      if (( domflags(      0,jbase+1,      1)
     &     +domflags(      0,jbase  ,      0)
     &     +domflags(      1,jbase+1,      0)) .eq. 3)
     .      domflags(      0,jbase+1,      0)   =   1
 
      if (( domflags(ibase+1,jbase+1,      1)
     &     +domflags(ibase+1,jbase  ,      0)
     &     +domflags(ibase  ,jbase+1,      0)) .eq. 3)
     .      domflags(ibase+1,jbase+1,      0)   =   1
 
      if (( domflags(      0,      0,kbase  )
     &     +domflags(      0,      1,kbase+1)
     &     +domflags(      1,      0,kbase+1)) .eq. 3)
     .      domflags(      0,      0,kbase+1)   =   1
 
      if (( domflags(ibase+1,      0,kbase  )
     &     +domflags(ibase+1,      1,kbase+1)
     &     +domflags(ibase  ,      0,kbase+1)) .eq. 3)
     .      domflags(ibase+1,      0,kbase+1)   =   1
 
      if (( domflags(      0,jbase+1,kbase  )
     &     +domflags(      0,jbase  ,kbase+1)
     &     +domflags(      1,jbase+1,kbase+1)) .eq. 3)
     .      domflags(      0,jbase+1,kbase+1)   =   1
 
      if (( domflags(ibase+1,jbase+1,kbase  )
     &     +domflags(ibase+1,jbase  ,kbase+1)
     &     +domflags(ibase  ,jbase+1,kbase+1)) .eq. 3)
     .      domflags(ibase+1,jbase+1,kbase+1)   =   1

      if (dprint) then
         write(outunit,*)" from domprep: domflags at level  ", lbase
         do 40 kk = 1, kbase
         k = kbase + 1 - kk
         write(outunit,*) 'plane k = ',k
         do 40 jj = 1, jbase
         j = jbase + 1 - jj
         write(outunit,100)(domflags(i,j,k),i=1,ibase)
 100     format(80i1)
 40      continue
      endif

      return
      end
