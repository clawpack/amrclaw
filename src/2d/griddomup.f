c
c ----------------------------------------------------
c
      subroutine griddomup(iflags,iflags2,ilo,ihi,jlo,jhi,
     .                     mbuff,lev,ilofine,ihifine,jlofine,jhifine)

      use amr_module
      implicit double precision (a-h, o-z)

      integer*1  iflags (ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)
      integer*1  iflags2(ilofine-mbuff:ihifine+mbuff,
     .                   jlofine-mbuff:jhifine+mbuff)

c
c ::::::::::::::::::::::::::: DOMUP :::::::::::::::::::::
c 
c  domain flags for THIS GRID  are in iflags. copy into iflags2, which is
c  at a different level, dimensions differently
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (dprint) then
         write(outunit,*)" from griddomup: flags (before expansion,",
     .                   " with buff cells)"
         do 5 j=jhi+mbuff,jlo-mbuff,-1
         write(outunit,100)(iflags(i,j),i=ilo-mbuff,ihi+mbuff)
 5       continue
      endif
c
      lratiox = intratx(lev)
      lratioy = intraty(lev)

      do 10 j = jlofine-mbuff,jhifine+mbuff
      do 10 i = ilofine-mbuff,ihifine+mbuff
         iflags2(i,j) = 0
 10   continue

c
c careful with buffer - cant just take coarse grid buffer and refine it
c since have same size buffer on finer grid
      do 20 j = jlo,jhi
      do 20 i = ilo,ihi
          ifine = i * lratiox - 1  ! subtract 1 so can add in next loop
          jfine = j * lratioy - 1
          do 15 mj = 1, lratiox
          do 15 mi = 1, lratioy
            iflags2(ifine+mi,jfine+mj) = iflags(i,j)  
 15       continue
 20       continue
c
c need to be careful due to possibly odd buffer size
c cant just take coarse cell and do all fine fine cells within, as above loop
c  
c    handle left and right buffer zones
      do 25 j = jlofine-mbuff, jhifine+mbuff
      do 23 i = ihifine+1, ihifine+mbuff
c       get coarse indices for fine pt. i,j (in buffer zone)
        ic = i/lratiox
        jc = j/lratioy
        iflags2(i,j) = iflags(ic,jc)
 23   continue
      do 24 i = ilofine-mbuff, ilofine-1
c       get coarse indices for fine pt. i,j (in buffer zone)
        ic = i/lratiox
        jc = j/lratioy
        iflags2(i,j) = iflags(ic,jc)
 24   continue
 25   continue

c    handle top and bottom buffer zones
      do 33 i = ilofine, ihifine
      do 35 j = jlofine-mbuff, jlofine-1
c       get coarse indices for fine pt. i,j (in buffer zone)
        ic = i/lratiox
        jc = j/lratioy
        iflags2(i,j) = iflags(ic,jc)
 35    continue
      do 34 j = jhifine+1, jhifine+mbuff
c       get coarse indices for fine pt. i,j (in buffer zone)
        ic = i/lratiox
        jc = j/lratioy
        iflags2(i,j) = iflags(ic,jc)
 34   continue
 33   continue
c
c  take care of periodicity again or if border of domain touches a 
c  physical boundary then set domain in ghost cell as well
c
C THIS DOESNT WORK - NEED TO DO SOMETHING ELSE FOR PERIODIC SINCE NOT WHOEL
C DOMAIN NOW ONLY A GRID

c      if (xperdom .or. yperdom) then
       if (xperdom .or. .not. xperdom) then
          write(*,*) "IN GRIDDOMUP: NOT YET SET FOR PERIODIC DOMAINS "
          go to 90
      endif

      if (xperdom) then
         do 37 j = 0, jsize+1
           iflags2(0,j)       = iflags2(isize,j)
           iflags2(isize+1,j) = iflags2(1,j)
 37      continue
       else
       do 55 j = 1, jsize
         if (iflags2(1,j) .eq. 1) iflags2(0,j) = 1
         if (iflags2(isize,j) .eq. 1) iflags2(isize+1,j) = 1
 55    continue
      endif
      if (yperdom) then
         do 45 i = 0, isize+1
           iflags2(i,0)       = iflags2(i,jsize)
           iflags2(i,jsize+1) = iflags2(i,1)
 45      continue
       else if (spheredom) then
         do 46 i = 0, isize+1
           iflags2(i,0)       = iflags2(isize+1-i,1)
           iflags2(i,jsize+1) = iflags2(isize+1-i,jsize)
 46      continue

       else
         do 65 i = 1, isize
           if (iflags2(i,1) .eq. 1) iflags2(i,0) = 1
           if (iflags2(i,jsize) .eq. 1) iflags2(i,jsize+1) = 1
 65      continue
      endif

c
c the 4 corners
c
        if (iflags2(0,1)+iflags2(1,0) .eq. 2) iflags2(0,0)=1
        if (iflags2(isize,0)+iflags2(isize+1,1) .eq. 2)
     .          iflags2(isize+1,0)=1
        if (iflags2(isize,jsize+1)+iflags2(isize+1,jsize) .eq. 2)
     .          iflags2(isize+1,jsize+1)=1
        if (iflags2(0,jsize)+iflags2(1,jsize+1) .eq. 2)
     .          iflags2(0,jsize+1)=1


 90   continue
      if (dprint) then
         write(outunit,*)"from griddomup: flags (after ref 1 level up,",
     .                   "with buff cells)"
         do 70 j = jlofine-mbuff,jhifine+mbuff,-1
         write(outunit,100)(iflags2(i,j),i=ilofine-mbuff,ihifine+mbuff)
 100     format(80i1)
 70      continue
      endif

      return
      end
