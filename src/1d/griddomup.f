c
c ----------------------------------------------------
c
      subroutine griddomup(iflags,iflags2,ilo,ihi,
     .                     mbuff,lev,ilofine,ihifine)

      use amr_module
      implicit double precision (a-h, o-z)

      integer*1  iflags (ilo-mbuff:ihi+mbuff)
      integer*1  iflags2(ilofine-mbuff:ihifine+mbuff)

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
         write(outunit,100)(iflags(i),i=ilo-mbuff,ihi+mbuff)
      endif
c
      lratiox = intratx(lev)

      do 10 i = ilofine-mbuff,ihifine+mbuff
         iflags2(i) = 0
 10   continue

c
c careful with buffer - cant just take coarse grid buffer and refine it
c since have same size buffer on finer grid
      do 20 i = ilo,ihi
          ifine = i * lratiox - 1  ! subtract 1 so can add in next loop
          do 15 mi = 1, lratiox
            iset = min(ifine+mi,ihifine+mbuff)  ! might as well include buffer, though done later
            iset = max(iset,ilofine-mbuff)      ! but so expensive
            iflags2(iset) = iflags(i)
 15       continue
 20       continue
c
c need to be careful due to possibly odd buffer size
c cant just take coarse cell and do all fine fine cells within, as above loop
c  
c    handle left and right buffer zones
      do 23 i = ihifine+1, ihifine+mbuff
c       get coarse indices for fine pt. i (in buffer zone)
        ic = i/lratiox
        iflags2(i) = iflags(ic)
 23   continue
      do 24 i = ilofine-mbuff, ilofine-1
c       get coarse indices for fine pt. i (in buffer zone)
        ic = i/lratiox
        iflags2(i) = iflags(ic)
 24   continue
c
c do not need to do something special for periodicity. already taken into account when
c setting enlarged grid with buffer zone at level lbase
c
 90   continue
      if (dprint) then
         write(outunit,*)"from griddomup: flags (after ref 1 level up,",
     .                   "with buff cells)"
         write(outunit,100)(iflags2(i),i=ilofine-mbuff,ihifine+mbuff)
 100     format(80i1)
      endif

      return
      end
