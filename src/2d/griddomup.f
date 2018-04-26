c
!> **iflags** described flagged cells in a rectangular region
!! described by **ilo**, **ihi**, **jlo**, **jhi** in level **lev**
!! index space
!! This subroutine projects **iflags** to **iflag**, which has flagging
!! information in a rectangular region described by **ilofine**,
!! **ihifine**, **jlofine**, **jhifine** in level **lev**+1 index space
!!
!!
!! **input**:
!! * iflags
!! **output**:
!! * iflags2
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
          do 15 mj = 1, lratioy
          do 15 mi = 1, lratiox
            iset = min(ifine+mi,ihifine+mbuff)  ! might as well include buffer, though done later
            jset = min(jfine+mj,jhifine+mbuff)  ! needed since grids dont align over many levels
            iset = max(iset,ilofine-mbuff)      ! but so expensive
            jset = max(jset,jlofine-mbuff)
            iflags2(iset,jset) = iflags(i,j)  
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
c do not need to do something special for periodicity. already taken into account when
c setting enlarged grid woth buffer zone at level lbase
c
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
