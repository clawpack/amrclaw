c
c :::::::::::::::::::::::::  GRIDDOMSHRINK ::::::::::::::::::::::::::::
c
!> Shrink domain flags one cell for allowable properly nested domain
!! This is needed even for lcheck = lbase. More shrinking needed
!! for finer levels.
!! Flags starts in iflags2, should end in iflags array
!!
!! The output **iflags** has flagged zone one cell smaller (shrinked)
!!than input **iflags2.
!! 
!! The flags and indices, **ilo**, **ihi**, **jlo** and **jhi** are all with
!! respect to level **level** index space
!!
!! **input**:
!! * iflags2
!! **output**:
!! * iflags
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c ----------------------------------------------------
c
      subroutine griddomshrink(iflags2,ilo,ihi,jlo,jhi,mbuff,iflags,
     .                         level)

      use amr_module
      implicit double precision (a-h, o-z)


      integer*1  iflags (ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)
      integer*1  iflags2(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)



      if (dprint) then
        write(outunit,*)" from griddomshrink: on entry, iflags2"
        do 10 j = jhi+mbuff,jlo-mbuff,-1
        write(outunit,100)(iflags2(i,j),i=ilo-mbuff,ihi+mbuff)
 100    format(80i1)
 10     continue
      endif

c     NB this untagging alg. includes corner cells in determining proper 
c     nesting.  not always nec., or always done
      do 41 j = jlo-mbuff+1,jhi+mbuff-1
      do 40 i = ilo-mbuff+1,ihi+mbuff-1
         iflags(i,j) = iflags2(i,j)
         if (iflags2(i  ,j  ) .le. 0 .or.
     1       iflags2(i+1,j  ) .le. 0 .or. iflags2(i-1,j  ) .le. 0 .or. 
     2       iflags2(i+1,j+1) .le. 0 .or. iflags2(i-1,j+1) .le. 0 .or. 
     3       iflags2(i  ,j-1) .le. 0 .or. iflags2(i  ,j+1) .le. 0 .or.
     4       iflags2(i+1,j-1) .le. 0 .or. iflags2(i-1,j-1) .le. 0) then
                 iflags(i,j) = 0
          endif
          iflags(ilo-mbuff,j) = 0   ! set last border to 0 instead of leaving uninitialized
          iflags(ihi+mbuff,j) = 0
 40   continue
 41   continue
      do i = ilo-mbuff,ihi+mbuff   ! finish zeroing out first and last col
         iflags(i,jlo-mbuff) = 0
         iflags(i,jhi+mbuff) = 0
      end do

c  dont need to handle periodicity here.  Setting of initial grid included enough room to shrink 1
c  for proper nesting.  If expand up then will need to add periodic domain flagging

c
c if border of domain touches a physical boundary then set domain in
c ghost cell as well
c
      call setPhysBndryFlags(iflags,ilo,ihi,jlo,jhi,mbuff,level)

 99   if (dprint) then
        write(outunit,*)" from griddomshrink: on exit, iflags"
        do 70 j = jhi+mbuff-1, jlo-mbuff+1, -1
        write(outunit,101)(iflags(i,j),i=ilo-mbuff+1,ihi+mbuff-1)
 101    format(80i1)
 70     continue
      endif

      return
      end
