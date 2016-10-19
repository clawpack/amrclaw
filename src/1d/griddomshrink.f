c
c ----------------------------------------------------
c
      subroutine griddomshrink(iflags2,ilo,ihi,mbuff,iflags,
     .                         level)

      use amr_module
      implicit double precision (a-h, o-z)


      integer*1  iflags (ilo-mbuff:ihi+mbuff)
      integer*1  iflags2(ilo-mbuff:ihi+mbuff)
c
c :::::::::::::::::::::::::  GRIDDOMSHRINK ::::::::::::::::::::::::::::
c
c  shrink domain flags one cell for allowable properly nested domain
c  This is needed even for lcheck = lbase. More shrinking needed
c  for finer levels.
c  flags starts in iflags2, should end in iflags array
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (dprint) then
        write(outunit,*)" from griddomshrink: on entry, iflags2"
        write(outunit,100)(iflags2(i),i=ilo-mbuff,ihi+mbuff)
 100    format(80i1)
      endif

c     NB this untagging alg. includes corner cells in determining proper 
c     nesting.  not always nec., or always done
      do 40 i = ilo-mbuff+1,ihi+mbuff-1
         iflags(i) = iflags2(i)
         if (iflags2(i  ) .le. 0 .or.
     1       iflags2(i+1) .le. 0 .or. iflags2(i-1) .le. 0) then
                 iflags(i) = 0
          endif
          iflags(ilo-mbuff) = 0   ! set last border to 0 instead of leaving uninitialized
          iflags(ihi+mbuff) = 0
 40   continue

c  dont need to handle periodicity here.  Setting of initial grid included enough room to shrink 1
c  for proper nesting.  If expand up then will need to add periodic domain flagging

c
c if border of domain touches a physical boundary then set domain in
c ghost cell as well
c
      call setPhysBndryFlags(iflags,ilo,ihi,mbuff,level)

 99   if (dprint) then
        write(outunit,*)" from griddomshrink: on exit, iflags"
        write(outunit,101)(iflags(i),i=ilo-mbuff+1,ihi+mbuff-1)
 101    format(80i1)
      endif

      return
      end
