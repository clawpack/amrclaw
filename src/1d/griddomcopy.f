c
c ----------------------------------------------------
c
      subroutine griddomcopy(i1, i2, ilo,ihi,mbuff)

      use amr_module
      implicit double precision (a-h, o-z)


      integer*1  i2(ilo-mbuff:ihi+mbuff)
      integer*1  i1(ilo-mbuff:ihi+mbuff)

c
c ::::::::::::::::::::::::::: GRIDDOMCOPY :::::::::::::::::::::
c 
c  griddomain flags need to be in different place
c  (can f90 do this in one statement, without sub call?)
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::


      do 10 i = ilo-mbuff,ihi+mbuff
         i1(i) = i2(i)
 10   continue
c
      return
      end
