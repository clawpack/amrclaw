c
c ----------------------------------------------------
c
      subroutine griddomcopy(i1, i2, ilo,ihi,jlo,jhi,mbuff)

      use amr_module
      implicit double precision (a-h, o-z)


      integer*1  i2(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)
      integer*1  i1(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)

c
c ::::::::::::::::::::::::::: GRIDDOMCOPY :::::::::::::::::::::
c 
c  griddomain flags need to be in different place
c  (can f90 do this in one statement, without sub call?)
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::


      do 10 j = jlo-mbuff,jhi+mbuff
      do 10 i = ilo-mbuff,ihi+mbuff
         i1(i,j) = i2(i,j)
 10   continue
c
      return
      end
