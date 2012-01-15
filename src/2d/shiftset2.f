c
c ----------------------------------------------------------
c
       subroutine shiftset2(rectflags,ilo,ihi,jlo,jhi,mbuff)

       use amr_module
       implicit double precision (a-h, o-z)
       dimension rectflags(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)
       dimension copyflags(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)


c :::::::::::::::::::::: CSHIFT :::::::::::::::::::::::::::::::
c shift by + or - 1 in either direction to do 
c  bitwise calculus for proper nesting, buffering, etc.
c similar to cshift on CM machine.
c includes periodic buffering as well.
c
c NEWER VERSION: DOES ALL DIRS AT SAME TIME
c rectflags array has been augmented by enough border
c           cells to do buffering in place in the grid
c later will look to see if flagged pts are properly nested
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


       do j = jlo-mbuff, jhi+mbuff
       do i = ilo-mbuff, ihi+mbuff
          copyflags(i,j) = 0
       end do
       end do

c  note: not looking at ghost cells, only real cells should be flagged
c  but can buffer the flags into the ghost zone
       do j = jlo, jhi
       do i = ilo, ihi
          rflag = rectflags(i,j)
          if (rflag .gt. 0)  then
c            cell is flagged, buffer in all dirs by ibuff
c            (note this is not nec same as mbuff)
c            use second array to avoid propagation
             mlo = i - ibuff
             mhi = i + ibuff
             klo = j - ibuff
             khi = j + ibuff
             do k = klo, khi 
             do m = mlo, mhi
                copyflags(m,k) = rflag  ! copy the flag (doesnt distinguish buffer flag from orig flag)
             end do
             end do
          endif

       end do
       end do


c   copy back. need flags in original array, not temp scratch array

       do 60 j = jlo-mbuff, jhi+mbuff
       do 60 i = ilo-mbuff, ihi+mbuff
         rectflags(i,j) = copyflags(i,j)
 60    continue


       return
       end
