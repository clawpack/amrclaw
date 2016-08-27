c
c ----------------------------------------------------------
c
       subroutine shiftset(rectflags,ilo,ihi,mbuff)

       use amr_module
       implicit double precision (a-h, o-z)
       dimension rectflags(ilo-mbuff:ihi+mbuff)
       dimension copyflags(ilo-mbuff:ihi+mbuff)


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


       do i = ilo-mbuff, ihi+mbuff
          copyflags(i) = 0
       end do

c  note: not looking at ghost cells, only real cells should be flagged
c  but can buffer the flags into the ghost zone
       do i = ilo, ihi
          rflag = rectflags(i)
          if (rflag .gt. 0)  then
c            cell is flagged, buffer in all dirs by ibuff
c            (note this is not nec same as mbuff)
c            use second array to avoid propagation
             mlo = i - ibuff
             mhi = i + ibuff
             do m = mlo, mhi
                copyflags(m) = rflag  ! copy the flag (doesnt distinguish buffer flag from orig flag)
             end do
          endif

       end do


c   copy back. need flags in original array, not temp scratch array

       do 60 i = ilo-mbuff, ihi+mbuff
         rectflags(i) = copyflags(i)
 60    continue


       return
       end
