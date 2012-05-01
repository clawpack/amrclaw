c
c -----------------------------------------------------------
c
       subroutine setNeg(rectflags,ilo,ihi,jlo,jhi,mbuff,ico)

       use amr_module
       implicit double precision (a-h, o-z)

       dimension rectflags(ilo-mbuff:ihi+mbuff, jlo-mbuff:jhi+mbuff)

c ****************************************************************
c  setNeg = set any flagged point in buffer region (exterior to grid proper,
c            this includes border of cells just inside grid)
c           to negative value. If it turns out to be properly nested
c           it will be reset to positive
c ****************************************************************

       ico = 0

c remember to handle diagonals
 
       do j = jlo-mbuff, jhi+mbuff
          do  i = ilo-mbuff, ilo
             if (rectflags(i,j) .gt. 0.) then
               rectflags(i,j) = -rectflags(i,j)    ! if not flagged stays 0
               ico = ico + 1
             endif
          end do
          do  i = ihi, ihi+mbuff
             if (rectflags(i,j) .gt. 0.) then
                rectflags(i,j) = -rectflags(i,j)    ! if not flagged stays 0
                ico = ico + 1
             endif
          end do
       end do

c
c   next do top and bottom; corners already handled
       do i = ilo, ihi
          do j = jlo-mbuff, jlo
             if (rectflags(i,j) .gt. 0.) then
                rectflags(i,j) = -rectflags(i,j)
                ico = ico + 1
             endif
          end do
          do j = jhi, jhi+mbuff
             if (rectflags(i,j) .gt. 0.) then
                rectflags(i,j) = -rectflags(i,j)
                ico = ico + 1
             endif
          end do
       end do

       return
       end
