
       subroutine resize_bndryList()
       
          use amr_module
          implicit double precision (a-h,o-z)
          integer, allocatable, target, dimension(:,:)  :: new_bndList

          ! get new bndry space
          new_bndListSize = 1.5*bndListSize
          print *,"Expanding size of boundary list from ",bndListSize,
     .            " to ",new_bndListSize
          allocate(new_bndList(new_bndListSize,2),STAT=istatus)
          if (istatus > 0) then
             write(*,*)" could not get new bndry list space"
             stop
          endif
          new_bndList(1:bndListSize,1:2) = bndList
          call move_alloc(new_bndList,bndList)

          ! thread new bndry list space
          do i = bndListSize+1, new_bndListSize
              bndList(i,nextfree) = i+1
          end do

          bndList(new_bndListSize,nextfree) = null
          ndfree_bnd = bndListSize+1
          bndListSize = new_bndListSize

          end
