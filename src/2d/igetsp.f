c
!> Allocate contiguous space of length nword in main storage array
!! alloc. 
!!
!! Return pointer to this storage
c ----------------------------------------------------------
c
      function igetsp (nwords)
c
      use amr_module
      implicit double precision (a-h,o-z)

c
c ::::::::::::::::::::::::::: IGETSP ::::::::::::::::::::::::::::
c
c  allocate contiguous space of length nword in main storage array
c  alloc. user code (or pointer to the owner of this storage)
c  is  mptr.  lenf = current length of lfree list.
c
c ::::::::::::::::::::::::::: IGETSP ::::::::::::::::::::::::::::
c

!$OMP CRITICAL (MemMgmt)

c  find first fit from free space list
c
 10   continue
      itake = 0
      do 20 i = 1, lenf
         if (lfree(i,2) .lt. nwords) go to 20
         itake = i
         go to 25
 20   continue
      go to 900
c
c  anything left?
c
 25   left = lfree(itake,2) - nwords
      igetsp = lfree(itake,1)
      iendtake = lfree(itake,1) + nwords
      if (lendim .lt. iendtake) lendim = iendtake
c
c  the following code which is ignored for now adds the new
      if (left .le. 0) go to 30
      lfree(itake,2) = left
      lfree(itake,1) = iendtake
      go to 99
c
c  item is totally removed.  move next items in list up one.
c
 30   lenf = lenf - 1
      do 40 i = itake, lenf
      lfree(i,1) = lfree(i+1,1)
 40   lfree(i,2) = lfree(i+1,2)
      go to 99
c
 900  write(outunit,901) nwords
      write(*,901)       nwords
 901  format('  require ',i10,' words - either none left or not big',
     1         '  enough space')
      write(outunit,902) ((lfree(i,j),j=1,2),i=1,lenf)
      write(*,902)       ((lfree(i,j),j=1,2),i=1,lenf)
 902  format(' free list: ',//,2x,50(i10,4x,i10,/,2x))
      
      ! Dynamic memory adjustment
      ! Attempt to allocate new memory
      factor = 2.0d0
      istatus = 1
      old_memsize = memsize
      do while (istatus > 0)
          factor = 0.5d0 * factor
          if (factor < 0.1d0) then
            write(*,*)'Memory allocation failed'
            write(outunit,*)'Memory allocation failed'
            write(*,*)'Current memory size ',memsize
            write(outunit,*)'Current memory size ',memsize
            stop
          endif
c
c         # test to make sure not overflowing integer indexing range.
c         # with integer*4, maximum integer is roughly 2**30, so indexing
c         # doesn't work if array is longer than this!

          max_size = huge(memsize)  ! maximum integer based on type(memsize)
          iprev_max_size = max_size / ceiling(1.d0+factor)
          ! If memsize is bigger than iprev_max_size then 
          ! memsize * ceiling(1.d0+factor) will be too big.
          ! Need to test this way to avoid possible integer overflow in test.
          if (memsize > iprev_max_size ) then
            write(*,*)"Requested memory larger than can be "
            write(*,*)" indexed by default integer size. "
            write(*,*)" Suggest compiling with integer*8 option"
            write(*,*)" Asking for less memory: factor=",factor
            write(outunit,*)"Requested memory larger than can be "
            write(outunit,*)" indexed with default integer size. "
            write(outunit,*)" Suggest compiling with integer*8 option"
            write(outunit,*)" Asking for less memory for now:",factor
            cycle
          else
            new_size = ceiling((1.d0+factor) * memsize)
            iadd_size = ceiling(factor * memsize)
            call resize_storage(new_size,istatus) ! success returns istatus 0
          endif
      enddo
          
      if (lfree(lenf-1,1) + lfree(lenf-1,2) - 1 .eq. old_memsize) then
          ! Merge new block with last free block on list, adjust sentinel to
          ! reflect new memory size
          lfree(lenf-1,2) = iadd_size + lfree(lenf-1,2)
          lfree(lenf,1) = new_size + 2
      else
          ! New free block added to end, make new sentinel 
          lfree(lenf,1) = old_memsize+1
          lfree(lenf,2) = iadd_size
          lfree(lenf+1,1) = new_size+2
          lfree(lenf+1,2) = 0
          lenf = lenf + 1
      endif
      go to 10

 99   lentot = lentot + nwords
      if (lenmax .lt. lentot) lenmax = lentot
      if (sprint) write(outunit,100) nwords, igetsp, lentot, lenmax
 100  format('  allocating ',i8,' words in location ',i8,
     1       ' lentot, lenmax ', 2i10)

!$OMP END CRITICAL (MemMgmt)

      return
      end
