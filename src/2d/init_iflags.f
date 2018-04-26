c
c --------------------------------------------------------------------------
c
       subroutine init_iflags(iflags,isize,jsize)
c
c      # Need this routine to initialize since iflags is part of double 
c      # precision alloc array but is used as an integer(kind=1) to save storage.

       implicit double precision (a-h,o-z)

       integer(kind=1)  iflags (0:isize+1,0:jsize+1)

         do j = 1, jsize
         do i = 1, isize
           iflags(i,j) = 0
           enddo
         enddo

       return
       end
