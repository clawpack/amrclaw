subroutine setprob
      implicit double precision (a-h,o-z)
      character*25 fname
      common /cparam/ ubar,vbar,wbar
!
!     # Set the velocity for scalar advection
!     # These values are passed to setaux in a common block
!

      iunit = 7
      fname = 'setprob.data'
      call opendatafile(iunit, fname)
                
      read(7,*) ubar
      read(7,*) vbar
      read(7,*) wbar

end subroutine setprob

