c     ==================
      subroutine setprob
c     ==================

      implicit double precision (a-h,o-z)
      character*12 fname
      common /cparam/  gamma, gamma1
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)

c
       read(7,*) gamma
       read(7,*) gamma1

      return
      end
