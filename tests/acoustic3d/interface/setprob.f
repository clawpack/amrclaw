c     ==================
      subroutine setprob
c     ==================

      implicit double precision (a-h,o-z)
      character*12 fname
      common /comaux/ z1,c1,z2,c2
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)

c
       read(7,*) z1
       read(7,*) c1
       read(7,*) z2
       read(7,*) c2

      return
      end
