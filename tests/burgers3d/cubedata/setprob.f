c     ==================
      subroutine setprob
c     ==================

      implicit double precision (a-h,o-z)
      common /comrp/ coeff(3)
c
c     # for 3d Burgers' equation
c     #    q_t  +  u*(.5*q^2)_x + v*(.5*q^2)_y + w*(.5*q^2)_z = 0
c     # where u,v,w are a given scalars, stored in the vector coeff

      open(unit=7,file='setprob.data',status='old',form='formatted')

      read(7,*) coeff(1)
      read(7,*) coeff(2)
      read(7,*) coeff(3)

      return
      end
