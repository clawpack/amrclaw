c     =====================================================
      subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,q,maux,aux)
c     =====================================================
c
c
      implicit double precision (a-h,o-z)
      dimension q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
      common /cqinit/ A1,beta1,x1,y1, A2,beta2,x2,y2

      do i=1,mx
        xc = xlower + (i-0.5d0)*dx
        do j=1,my
            yc = ylower + (j-0.5d0)*dy
            call mapc2p(xc,yc,xp,yp)
            q(1,i,j) = A1*exp(-beta1*((xp - x1)**2 + (yp - y1)**2)) 
     &                  +A2*exp(-beta2*((xp - x2)**2 + (yp - y2)**2))
        enddo
      enddo

      return
      end
