c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  ilr,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for the Euler equations
c     #  with a tracer variable.
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c     # Uses Roe averages and other quantities which were 
c     # computed in rpn2eu and stored in the common block comroe.
c
      dimension     ql(meqn, 1-mbc:maxm+mbc)
      dimension     qr(meqn, 1-mbc:maxm+mbc)
      dimension   asdq(meqn, 1-mbc:maxm+mbc)
      dimension bmasdq(meqn, 1-mbc:maxm+mbc)
      dimension bpasdq(meqn, 1-mbc:maxm+mbc)
c
      common /cparam/  gamma,gamma1
      dimension waveb(5,4),sb(4)
	  
c
      if (ixy.eq.1) then
          mu = 2
          mv = 3
        else
          mu = 3
          mv = 2
        endif
c
         do 20 i = 2-mbc, mx+mbc
C 			Calculate Roe averages         
			rhsqrtl = dsqrt(qr(1,i-1))
            rhsqrtr = dsqrt(ql(1,i))
            pl = gamma1*(qr(4,i-1) - 0.5d0*(qr(2,i-1)**2 +
     &        qr(3,i-1)**2)/qr(1,i-1))
            pr = gamma1*(ql(4,i) - 0.5d0*(ql(2,i)**2 +
     &        ql(3,i)**2)/ql(1,i))
            rhsq2 = rhsqrtl + rhsqrtr
            u = (qr(mu,i-1)/rhsqrtl + ql(mu,i)/rhsqrtr) / rhsq2
            v = (qr(mv,i-1)/rhsqrtl + ql(mv,i)/rhsqrtr) / rhsq2
            enth = (((qr(4,i-1)+pl)/rhsqrtl
     &             + (ql(4,i)+pr)/rhsqrtr)) / rhsq2
            u2v2 = u**2 + v**2
            a2 = gamma1*(enth - .5d0*u2v2)
            a = dsqrt(a2)
            g1a2 = gamma1 / a2
            euv = enth - u2v2
			 
            a3 = g1a2 * (euv*asdq(1,i) 
     &             + u*asdq(mu,i) + v*asdq(mv,i) - asdq(4,i))
            a2 = asdq(mu,i) - u*asdq(1,i)
            a4 = (asdq(mv,i) + (a-v)*asdq(1,i) - a*a3)
     &              / (2.d0*a)
            a1 = asdq(1,i) - a3 - a4
c
            waveb(1,1) = a1
            waveb(mu,1) = a1*u
            waveb(mv,1) = a1*(v-a)
            waveb(4,1) = a1*(enth - v*a)
            waveb(5,1) = 0.d0
            sb(1) = v - a
c
            waveb(1,2) = a3
            waveb(mu,2) = a3*u + a2
            waveb(mv,2) = a3*v
            waveb(4,2) = a3*0.5d0*u2v2 + a2*u
            waveb(5,2) = 0.d0
            sb(2) = v
c
            waveb(1,3) = a4
            waveb(mu,3) = a4*u
            waveb(mv,3) = a4*(v+a)
            waveb(4,3) = a4*(enth+v*a)
            waveb(5,3) = 0.d0
            sb(3) = v + a
c
            waveb(1,4) = 0.d0
            waveb(mu,4) = 0.d0
            waveb(mv,4) = 0.d0
            waveb(4,4) = 0.d0
            waveb(5,4) = asdq(5,i)
            sb(4) = v
c
c           # compute the flux differences bmasdq and bpasdq
c
            do 10 m=1,meqn
               bmasdq(m,i) = 0.d0
               bpasdq(m,i) = 0.d0
               do 10 mw=1,4
                  bmasdq(m,i) = bmasdq(m,i) 
     &                         + dmin1(sb(mw), 0.d0) * waveb(m,mw)
                  bpasdq(m,i) = bpasdq(m,i)
     &                         + dmax1(sb(mw), 0.d0) * waveb(m,mw)
   10             continue
c                 
   20       continue
c
      return
      end
