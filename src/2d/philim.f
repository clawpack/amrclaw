c
c
c     =====================================================
      real(CLAW_REAL) function philim(a,b,meth)
c     =====================================================
      implicit real(CLAW_REAL) (a-h,o-z)
c
c     # Compute a limiter based on wave strengths a and b.
c     # meth determines what limiter is used.
c     # a is assumed to be nonzero.
c
c     # NOTE: This routine is obsolete.  Instead of using limiter.f,
c     # which calls philim.f for every wave, it is more efficient to 
c     # use inlinelimiter.f, which eliminates all these function calls
c     # to philim.  If you wish to change the limiter function and are
c     # using inlinelimiter.f, the formulas must be changed in that routine.
c
      r = b/a
      go to (10,20,30,40,50) meth

c
   10 continue
c     --------
c     # minmod
c     --------
      philim = max(0.d0, min(1.d0, r))
      return
c
   20 continue
c     ----------
c     # superbee
c     ----------
      philim = max(0.d0, min(1.d0, 2.d0*r), min(2.d0, r))
      return
c
   30 continue
c     ----------
c     # van Leer
c     ----------
      philim = (r + abs(r)) / (1.d0 + abs(r))
      return
c
   40 continue
c     ------------------------------
c     # monotinized centered 
c     ------------------------------
      c = (1.d0 + r)/2.d0
      philim = max(0.d0, min(c, 2.d0, 2.d0*r))
      return
c
   50 continue
c     ------------------------------
c     # Beam-Warming
c     ------------------------------
      philim = r

      return
      end
