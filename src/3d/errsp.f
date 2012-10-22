c
c ----------------------------------------------------------
c
      subroutine errsp(ect,sperr,mitot,mjtot,mktot,nvar,mptr,ng,
     &                 eprint, outunit)
c
c **********************
c ****** OBSOLETE ******  No longer used in amrclaw.
c **********************
c
c ::::::::::::::::::::: ERRSP ::::::::::::::::::::::::::::::::::
c estimate spatial only component of the error
c rect   = grid values including ghost cells
c sperr  = user computed fn. (density gradient here)
c           (if sperr(i,j,k) > tolsp cell will be flagged)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit double precision (a-h, o-z)

c      include  "call.i"

      dimension   rect(mitot,mjtot,mktot,nvar)
      dimension  sperr(mitot,mjtot,mktot)

c
c     # This routine is obsolete and no longer used.
c
      write(6,*) '*** The new version of AMRCLAW requires some changes'
      write(6,*) 'to Makefile:'
      write(6,*) '   Remove errsp.f (no longer used)'
      write(6,*) '   Add the following library routines to Makefile:'
      write(6,*) '      $(CLAW)/amrclaw/3d/lib/bufnst.f'
      write(6,*) '      $(CLAW)/amrclaw/3d/lib/spest.f'
      write(6,*) '      $(CLAW)/amrclaw/3d/lib/flag2refine.f'
      write(6,*) '      $(CLAW)/amrclaw/3d/lib/allowflag.f'
      write(6,*) 'See the documentation in flag2refine.f and '
      write(6,*) '    allowflag.f for more information.'


      stop

      end
