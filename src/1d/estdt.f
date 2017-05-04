c
c-----------------------------------------------------------------------
c
      subroutine estdt(val,mitot,nvar,dx,dt,nghost,aux,
     &                 naux,cfl)
c
c :::::::::::::::::::::::: ESTDT :::::::::::::::::::::::::::;
c  estimate the initial time step for the given values
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

       implicit double precision (a-h, o-z)
       dimension val(nvar,mitot)
       dimension aux(naux,mitot)
c
c
       return
       end
