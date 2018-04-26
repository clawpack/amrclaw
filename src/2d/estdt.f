c
c-----------------------------------------------------------------------
c
      subroutine estdt(val,mitot,mjtot,nvar,dx,dy,dt,nghost,aux,
     &                 naux,cfl)
c
c :::::::::::::::::::::::: ESTDT :::::::::::::::::::::::::::;
c  estimate the initial time step for the given values
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

       implicit double precision (a-h, o-z)
       dimension val(nvar,mitot,mjtot)
       dimension aux(naux,mitot,mjtot)
c
c
       return
       end
