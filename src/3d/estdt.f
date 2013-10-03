c
c-----------------------------------------------------------------------
c
       subroutine estdt(val,mitot,mjtot,mktot,nvar,dx,dy,dz,dt,nghost,
     &                  aux,naux, cfl)
c
c :::::::::::::::::::::::: ESTDT :::::::::::::::::::::::::::;
c  estimate the initial time step for the given values
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

       implicit double precision (a-h, o-z)

       dimension val(nvar,mitot,mjtot,mktot)
       dimension aux(naux,mitot,mjtot,mktot)
c
c
       return
       end
