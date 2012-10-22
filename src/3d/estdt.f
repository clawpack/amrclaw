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

       dimension val(mitot,mjtot,mktot,nvar)
       dimension aux(mitot,mjtot,mktot,naux)
c
c
       return
       end
