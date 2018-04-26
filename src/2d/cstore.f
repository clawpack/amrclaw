c
c ------------------------------------------------------------
c
       subroutine cstore(qc,nrow,ncol,nvar,qc1d,lenbc,naux,auxc,auxc1d)

       implicit double precision (a-h, o-z)

       dimension qc(nvar,nrow,ncol)
       dimension qc1d(nvar,lenbc)
       dimension auxc(naux,nrow,ncol)
       dimension auxc1d(naux,lenbc)
c
!>      store coarse perimeter worth of solution into 1d array.
c      go around fine grid in following order
c                2
c           __________
c        1 |          | 3
c           __________
c               4
c
c  save first interior cell of enlarged grid corresponding to
c  fine grid bordering cell. note that since fine grid is smaller,
c  the cell is one in. coarse (temporary) grid has no ghost cells

c      side 1
       index = 0
       do 10 j = 2, ncol-1
         index = index + 1
         do 5 ivar = 1, nvar
 5         qc1d(ivar,index) = qc(ivar,1,j)
         do 6 iaux = 1, naux
 6         auxc1d(iaux,index) = auxc(iaux,1,j)
 10    continue

c      side 2
       do 20 i = 2, nrow-1
         index = index + 1
         do 15 ivar = 1, nvar
 15        qc1d(ivar,index) = qc(ivar,i,ncol)
         do 16 iaux = 1, naux
 16        auxc1d(iaux,index) = auxc(iaux,i,ncol)
 20    continue

c      side 3
       do 30 j = 2, ncol-1
         index = index + 1
         do 25 ivar = 1, nvar
 25        qc1d(ivar,index) = qc(ivar,nrow,j)
         do 26 iaux = 1, naux
 26        auxc1d(iaux,index) = auxc(iaux,nrow,j)
 30    continue

c      side 4
       do 40 i = 2, nrow-1
         index = index + 1
         do 35 ivar = 1, nvar
 35        qc1d(ivar,index) = qc(ivar,i,1)
         do 36 iaux = 1, naux
 36        auxc1d(iaux,index) = auxc(iaux,i,1)
 40    continue

       return
       end
