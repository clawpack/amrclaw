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
           qc1d(ivar,index) = qc(ivar,1,j)
 5       end do
         do 6 iaux = 1, naux
           auxc1d(iaux,index) = auxc(iaux,1,j)
 6       end do
 10    continue

c      side 2
       do 20 i = 2, nrow-1
         index = index + 1
         do 15 ivar = 1, nvar
           qc1d(ivar,index) = qc(ivar,i,ncol)
 15      end do
         do 16 iaux = 1, naux
           auxc1d(iaux,index) = auxc(iaux,i,ncol)
 16      end do
 20    continue

c      side 3
       do 30 j = 2, ncol-1
         index = index + 1
         do 25 ivar = 1, nvar
           qc1d(ivar,index) = qc(ivar,nrow,j)
 25      end do
         do 26 iaux = 1, naux
           auxc1d(iaux,index) = auxc(iaux,nrow,j)
 26      end do
 30    continue

c      side 4
       do 40 i = 2, nrow-1
         index = index + 1
         do 35 ivar = 1, nvar
           qc1d(ivar,index) = qc(ivar,i,1)
 35        end do 
         do 36 iaux = 1, naux
           auxc1d(iaux,index) = auxc(iaux,i,1)
 36      end do
 40    continue

       return
       end
