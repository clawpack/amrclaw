c
c ------------------------------------------------------------
c
       subroutine cstore(qc,nrow,ncol,nfil,nvar,
     &                   qc1d,lenbc,naux,auxc,auxc1d)

       implicit double precision (a-h, o-z)

       dimension qc(nrow,ncol,nfil,nvar)
       dimension qc1d(lenbc,nvar)
       dimension auxc(nrow,ncol,nfil,naux)
       dimension auxc1d(lenbc,naux)
c
c      store coarse perimeter worth of solution into 1d array.
c      go around fine grid in following order
c                2
c           __________
c        1 |          | 3
c           __________
c               4
c
c  then side 5 and then side 6.
c  In 3D, side 1 is x constant,
c              2 is y constant,
c              3 is x constant,
c              4 is y constant,
c              5 is z constant, (bottom)
c              6 is z constant. (top)
c
c  save first interior cell of enlarged grid corresponding to
c  fine grid bordering cell. note that since fine grid is smaller,
c  the cell is one in. coarse (temporary) grid has no ghost cells

c      side 1
       index = 0
       do 10 k = 2, nfil-1
       do 10 j = 2, ncol-1
	 index = index + 1
	 do 5 ivar = 1, nvar
 5	   qc1d(index,ivar) = qc(1,j,k,ivar)
	 do 6 iaux = 1, naux
 6	   auxc1d(index,iaux) = auxc(1,j,k,iaux)
 10    continue

c      side 2
       do 20 k = 2, nfil-1
       do 20 i = 2, nrow-1
	 index = index + 1
	 do 15 ivar = 1, nvar
 15	   qc1d(index,ivar) = qc(i,ncol,k,ivar)
	 do 16 iaux = 1, naux
 16	   auxc1d(index,iaux) = auxc(i,ncol,k,iaux)
 20    continue

c      side 3
       do 30 k = 2, nfil-1
       do 30 j = 2, ncol-1
	 index = index + 1
	 do 25 ivar = 1, nvar
 25	   qc1d(index,ivar) = qc(nrow,j,k,ivar)
	 do 26 iaux = 1, naux
 26	   auxc1d(index,iaux) = auxc(nrow,j,k,iaux)
 30    continue

c      side 4
       do 40 k = 2, nfil-1
       do 40 i = 2, nrow-1
	 index = index + 1
	 do 35 ivar = 1, nvar
 35	   qc1d(index,ivar) = qc(i,1,k,ivar)
	 do 36 iaux = 1, naux
 36	   auxc1d(index,iaux) = auxc(i,1,k,iaux)
 40    continue

c      side 5
       do 50 j = 2, ncol-1
       do 50 i = 2, nrow-1
         index = index + 1
         do 45 ivar = 1, nvar
 45        qc1d(index,ivar) = qc(i,j,1,ivar)
         do 46 iaux = 1, naux
 46        auxc1d(index,iaux) = auxc(i,j,1,iaux)
 50    continue

c      side 6
       do 60 j = 2, ncol-1
       do 60 i = 2, nrow-1
         index = index + 1
         do 55 ivar = 1, nvar
 55        qc1d(index,ivar) = qc(i,j,nfil,ivar)
         do 56 iaux = 1, naux
 56        auxc1d(index,iaux) = auxc(i,j,nfil,iaux)
 60    continue

       return
       end
