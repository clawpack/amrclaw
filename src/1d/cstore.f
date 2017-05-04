c
c ------------------------------------------------------------
c
       subroutine cstore(qc,nrow,nvar,qc1d,lenbc,naux,auxc,auxc1d)

       implicit double precision (a-h, o-z)

       dimension qc(nvar,nrow)
       dimension qc1d(nvar,lenbc)
       dimension auxc(naux,nrow)
       dimension auxc1d(naux,lenbc)
c
c      store coarse perimeter worth of solution into 1d array.
c      go around fine grid in following order
c        1 |          | 2
c
c  save first interior cell of enlarged grid corresponding to
c  fine grid bordering cell. note that since fine grid is smaller,
c  the cell is one in. coarse (temporary) grid has no ghost cells

c      side 1
       index = 0

        index = index + 1
        do 5 ivar = 1, nvar
 5         qc1d(ivar,index) = qc(ivar,1)
        do 6 iaux = 1, naux
 6         auxc1d(iaux,index) = auxc(iaux,1)

c      side 2
        index = index + 1
        do 25 ivar = 1, nvar
 25        qc1d(ivar,index) = qc(ivar,nrow)
        do 26 iaux = 1, naux
 26        auxc1d(iaux,index) = auxc(iaux,nrow)

       return
       end
