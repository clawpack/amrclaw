c
c --------------------------------------------------------------------
c
       subroutine signs(badpts,npts,iscr,idim,ist,iend,
     &                  ilo,ihi)
c
       implicit double precision (a-h,o-z)
       dimension badpts(1,npts)
       dimension iscr(idim)
c
c :::::::::::::::::::::::::::: SIGNS ::::::::::::::::::::::::::::::
c  compute signatures = number of flagged cells in each row/column.
c  also return first and last nonzero row/column, so don't have
c  to waste time over entire region.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
       ilo= 1
       ihi= idim
       do 10 i = 1, idim
 10       iscr(i) = 0
c
c count all flagged points in a given row in one pass through
c the points, i.e. a bin count
c
       do 20 ipt = ist, iend
          iloc  =  badpts(1,ipt)+1.1
          iscr(iloc) = iscr(iloc)+1
 20    continue
c
       do 30 ipt = 1, idim
         if (iscr(ipt) .ne. 0) then
            ilo = ipt
            go to 40
         endif
 30    continue
 40    do 50 ipt = 1, idim
          if (iscr(idim+1-ipt) .ne. 0) then
             ihi = idim+1-ipt
             go to 99
          endif
 50    continue

 99    return
       end
