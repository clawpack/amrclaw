c
c -----------------------------------------------------------
c
      subroutine findcut(icl,iscr,idim,index,iside,
     1                   ilo,ihi)
c
c ::::::::::::::::::::: FINDCUT ::::::::::::::::::::::::::::;
c   find best place to split the 1D array of flagged points
c   either split at a hole, or use signatures to find
c   zero crossing of laplacian.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension iscr(idim)

c Modified 6/02:
c Include call.i to get def's of horizontal.
c      integer   horizontal
c      parameter(horizontal = 1)

      parameter(ithres = 2)
      parameter(minoff = 2)
c
c  look for holes in horizontal direction
c
       do 10 i = ilo, ihi
          if (iscr(i) .eq. 0) then
             index = i
             iside = horizontal
             return
          endif
 10    continue

c
c  no holes - find 2nd derivative of signatures for best cut.
c  overwrite signature arrays. don't make cuts less than minoff
c  from boundary
c
      ipre = iscr(ilo)
      do 50 i = ilo+1, ihi-1
         icur = iscr(i)
         iscr(i) = iscr(i+1)-2*icur+ipre
         ipre = icur
 50   continue

      locmaxi = 0
      indexi  = 0
      imid    = (ilo + ihi) / 2
      do 60 i = ilo+minoff, ihi-minoff+1
           itemp1 = iscr(i-1)
           itemp2 = iscr(i)
           locdif = iabs(itemp1-itemp2)
           if (itemp1*itemp2.lt.0) then
                if (locdif .gt. locmaxi) then
                 locmaxi = locdif
                 indexi = i
                else if (locdif .eq. locmaxi) then
                    if (iabs(i-imid).lt.iabs(indexi-imid)) indexi = i
                endif
           endif
 60   continue
c
c      ::::: set index for splitting
c
       index  = indexi
       iside  = horizontal
       locmax = locmaxi


c      ::::: if inflection pt. not over the threshold, signal
c      ::::: with index 0 (out of range)
       if (locmax .lt. ithres) index = 0

       return
       end
