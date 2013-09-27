c
c -----------------------------------------------------------
c
      subroutine findcut(icl,iscr,jscr,kscr,idim,jdim,kdim,index,iside,
     1                   ilo,ihi,jlo,jhi,klo,khi)
c
c ::::::::::::::::::::: FINDCUT ::::::::::::::::::::::::::::;
c   find best place to split the 3D array of flagged points
c   either split at a hole, or use signatures to find
c   zero crossing of laplacian.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      use amr_module
      implicit double precision (a-h,o-z)


      parameter(ithres = 2)
      parameter(minoff = 2)

      dimension iscr(idim), jscr(jdim), kscr(kdim)
      dimension locval(4,3)

c
c  look for holes first in horizontal then vertical then transverse direction
c
       do 10 i = ilo, ihi
          if (iscr(i) .eq. 0) then
             index = i
             iside = iplane
             return
          endif
 10    continue

       do 20 j = jlo, jhi
          if (jscr(j) .eq. 0) then
              index = j
              iside = jplane
              return
          endif
 20    continue

       do 30 k = klo, khi
          if (kscr(k) .eq. 0) then
              index = k
              iside = kplane
              return
          endif
 30    continue

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


       jpre   = jscr(jlo)
       do 130 j = jlo+1, jhi-1
           jcur = jscr(j)
           jscr(j) = jscr(j+1)-2*jcur+jpre
           jpre = jcur
 130   continue

       locmaxj = 0
       indexj  = 0
       jmid    = (jlo + jhi) / 2
       do 160 j = jlo+minoff, jhi-minoff+1
           jtemp1 = jscr(j-1)
           jtemp2 = jscr(j)
           locdif = iabs(jtemp1-jtemp2)
           if (jtemp1*jtemp2.lt.0) then
               if (locdif .gt. locmaxj) then
                  locmaxj = locdif
                  indexj = j
                else if (locdif .eq. locmaxj) then
                       if (iabs(j-jmid).lt.iabs(indexj-jmid)) indexj = j
               endif
           endif
 160   continue

       kpre   = kscr(klo)
       do 170 k = klo+1, khi-1
           kcur = kscr(k)
           kscr(k) = kscr(k+1)-2*kcur+kpre
           kpre = kcur
 170   continue
 
       locmaxk = 0
       indexk  = 0
       kmid    = (klo + khi) / 2
       do 180 k = klo+minoff, khi-minoff+1
           ktemp1 = kscr(k-1)
           ktemp2 = kscr(k)
           locdif = iabs(ktemp1-ktemp2)
           if (ktemp1*ktemp2.lt.0) then
               if (locdif .gt. locmaxk) then
                  locmaxk = locdif
                  indexk = k
                else if (locdif .eq. locmaxk) then
                       if (iabs(k-kmid).lt.iabs(indexk-kmid)) indexk = k
               endif
           endif
 180   continue
c
c      ::::: choose max dif for splitting
c      Select the side with the largest locmax.
c      If two are equal and greater than the third,
c      or if all are equal, among the greatest select
c      the one whose locmax is closest to its mid.
c      If two are equal and greater than the third,
c      or if all are equal, among the greatest select
c      the first one arbitrarily.
c
       locval(1,1) = locmaxi
       locval(2,1) = iplane
       locval(3,1) = indexi
       locval(4,1) = imid
       locval(1,2) = locmaxj
       locval(2,2) = jplane
       locval(3,2) = indexj
       locval(4,2) = jmid
       locval(1,3) = locmaxk
       locval(2,3) = kplane
       locval(3,3) = indexk
       locval(4,3) = kmid
 
c      sort in descending order on the value of locmax [locval(1, )]
       do 195 ione=1,2
       do 195 itwo=ione+1,3
          if (locval(1,ione) .lt. locval(1,itwo)) then
             do 190 ii=1,4
                loctmp          = locval(ii,ione)
                locval(ii,ione) = locval(ii,itwo)
                locval(ii,itwo) = loctmp
 190         continue
          end if
 195   continue
       if     (locval(1,1) .gt. locval(1,2)) then
c        locval(1,1) is uniquely greatest
         icount = 1
       elseif (locval(1,2) .gt. locval(1,3)) then
c        locval(1,1) = locval(1,2) but > locval(1,3)
         icount = 2
       else
c        locval(1,1) = locval(1,2) = locval(1,3)
         icount = 3
       end if
 
       if (icount .ne. 1) then
c         sort the first icount entries in ascending order on the abs.val. of
c         the difference between index and mid [locval(3, ) and locval(4, )]
          do 205 ione=1,icount-1
          do 205 itwo=ione+1,icount
             if ((iabs(locval(3,ione)-locval(4,ione))) .gt.
     &           (iabs(locval(3,itwo)-locval(4,itwo)))) then
                do 200 ii=1,4
                   loctmp          = locval(ii,ione)
                   locval(ii,ione) = locval(ii,itwo)
                   locval(ii,itwo) = loctmp
 200            continue
             end if
 205      continue
       end if
 
c      select the first entry
       index  = locval(3,1)
       iside  = locval(2,1)
       locmax = locval(1,1)
 
c      ::::: if inflection pt. not over the threshold, signal
c      ::::: with index 0 (out of range)
       if (locmax .lt. ithres) index = 0
 
       return
       end
