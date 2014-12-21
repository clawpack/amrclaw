c
c ---------------------------------------------------------
c
      subroutine smartbis(badpts,npts,cutoff,numptc,nclust,
     1                    lbase,intcorn,idim,jdim,kdim)
c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension     badpts(3,npts),intcorn(nsize,maxcl)
      dimension     iscr(idim), jscr(jdim), kscr(kdim)
      integer       nclust, numptc(maxcl)
      parameter     (usemin=.4)
c
c :::::::::::::::::::::::::::: SMARTBIS :::::::::::::::::::::::::;
c smart bisect rectangles until cutoff reached for each.
c replaced old bisection routine that cut all grids in half.
c now look for good place to do the cut, based on holes or signatures.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
c     ## initially all points in 1 cluster
      nclust      = 1
      numptc(1)   = npts

      if (gprint) write(outunit,100) nclust
 100  format(' starting smart bisection with ',i5,' clusters')
c
      icl         = 1
      ist         = 1
      iend        = numptc(icl)
c
 10   call moment(intcorn(1,icl),badpts(1,ist),numptc(icl),usenew)
      if (gprint) write(outunit,101) icl,numptc(icl),usenew
 101  format('   testing cluster ',i4,' with ',i9,' pts. use ',e12.4)
c
      if (usenew .lt. cutoff) go to 20
c
c  this cluster ok - on to next
c
      if (.not. gprint) go to 15
         write(outunit,102) icl,numptc(icl),usenew
 102     format('     accepting smart bisected cluster',i4,' with ',i9,
     1          ' pts. use = ',e10.3)
 15   icl   = icl + 1
      if (icl .gt. nclust) go to 200
      ist   = iend + 1
      iend  = ist + numptc(icl) - 1
      go to 10
c
c  smart bisect rectangle (and its cluster) in best location
c

 20   if (nclust .lt. maxcl) go to 25
          write(outunit,900) maxcl
          write(*      ,900) maxcl
 900      format('  too many clusters:  > ',i5)
          stop
 25   continue
c
c smart bisection computes signatures, finds best cut and splits there
c
      call signs(badpts,npts,iscr,jscr,kscr,idim,jdim,kdim,
     &           ist,iend,ilo,ihi,jlo,jhi,klo,khi)
      call findcut(icl,iscr,jscr,kscr,idim,jdim,kdim,index,iside,
     &             ilo,ihi,jlo,jhi,klo,khi)
      if (index .eq. 0) then
         icl = icl + 1
         if (icl .gt. nclust) go to 200
         ist = iend + 1
         iend = ist + numptc(icl) - 1
         go to 10
      endif
c
      if     (iside .eq. kplane) then
         fmid = (index-.5)
         idir = 3
      elseif (iside .eq. jplane) then
         fmid = (index-.5)
         idir = 2
      else
         fmid = (index-.5)
         idir = 1
      endif
c
      itop = ist - 1
      ibot = iend + 1
      i    = ist

c 2/28/02 : Here, the logic starts differing from the 2d version.
   50 if  (badpts(idir,i) .lt. fmid) then
c         point in bottom half. switch with a bottom point
c         which does not belong in the bottom half.
   60     ibot           = ibot - 1
          if  (badpts(idir,ibot) .lt. fmid) then
c             this point should stay in the bottom.
c             do not switch it with point i.
              if  (ibot .eq. i) then
c                 ibot and i point to the same badpt.
c                 the subdivision is complete.
                  go to 80
              else
c                 decrement the bottom pointer
c                 Remember: in order to get to this part of the code,
c                 i must point to a badpt which belongs in the bottom
c                 part of the list. Hence, no array bounds test is
c                 necessary here because ibot cannot decrease below i.
                  go to 60
              end if
          else
c             do the switch in each coordinate direction
              do 61 ndim=1,3
              temp              = badpts(ndim,ibot)
              badpts(ndim,ibot) = badpts(ndim,i)
              badpts(ndim,i)    = temp
   61         continue
              if (itop+1 .lt. ibot) go to 50
          end if
      else
c         point in top half. let it stay, increment counter.
c         itop always points to a badpt in the top half (ie, .ge. fmid)
          itop = itop + 1
          if  (itop+1 .ge. ibot) then
              go to 80
          else
              i = i + 1
              go to 50
          end if
      end if
c
c done smartbisecting icl-th clusters. adjust counts, repeat bisect stage .
c
 80   numptc(icl) = itop - ist + 1
      ibump       = icl + 1
c
c  bump down remaining clusters to make room for the new half of one.
c
      if (ibump .gt. nclust) go to 120
      do 90 ico         = ibump, nclust
      nmove             = nclust - ico + ibump
 90   numptc(nmove + 1) = numptc(nmove)

 120  numptc(ibump)     = iend - ibot + 1
      nclust            = nclust + 1
      iend              = itop
c
c     other half of the cluster has been inserted into cluster list.
c     icl remains the same - need to redo it.
      go to 10
c
c  done: there are nclust clusters.
c
 200  continue
c
      return
      end
