c
!> Smart bisect rectangles until cutoff reached for each.
!! replaced old bisection routine that cut all grids in half.
!! now look for good place to do the cut, based on holes or signatures.
!!
!! Reference: Berger, M.J., Rigoutsos, I., 1991. An Algorithm for Point Clustering and Grid Generation. IEEE Trans. Syst. Man Cybern. 21, 1278–1286
!!
!! \param[in,out] badpts x,y centered coordinate of all flagged cells.
!! At output, cells for each cluster are stored consecutively in this
!! array.
!! \param[in] npts number of flagged cells
!! \param[in] cutoff required minimum efficiency of each cluster grid
!! \param[out] numptc number of cells in each cluster (grid)
!! \param[out] nclust number of clusters
!! \param[in] lbase base level of current refinement
!! \param[out] intcorn Some information of each generated cluster (grid)
!! \param[in] idim region size (in number of cells) for current AMR level in *i* direction
!! \param[out] idim region size (in number of cells) for current AMR level in *i* direction
!!
!!
!!
c
c ---------------------------------------------------------
c
      subroutine smartbis(badpts,npts,cutoff,numptc,nclust,
     1                    lbase,intcorn,idim,jdim)
c     1                    lbase,intcorn,iscr,jscr,idim,jdim)
c
c  iscr, jscr now stackbased, no need for use of alloc
c
      use amr_module
      implicit double precision (a-h,o-z)

      dimension     badpts(2,npts),intcorn(nsize,maxcl)
c
c  iscr, jscr now stackbased, no need for use of alloc
      dimension     iscr(idim), jscr(jdim)
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
      icl         = 1 ! cluster iterator 
      ist         = 1 ! pt iterator in badpts
      iend        = numptc(icl) ! pt iterator in badpts
c
 10   call moment(intcorn(1,icl),badpts(1,ist),numptc(icl),usenew)
      if (gprint) write(outunit,101) icl,numptc(icl),usenew
 101  format('   testing cluster ',i8,' with ',i9,' pts. use ',e12.4)
c
      if (usenew .lt. cutoff) go to 20
c
c  this cluster ok - on to next
c
      if (.not. gprint) go to 15
         write(outunit,102) icl,numptc(icl),usenew
 102     format('     accepting smart bisected cluster',i4,' with ',i5,
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
      call signs(badpts,npts,iscr,jscr,idim,jdim,
     &           ist,iend,ilo,ihi,jlo,jhi)
      call findcut(icl,iscr,jscr,idim,jdim,index,iside,
     &             ilo,ihi,jlo,jhi)
      if (index .eq. 0) then

c         if (usenew .gt. usemin) then
c           icl = icl + 1
c           if (icl .gt. nclust) go to 200
c           ist = iend + 1
c           iend = ist + numptc(icl) - 1
c           go to 10
c         else
c c         bisect in long direction
c           if (ihi-ilo .gt. jhi-jlo) then
c              iside = horizontal
c              index = (ilo + ihi)/2
c           else
c              iside = vertical
c              index = (jlo + jhi)/2
c           endif
c          endif

c 2/28/02 : 3d version uses this branch only;  no 'if' statement.
         icl = icl + 1
         if (icl .gt. nclust) go to 200
         ist = iend + 1
         iend = ist + numptc(icl) - 1
         go to 10
      endif
c
      if (iside .eq. vertical) then
c        fmid = (index-.5)*hy
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
 50   if (badpts(idir,i) .lt. fmid) go to 60
c
c  point in top half. let it stay, increment counter
c
        itop = itop + 1
        if (itop+1 .ge. ibot) go to 80
             i = i + 1
             go to 50
c
c  point in bottom half. switch with a bottom point that's not yet
c  checked, and increment bot. pointer
c
 60    ibot           = ibot - 1
       temp           = badpts(1,ibot)
       badpts(1,ibot) = badpts(1,i)
       badpts(1,i)    = temp
       temp           = badpts(2,ibot)
       badpts(2,ibot) = badpts(2,i)
       badpts(2,i)    = temp
       if (itop+1 .lt. ibot) go to 50
c
c done smartbisecting icl'th clusters. adjust counts, repeat bisect stage .
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
