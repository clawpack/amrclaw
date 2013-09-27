c
c ----------------------------------------------------------
c
      subroutine prepc(level,nvar)
c
      use amr_module
      implicit double precision (a-h,o-z)

      parameter(numbcs=6)
c
c :::::::::::::::::::: PREPC ::::::::::::::::::::::::::::::::::::::
c
c this routine called because regridding just changed the fine grids.
c modify coarse grid boundary lists to store fluxes in appropriate
c fine grids lists.
c assume new fine grids have node(cfluxptr) initialized to point to null
c
c  first compute max. possible number of list cells. allocate
c  initially so that one pass through is enough.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      maxsp  = 0
      mkid   = lstart(level+1)
 10   if (mkid .eq. 0) go to 20
         ikeep  = (node(ndihi,mkid)-node(ndilo,mkid)+1)/intratx(level)
         jkeep  = (node(ndjhi,mkid)-node(ndjlo,mkid)+1)/intraty(level)
         kkeep  = (node(ndkhi,mkid)-node(ndklo,mkid)+1)/intratz(level)
         maxsp  = maxsp + 2*(ikeep*jkeep + jkeep*kkeep + kkeep*ikeep)
      mkid = node(levelptr,mkid)
      go to 10
 20   listsp(level) = maxsp
      if (maxsp .eq. 0) go to 99
c
      hxpar   = hxposs(level)
      hypar   = hyposs(level)
      hzpar   = hzposs(level)
      hxkid   = hxposs(level+1)
      hykid   = hyposs(level+1)
      hzkid   = hzposs(level+1)
      imax    = iregsz(level) - 1
      jmax    = jregsz(level) - 1
      kmax    = kregsz(level) - 1

      mpar = lstart(level)
 30   if (mpar .eq. 0) go to 99
c
       ispot   = 0
       ilo     = node(ndilo,mpar)
       jlo     = node(ndjlo,mpar)
       klo     = node(ndklo,mpar)
       ihi     = node(ndihi,mpar)
       jhi     = node(ndjhi,mpar)
       khi     = node(ndkhi,mpar)
       locbc   = igetsp(numbcs*maxsp)
c      #  initialize list space to 0 (0 terminator indicates end of bc list)
       do 35 i = 1,numbcs*maxsp
 35      alloc(locbc+i-1) = 0.d0
       node(cfluxptr,mpar) = locbc
c
       mkid = lstart(level+1)
 40    if (mkid .eq. 0) go to 60

          iclo = node(ndilo,mkid)/intratx(level)
          jclo = node(ndjlo,mkid)/intraty(level)
          kclo = node(ndklo,mkid)/intratz(level)
          ichi = node(ndihi,mkid)/intratx(level)
          jchi = node(ndjhi,mkid)/intraty(level)
          kchi = node(ndkhi,mkid)/intratz(level)

          iplo = max(ilo,iclo)
          jplo = max(jlo,jclo)
          kplo = max(klo,kclo)
          iphi = min(ihi,ichi)
          jphi = min(jhi,jchi)
          kphi = min(khi,kchi)

c   regular intersections
          if ((iplo .le. iphi+1 .and. jplo .le. jphi+1) .and.
     &                               (kplo .le. kphi+1))
     1          call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,klo,khi,iclo,ichi,jclo,jchi,kclo,kchi,
     3          nghost)
 
c   for fine grids touching periodic boundary on right
          if  (xperdom .and. ilo .eq. 0 .and. ichi .eq. imax)
     1        call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,klo,khi,
     &          iclo-iregsz(level),ichi-iregsz(level),
     3          jclo,jchi,kclo,kchi,nghost)
 
c   for fine grids touching periodic boundary on left
          if  (xperdom .and. iclo .eq. 0 .and. ihi .eq. imax)
     1        call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,klo,khi,
     &          iclo+iregsz(level),ichi+iregsz(level),
     3          jclo,jchi,kclo,kchi,nghost)
 
c   for fine grids touching periodic boundary on rear
          if  (yperdom .and. jlo .eq. 0 .and. jchi .eq. jmax)
     1          call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,klo,khi,iclo,ichi,
     3          jclo-jregsz(level),jchi-jregsz(level),kclo,kchi,nghost)
 
c   for fine grids touching periodic boundary on front
          if  (yperdom .and. jclo .eq. 0 .and. jhi .eq. jmax)
     1        call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,klo,khi,iclo,ichi,
     3          jclo+jregsz(level),jchi+jregsz(level),kclo,kchi,nghost)
 
c   for fine grids touching periodic boundary on top
          if  (zperdom .and. klo .eq. 0 .and. kchi .eq. kmax)
     1          call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,klo,khi,iclo,ichi,
     3          jclo,jchi,kclo-kregsz(level),kchi-kregsz(level),nghost)

c   for fine grids touching periodic boundary on bottom
          if  (zperdom .and. kclo .eq. 0 .and. khi .eq. kmax)
     1        call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,klo,khi,iclo,ichi,
     3          jclo,jchi,kclo+kregsz(level),kchi+kregsz(level),nghost)

 50     mkid = node(levelptr,mkid)
        go to 40
c
c  done with subgrid cycle. if no cells would need fixing, all done
c  else cycle through again to set up list with info. for bc processing
c
 60     continue
c
c  for now, leave unused space allocated to the grid. alternative is
c  to return (maxsp-ispot) amount starting at loc. node(cfluxptr,mpar)+ispot.
c
       mpar = node(levelptr,mpar)
       go to 30
c
 99    return
       end
