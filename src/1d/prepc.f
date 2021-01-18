c
c ----------------------------------------------------------
c
      subroutine prepc(level,nvar)
c
      use amr_module
      implicit double precision (a-h,o-z)

      parameter(numbcs=4)
c
c :::::::::::::::::::: PREPC ::::::::::::::::::::::::::::::::::::::
c
c this routine called because regridding just changed the fine grids.
c modify coarse grid boundary lists to store fluxes in appropriate
c fine grids lists.
c assume new fine grids have node(cfluxptr) initialized to null
c
c  first compute max. possible number of list cells. allocate
c  initially so that one pass through is enough.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      maxsp  = 0
      mkid   = lstart(level+1)

 10   if (mkid .eq. 0) go to 20
c        Adding twice perimeter (which is just one cell
c         on each side)
         maxsp  = maxsp + 2
      mkid = node(levelptr,mkid)
      go to 10
 20   listsp(level) = maxsp
      if (maxsp .eq. 0) go to 99
c
      hxpar   = hxposs(level)
      hxkid   = hxposs(level+1)
      imax    = iregsz(level) - 1

      mpar = lstart(level)
 30   if (mpar .eq. 0) go to 99
c
       ispot   = 0
       ilo     = node(ndilo,mpar)
       ihi     = node(ndihi,mpar)
       locbc   = igetsp(numbcs*maxsp)
c      #  initialize list to 0 (0 terminator indicates end of bc list)
       do 35 i = 1,numbcs*maxsp
 35      alloc(locbc+i-1) = 0.d0
       node(cfluxptr,mpar) = locbc
c
       mkid = lstart(level+1)
 40    if (mkid .eq. 0) go to 60

          iclo = node(ndilo,mkid)/intratx(level)
          ichi = node(ndihi,mkid)/intratx(level)

          iplo = max(ilo,iclo)
          iphi = min(ihi,ichi)

c   regular intersections (will check in setuse that no duplication)
c   this first call is only interior interfaces. 

          if (iplo .le. iphi+1) then
               kflag = 1 ! interior stuff, no mappings

                call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,iclo,ichi,kflag)

          endif

c   for fine grids touching periodic boundary on right
          if  (xperdom .and. ilo .eq. 0 .and. ichi .eq. imax) then
              kflag = 1 ! periodic in x
              call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,iclo-iregsz(level),ichi-iregsz(level),
     3          kflag)
           endif

c   for fine grids touching periodic boundary on left
          if  (xperdom .and. iclo .eq. 0 .and. ihi .eq. imax) then
              kflag = 1
              call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,iclo+iregsz(level),ichi+iregsz(level),
     3          kflag)
          endif

 50     mkid = node(levelptr,mkid)
        go to 40
c
c  done with subgrid cycle. if no cells would need fixing, all done
c  else cycle through again to set up list with info. for bc processing
c
 60     continue
c
c  for now, leave unused space allocated to the grid. alternative is to
c  return (maxsp-ispot) amt starting at loc node(cfluxptr,mpar)+ispot.
c
       mpar = node(levelptr,mpar)
       go to 30
c
 99    return
       end
