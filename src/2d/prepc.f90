!
!> This routine is called because regridding just changed the fine grids.
!! It modifies boundary list of each level **level** grid 
!! such that each grid knows along its boundary, what fine cells 
!! its cell borders and where to get saved flux for conservation fix-up
!! 
!! \param[in] level boudnary lists of grids on this level get updated
!! \param[in] nvar number of equations for the system
! ----------------------------------------------------------
!
subroutine prepc(level,nvar)
!
      use amr_module
#ifdef CUDA
        use cuda_module, only: device_id
        use memory_module, only: gpu_allocate, cpu_allocate_pinned
#endif
      implicit double precision (a-h,o-z)

!
! :::::::::::::::::::: PREPC ::::::::::::::::::::::::::::::::::::::
!
! this routine called because regridding just changed the fine grids.
! modify coarse grid boundary lists to store fluxes in appropriate
! fine grids lists.
! assume new fine grids have node(cfluxptr) initialized to null
!
!  first compute max. possible number of list cells. allocate
!  initially so that one pass through is enough.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      maxsp  = 0
      mkid   = lstart(level+1)
 10   if (mkid .eq. 0) go to 20
         ikeep  = (node(ndihi,mkid)-node(ndilo,mkid)+1)/intratx(level)
         jkeep  = (node(ndjhi,mkid)-node(ndjlo,mkid)+1)/intraty(level)
         maxsp  = maxsp + 2*(ikeep+jkeep)
      mkid = node(levelptr,mkid)
      go to 10
 20   listsp(level) = maxsp
      if (maxsp .eq. 0) go to 99
!
      hxpar   = hxposs(level)
      hypar   = hyposs(level)
      hxkid   = hxposs(level+1)
      hykid   = hyposs(level+1)
      imax    = iregsz(level) - 1
      jmax    = jregsz(level) - 1

      ! The code from below to the end of this file does the following:
      ! for mpar in (all level "level" grid):
      !     for mkid (all level "level"+1 grid):
      !         if mkid is encompassed by mpar:
      !             setuse(mkid)
      !         endif
      !     endfor
      ! endfor
      mpar = lstart(level)

 30   if (mpar .eq. 0) go to 99
!
       ispot   = 0
       ilo     = node(ndilo,mpar)
       jlo     = node(ndjlo,mpar)
       ihi     = node(ndihi,mpar)
       jhi     = node(ndjhi,mpar)
       locbc   = igetsp(5*maxsp)
!      initialize list to 0 (0 terminator indicates end of bc list)
       do 35 i = 1,5*maxsp
 35      alloc(locbc+i-1) = 0.d0
       node(cfluxptr,mpar) = locbc
#ifdef NOTHING
        call cpu_allocate_pinned(cflux(mpar)%ptr, 1, 5, 1, maxsp)
        call gpu_allocate(cflux_d(mpar)%ptr, device_id, 1, 5, 1, maxsp)
        do jj = 1,maxsp
            do ii = 1, 5
                cflux(mpar)%ptr(ii,jj) = 0
                ! cflux_d(mpar)%ptr(ii,jj) = 0
            enddo
        enddo
        ! cflux(mpar)%ptr = 0
        ! cflux_d(mpar)%ptr = 0
#endif
!
       mkid = lstart(level+1)
 40    if (mkid .eq. 0) go to 60

          ! ilo is smallest index i of grid mpar on level "level"
          ! iclo tells horizontal index of the cell on level "level" that 
          ! encompasses the left-most cell of grid mkid 
          ! both ilo and iclo is based on cell index coordinate system
          ! on level "level"
          ! By comparing the value of ilo and iclo, we can tell their
          ! relative location. For instance, if ilo < iclo, then the left
          ! border of grid mpar is to the left of left border of grid
          ! mkid. By comparing all four borders, we can tell if grid mkid is
          ! encompassed by grid mpar.
          iclo = node(ndilo,mkid)/intratx(level)
          jclo = node(ndjlo,mkid)/intraty(level)
          ichi = node(ndihi,mkid)/intratx(level)
          jchi = node(ndjhi,mkid)/intraty(level)

          iplo = max(ilo,iclo)
          jplo = max(jlo,jclo)
          iphi = min(ihi,ichi)
          jphi = min(jhi,jchi)

!   regular intersections (will check in setuse that no duplication)
!   this first call is only interior interfaces. 

          ! check if this level "level"+1 grid, mkid, is encompassed by the
          ! level "level" grid, mpar.
          ! If so, write info along the border of this grid, mkid, to
          ! listbc
          if (iplo .le. iphi+1 .and. jplo .le. jphi+1) then
               kflag = 1 ! interior stuff, no mappings
                call setuse(alloc(locbc),maxsp,ispot,mkid, &
                ilo,ihi,jlo,jhi,iclo,ichi,jclo,jchi,kflag)
          endif

!   for fine grids touching periodic boundary on right
          if  (xperdom .and. ilo .eq. 0 .and. ichi .eq. imax) then
              kflag = 1 ! periodic in x
              call setuse(alloc(locbc),maxsp,ispot,mkid,&
                ilo,ihi,jlo,jhi,iclo-iregsz(level),ichi-iregsz(level),&
                jclo,jchi,kflag)
           endif

!   for fine grids touching periodic boundary on left
          if  (xperdom .and. iclo .eq. 0 .and. ihi .eq. imax) then
              kflag = 1
              call setuse(alloc(locbc),maxsp,ispot,mkid,&
                ilo,ihi,jlo,jhi,iclo+iregsz(level),ichi+iregsz(level),&
                jclo,jchi,kflag)
          endif

!   for fine grids touching periodic boundary on top
          if  (yperdom .and. jlo .eq. 0 .and. jchi .eq. jmax) then
                kflag = 1
                call setuse(alloc(locbc),maxsp,ispot,mkid,&
                ilo,ihi,jlo,jhi,iclo,ichi,&
                jclo-jregsz(level),jchi-jregsz(level),kflag)
          endif

!   for fine grids touching periodic boundary on bottom
          if  (yperdom .and. jclo .eq. 0 .and. jhi .eq. jmax)  then
              kflag = 1
              call setuse(alloc(locbc),maxsp,ispot,mkid, &
                ilo,ihi,jlo,jhi,iclo,ichi, &
                jclo+jregsz(level),jchi+jregsz(level),kflag)
          endif

!   for fine grids touching boundary on top in spherically mapped case
!   and coarse grid touches top too. see if (mapped) x extent overlap.
          if  (spheredom .and. jhi .eq. jmax .and. jchi .eq. jmax) then
               kflag = 2
!              write(dbugunit,*)" for coarse grid ",mpar
               iwrap2 = iregsz(level) - iclo - 1  !higher mapped index
               iwrap1 = iregsz(level) - ichi - 1  !lower mapped index
               if (max(ilo,iwrap1) .le. min(ihi,iwrap2)) then
                  call setuse(alloc(locbc),maxsp,ispot,mkid,&
                              ilo,ihi,jlo,jhi,iclo,ichi,&
                              jclo,jchi,kflag)
               endif
          endif

!   fine grids touching boundary on bottom for spherically mapped case
!   coarse grid touches bottom too. see if (mapped) x extents overlap
          if  (spheredom .and. jclo .eq. 0 .and. jlo .eq. 0) then
               kflag = 3
               iwrap2 = iregsz(level) - iclo - 1  !higher mapped index
               iwrap1 = iregsz(level) - ichi - 1  !lower mapped index
               if (max(ilo,iwrap1) .le. min(ihi,iwrap2)) then
                  call setuse(alloc(locbc),maxsp,ispot,mkid,&
                              ilo,ihi,jlo,jhi,iclo,ichi,&
                              jclo,jchi,kflag)
               endif
          endif

 50     mkid = node(levelptr,mkid)
        go to 40
!
!  done with subgrid cycle. if no cells would need fixing, all done
!  else cycle through again to set up list with info. for bc processing
!
 60     continue
!
!  for now, leave unused space allocated to the grid. alternative is to
!  return (maxsp-ispot) amt starting at loc node(cfluxptr)+ispot.
!
#ifdef CUDA
        ! call alloc_to_int(cflux(mpar)%ptr, alloc(locbc), maxsp)
#endif
       mpar = node(levelptr,mpar)
       go to 30
!
 99    return
end subroutine prepc

subroutine alloc_to_int(dst_int, src_real, maxsp)
      implicit double precision (a-h,o-z)
      dimension src_real(5,maxsp)
      integer, intent(inout) :: dst_int(5, maxsp)
      integer :: i,j
      do j = 1, maxsp
          do i = 1,5
              dst_int(i,j) = src_real(i,j)
          enddo
      enddo
end subroutine alloc_to_int
