c
!> For new fine grids (grids on level **level**), allocate space for 
!! saving fluxes at boundary (even if it's physical boundary) after 
!! each integration step.
!! These fluxes are for future conservative fixing of coarse grids
!! (level **level**-1 grids).
!! The address of this space is stored in node(ffluxptr, mkid) for grid
!! mkid. 
!!
!! Note that if the refinment ratio is r, fluxes from every r cells 
!! on grid **mkid** are saved to one slot in this space, since every r
!! cells on grid **mkid** border one cell on level **level**-1 grids.
!! \param[in] level boudnary lists of grids on this level get updated
!! \param[in] nvar number of equations for the system
!! \param[in] naux number of auxiliary variables for the system
c ----------------------------------------------------------
c
      subroutine prepf(level,nvar,naux)
c
      use amr_module
      implicit double precision (a-h,o-z)

c
c ::::::::::::::::::::::::::: PREPF :::::::::::::::::::::::::::::
c
c prepare new fine grids to save fluxes after each integration step
c for future conservative fixing of coarse grids
c save all boundary fluxes of fine grid (even if phys. bndry.) -
c but only save space for every intrat. (remember - 4 fluxes).
c
c :::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::::::
c
      mkid = lstart(level)
 10   if (mkid .eq. 0) go to 99
          nx    = node(ndihi,mkid)-node(ndilo,mkid) + 1
          ny    = node(ndjhi,mkid)-node(ndjlo,mkid) + 1
          ikeep = nx/intratx(level-1)
          jkeep = ny/intraty(level-1)
          lenbc = 2*(ikeep+jkeep)
c
c         get twice the storage, one for plus or minus fluxes, the
c         other for coarse solution for wave fixing. also for aux vars.
c
          node(ffluxptr,mkid) = igetsp(2*nvar*lenbc+naux*lenbc)
          ist   = node(ffluxptr,mkid)

          do 20 i = 1, 2*nvar*lenbc + naux*lenbc
             alloc(ist+i-1) = 0.d0
 20       continue
          mkid = node(levelptr,mkid)
          go to 10

 99    return
       end
