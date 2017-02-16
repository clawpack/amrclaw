!
! :::::::::::::::::::::::::::::: FILVAL ::::::::::::::::::::::::::
!
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! ------------------------------------------------------------------
!
subroutine filval(val, mitot, dx, level, time,  mic, &
                  xleft, xright, nvar, mptr, ilo, ihi, &
                  aux, naux)

    use amr_module, only: xlower, intratx, nghost, xperdom
    use amr_module, only: xupper, alloc
    use amr_module, only: outunit, NEEDS_TO_BE_SET, mcapa
    use amr_module, only: iregsz
    
    !for setaux timing
    use amr_module, only: timeSetaux, timeSetauxCPU

    implicit none

    ! Input
    integer, intent(in) :: mitot, level, mic, nvar, mptr, ilo, ihi
    integer, intent(in) :: naux
    real(kind=8), intent(in) :: dx, time, xleft, xright

    ! Output
    real(kind=8), intent(in out) :: val(nvar,mitot), aux(naux,mitot)

    ! Local storage
    integer :: refinement_ratio_x, iclo, ichi, ng
    integer :: ivar, i, ico, ifine, nx
    real(kind=8) :: valc(nvar,mic), auxc(naux,mic)
    real(kind=8) :: dx_coarse, xl, xr, area
    real(kind=8) :: s1m, s1p, slopex, xoff
    real(kind=8) :: fliparray((mitot)*nghost*(nvar+naux))
    real(kind=8) :: setflags(mitot),maxauxdif,aux2(naux,mitot)
    integer :: mjb
    logical :: sticksoutxfine,sticksoutxcrse
    
    !for setaux timing
    integer :: clock_start, clock_finish, clock_rate
    real(kind=8) :: cpu_start, cpu_finish


    ! External function definitions
    real(kind=8) :: get_max_speed

    refinement_ratio_x = intratx(level-1)
    dx_coarse  = dx * refinement_ratio_x
    xl      = xleft  - dx_coarse
    xr      = xright + dx_coarse

    ! if topo not yet final then aux is set outside filval (in gfixup)
    ! and so aux has real data already, (ie dont overwrite here)

    ! set integer indices for coarser patch enlarged by 1 cell
    ! (can stick out of domain). proper nesting will insure this one
    ! call is sufficient.
    iclo   = ilo / refinement_ratio_x - 1
    ichi   = (ihi + 1) / refinement_ratio_x - 1 + 1
    ng     = 0

    sticksoutxfine = ( (ilo .lt. 0) .or. (ihi .ge. iregsz(level)))
    sticksoutxcrse = ((iclo .lt. 0) .or. (ichi .ge. iregsz(level-1)))

    if (naux == 0) then
        if (xperdom .and. sticksoutxcrse) then
            call preintcopy(valc,mic,nvar,iclo,ichi,level-1,fliparray)
        else
            call intcopy(valc,mic,nvar,iclo,ichi,level-1,1)
        endif
    else  
        ! intersect grids and copy all (soln and aux)
        auxc(1,:) = NEEDS_TO_BE_SET
        if (xperdom .and. sticksoutxcrse) then
            call preicall(valc,auxc,mic,nvar,naux,iclo,ichi, &
                          level-1,fliparray)
        else
            call icall(valc,auxc,mic,nvar,naux,iclo,ichi,level-1,1)
        endif
!!$        do i = 1, mic
!!$          if (auxc(1,:) == NEEDS_TO_BE_SET) then
!!$             write(*,*)" *** coarsenened new fine grid not completely set from previously"  &
!!$                       "existing coarse grids ***"
!!$             stop
!!$          endif
!!$        end do
           ! no ghost cells on coarse enlarged patch. set any remaining
           ! vals. should only be bcs that stick out of domain.
           call system_clock(clock_start, clock_rate)
           call cpu_time(cpu_start)
           call setaux(ng,mic,xl,dx_coarse,naux,auxc)
           call system_clock(clock_finish, clock_rate)
           call cpu_time(cpu_finish)
           timeSetaux = timeSetaux + clock_finish - clock_start
           timeSetauxCPU = timeSetauxCPU + cpu_finish - cpu_start
    endif

    call bc1amr(valc,auxc,mic,nvar,naux,dx_coarse,level-1,   &
                time,xl,xr)

!  NOTE change in order of code.  Since the interp from coarse to fine needs the aux
!       arrays set already, the fine copy is done first, to set up the aux arrays.
!       we can do this since we have the flag array to test where to overwrite.

!  SO this is no longer overwriting but setting for the first time.
! overwrite interpolated values with fine grid values, if available.
    nx = mitot - 2*nghost

    if (naux .gt. 0) then 
!       ## NEEDS_TO_BE_SET is signal that aux array not set.
!       ## after calling icall to copy aux from other grids
!       ## any remaining NEEDS_TO_BE_SET signals will be set in setaux.
!       ## it also signals where soln was copied, so it wont be
!       ## overwritten with coarse grid interpolation
        aux(1,:) = NEEDS_TO_BE_SET  ! indicates fine cells not yet set.

        if (xperdom.and.sticksoutxfine) then
            call preicall(val,aux,mitot,nvar,naux,ilo-nghost,ihi+nghost, &
                          level,fliparray)
        else
            call icall(val,aux,mitot,nvar,naux,ilo-nghost,ihi+nghost,  &
                      level,1)
        endif
        setflags = aux(1,:)   ! save since will overwrite in setaux when setting all aux vals
           ! need this so we know where to use coarse grid to set fine solution w/o overwriting
           ! set remaining aux vals not set by copying from prev existing grids
        call system_clock(clock_start, clock_rate)
        call cpu_time(cpu_start)
        call setaux(nghost,nx,xleft,dx,naux,aux)
        call system_clock(clock_finish, clock_rate)
        call cpu_time(cpu_finish)
        timeSetaux = timeSetaux + clock_finish - clock_start
        timeSetauxCPU = timeSetauxCPU + cpu_finish - cpu_start
    else ! either no aux exists, or cant reuse yet  
         ! so only call intcopy (which copies soln) and not icall.
         ! in this case flag q(1,:) to NEEDS_TO_BE_SET flag so wont be overwritten
         ! by coarse grid interp.  this is needed due to reversing order of
         ! work - first copy from fine grids, then interpolate from coarse grids
        val(1,:) = NEEDS_TO_BE_SET
        if (xperdom.and.sticksoutxfine) then
            call preintcopy(val,mitot,nvar,ilo-nghost,ihi+nghost,     &
                            level,fliparray)
        else
            call intcopy(val,mitot,nvar,ilo-nghost,ihi+nghost,  &
                         level,1)
        endif
        setflags = val(1,:)  ! remaining flags signals need to set
    endif

   
    ! Prepare slopes - use min-mod limiters

     do i=2, mic-1
       do ivar = 1, nvar
 
            s1p = valc(ivar,i+1) - valc(ivar,i)
            s1m = valc(ivar,i)   - valc(ivar,i-1)
            slopex = min(abs(s1p), abs(s1m)) &
                         * sign(1.d0,valc(ivar,i+1) - valc(ivar,i-1))
            ! if there's a sign change, set slope to 0.
            if ( s1m*s1p <=  0.d0) slopex = 0.d0


            do ico = 1,refinement_ratio_x
              xoff = (real(ico,kind=8) - 0.5d0) / refinement_ratio_x - 0.5d0
              ifine = (i-2) * refinement_ratio_x + nghost + ico

              if (setflags(ifine) .eq. NEEDS_TO_BE_SET) then
                val(ivar,ifine) = valc(ivar,i) + xoff*slopex
              endif

            end do

       enddo !end of ivar loop
     enddo !end of coarse i loop

    ! adjust to conserve kappa*q, but only where coarse grid was interpolated
    ! so now need to pass setflags to this subr.
    if (mcapa .ne. 0) then  
        call fixcapaq(val,aux,mitot,valc,auxc,mic,nvar,naux,level-1,setflags)
    endif
 
!!$! CHECK BY CALLING SETAUX AND SETTING ALL, THEN DIFFING
!!$   mjb = 0
!!$   if (naux .gt. 0  .and. mjb .eq. 1) then
!!$      aux2(1,:,:) = NEEDS_TO_BE_SET   ! indicates fine cells not yet set
!!$      call setaux(nghost,nx,ny,xleft,ybot,dx,dy,naux,aux2)
!!$      maxauxdif = 1.d-13
!!$      do i = 1, mitot
!!$      do j = 1, mjtot
!!$         if (abs(aux(1,i,j)-aux2(1,i,j)) .gt. maxauxdif) then
!!$            maxauxdif = abs(aux(1,i,j)-aux2(1,i,j))
!!$            write(*,444)i,j,aux(1,i,j),aux2(1,i,j),maxauxdif
!!$444         format("i,j = ",2i4," auxs ",2e15.7," maxauxdif ",e12.5)
!!$         endif
!!$      end do
!!$      end do
!!$      if (maxauxdif .gt. 2.d-13) then
!!$         write(*,*)" maxauxdif = ",maxauxdif," with mitot,mjtot ",mitot,mjtot, &
!!$              " on grid ",mptr," level ",level
!!$      endif
!!$   endif

end subroutine filval

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine dumpaux(aux,naux,mitot)
   implicit none
   real(kind=8) :: aux(naux,mitot)
   integer :: naux,mitot,i,iaux

   do i = 1, mitot 
      write(*,444) i,(aux(iaux,i),iaux=1,naux)
 444  format(2i4,5e12.5)
   end do

end subroutine dumpaux
