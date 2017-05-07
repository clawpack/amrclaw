!> Fill grid **mptr** on level **level** by copying values from OLD level **level**
!! grids if available, otherwise by interpolating values from coarser grids.
subroutine filval(val, mitot, mjtot, dx, dy, level, time,  mic,          &
                  mjc, xleft, xright, ybot, ytop, nvar, mptr, ilo, ihi,  &
                  jlo, jhi, aux, naux)

    use amr_module, only: xlower, ylower, intratx, intraty, nghost, xperdom
    use amr_module, only: yperdom, spheredom, xupper, yupper, alloc
    use amr_module, only: outunit, NEEDS_TO_BE_SET, mcapa
    use amr_module, only: newstl, iregsz, jregsz
    
    implicit none

    ! Input
    integer, intent(in) :: mitot, mjtot, level, mic, mjc, nvar, mptr, ilo, ihi
    integer, intent(in) :: jlo, jhi, naux
    real(kind=8), intent(in) :: dx, dy, time, xleft, xright, ybot, ytop

    ! Output
    real(kind=8), intent(in out) :: val(nvar,mitot,mjtot), aux(naux,mitot,mjtot)

    ! Local storage
    integer :: refinement_ratio_x, refinement_ratio_y, iclo, jclo, ichi, jchi, ng
    integer :: ivar, i, j, ico, jco, ifine, jfine, nx, ny
    real(kind=8) :: valc(nvar,mic,mjc), auxc(naux,mic,mjc)
    real(kind=8) :: dx_coarse, dy_coarse, xl, xr, yb, yt, area
    real(kind=8) :: s1m, s1p, slopex, slopey, xoff, yoff
    real(kind=8) :: fliparray((mitot+mjtot)*nghost*(nvar+naux))
    real(kind=8) :: setflags(mitot,mjtot),maxauxdif,aux2(naux,mitot,mjtot)
    integer :: mjb
    integer :: jm, im, nm 
    logical :: sticksoutxfine, sticksoutyfine,sticksoutxcrse,sticksoutycrse
    
    !for setaux timing
    integer :: clock_start, clock_finish, clock_rate
    real(kind=8) :: cpu_start, cpu_finish


    ! External function definitions
    real(kind=8) :: get_max_speed


    refinement_ratio_x = intratx(level-1)
    refinement_ratio_y = intraty(level-1)
    dx_coarse  = dx * refinement_ratio_x
    dy_coarse  = dy * refinement_ratio_y
    xl      = xleft  - dx_coarse
    xr      = xright + dx_coarse
    yb      = ybot   - dy_coarse
    yt      = ytop   + dy_coarse

    ! if topo not yet final then aux is set outside filval (in gfixup)
    ! and so aux has real data already, (ie dont overwrite here)

    ! set integer indices for coarser patch enlarged by 1 cell
    ! (can stick out of domain). proper nesting will insure this one
    ! call is sufficient.
    iclo   = ilo / refinement_ratio_x - 1
    jclo   = jlo / refinement_ratio_y - 1
    ichi   = (ihi + 1) / refinement_ratio_x - 1 + 1
    jchi   = (jhi + 1) / refinement_ratio_y - 1 + 1
    ng     = 0

    sticksoutxfine = ( (ilo .lt. 0) .or. (ihi .ge. iregsz(level)))
    sticksoutyfine = ( (jlo .lt. 0) .or. (jhi .ge. jregsz(level)))
    sticksoutxcrse = ((iclo .lt. 0) .or. (ichi .ge. iregsz(level-1)))
    sticksoutycrse = ((jclo .lt. 0) .or. (jchi .ge. jregsz(level-1)))

    if (naux == 0) then
        if ((xperdom .and. sticksoutxcrse).or. (yperdom.and.sticksoutycrse)  .or. spheredom) then
            call preintcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,level-1,fliparray)
        else
            call intcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,level-1,1,1)
        endif
    else  
        ! intersect grids and copy all (soln and aux)
        auxc(1,:,:) = NEEDS_TO_BE_SET
        if ((xperdom .and.sticksoutxcrse) .or. (yperdom.and.sticksoutycrse) .or. spheredom) then
            call preicall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi, &
                          level-1,fliparray)
        else
            call icall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi,level-1,1,1)
        endif
!!$        do i = 1, mic
!!$        do j = 1, mjc
!!$          if (auxc(1,:,:) == NEEDS_TO_BE_SET) then
!!$             write(*,*)" *** coarsenened new fine grid not completely set from previously"  &
!!$                       "existing coarse grids ***"
!!$             stop
!!$          endif
!!$        end do
!!$        end do
           ! no ghost cells on coarse enlarged patch. set any remaining
           ! vals. should only be bcs that stick out of domain.
           call setaux(ng,mic,mjc,xl,yb,dx_coarse,dy_coarse,naux,auxc)
    endif

    call bc2amr(valc,auxc,mic,mjc,nvar,naux,dx_coarse,dy_coarse,level-1,time,xl,xr,yb, &
                yt,xlower,ylower,xupper,yupper,xperdom,yperdom,spheredom)



!  NOTE change in order of code.  Since the interp from coarse to fine needs the aux
!       arrays set already, the fine copy is done first, to set up the aux arrays.
!       we can do this since we have the flag array to test where to overwrite.

!  SO this is no longer overwriting but setting for the first time.
! overwrite interpolated values with fine grid values, if available.
    nx = mitot - 2*nghost
    ny = mjtot - 2*nghost

    if (naux .gt. 0) then 
!       ## NEEDS_TO_BE_SET is signal that aux array not set.
!       ## after calling icall to copy aux from other grids
!       ## any remaining NEEDS_TO_BE_SET signals will be set in setaux.
!       ## it also signals where soln was copied, so it wont be
!       ## overwritten with coarse grid interpolation
        aux(1,:,:) = NEEDS_TO_BE_SET  ! indicates fine cells not yet set. 

        if ((xperdom.and.sticksoutxfine) .or. (yperdom.and.sticksoutyfine) .or. spheredom) then
            call preicall(val,aux,mitot,mjtot,nvar,naux,ilo-nghost,ihi+nghost, & 
                          jlo-nghost,jhi+nghost,level,fliparray)
        else
            call icall(val,aux,mitot,mjtot,nvar,naux,ilo-nghost,ihi+nghost,  &
                      jlo-nghost,jhi+nghost,level,1,1)   
        endif
        setflags = aux(1,:,:)   ! save since will overwrite in setaux when setting all aux vals
           ! need this so we know where to use coarse grid to set fine solution w/o overwriting
           ! set remaining aux vals not set by copying from prev existing grids
        call setaux(nghost,nx,ny,xleft,ybot,dx,dy,naux,aux)
    else ! either no aux exists, or cant reuse yet  
         ! so only call intcopy (which copies soln) and not icall.
         ! in this case flag q(1,:) to NEEDS_TO_BE_SET flag so wont be overwritten
         ! by coarse grid interp.  this is needed due to reversing order of
         ! work - first copy from fine grids, then interpolate from coarse grids
        val(1,:,:) = NEEDS_TO_BE_SET
        if ((xperdom.and.sticksoutxfine) .or. (yperdom.and.sticksoutyfine) .or. spheredom) then
            call preintcopy(val,mitot,mjtot,nvar,ilo-nghost,ihi+nghost,     &
                            jlo-nghost,jhi+nghost,level,fliparray)
        else
            call intcopy(val,mitot,mjtot,nvar,ilo-nghost,ihi+nghost,  &
                         jlo-nghost,jhi+nghost,level,1,1)   
        endif
        setflags = val(1,:,:)  ! remaining flags signals need to set
    endif

   
    ! Prepare slopes - use min-mod limiters

    do j=2, mjc-1
     do i=2, mic-1
       do ivar = 1, nvar
 
            s1p = valc(ivar,i+1,j) - valc(ivar,i,j)
            s1m = valc(ivar,i,j)   - valc(ivar,i-1,j)
            slopex = min(abs(s1p), abs(s1m)) &
                         * sign(1.d0,valc(ivar,i+1,j) - valc(ivar,i-1,j))
            ! if there's a sign change, set slope to 0.
            if ( s1m*s1p <=  0.d0) slopex = 0.d0

            s1p = valc(ivar,i,j+1) - valc(ivar,i,j)
            s1m = valc(ivar,i,j)   - valc(ivar,i,j-1)
            slopey = min(abs(s1p), abs(s1m))  &
                         * sign(1.0d0, valc(ivar,i,j+1) - valc(ivar,i,j-1))
            if ( s1m*s1p <= 0.d0)  slopey = 0.d0

            ! Interpolate from coarse cells to fine grid to find depth
             do jco = 1,refinement_ratio_y
               yoff = (real(jco,kind=8) - 0.5d0) / refinement_ratio_y - 0.5d0
               jfine = (j-2) * refinement_ratio_y + nghost + jco

               do ico = 1,refinement_ratio_x
                 xoff = (real(ico,kind=8) - 0.5d0) / refinement_ratio_x - 0.5d0
                 ifine = (i-2) * refinement_ratio_x + nghost + ico

                 if (setflags(ifine,jfine) .eq. NEEDS_TO_BE_SET) then
                    val(ivar,ifine,jfine) = valc(ivar,i,j) + xoff*slopex + yoff*slopey
                 endif

               end do
            end do

       enddo !end of ivar loop
     enddo !end of coarse i loop
    enddo !end of coarse j loop

    ! adjust to conserve kappa*q, but only where coarse grid was interpolated
    ! so now need to pass setflags to this subr.
    if (mcapa .ne. 0) then  
        call fixcapaq(val,aux,mitot,mjtot,valc,auxc,mic,mjc,nvar,naux,level-1,setflags)
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

subroutine dumpaux(aux,naux,mitot,mjtot)
   implicit none
   real(kind=8) :: aux(naux,mitot,mjtot)
   integer :: naux,mitot,mjtot,i,j,iaux

   do j = 1, mjtot 
   do i = 1, mitot 
      write(*,444) i,j,(aux(iaux,i,j),iaux=1,naux)
 444  format(2i4,5e12.5)
   end do
   end do

end subroutine dumpaux
