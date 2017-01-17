! -------------------------------------------------------------------
subroutine flag2refine(mx,my,mz,mbc,meqn,maux,xlower,ylower, &
                        zlower,dx,dy,dz,t, &
                        level,tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)
! -------------------------------------------------------------------

! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! The logical function allowflag(x,y,t,level) is called to check whether
! further refinement at this level is allowed at this particular location
! and time.  The default library version of this routine returns .true.
! for all arguments.  Copy that routine to the application directory and
! modify it if needed to restrict the region where refinement is allowed.
!
! First, each point is checked against the min_level and max_level
! requirements of any regions present.  If no changes need to be made,
! the infinity norm of the stress tensor is checked against the user
! specified tolsp value.  This function assumes the first 6 components of
! q are the 6 stress tensor components.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
    use regions_module, only: regions, num_regions

    implicit none

    integer, intent(in) :: mx, my, mz, mbc, meqn, maux, level
    real (kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real (kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real (kind=8), intent(inout) :: amrflags(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    logical     allowflag
    external    allowflag
    real (kind=8), intent(in) :: DOFLAG, DONTFLAG, xlower, ylower, zlower, dx, dy, dz, tolsp, t

!    real (kind=8) :: xcell, ycell, zcell, max_stress
    real(kind=8) :: xcell, ycell, zcell, dq(meqn), dqi(meqn), dqj(meqn), dqk(meqn)
    integer :: i, j, k, m, min_level, max_level
    integer :: infinity = 1e3

!   # loop over interior points on this grid:
    do k = 1,mz
        zcell = zlower + (k-0.5d0)*dz
        do j = 1,my
            ycell = ylower + (j-0.5d0)*dy
            do i = 1,mx
                xcell = xlower + (i-0.5d0)*dx

!               # obtain the overall min and max levels from any regions intersecting the cell
                min_level = 0
                max_level = 0
                do m =1,num_regions
                    if (regions(m)%t_low .le. t .and. t .le. regions(m)%t_hi .and. &
                        regions(m)%x_low .le. xcell + 0.5d0*dx .and. xcell - 0.5d0*dx .le. regions(m)%x_hi .and. &
                        regions(m)%y_low .le. ycell + 0.5d0*dy .and. ycell - 0.5d0*dy .le. regions(m)%y_hi .and. &
                        regions(m)%z_low .le. zcell + 0.5d0*dz .and. zcell - 0.5d0*dz .le. regions(m)%z_hi) then
                        min_level = max(min_level, regions(m)%min_level)
                        max_level = max(max_level, regions(m)%max_level)
                    end if
                end do

!               # if the cell intersects any region, make sure that cell is refined as specified
!               # if nothing needs to be changed, use specified tolerance
                if (min_level > 0 .and. level < min_level) then
                    amrflags(i,j,k) = DOFLAG
                else if (max_level > 0 .and. max_level <= level) then
                    amrflags(i,j,k) = DONTFLAG
                else if (allowflag(xcell,ycell,zcell,t,level)) then
!                    max_stress = 0.d0
                    dq = 0.d0
                    dqi = abs(q(:,i+1,j,k) - q(:,i-1,j,k))
                    dqj = abs(q(:,i,j+1,k) - q(:,i,j-1,k))
                    dqk = abs(q(:,i,j,k+1) - q(:,i,j,k-1))
                    dq = max(dq,dqi,dqj,dqk)

                    do m = 1,meqn
!                        max_stress = max(max_stress, dabs(q(m,i,j,k)))
                        if (dq(m) > tolsp) then
                             amrflags(i,j,k) = DOFLAG
                        else
                             amrflags(i,j,k) = DONTFLAG
                        end if
                    end do
                else
                    amrflags(i,j,k) = DONTFLAG
                end if

            end do
        end do
    end do

    return
end
