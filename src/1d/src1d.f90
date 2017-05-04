subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)

    ! This routine should be a simplified version of src1
    ! which applies source terms for a 1-d slice of data along the
    ! edge of a grid.  This is called only from qad where the conservative
    ! fix-up is applied and is used to apply source terms over partial
    ! time steps to the coarse grid cell values used in solving Riemann 
    ! problems at the interface between coarse and fine grids.
 
    ! If the source terms depend only on q, it should be easy to 
    ! adapt src2 to create this routine, just loop over 1:mx1d.
    ! If the source terms are more complicated, it may not be easy.
 
    ! The code may work fine without applying source terms in this
    ! context, so using this dummy routine might be successful even when
    ! source terms are present. 
 
    ! This default version does nothing.

    ! For the 1d code, does src1 and src1d need to be different files?
 
    implicit none
    integer, intent(in) :: meqn,mbc,mx1d,maux
    real(kind=8), intent(in) :: t, dt
    real(kind=8), intent(in) ::  aux1d(maux,mx1d)
    real(kind=8), intent(inout) ::  q1d(meqn,mx1d)

end subroutine src1d
