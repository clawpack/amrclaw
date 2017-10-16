module memory_module

    use iso_c_binding

    implicit none

    integer, parameter, private :: szr = 8

    interface gpu_allocate
        module procedure gpu_allocate_r1
        module procedure gpu_allocate_r2
        module procedure gpu_allocate_r3
    end interface
    interface gpu_deallocate
        module procedure gpu_deallocate_r1
        module procedure gpu_deallocate_r2
        module procedure gpu_deallocate_r3
    end interface


contains
    subroutine gpu_allocate_r1(a, dev_id, lo1, hi1)
        use cudafor
        real(kind=8), pointer, device, intent(inout) :: a(:)
        integer, intent(in) :: lo1, hi1
        integer, intent(in) :: dev_id
        integer :: n1
        integer (kind=c_size_t) :: sz
        type(c_devptr) :: cp
        real(kind=8), pointer, device :: fp(:)
        integer :: istat
        n1 = max(hi1-lo1+1, 1)
        sz = int(n1,c_size_t)
        istat = cudaSetDevice(dev_id)
        istat = cudaMalloc(cp, szr*sz)
        call c_f_pointer(cp, fp, shape=(/n1/))
        call shift_bound_d1(fp, lo1, a)
    contains
        subroutine shift_bound_d1(fp, lo1, a)
            integer, intent(in) :: lo1
            real(kind=8), target,  device, intent(in) :: fp(lo1:)
            real(kind=8), pointer, device, intent(inout) :: a(:)
            a => fp
        end subroutine shift_bound_d1
    end subroutine gpu_allocate_r1

    subroutine gpu_allocate_r2(a, dev_id, lo1, hi1, lo2, hi2)
        use cudafor
        real(kind=8), pointer, device, intent(inout) :: a(:,:)
        integer, intent(in) :: lo1, hi1, lo2, hi2
        integer, intent(in) :: dev_id
        integer :: n1, n2
        integer (kind=c_size_t) :: sz
        type(c_devptr) :: cp
        real(kind=8), pointer, device :: fp(:,:)
        integer :: istat
        n1 = max(hi1-lo1+1, 1)
        n2 = max(hi2-lo2+1, 1)
        sz = int(n1,c_size_t) * int(n2,c_size_t)
        istat = cudaSetDevice(dev_id)
        istat = cudaMalloc(cp, szr*sz)
        call c_f_pointer(cp, fp, shape=(/n1,n2/))
        call shift_bound_d2(fp, lo1, lo2, a)
    contains
        subroutine shift_bound_d2(fp, lo1, lo2, a)
            integer, intent(in) :: lo1, lo2
            real(kind=8), target,  device, intent(in) :: fp(lo1:,lo2:)
            real(kind=8), pointer, device, intent(inout) :: a(:,:)
            a => fp
        end subroutine shift_bound_d2
    end subroutine gpu_allocate_r2

    subroutine gpu_allocate_r3(a, dev_id, lo1, hi1, lo2, hi2, lo3, hi3)
        use cudafor
        real(kind=8), pointer, device, intent(inout) :: a(:,:,:)
        integer, intent(in) :: lo1, hi1, lo2, hi2, lo3, hi3
        integer, intent(in) :: dev_id
        integer :: n1, n2, n3
        integer (kind=c_size_t) :: sz
        type(c_devptr) :: cp
        real(kind=8), pointer, device :: fp(:,:,:)
        integer :: istat
        n1 = max(hi1-lo1+1, 1)
        n2 = max(hi2-lo2+1, 1)
        n3 = max(hi3-lo3+1, 1)
        sz = int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t)
        istat = cudaSetDevice(dev_id)
        istat = cudaMalloc(cp, szr*sz)
        call c_f_pointer(cp, fp, shape=(/n1,n2,n3/))
        call shift_bound_d3(fp, lo1, lo2, lo3, a)
    contains
        subroutine shift_bound_d3(fp, lo1, lo2, lo3, a)
            integer, intent(in) :: lo1, lo2, lo3
            real(kind=8), target,  device, intent(in) :: fp(lo1:,lo2:,lo3:)
            real(kind=8), pointer, device, intent(inout) :: a(:,:,:)
            a => fp
        end subroutine shift_bound_d3
    end subroutine gpu_allocate_r3

    subroutine gpu_deallocate_r1(a)
        use cudafor
        real(kind=8), pointer, device, intent(inout) :: a(:)
        integer :: istat
        istat = cudaFree(a)
        if (istat .ne. cudaSuccess) then
            write(*,*) 'CUDA  error:', cudaGetErrorString(cudaGetLastError())
            stop
        endif
    end subroutine gpu_deallocate_r1

    subroutine gpu_deallocate_r2(a)
        use cudafor
        real(kind=8), pointer, device, intent(inout) :: a(:,:)
        integer :: istat
        istat = cudaFree(a)
        if (istat .ne. cudaSuccess) then
            write(*,*) 'CUDA  error:', cudaGetErrorString(cudaGetLastError())
            stop
        endif
    end subroutine gpu_deallocate_r2

    subroutine gpu_deallocate_r3(a)
        use cudafor
        real(kind=8), pointer, device, intent(inout) :: a(:,:,:)
        integer :: istat
        istat = cudaFree(a)
        if (istat .ne. cudaSuccess) then
            write(*,*) 'CUDA  error:', cudaGetErrorString(cudaGetLastError())
            stop
        endif
    end subroutine gpu_deallocate_r3



end module memory_module
