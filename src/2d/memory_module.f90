module memory_module

    use iso_c_binding

    implicit none

#if (CLAW_REAL == 8)
    integer, parameter, private :: szr = 8
#else
    integer, parameter, private :: szr = 4
#endif
    integer (kind=c_size_t), parameter, private :: szi = 4_c_size_t


    interface gpu_allocate
        module procedure gpu_allocate_r1
        module procedure gpu_allocate_r2
        module procedure gpu_allocate_r3
        module procedure gpu_allocate_r4
        module procedure gpu_allocate_i1
        module procedure gpu_allocate_i2
        module procedure gpu_allocate_i3
    end interface
    interface gpu_deallocate
        module procedure gpu_deallocate_r1
        module procedure gpu_deallocate_r2
        module procedure gpu_deallocate_r3
        module procedure gpu_deallocate_r4
        module procedure gpu_deallocate_i1
        module procedure gpu_deallocate_i2
        module procedure gpu_deallocate_i3
    end interface

    interface cpu_allocate_pinned
        module procedure cpu_allocate_pinned_r1
        module procedure cpu_allocate_pinned_r2
        module procedure cpu_allocate_pinned_r3
        module procedure cpu_allocate_pinned_r4
        module procedure cpu_allocate_pinned_i1
        module procedure cpu_allocate_pinned_i2
        module procedure cpu_allocate_pinned_i3
    end interface cpu_allocate_pinned

    interface cpu_deallocate_pinned
        module procedure cpu_deallocate_pinned_r1
        module procedure cpu_deallocate_pinned_r2
        module procedure cpu_deallocate_pinned_r3
        module procedure cpu_deallocate_pinned_r4
        module procedure cpu_deallocate_pinned_i1
        module procedure cpu_deallocate_pinned_i2
        module procedure cpu_deallocate_pinned_i3
    end interface cpu_deallocate_pinned

    interface 
        subroutine clawpack_mempool_init() bind(C, name='clawpack_mempool_init')
            use, intrinsic :: iso_c_binding
            implicit none
        end subroutine

        function clawpack_mempool_alloc_gpu (nbytes, dev_id) result(p) bind(C, name='clawpack_mempool_alloc_gpu')
            use cudafor
            use, intrinsic :: iso_c_binding
            type(c_devptr) :: p
            integer(kind=c_size_t), intent(in), value :: nbytes
            integer(kind=c_int), intent(in), value :: dev_id
        end function clawpack_mempool_alloc_gpu

        subroutine clawpack_mempool_free_gpu (p, dev_id) bind(C, name='clawpack_mempool_free_gpu')
            use cudafor
            use, intrinsic :: iso_c_binding
            type(c_devptr), value :: p
            integer(kind=c_int), intent(in), value :: dev_id
        end subroutine clawpack_mempool_free_gpu

        function clawpack_mempool_alloc_pinned (nbytes) result(p) bind(C, name='clawpack_mempool_alloc_pinned')
            use, intrinsic :: iso_c_binding
            type(c_ptr) :: p
            integer(kind=c_size_t), intent(in), value :: nbytes
        end function clawpack_mempool_alloc_pinned

        subroutine clawpack_mempool_free_pinned (p) bind(C, name='clawpack_mempool_free_pinned')
            use, intrinsic :: iso_c_binding
            type(c_ptr), value :: p
        end subroutine clawpack_mempool_free_pinned

        subroutine clawpack_real_array_init (p, n) bind(C, name='clawpack_real_array_init')
            use, intrinsic :: iso_c_binding
            type(c_ptr), value :: p
            integer(kind=c_size_t), intent(in), value :: n
        end subroutine clawpack_real_array_init

    end interface

contains
    subroutine gpu_allocate_r1(a, dev_id, lo1, hi1)
        use cudafor
        real(CLAW_REAL), pointer, device, intent(inout) :: a(:)
        integer, intent(in) :: lo1, hi1
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: n1
        integer (kind=c_size_t) :: sz
        type(c_devptr) :: cp
        real(CLAW_REAL), pointer, device :: fp(:)
        integer :: istat

        n1 = max(hi1-lo1+1, 1)
        sz = int(n1,c_size_t)
        dev_id_c = int(dev_id,c_int)
        cp = clawpack_mempool_alloc_gpu(szr*sz, dev_id_c)

        call c_f_pointer(cp, fp, shape=(/n1/))
        call shift_bound_d1(fp, lo1, a)
    contains
        subroutine shift_bound_d1(fp, lo1, a)
            integer, intent(in) :: lo1
            real(CLAW_REAL), target,  device, intent(in) :: fp(lo1:)
            real(CLAW_REAL), pointer, device, intent(inout) :: a(:)
            a => fp
        end subroutine shift_bound_d1
    end subroutine gpu_allocate_r1

    subroutine gpu_allocate_r2(a, dev_id, lo1, hi1, lo2, hi2)
        use cudafor
        real(CLAW_REAL), pointer, device, intent(inout) :: a(:,:)
        integer, intent(in) :: lo1, hi1, lo2, hi2
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: n1, n2
        integer (kind=c_size_t) :: sz
        type(c_devptr) :: cp
        real(CLAW_REAL), pointer, device :: fp(:,:)
        integer :: istat

        n1 = max(hi1-lo1+1, 1)
        n2 = max(hi2-lo2+1, 1)
        sz = int(n1,c_size_t) * int(n2,c_size_t)
        dev_id_c = int(dev_id,c_int)
        cp = clawpack_mempool_alloc_gpu(szr*sz, dev_id_c)

        call c_f_pointer(cp, fp, shape=(/n1,n2/))
        call shift_bound_d2(fp, lo1, lo2, a)
    contains
        subroutine shift_bound_d2(fp, lo1, lo2, a)
            integer, intent(in) :: lo1, lo2
            real(CLAW_REAL), target,  device, intent(in) :: fp(lo1:,lo2:)
            real(CLAW_REAL), pointer, device, intent(inout) :: a(:,:)
            a => fp
        end subroutine shift_bound_d2
    end subroutine gpu_allocate_r2

    subroutine gpu_allocate_r3(a, dev_id, lo1, hi1, lo2, hi2, lo3, hi3)
        use cudafor
        real(CLAW_REAL), pointer, device, intent(inout) :: a(:,:,:)
        integer, intent(in) :: lo1, hi1, lo2, hi2, lo3, hi3
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: n1, n2, n3
        integer (kind=c_size_t) :: sz
        type(c_devptr) :: cp
        real(CLAW_REAL), pointer, device :: fp(:,:,:)
        integer :: istat

        n1 = max(hi1-lo1+1, 1)
        n2 = max(hi2-lo2+1, 1)
        n3 = max(hi3-lo3+1, 1)
        sz = int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t)
        dev_id_c = int(dev_id,c_int)
        cp = clawpack_mempool_alloc_gpu(szr*sz, dev_id_c)

        call c_f_pointer(cp, fp, shape=(/n1,n2,n3/))
        call shift_bound_d3(fp, lo1, lo2, lo3, a)
    contains
        subroutine shift_bound_d3(fp, lo1, lo2, lo3, a)
            integer, intent(in) :: lo1, lo2, lo3
            real(CLAW_REAL), target,  device, intent(in) :: fp(lo1:,lo2:,lo3:)
            real(CLAW_REAL), pointer, device, intent(inout) :: a(:,:,:)
            a => fp
        end subroutine shift_bound_d3
    end subroutine gpu_allocate_r3

    subroutine gpu_allocate_r4(a, dev_id, lo1, hi1, lo2, hi2, lo3, hi3, lo4, hi4)
        use cudafor
        real(CLAW_REAL), pointer, device, intent(inout) :: a(:,:,:,:)
        integer, intent(in) :: lo1, hi1, lo2, hi2, lo3, hi3, lo4, hi4
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: n1, n2, n3, n4
        integer (kind=c_size_t) :: sz
        type(c_devptr) :: cp
        real(CLAW_REAL), pointer, device :: fp(:,:,:,:)
        integer :: istat

        n1 = max(hi1-lo1+1, 1)
        n2 = max(hi2-lo2+1, 1)
        n3 = max(hi3-lo3+1, 1)
        n4 = max(hi4-lo4+1, 1)
        sz = int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t) * int(n4,c_size_t)
        dev_id_c = int(dev_id,c_int)
        cp = clawpack_mempool_alloc_gpu(szr*sz, dev_id_c)

        call c_f_pointer(cp, fp, shape=(/n1,n2,n3,n4/))
        call shift_bound_d4(fp, lo1, lo2, lo3, lo4, a)
    contains
        subroutine shift_bound_d4(fp, lo1, lo2, lo3, lo4, a)
            integer, intent(in) :: lo1, lo2, lo3, lo4
            real(CLAW_REAL), target,  device, intent(in) :: fp(lo1:,lo2:,lo3:,lo4:)
            real(CLAW_REAL), pointer, device, intent(inout) :: a(:,:,:,:)
            a => fp
        end subroutine shift_bound_d4
    end subroutine gpu_allocate_r4

    subroutine gpu_allocate_i1(a, dev_id, lo1, hi1)
        use cudafor
        integer, pointer, device, intent(inout) :: a(:)
        integer, intent(in) :: lo1, hi1
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: n1
        integer (kind=c_size_t) :: sz
        type(c_devptr) :: cp
        integer, device, pointer :: fp(:)
        n1 = max(hi1-lo1+1, 1)
        sz = szi * int(n1,c_size_t)
        dev_id_c = int(dev_id,c_int)
        cp = clawpack_mempool_alloc_gpu(sz, dev_id_c)
        call c_f_pointer(cp, fp, shape=(/n1/))
        call shift_bound_i1(fp, lo1, a)
    contains
        subroutine shift_bound_i1(fp, lo1, a)
            integer, intent(in) :: lo1
            integer, target, device, intent(in) :: fp(lo1:)
            integer, pointer, device, intent(inout) :: a(:)
            a => fp
        end subroutine shift_bound_i1
    end subroutine gpu_allocate_i1

    subroutine gpu_allocate_i2(a, dev_id, lo1, hi1, lo2, hi2)
        use cudafor
        integer, pointer, device, intent(inout) :: a(:,:)
        integer, intent(in) :: lo1, hi1, lo2, hi2
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: n1, n2
        integer (kind=c_size_t) :: sz
        type(c_devptr) :: cp
        integer, device, pointer :: fp(:,:)
        n1 = max(hi1-lo1+1, 1)
        n2 = max(hi2-lo2+1, 1)
        sz = szi * int(n1,c_size_t) * int(n2,c_size_t)
        dev_id_c = int(dev_id,c_int)
        cp = clawpack_mempool_alloc_gpu(sz, dev_id_c)
        call c_f_pointer(cp, fp, shape=(/n1,n2/))
        call shift_bound_i2(fp, lo1, lo2, a)
    contains
        subroutine shift_bound_i2(fp, lo1, lo2, a)
            integer, intent(in) :: lo1, lo2
            integer, target, device, intent(in) :: fp(lo1:, lo2:)
            integer, pointer, device, intent(inout) :: a(:,:)
            a => fp
        end subroutine shift_bound_i2
    end subroutine gpu_allocate_i2

    subroutine gpu_allocate_i3(a, dev_id, lo1, hi1, lo2, hi2, lo3, hi3)
        use cudafor
        integer, pointer, device, intent(inout) :: a(:,:,:)
        integer, intent(in) :: lo1, hi1, lo2, hi2, lo3, hi3
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: n1, n2, n3
        integer (kind=c_size_t) :: sz
        type(c_devptr) :: cp
        integer, device, pointer :: fp(:,:,:)
        n1 = max(hi1-lo1+1, 1)
        n2 = max(hi2-lo2+1, 1)
        n3 = max(hi3-lo3+1, 1)
        sz = szi * int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t)
        dev_id_c = int(dev_id,c_int)
        cp = clawpack_mempool_alloc_gpu(sz, dev_id_c)
        call c_f_pointer(cp, fp, shape=(/n1,n2,n3/))
        call shift_bound_i3(fp, lo1, lo2, lo3, a)
    contains
        subroutine shift_bound_i3(fp, lo1, lo2, lo3, a)
            integer, intent(in) :: lo1, lo2, lo3
            integer, target, device, intent(in) :: fp(lo1:, lo2:, lo3:)
            integer, pointer, device, intent(inout) :: a(:,:,:)
            a => fp
        end subroutine shift_bound_i3
    end subroutine gpu_allocate_i3


    subroutine gpu_deallocate_r1(a, dev_id)
        use cudafor
        real(CLAW_REAL), pointer, device, intent(inout) :: a(:)
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: lo(1)
        type(c_devptr) :: cp
        lo = lbound(a)
        cp = c_devloc(a(lo(1)))
        dev_id_c = int(dev_id,c_int)
        call clawpack_mempool_free_gpu(cp, dev_id_c)
        a => Null()
    end subroutine gpu_deallocate_r1

    subroutine gpu_deallocate_r2(a, dev_id)
        use cudafor
        real(CLAW_REAL), pointer, device, intent(inout) :: a(:,:)
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: lo(2)
        type(c_devptr) :: cp
        lo = lbound(a)
        cp = c_devloc(a(lo(1),lo(2)))
        dev_id_c = int(dev_id,c_int)
        call clawpack_mempool_free_gpu(cp, dev_id_c)
        a => Null()
    end subroutine gpu_deallocate_r2

    subroutine gpu_deallocate_r3(a, dev_id)
        use cudafor
        real(CLAW_REAL), pointer, device, intent(inout) :: a(:,:,:)
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: lo(3)
        type(c_devptr) :: cp
        lo = lbound(a)
        cp = c_devloc(a(lo(1),lo(2),lo(3)))
        dev_id_c = int(dev_id,c_int)
        call clawpack_mempool_free_gpu(cp, dev_id_c)
        a => Null()
    end subroutine gpu_deallocate_r3

    subroutine gpu_deallocate_r4(a, dev_id)
        use cudafor
        real(CLAW_REAL), pointer, device, intent(inout) :: a(:,:,:,:)
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: lo(4)
        type(c_devptr) :: cp
        lo = lbound(a)
        cp = c_devloc(a(lo(1),lo(2),lo(3),lo(4)))
        dev_id_c = int(dev_id,c_int)
        call clawpack_mempool_free_gpu(cp, dev_id_c)
        a => Null()
    end subroutine gpu_deallocate_r4

    subroutine gpu_deallocate_i1(a, dev_id)
        use cudafor
        integer, pointer, device, intent(inout) :: a(:)
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
        integer :: lo(1)
        type(c_devptr) :: cp
        lo = lbound(a)
        cp = c_devloc(a(lo(1)))
        dev_id_c = int(dev_id,c_int)
        call clawpack_mempool_free_gpu(cp, dev_id_c)
        a => Null()
    end subroutine gpu_deallocate_i1

    subroutine gpu_deallocate_i2(a, dev_id)
	use cudafor
	integer, pointer, device, intent(inout) :: a(:,:)
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
	integer :: lo(2)
	type(c_devptr) :: cp
	lo = lbound(a)
	cp = c_devloc(a(lo(1),lo(2)))
        dev_id_c = int(dev_id,c_int)
	call clawpack_mempool_free_gpu(cp, dev_id_c)
	a => Null()
    end subroutine gpu_deallocate_i2

    subroutine gpu_deallocate_i3(a, dev_id)
	use cudafor
	integer, pointer, device, intent(inout) :: a(:,:,:)
        integer, intent(in) :: dev_id
        integer (kind=c_int) :: dev_id_c
	integer :: lo(3)
	type(c_devptr) :: cp
	lo = lbound(a)
	cp = c_devloc(a(lo(1),lo(2),lo(3)))
        dev_id_c = int(dev_id,c_int)
	call clawpack_mempool_free_gpu(cp, dev_id_c)
	a => Null()
    end subroutine gpu_deallocate_i3


! ####################### allocating real and int arrays on CPU ################## !

    subroutine cpu_allocate_pinned_r1(a, lo1, hi1)
        real(CLAW_REAL), pointer, intent(inout) :: a(:)
        integer, intent(in) :: lo1, hi1
        integer :: n1
        integer (kind=c_size_t) :: sz
        type(c_ptr) :: cp
        real(CLAW_REAL), pointer :: fp(:)
        n1 = max(hi1-lo1+1, 1)
        sz = int(n1,c_size_t)
        cp = clawpack_mempool_alloc_pinned(szr*sz)
        call clawpack_real_array_init(cp, sz)
        call c_f_pointer(cp, fp, shape=(/n1/))
        call shift_bound_d1(fp, lo1, a)
    contains
        subroutine shift_bound_d1(fp, lo1, a)
            integer, intent(in) :: lo1
            real(CLAW_REAL), target, intent(in) :: fp(lo1:)
            real(CLAW_REAL), pointer, intent(inout) :: a(:)
            a => fp
        end subroutine shift_bound_d1
    end subroutine cpu_allocate_pinned_r1

    subroutine cpu_allocate_pinned_r2(a, lo1, hi1, lo2, hi2)
        real(CLAW_REAL), pointer, intent(inout) :: a(:,:)
        integer, intent(in) :: lo1, hi1, lo2, hi2
        integer :: n1, n2
        integer (kind=c_size_t) :: sz
        type(c_ptr) :: cp
        real(CLAW_REAL), pointer :: fp(:,:)
        n1 = max(hi1-lo1+1, 1)
        n2 = max(hi2-lo2+1, 1)
        sz = int(n1,c_size_t) * int(n2,c_size_t) 
        cp = clawpack_mempool_alloc_pinned(szr*sz)
        call clawpack_real_array_init(cp, sz)
        call c_f_pointer(cp, fp, shape=(/n1,n2/))
        call shift_bound_d2(fp, lo1, lo2, a)
    contains
        subroutine shift_bound_d2(fp, lo1, lo2, a)
            integer, intent(in) :: lo1, lo2
            real(CLAW_REAL), target, intent(in) :: fp(lo1:,lo2:)
            real(CLAW_REAL), pointer, intent(inout) :: a(:,:)
            a => fp
        end subroutine shift_bound_d2
    end subroutine cpu_allocate_pinned_r2

    subroutine cpu_allocate_pinned_r3(a, lo1, hi1, lo2, hi2, lo3, hi3)
        real(CLAW_REAL), pointer, intent(inout) :: a(:,:,:)
        integer, intent(in) :: lo1, hi1, lo2, hi2, lo3, hi3
        integer :: n1, n2, n3
        integer (kind=c_size_t) :: sz
        type(c_ptr) :: cp
        real(CLAW_REAL), pointer :: fp(:,:,:)
        n1 = max(hi1-lo1+1, 1)
        n2 = max(hi2-lo2+1, 1)
        n3 = max(hi3-lo3+1, 1)
        sz = int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t)
        cp = clawpack_mempool_alloc_pinned(szr*sz)
        call clawpack_real_array_init(cp, sz)
        call c_f_pointer(cp, fp, shape=(/n1,n2,n3/))
        call shift_bound_d3(fp, lo1, lo2, lo3, a)
    contains
        subroutine shift_bound_d3(fp, lo1, lo2, lo3, a)
            integer, intent(in) :: lo1, lo2, lo3
            real(CLAW_REAL), target, intent(in) :: fp(lo1:,lo2:,lo3:)
            real(CLAW_REAL), pointer, intent(inout) :: a(:,:,:)
            a => fp
        end subroutine shift_bound_d3
    end subroutine cpu_allocate_pinned_r3

    subroutine cpu_allocate_pinned_r4(a, lo1, hi1, lo2, hi2, lo3, hi3, lo4, hi4)
        real(CLAW_REAL), pointer, intent(inout) :: a(:,:,:,:)
        integer, intent(in) :: lo1, hi1, lo2, hi2, lo3, hi3, lo4, hi4
        integer :: n1, n2, n3, n4
        integer (kind=c_size_t) :: sz
        type(c_ptr) :: cp
        real(CLAW_REAL), pointer :: fp(:,:,:,:)
        n1 = max(hi1-lo1+1, 1)
        n2 = max(hi2-lo2+1, 1)
        n3 = max(hi3-lo3+1, 1)
        n4 = max(hi4-lo4+1, 1)
        sz = int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t) * int(n4,c_size_t)
        cp = clawpack_mempool_alloc_pinned(szr*sz)
        call clawpack_real_array_init(cp, sz)
        call c_f_pointer(cp, fp, shape=(/n1,n2,n3,n4/))
        call shift_bound_d4(fp, lo1, lo2, lo3, lo4, a)
    contains
        subroutine shift_bound_d4(fp, lo1, lo2, lo3, lo4, a)
            integer, intent(in) :: lo1, lo2, lo3, lo4
            real(CLAW_REAL), target, intent(in) :: fp(lo1:,lo2:,lo3:,lo4:)
            real(CLAW_REAL), pointer, intent(inout) :: a(:,:,:,:)
            a => fp
        end subroutine shift_bound_d4
    end subroutine cpu_allocate_pinned_r4

    subroutine cpu_allocate_pinned_i1(a, lo1, hi1)
        integer, pointer, intent(inout) :: a(:)
        integer, intent(in) :: lo1, hi1
        integer :: n1
        integer (kind=c_size_t) :: sz
        type(c_ptr) :: cp
        integer, pointer :: fp(:)
        n1 = max(hi1-lo1+1, 1)
        sz = szi * int(n1,c_size_t)
        cp = clawpack_mempool_alloc_pinned(sz)
        call c_f_pointer(cp, fp, shape=(/n1/))
        call shift_bound_i1(fp, lo1, a)
    contains
        subroutine shift_bound_i1(fp, lo1, a)
            integer, intent(in) :: lo1
            integer, target, intent(in) :: fp(lo1:)
            integer, pointer, intent(inout) :: a(:)
            a => fp
        end subroutine shift_bound_i1
    end subroutine cpu_allocate_pinned_i1


    subroutine cpu_allocate_pinned_i2(a, lo1, hi1, lo2, hi2)
	integer, pointer, intent(inout) :: a(:,:)
	integer, intent(in) :: lo1, hi1, lo2, hi2
	integer :: n1, n2
	integer (kind=c_size_t) :: sz
	type(c_ptr) :: cp
	integer, pointer :: fp(:,:)
	n1 = max(hi1-lo1+1, 1)
	n2 = max(hi2-lo2+1, 1)
	sz = szi * int(n1,c_size_t) * int(n2,c_size_t)
	cp = clawpack_mempool_alloc_pinned(sz)
	call c_f_pointer(cp, fp, shape=(/n1,n2/))
	call shift_bound_i2(fp, lo1, lo2, a)
    contains
	subroutine shift_bound_i2(fp, lo1, lo2, a)
	    integer, intent(in) :: lo1, lo2
	    integer, target, intent(in) :: fp(lo1:,lo2:)
	    integer, pointer, intent(inout) :: a(:,:)
	    a => fp
	end subroutine shift_bound_i2
    end subroutine cpu_allocate_pinned_i2

    subroutine cpu_allocate_pinned_i3(a, lo1, hi1, lo2, hi2, lo3, hi3)
	integer, pointer, intent(inout) :: a(:,:,:)
	integer, intent(in) :: lo1, hi1, lo2, hi2, lo3, hi3
	integer :: n1, n2, n3
	integer (kind=c_size_t) :: sz
	type(c_ptr) :: cp
	integer, pointer :: fp(:,:,:)
	n1 = max(hi1-lo1+1, 1)
	n2 = max(hi2-lo2+1, 1)
	n3 = max(hi3-lo3+1, 1)
	sz = szi * int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t)
	cp = clawpack_mempool_alloc_pinned(sz)
	call c_f_pointer(cp, fp, shape=(/n1,n2,n3/))
	call shift_bound_i3(fp, lo1, lo2, lo3, a)
    contains
	subroutine shift_bound_i3(fp, lo1, lo2, lo3, a)
	    integer, intent(in) :: lo1, lo2, lo3
	    integer, target, intent(in) :: fp(lo1:,lo2:,lo3:)
	    integer, pointer, intent(inout) :: a(:,:,:)
	    a => fp
	end subroutine shift_bound_i3
    end subroutine cpu_allocate_pinned_i3

    subroutine cpu_deallocate_pinned_r1(a)
        real(CLAW_REAL), pointer, intent(inout) :: a(:)
        integer :: lo(1)
        type(c_ptr) :: cp
        lo = lbound(a)
        cp = c_loc(a(lo(1)))
        call clawpack_mempool_free_pinned(cp)
        a => Null()
    end subroutine cpu_deallocate_pinned_r1

    subroutine cpu_deallocate_pinned_r2(a)
        real(CLAW_REAL), pointer, intent(inout) :: a(:,:)
        integer :: lo(2)
        type(c_ptr) :: cp
        lo = lbound(a)
        cp = c_loc(a(lo(1),lo(2)))
        call clawpack_mempool_free_pinned(cp)
        a => Null()
    end subroutine cpu_deallocate_pinned_r2

    subroutine cpu_deallocate_pinned_r3(a)
        real(CLAW_REAL), pointer, intent(inout) :: a(:,:,:)
        integer :: lo(3)
        type(c_ptr) :: cp
        lo = lbound(a)
        cp = c_loc(a(lo(1),lo(2),lo(3)))
        call clawpack_mempool_free_pinned(cp)
        a => Null()
    end subroutine cpu_deallocate_pinned_r3

    subroutine cpu_deallocate_pinned_r4(a)
        real(CLAW_REAL), pointer, intent(inout) :: a(:,:,:,:)
        integer :: lo(4)
        type(c_ptr) :: cp
        lo = lbound(a)
        cp = c_loc(a(lo(1),lo(2),lo(3),lo(4)))
        call clawpack_mempool_free_pinned(cp)
        a => Null()
    end subroutine cpu_deallocate_pinned_r4

    subroutine cpu_deallocate_pinned_i1(a)
        integer, pointer, intent(inout) :: a(:)
        integer :: lo(1)
        type(c_ptr) :: cp
        lo = lbound(a)
        cp = c_loc(a(lo(1)))
        call clawpack_mempool_free_pinned(cp)
        a => Null()
    end subroutine cpu_deallocate_pinned_i1

    subroutine cpu_deallocate_pinned_i2(a)
	integer, pointer, intent(inout) :: a(:,:)
	integer :: lo(2)
	type(c_ptr) :: cp
	lo = lbound(a)
	cp = c_loc(a(lo(1),lo(2)))
	call clawpack_mempool_free_pinned(cp)
	a => Null()
    end subroutine cpu_deallocate_pinned_i2

    subroutine cpu_deallocate_pinned_i3(a)
	integer, pointer, intent(inout) :: a(:,:,:)
	integer :: lo(3)
	type(c_ptr) :: cp
	lo = lbound(a)
	cp = c_loc(a(lo(1),lo(2),lo(3)))
	call clawpack_mempool_free_pinned(cp)
	a => Null()
    end subroutine cpu_deallocate_pinned_i3



end module memory_module
