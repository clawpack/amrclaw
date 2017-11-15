module cuda_module

    use cudafor, only: cuda_stream_kind, cudaEvent, dim3
    use amr_module, only: inunit
    use memory_module, only: clawpack_mempool_init

    implicit none
    save

    integer, parameter :: max_cuda_streams = 20
    integer, parameter :: max_num_devices = 16
    ! Note that CUDA enumerates devices starting from 0
    integer(kind=cuda_stream_kind) :: cuda_streams(max_cuda_streams, 0:max_num_devices-1)

    integer :: num_devices
    integer :: device_id ! which GPU is used

    ! For timing
    ! use 0-index array
    integer, parameter :: max_cuda_timer =  100
    real(kind=8) :: elapsed_time(max_cuda_timer)
    character(len=20) :: timer_name(max_cuda_timer)
    logical :: timer_initialized(max_cuda_timer)
    integer :: n_calls(max_cuda_timer)
    type(cudaEvent) :: event_start(max_cuda_timer)
    type(cudaEvent) :: event_stop(max_cuda_timer)
    integer :: n_timer

    ! current number of cuda threads and cuda blocks
    type(dim3) :: numBlocks, numThreads

    character(len=*), parameter :: gpufile = 'gpu.data'

    interface toString
        module procedure toString1
        module procedure toString2
    end interface toString

    interface write_grid
        module procedure write_grid2
    end interface write_grid

#ifdef CUDA
    type grid2d
        sequence ! force the derived type to be stored contiguously
        double precision, dimension(:,:,:), pointer, contiguous :: dataptr=>null()
    end type grid2d
    type grid2d_device
        sequence ! force the derived type to be stored contiguously
        double precision, dimension(:,:,:), pointer, contiguous, device :: dataptr=>null()
    end type grid2d_device
#endif

contains

    subroutine initialize_cuda() 
        use cudafor, only: cudaStreamCreate, cudaGetDeviceCount, cudaSetDevice &
            , cudaGetDeviceCount, cudaGetDeviceProperties
        use cudafor, only: cudaDeviceProp
        implicit none
        integer :: i, p, cudaResult
        type(cudaDeviceProp) :: prop


        print *, 'Initialize CUDA environment ...'

        ! check how many devices are available
        cudaResult = cudaGetDeviceCount(num_devices) 
        print *, 'Found ', num_devices, ' GPUs.'
        if (num_devices > max_num_devices) then
            print *, "Found too many devices. Can't handle this so far. "
            stop
        endif

        do p = 0, num_devices-1
            cudaResult = cudaGetDeviceProperties(prop, p)
            write(*,"(' Device Number: ',i0)") p
            write(*,"('   Device name: ',a)") trim(prop%name)
            write(*,"('   Memory Clock Rate (KHz): ', i0)") &
              prop%memoryClockRate
            write(*,"('   Memory Bus Width (bits): ', i0)") &
              prop%memoryBusWidth
            write(*,"('   Peak Memory Bandwidth (GB/s): ', f6.2)") &
              2.0*prop%memoryClockRate*(prop%memoryBusWidth/8)/10.0**6
            write(*,*)        
        enddo

        do p = 0, num_devices-1
            cudaResult = cudaSetDevice(p)
            call check_cuda_error(cudaResult)
            do i = 1, max_cuda_streams
               cudaResult = cudaStreamCreate(cuda_streams(i,p))
                call check_cuda_error(cudaResult)
            enddo
        enddo

        call opendatafile(inunit,gpufile)
        read(inunit,*) device_id    ! which gpu will be used
        cudaResult = cudaSetDevice(device_id)
        call check_cuda_error(cudaResult)
        print *, 'Use the GPU with id: ',device_id

        ! Initialize memory pool
        call clawpack_mempool_init()

        n_timer = 0

    end subroutine initialize_cuda

    subroutine check_cuda_error(istat) 
        use cudafor
        integer :: istat
        if (istat /= cudaSuccess) then
            write(*,*) 'CUDA  error:', cudaGetErrorString(cudaGetLastError())
            stop
        endif
    end subroutine check_cuda_error

    integer function stream_from_index(idx)
    	implicit none
        integer :: idx

        ! note that available streams are indexed from 1 to 100
        ! reserve the stream 1 to 10 for special purposes
        if (idx < 0 .and. idx >= -10) then
            stream_from_index = -idx
        else
            ! stream_from_index below ranges from 11 to max_cuda-streams
            stream_from_index = MOD(idx, max_cuda_streams-10) + 11
        endif

    end function stream_from_index

    ! return a cuda stream based on an index and device id
    function get_cuda_stream(idx, dev_id) result(s)

        implicit none

        integer :: idx, dev_id
        integer(kind=cuda_stream_kind) :: s

        s = cuda_streams(stream_from_index(idx), dev_id)

    end function get_cuda_stream

    subroutine threads_and_blocks(lo, hi, numBlocks, numThreads)
	use cudafor, only: dim3
	implicit none

	integer, intent(in)       :: lo(SPACEDIM), hi(SPACEDIM)
	type(dim3), intent(inout) :: numBlocks, numThreads

	integer :: tile_size(SPACEDIM)

	tile_size = hi - lo + 1

	if (SPACEDIM .eq. 1) then

	    numThreads % x = 256
	    numThreads % y = 1
	    numThreads % z = 1

	    numBlocks % x = (tile_size(1) + numThreads % x - 1) / numThreads % x
	    numBlocks % y = 1
	    numBlocks % z = 1

	else if (SPACEDIM .eq. 2) then

	    numThreads % x = 16
	    numThreads % y = 16
	    numThreads % z = 1

	    numBlocks % x = (tile_size(1) + numThreads % x - 1) / numThreads % x
	    numBlocks % y = (tile_size(2) + numThreads % y - 1) / numThreads % y
	    numBlocks % z = 1

	else

	    numThreads % x = 8
	    numThreads % y = 8
	    numThreads % z = 8

	    numBlocks % x = (tile_size(1) + numThreads % x - 1) / numThreads % x
	    numBlocks % y = (tile_size(2) + numThreads % y - 1) / numThreads % y
	    numBlocks % z = (tile_size(3) + numThreads % z - 1) / numThreads % z

	endif
    end subroutine threads_and_blocks

    subroutine copy_cpu_to_gpu_async(p_d, p_h, sz, idx, dev_id)

        use cudafor
        use iso_c_binding, only: c_ptr, c_size_t

        implicit none

        type(c_devptr), value :: p_d
        type(c_ptr), value :: p_h
        integer(c_size_t) :: sz
        integer :: idx, dev_id

        integer :: s
        integer :: cudaResult

        cudaResult = cudaSetDevice(dev_id)
        call check_cuda_error(cudaResult)
        s = stream_from_index(idx)
        cudaResult = cudaMemcpyAsync(p_d, p_h, sz, cudaMemcpyHostToDevice, cuda_streams(s,dev_id))

    end subroutine copy_cpu_to_gpu_async



    subroutine copy_gpu_to_cpu_async(p_h, p_d, sz, idx, dev_id)

        use cudafor
        use iso_c_binding, only: c_ptr, c_size_t

        implicit none

        type(c_ptr), value :: p_h
        type(c_devptr), value :: p_d
        integer(c_size_t) :: sz
        integer :: idx, dev_id

        integer :: s
        integer :: cudaResult

        cudaResult = cudaSetDevice(dev_id)
        call check_cuda_error(cudaResult)
        s = stream_from_index(idx)
        cudaResult = cudaMemcpyAsync(p_h, p_d, sz, cudaMemcpyDeviceToHost, cuda_streams(s,dev_id))

    end subroutine copy_gpu_to_cpu_async

    subroutine copy_cpu_to_gpu(p_d, p_h, sz, dev_id)

        use cudafor
        use iso_c_binding, only: c_ptr, c_size_t

        implicit none

        type(c_devptr), value :: p_d
        type(c_ptr), value :: p_h
        integer(c_size_t) :: sz
        integer :: dev_id
        integer :: cudaResult


        cudaResult = cudaSetDevice(dev_id)
        call check_cuda_error(cudaResult)
        cudaResult = cudaMemcpy(p_d, p_h, sz, cudaMemcpyHostToDevice)

    end subroutine copy_cpu_to_gpu



    subroutine copy_gpu_to_cpu(p_h, p_d, sz, dev_id)

        use cudafor
        use iso_c_binding, only: c_ptr, c_size_t

        implicit none

        type(c_ptr), value :: p_h
        type(c_devptr), value :: p_d
        integer(c_size_t) :: sz
        integer :: dev_id
        integer :: cudaResult

        cudaResult = cudaSetDevice(dev_id)
        call check_cuda_error(cudaResult)
        cudaResult = cudaMemcpy(p_h, p_d, sz, cudaMemcpyDeviceToHost)

    end subroutine copy_gpu_to_cpu


    subroutine wait_for_all_gpu_tasks(dev_id)

        use cudafor

        implicit none

        integer :: cudaResult
        integer :: dev_id

        cudaResult = cudaSetDevice(dev_id)

#ifdef DEBUG
        call check_cuda_error(cudaResult)
#endif

        cudaResult = cudaDeviceSynchronize()

#ifdef DEBUG
        call check_cuda_error(cudaResult)
#endif

    end subroutine wait_for_all_gpu_tasks

    subroutine wait_for_stream(idx, dev_id) 

        use cudafor, only: cudaSetDevice, cudaStreamSynchronize

        implicit none

        integer :: cudaResult
        integer :: s, idx, dev_id

        s = stream_from_index(idx)
        cudaResult = cudaSetDevice(dev_id)
        cudaResult = cudaStreamSynchronize(cuda_streams(s, dev_id))

    end subroutine wait_for_stream

    ! put a timer with name t_name and ID id
    ! id should be in the range 1:max_cuda_timer
    subroutine timer_take(t_name, id)
        implicit none
        integer, intent(in) :: id
        character(len=10), intent(in) :: t_name
        if (timer_initialized(id) .eq. .false.) then
            timer_name(id) = t_name
            elapsed_time(id) = 0.0
            n_calls(id) = 0
            timer_initialized(id) = .true.
        endif
    end subroutine timer_take

    subroutine timer_start(id)
        use cudafor, only: cudaEventCreate, cudaEventRecord
        implicit none
        integer :: id, cuda_result
        cuda_result = cudaEventCreate(event_start(id))
        cuda_result = cudaEventCreate(event_stop(id))
        cuda_result = cudaEventRecord(event_start(id), 0)
        n_calls(id) = n_calls(id) + 1
    end subroutine timer_start

    subroutine timer_stop(id)
        use cudafor, only: cudaEventRecord, cudaEventDestroy, & 
            cudaEventSynchronize, cudaEventElapsedTime
        implicit none
        integer :: id, cuda_result
        real :: local_time
        cuda_result = cudaEventRecord(event_stop(id), 0)
        cuda_result = cudaEventSynchronize(event_stop(id))
        cuda_result = cudaEventElapsedTime(local_time, event_start(id), event_stop(id))
        cuda_result = cudaEventDestroy(event_start(id))
        cuda_result = cudaEventDestroy(event_stop(id))
        elapsed_time(id) = elapsed_time(id) + local_time/1000
    end subroutine timer_stop

    ! returned time is in second
    subroutine get_cuda_time(id, time) 
        implicit none
        integer, intent(in) :: id
        real(kind=8), intent(out) :: time
        time = elapsed_time(id)
    end subroutine get_cuda_time

    subroutine get_cuda_num_calls(id, n) 
        implicit none
        integer, intent(in) :: id
        integer, intent(out) :: n
        n= n_calls(id)
    end subroutine get_cuda_num_calls


    subroutine finalize_cuda() 
    end subroutine finalize_cuda

    function get_num_devices_used() result(nDev) bind(C, name='get_num_devices_used')
        use iso_c_binding, only: c_int
        implicit none
        integer(c_int) :: nDev
        nDev = int(num_devices, c_int)
    end function get_num_devices_used

    ! Do max reduction on array a_x 
    ! and put the maximum value in max_array(blockIdx%x, blockIdx%y) 
    attributes(device) &
    subroutine max_reduce_device_2d(a_s, mx, my, max_array, rmx, rmy)
	use cudafor

	implicit none

        ! we assume index of the entire grid is (1:mx, 1:my)
        ! TODO: add functionality to handle (lox:hix, loy:hiy)
	integer, value, intent(in) :: mx, my, rmx, rmy
	double precision, intent(out) :: max_array(rmx, rmy)

	! local
	integer :: i, j, tidx, tidy, sy
	double precision :: a_s(blockDim%x, blockDim%y)

	tidx = threadIdx%x
	tidy = threadIdx%y
	i = (blockIdx%x-1) * blockDim%x + threadIdx%x
	j = (blockIdx%y-1) * blockDim%y + threadIdx%y

	if (i > mx .or. j > my) return

	if (blockIdx%x * blockDim%x >= mx .or. &
	    blockIdx%y * blockDim%y >= my ) then ! Part of this block is outside the grid
	    sy = blockDim%y/2
	    ! I found that unrolling this actually gives worse performance
	    do while (sy > 0)
		if (tidy <= sy .and. (j+sy) <= my) then
		    a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx,tidy+sy))
		endif
		call syncthreads()
		sy = sy/2
	    enddo

	    ! reduce max of first row
	    if (tidy == 1) then
		if (i+8 <= mx) then
		    a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+8,tidy))
		    call syncthreads()
		endif
		if (i+4 <= mx) then
		    a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+4,tidy))
		    call syncthreads()
		endif
		if (i+2 <= mx) then
		    a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+2,tidy))
		    call syncthreads()
		endif
		if (i+1 <= mx) then
		    a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+1,tidy))
		    call syncthreads()
		endif
	    endif
	    if (tidx == 1 .and. tidy == 1) then 
		max_array(blockIdx%x,blockIdx%y) = a_s(1,1)
	    endif

	else ! The entire 16x16 block is mapped to a 16x16 subgrid
	    if (tidy <= 8) then
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx,tidy+8))
		call syncthreads()
	    endif

	    if (tidy <= 4) then
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx,tidy+4))
		call syncthreads()
	    endif

	    if (tidy <= 2) then
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx,tidy+2))
		call syncthreads()
	    endif

	    if (tidy <= 1) then
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx,tidy+1))
		call syncthreads()
	    endif

	    if (tidy == 1) then
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+8,tidy))
		call syncthreads()
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+4,tidy))
		call syncthreads()
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+2,tidy))
		call syncthreads()
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+1,tidy))
		call syncthreads()
	    endif


	    if (tidx == 1 .and. tidy == 1) then 
		max_array(blockIdx%x,blockIdx%y) = a_s(1,1)
	    endif
	endif

    end subroutine

    ! Do a max reduction on current block and write output to local_max
    attributes(device) &
    subroutine max_reduce_device_local_2d(a_s, mx, my, local_max)
	use cudafor

	implicit none

        ! we assume index of the entire grid is (1:mx, 1:my)
        ! TODO: add functionality to handle (lox:hix, loy:hiy)
	integer, value, intent(in) :: mx, my ! dimension of the entire grid
	double precision, intent(out) :: local_max

	! local
	integer :: i, j, tidx, tidy, sy
	double precision :: a_s(blockDim%x, blockDim%y)

	tidx = threadIdx%x
	tidy = threadIdx%y
	i = (blockIdx%x-1) * blockDim%x + threadIdx%x
	j = (blockIdx%y-1) * blockDim%y + threadIdx%y

	if (i > mx .or. j > my) return

	if (blockIdx%x * blockDim%x >= mx .or. &
	    blockIdx%y * blockDim%y >= my ) then ! Part of this block is outside the grid
	    sy = blockDim%y/2
	    ! I found that unrolling this actually gives worse performance
	    do while (sy > 0)
		if (tidy <= sy .and. (j+sy) <= my) then
		    a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx,tidy+sy))
		endif
		call syncthreads()
		sy = sy/2
	    enddo

	    ! reduce max of first row
	    if (tidy == 1) then
		if (i+8 <= mx) then
		    a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+8,tidy))
		    call syncthreads()
		endif
		if (i+4 <= mx) then
		    a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+4,tidy))
		    call syncthreads()
		endif
		if (i+2 <= mx) then
		    a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+2,tidy))
		    call syncthreads()
		endif
		if (i+1 <= mx) then
		    a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+1,tidy))
		    call syncthreads()
		endif
	    endif
	    if (tidx == 1 .and. tidy == 1) then 
		local_max = a_s(1,1)
	    endif

	else ! The entire 16x16 block is mapped to a 16x16 subgrid
	    if (tidy <= 8) then
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx,tidy+8))
		call syncthreads()
	    endif

	    if (tidy <= 4) then
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx,tidy+4))
		call syncthreads()
	    endif

	    if (tidy <= 2) then
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx,tidy+2))
		call syncthreads()
	    endif

	    if (tidy <= 1) then
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx,tidy+1))
		call syncthreads()
	    endif

	    if (tidy == 1) then
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+8,tidy))
		call syncthreads()
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+4,tidy))
		call syncthreads()
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+2,tidy))
		call syncthreads()
		a_s(tidx,tidy) = max(a_s(tidx,tidy) , a_s(tidx+1,tidy))
		call syncthreads()
	    endif


	    if (tidx == 1 .and. tidy == 1) then 
		local_max = a_s(1,1)
	    endif
	endif

    end subroutine


    ! Convert integer k to str
    function toString1(k) result(str)
        !   "Convert an integer to string."
        integer, intent(in) :: k
        character(len=20) :: str
        write (str, *) k
        str = adjustl(str)
    end function toString1

    ! Convert integer k to str
    ! Always output a tring of 'length' characters
    ! Fill empty space on the left of 'k' with '0's
    function toString2(k,length) result(str)
        !   "Convert an integer to string."
        integer, intent(in) :: k, length
        character(len=20) :: str
        character(len=8) :: fmt ! format descriptor

        fmt = '(I'//trim(toString1(length))//'.'//trim(toString1(length))//')' ! an integer of width 5 with zeros at the left
        write (str,fmt) k ! converting integer to string using a 'internal file'
        str = adjustl(str)
    end function toString2



    !> ################################################ !
    ! Write a 2d grid data, q, to file, fname
    ! Example:
    ! 
    ! do m = 1,meqn
    !     call write_grid( q(m,:,:), 1-mbc,mx+mbc, 1-mbc, my+mbc, 'q_'//trim(toString(m,3))//'.txt', frame_no)
    ! enddo
    !! ################################################ !
    subroutine write_grid2(q, lox, hix, loy, hiy, fname, iframe)
        integer, intent(in) :: lox, hix, loy, hiy, iframe
        real(kind=8), intent(in) :: q(lox:hix, loy:hiy)
        character(len=100), intent(in) :: fname
        integer :: i,j

        open (1, file=fname, position="append")
        write(1, *) "frame: ", iframe
        write(1, *) "i: ", lox, hix
        write(1, *) "j: ", loy, hiy
        do i = lox, hix
            do j = loy, hiy
                write(1, *) i, j, q(i,j)
            enddo
        enddo
        close(1)
    end subroutine write_grid2

    subroutine aos_to_soa_r2(soa, aos, nvar, xlo, xhi, ylo, yhi)
        implicit none
        integer, intent(in) :: nvar, xlo, xhi, ylo, yhi
        real(kind=8), intent(in) :: aos(1:nvar, xlo:xhi, ylo:yhi)
        real(kind=8), intent(inout) :: soa(xlo:xhi, ylo:yhi, 1:nvar)
        integer :: i,j,m

        do i = xlo, xhi
            do j = ylo, yhi
                do m = 1,nvar
                    soa(i,j,m) = aos(m,i,j)
                enddo
            enddo
        enddo
    end subroutine aos_to_soa_r2

    subroutine soa_to_aos_r2(aos, soa, nvar, xlo, xhi, ylo, yhi)
        implicit none
        integer, intent(in) :: nvar, xlo, xhi, ylo, yhi
        double precision, intent(inout) :: aos(1:nvar, xlo:xhi, ylo:yhi)
        double precision, intent(in) :: soa(xlo:xhi, ylo:yhi, 1:nvar)
        integer :: i,j,m
        do m = 1,nvar
            do j = ylo, yhi
                do i = xlo, xhi
                    aos(m,i,j) = soa(i,j,m) 
                enddo
            enddo
        enddo
    end subroutine soa_to_aos_r2

end module cuda_module
