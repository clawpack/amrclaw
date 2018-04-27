!> The module contains the definition of a "node descriptor" 
!! as well as other global variables used during the whole run.
!! 
!! The "node descriptor" is a tree data structure we use to keep track of 
!! all AMR grids.
!! Each grid in the grid hierarchy corresponds to a node in the tree. 
!! When a fine grid is nested in a coarse grid, its corresponding node in the tree 
!! is an offspring of the parent node corresponding to the coarse grid.
!! ![Tree data structure for grid management](./images/AMR_tree.png "amr_tree")
!!
!! Each node has a fixed number of items of information describing the grid,
!! stored and accessable as **node(property number, node number)** (which returns
!! a integer type as either data itself or a pointer to a memory location)
!! and **rnode(property number, node number)** (which returns a real type).
!! E.g. **node(nestlevel, 5)** stores AMR level of grid 5.
!! **rnode(cornxlo, 3)** stores x-coordinate of left border of grid 3.
!!
!! See explanation to each data members in this module 
!! for details of all properties of a node.
!!

#include "amr_macros.H"

module amr_module

    implicit none
       
    save
       
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::   data structure info.
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    integer, parameter :: rsize = 5
    integer, parameter :: nsize = 19

    !  :::::::   integer part of node descriptor
    !> node number (index) of next grid on the same level
    integer, parameter :: levelptr  = 1 

    !> temporary pointer
    integer, parameter :: tempptr   = 2

    ! TODO
    integer, parameter :: errptr    = 3

    !> AMR level of the grid
    integer, parameter :: nestlevel = 4

    !> Pointer to an 5 by **maxsp** array, which has boundary information 
    !! for this grid.
    !!
    !! ### About the Array
    !! The array is generally referred to as **listbc**. It stores information 
    !! about interface between grid **mptr** and all grids one
    !! level finer if they intersect. 
    !! These interfaces consist of multiple cell edges (or segments,
    !! whose length = cell size on grid **mptr**).
    !! Each column in **listbc** has five entries and stores information 
    !! associated with one segment. 
    !! If a global index for all these segments is denoted by **ispot**,
    !! the five entries for the \f$ ispot^{th}\f$ segment can be explained as below: 
    !! 
    !! - listbc(1,ispot) stores LOCAL (RALATIVE to left boundary of grid
    !! **mptr**) *i* index of the cell on grid **mptr** that border this 
    !! segment. 
    !! - listbc(2,ispot) stores LOCAL (RALATIVE to left boundary of grid
    !! **mptr**) *j* index of the cell on grid **mptr** that border this 
    !! segment. 
    !! - listbc(3,ispot) stores side number, which indicates which side 
    !! of this segment is a cell of grid **mptr** (coarse cell).
    !! If this segment is the left edge of a cell on grid **mptr**, this 
    !! number is 1. The number increases (2, 3, 4) as it goes around a coarse cell
    !! clockwise (top edge, right edge, bottom edge).
    !! - listbc(4,ispot) stores grid number of the finer grid that borders
    !! this segment
    !! - listbc(5,ispot) stores the position of this segment with respect to 
    !! the perimeter of this finer grid (grid listbc(4,ispot)). 
    !! The fine grid will save all its fluxes all around its
    !! perimeter. This entry tells where the coarse grid should
    !! take them from. Note that number of fluxes saved here is equal 
    !! to number of such segments (length = coarse cell size) needed to make
    !! up the perimeter.
#ifndef CUDA
    integer, parameter :: cfluxptr  = 5
#endif

    !> pointer to the address of memory storing fluxes in a layer around the grid, 
    !! to be used in conservation fixup near coarse-fine grid intersections.
    !! The memory allocated is 
    !! \f$ 2 \times nvar \times lenbc + naux \times lenbc \f$, 
    !! where \f$lenbc\f$ is 
    !! the number of its parent grid cells adjacent to (around) their boundary interfaces. 
    !! The coefficient 2 is for two spaces. 
    !! One is for plus or minus fluxes, and the other is for the coarse solution for wave fixing.
    !!
    !! From node(ffluxptr,mptr) to node(ffluxptr,mptr)+lenbc-1, 
    !! it stores value that should be added to the adjacent coarse cell
    !! during the synchronization for conservation fixup
    !!
    !! From node(ffluxptr,mptr)+lenbc to node(ffluxptr,mptr)+2*lenbc-1,
    !! it stores solution q on the adjacent coarse cell
#ifndef CUDA
    integer, parameter :: ffluxptr  = 6
#endif

    !> pointer to the address of memory storing the first copy of solution data on this grid,
    !! usually for storing new solution
    integer, parameter :: store1    = 7

    !> pointer to the address of memory storing the second copy of solution data on this grid, 
    !! usually for storing old solution
    integer, parameter :: store2    = 8

    !> global *i* index of left border of this grid
    integer, parameter :: ndilo     = 9

    !> global *i* index of right border of this grid
    integer, parameter :: ndihi     = 10

    !> global *j* index of lower border of this grid
    integer, parameter :: ndjlo     = 11

    !> global *j* index of upper border of this grid
    integer, parameter :: ndjhi     = 12

    !> pointer to the address of memory storing auxiliary data on this grid
    integer, parameter :: storeaux  = 13

    !> pointer to the address of memory storing flags for refinement on this grid
    integer, parameter :: storeflags  = 14

    !> number of flagged cells on this grid
    integer, parameter :: numflags  = 15

    !> domain flags, indexed within base level (lbase) index space
    integer, parameter :: domflags_base  = 16

    !> domain flags, indexed within level-of-this-grid level index space
    integer, parameter :: domflags2  = 17

    !> pointer (actually it's an index in the bndList array) to the first node 
    !! of a linked list, which stores all grids (on the same level) that border this grid
    integer, parameter :: bndListSt  = 18

    !> number of grids (on the same level) that border this grid
    integer, parameter :: bndListNum = 19

    ! :::::::  real part of node descriptor
    !> x-coordinate of the left border of this grid
    integer, parameter :: cornxlo  = 1
    !> y-coordinate of the lower border of this grid
    integer, parameter :: cornylo  = 2
    !> x-coordinate of the right border of this grid
    integer, parameter :: cornxhi  = 3
    !> y-coordinate of the upper border of this grid
    integer, parameter :: cornyhi  = 4
    !> current simulation time on this grid
    integer, parameter :: timemult = 5

    ! :::::::   for linking nodes
    integer, parameter :: nextfree = 2
    integer, parameter :: clawpack_null = 0
    integer, parameter :: nil  = 0

    integer, parameter :: gridNbor = 1 !use 1st col, 2nd col is nextfree - the link

    ! :::::::  for flagging points   
    ! TODO: can use one bit for this instead of real?
    ! needs no refine
    real(CLAW_REAL), parameter :: goodpt = 0.0
    ! needs refine
    real(CLAW_REAL), parameter :: badpt  = 2.0
    real(CLAW_REAL), parameter :: badpro = 3.0

    real(CLAW_REAL), parameter :: NEEDS_TO_BE_SET = 10.e33
    real(CLAW_REAL), parameter :: rinfinity = 10.e32
    integer, parameter :: iinfinity = 999999999
    integer, parameter :: horizontal = 1
    integer, parameter :: vertical = 2
    integer, parameter :: maxgr = 15000
    integer, parameter :: maxlv = 10

    !> maximum number of clusters (grids) on each grid level
    integer, parameter :: maxcl = 5000

    ! The max1d parameter should be changed if using OpenMP grid based 
    ! looping, usually set to max1d = 60
    integer, parameter :: max1d = 60

    integer, parameter :: maxvar = 10
    integer, parameter :: maxaux = 20
    integer, parameter :: maxwave = 10


    ! note use of sentinel in listStart
    integer :: listOfGrids(maxgr),listStart(0:maxlv+1)
    integer,parameter :: bndListSize = 8*maxgr
    integer :: bndList(bndListSize,2)  ! guess size, average # nbors 4? manage as linked list

    real(CLAW_REAL) hxposs(maxlv), hyposs(maxlv),possk(maxlv),rnode(rsize, maxgr) 



    real(CLAW_REAL) tol, tolsp
    integer ibuff,  mstart, ndfree, ndfree_bnd, &
        lfine, & !  level of the finest grid at current time
        node(nsize, maxgr), &
        icheck(maxlv),lstart(maxlv),newstl(maxlv), &
        listsp(maxlv),intratx(maxlv),intraty(maxlv), &
        kratio(maxlv), iregsz(maxlv),jregsz(maxlv), &
        iregst(maxlv),jregst(maxlv), &
        iregend(maxlv),jregend(maxlv), &
        numgrids(maxlv),numcells(maxlv), &
        iorder,&
        mxnest,& ! maximum allowed refined level set by the user
        kcheck
#ifdef CUDA
    
    ! ### integer array ###
    type cpu_2d_int_ptr_type
        integer, dimension(:,:), pointer, contiguous :: ptr=>null()
    end type cpu_2d_int_ptr_type

    type gpu_2d_int_ptr_type
        integer, dimension(:,:), pointer, contiguous, device :: ptr=>null()
    end type gpu_2d_int_ptr_type

    ! ### grid patch data array ###
    ! on CPU
    type cpu_1d_real_ptr_type
        real(CLAW_REAL), dimension(:), pointer, contiguous :: ptr=>null()
    end type cpu_1d_real_ptr_type

    type cpu_3d_real_ptr_type
        real(CLAW_REAL), dimension(:,:,:), pointer, contiguous :: ptr=>null()
    end type cpu_3d_real_ptr_type

    ! on GPU
    type gpu_1d_real_ptr_type
        real(CLAW_REAL), dimension(:), pointer, contiguous, device :: ptr=>null()
    end type gpu_1d_real_ptr_type

    type gpu_3d_real_ptr_type
        real(CLAW_REAL), dimension(:,:,:), pointer, contiguous, device :: ptr=>null()
    end type gpu_3d_real_ptr_type

    type gpu_4d_real_ptr_type
        real(CLAW_REAL), dimension(:,:,:,:), pointer, contiguous, device :: ptr=>null()
    end type gpu_4d_real_ptr_type


    ! fflux_hd(mptr)%ptr and fflux_dd(mptr)%ptr points to the same 
    ! memory location in GPU global memory.
    ! We never allocate fflux_dd(mptr)%ptr
    ! When we need fflux on device, we copy fflux_hd to fflux_dd
    ! by doing fflux_dd = fflux_hd and pass fflux_dd to CUDA kernels.
    type(cpu_1d_real_ptr_type), allocatable         :: fflux_hh(:)
    type(gpu_1d_real_ptr_type), allocatable         :: fflux_hd(:)
    type(gpu_1d_real_ptr_type), allocatable, device :: fflux_dd(:)

    ! We never allocate cflux_dd(mptr)%ptr
    ! When we need cflux on device, we copy cflux_hd to cflux_dd
    type(cpu_2d_int_ptr_type), allocatable         :: cflux_hh(:)
    type(gpu_2d_int_ptr_type), allocatable         :: cflux_hd(:)
    type(gpu_2d_int_ptr_type), allocatable, device :: cflux_dd(:)

    ! These are in SoA format, namely q(i,j,ivar)
    type(cpu_3d_real_ptr_type) :: grid_data(maxgr)

    type(gpu_3d_real_ptr_type) :: grid_data_d(maxgr)
    type(gpu_3d_real_ptr_type) :: fms_d(maxgr)
    type(gpu_3d_real_ptr_type) :: fps_d(maxgr)
    type(gpu_3d_real_ptr_type) :: gms_d(maxgr)
    type(gpu_3d_real_ptr_type) :: gps_d(maxgr)

    type(gpu_3d_real_ptr_type) :: sx_d(maxgr)
    type(gpu_3d_real_ptr_type) :: sy_d(maxgr)
    type(gpu_4d_real_ptr_type) :: wave_x_d(maxgr)
    type(gpu_4d_real_ptr_type) :: wave_y_d(maxgr)

#endif

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::  for alloc array/memory
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Static memory implementation
    ! parameter  (memsize = 10000000)
    ! common  /calloc/   alloc(memsize)

    ! Dynamic memory: 
    !real(CLAW_REAL), allocatable, target, dimension(:) :: storage
    !real(CLAW_REAL), pointer, dimension(:) :: alloc   ! old way, changed mjb Sept. 2014
    real(CLAW_REAL), allocatable, dimension(:) :: alloc    ! new way, use allocatable, not pointer
    integer memsize
       
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\
    ! :::::   for space management of alloc array
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    integer, parameter :: lfdim=5000
    integer lfree(lfdim,2),lenf

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::  domain description variables
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    logical xperdom, yperdom, spheredom
    real(CLAW_REAL) :: xupper, yupper, xlower, ylower
    integer :: nghost, mthbc(4)

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::  collect stats
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    real(CLAW_REAL)  rvoll(maxlv),evol,rvol,avenumgrids(maxlv)
    integer ::  iregridcount(maxlv), tvoll(maxlv)
    integer :: timeRegridding, timeUpdating, timeValout
    integer :: timeFlglvl,timeGrdfit2,timeGrdfit3,timeGrdfitAll
    integer :: timeBound,timeStepgrid
    integer :: timeFlagger, timeBufnst,timeTick
    real(CLAW_REAL) tvollCPU(maxlv), timeTickCPU
    real(CLAW_REAL) timeBoundCPU,timeStepgridCPU,timeRegriddingCPU
    real(CLAW_REAL) timeValoutCPU
    real(CLAW_REAL) timeUpdatingCPU

    integer lentot,lenmax,lendim

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::  method parameters
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    character(len=10), allocatable :: auxtype(:)
    integer  method(7), mwaves, mcapa, dimensional_split
    integer, allocatable :: mthlim(:)
    real(CLAW_REAL) cfl,cflmax,cflv1,cfl_level

    logical :: use_fwaves
    logical :: flag_richardson,flag_gradient
    integer :: verbosity_regrid

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::: Parameters and variables related to I/O and checkpointing
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    logical    printout,matlabout,ncarout

    ! variables for conservation checking:
    real(CLAW_REAL) tmass0

    ! variables for specifying output format
    integer :: output_style, nstop, nout, iout
    real(CLAW_REAL), allocatable :: tout(:)
    real(CLAW_REAL) :: t0, tfinal
    real(CLAW_REAL) :: tstart_thisrun  ! /= t0 in case of restart
    integer :: nq_components, naux_components, output_format
    integer, allocatable :: output_q_components(:)
    integer, allocatable :: output_aux_components(:)
    logical :: output_aux_onlyonce

    ! checkpointing:
    integer :: checkpt_style, nchkpt, checkpt_interval
    real(CLAW_REAL), allocatable :: tchk(:)

    integer :: matlabu

    !  USE UNITS NUMBERS < 89.
    ! 89 and + numthreads taken by gauge output
    integer, parameter :: parmunit = 12
    integer, parameter :: chkunit = 10
    integer, parameter :: inunit  = 5
    integer, parameter :: outunit = 66
    integer, parameter :: pltunit1 = 3
    integer, parameter :: rstunit = 9
    integer, parameter :: dbugunit = 11
    integer, parameter :: matunit = 70

    ! ::::  Debugging flags (verbose output)
    logical &
                dprint,     & !  domain flags output
                eprint,     & !  error estimation output
                edebug,     & !  even more error estimation output
                gprint,     & !  verbose grid generation (clustering,colating...)
                nprint,     & !  nestck reporting
                pprint,     & !  projec tagged pts.
                rprint,     & !  regridding -  summary of new grids
                sprint,     & !  space (memory) output
                tprint,     & !  tick (time stepping) reporting
                uprint        !  updating/upbnding reporting


    ! Restart file name:
    character(len=200) :: rstfile
    logical :: check_a

end module amr_module
