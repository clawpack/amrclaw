module amr_module

    implicit none
       
    save
       
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::   data structure info.
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    integer, parameter :: rsize = 7
    integer, parameter :: nsize = 21

    !  :::::::   integer part of node descriptor
    integer, parameter :: levelptr  = 1
    integer, parameter :: tempptr   = 2
    integer, parameter :: errptr    = 3
    integer, parameter :: nestlevel = 4
    integer, parameter :: cfluxptr  = 5
    integer, parameter :: ffluxptr  = 6
    integer, parameter :: store1    = 7
    integer, parameter :: store2    = 8
    integer, parameter :: ndilo     = 9
    integer, parameter :: ndihi     = 10
    integer, parameter :: ndjlo     = 11
    integer, parameter :: ndjhi     = 12
    integer, parameter :: ndklo     = 13
    integer, parameter :: ndkhi     = 14
    integer, parameter :: storeaux  = 15
    integer, parameter :: storeflags  = 16
    integer, parameter :: numflags  = 17
    integer, parameter :: domflags_base  = 18
    integer, parameter :: domflags2  = 19
    integer, parameter :: bndListSt  = 20
    integer, parameter :: bndListNum = 21

    ! :::::::  real part of node descriptor
    integer, parameter :: cornxlo  = 1
    integer, parameter :: cornylo  = 2
    integer, parameter :: cornzlo  = 3
    integer, parameter :: cornxhi  = 4
    integer, parameter :: cornyhi  = 5
    integer, parameter :: cornzhi  = 6
    integer, parameter :: timemult = 7

    ! :::::::   for linking nodes
    integer, parameter :: nextfree = 2
    integer, parameter :: null = 0
    integer, parameter :: nil  = 0

    integer, parameter :: gridNbor = 1  ! use first col, 2nd col is nextfee - the link

    ! :::::::  for flagging points   
    real(kind=8), parameter :: goodpt = 0.0
    real(kind=8), parameter :: badpt  = 2.0
    real(kind=8), parameter :: badpro = 3.0

    real(kind=8), parameter :: rinfinity = 10.e32
    integer, parameter :: iinfinity = 999999999
    integer, parameter :: iplane = 1
    integer, parameter :: jplane = 2
    integer, parameter :: kplane = 3
    !! integer, parameter :: maxgr = 5000 ! No longer fixed
    integer, parameter :: maxlv = 10
    integer, parameter :: maxcl = 5000

    ! The max1d parameter should be changed if using OpenMP grid based 
    ! looping, suggest setting max1d = 32 for 3d 
    integer, parameter :: max1d = 32 
    !integer, parameter :: max1d = 28 

    integer, parameter :: maxvar = 10
    integer, parameter :: maxaux = 20
    integer, parameter :: maxwave = 10
    integer, parameter :: maxout = 50  !until change amr to f90 and allocate

    ! put linked list of grids into array and save.
    ! order is coarsest level to finest. is redone after regridding
    ! and on restarting.  note use of sentinel in listStart

    !! CHANGED mjb 6/15/18 to make maxgr a variable
    !!integer :: listOfGrids(maxgr),listStart(0:maxlv+1)
    integer :: listStart(0:maxlv+1)
    !integer, parameter :: bndListSize = 8*maxgr
    !integer :: bndList(bndListSize,2) ! guess size, average # nbors 6? manage as linked list
    integer :: bndListSize
    integer,allocatable, dimension(:,:) :: bndList  ! new way is allocatable 


    !real(kind=8) hxposs(maxlv),hyposs(maxlv),hzposs(maxlv), &
    !             possk(maxlv),rnode(rsize, maxgr) 
    real(kind=8) hxposs(maxlv),hyposs(maxlv),hzposs(maxlv),possk(maxlv) 

    ! start of dynamic allocation for maxgr and associated arrays
    integer maxgr
    real(kind=8), allocatable, dimension(:,:) :: rnode    
    ! new way, use allocatable, not pointer
    integer, allocatable, dimension(:,:) :: node
    integer, allocatable, dimension(:) :: listOfGrids


    real(kind=8) tol, tolsp
    !integer ibuff, mstart,ndfree,ndfree_bnd,lfine,node(nsize, maxgr), &
    integer ibuff, mstart,ndfree,ndfree_bnd,lfine, &
            icheck(maxlv),lstart(maxlv),newstl(maxlv), &
            listsp(maxlv),intratx(maxlv),intraty(maxlv), &
            intratz(maxlv),kratio(maxlv), &
            iregsz(maxlv),jregsz(maxlv),kregsz(maxlv), &
            iregst(maxlv),jregst(maxlv),kregst(maxlv), &
            iregend(maxlv),jregend(maxlv),kregend(maxlv), &
            numgrids(maxlv),numcells(maxlv), &
            iorder,mxnest,kcheck

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::  for alloc array/memory
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Static memory implementation
    ! parameter  (memsize = 10000000)
    ! common  /calloc/   alloc(memsize)

    ! Dynamic memory: 
    !real(kind=8), allocatable, target, dimension(:) :: storage
    !real(kind=8), pointer, dimension(:) :: alloc
    real(kind=8), allocatable, dimension(:) :: alloc  ! new way, use alloctable, not pointer
    integer memsize
       
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\
    ! :::::   for space management of alloc array
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    integer, parameter :: lfdim=5000
    integer lfree(lfdim,2),lenf

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::  domain description variables
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    logical xperdom, yperdom ,zperdom
    real(kind=8) :: xupper,yupper,zupper,xlower,ylower,zlower
    integer :: nghost, mthbc(6)

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::  collect stats
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    real(kind=8)  rvoll(maxlv),evol,rvol,avenumgrids(maxlv)
    integer iregridcount(maxlv), tvoll(maxlv)
    integer lentot,lenmax,lendim
    integer timeRegridding, timeValout
    integer timeBound, timeStepgrid
    integer :: timeTick, tick_clock_start
    real(kind=8) tvollCPU(maxlv)
    real(kind=8) timeBoundCPU,timeStepgridCPU,timeRegriddingCPU
    real(kind=8) timeValoutCPU

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::  method parameters
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    character(len=10), allocatable :: auxtype(:)
    integer  method(7), mwaves, mcapa, dimensional_split
    integer, allocatable :: mthlim(:)
    real(kind=8) cfl,cflmax,cflv1,cfl_level

    logical :: use_fwaves
    logical :: flag_richardson,flag_gradient
    integer :: verbosity_regrid

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::: Parameters and variables related to I/O and checkpointing
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    logical    printout,matlabout,ncarout

    ! variables for conservation checking:
    real(kind=8) tstart,tmass0

    ! variables for specifying output format
    integer :: output_style, nstop, nout, iout
    real(kind=8), allocatable :: tout(:)
    real(kind=8) :: t0, tfinal
    real(kind=8) :: tstart_thisrun  ! /= t0 in case of restart
    integer :: nq_components, naux_components, output_format
    integer, allocatable :: output_q_components(:)
    integer, allocatable :: output_aux_components(:)
    logical :: output_aux_onlyonce

    ! checkpointing:
    integer :: checkpt_style, nchkpt, checkpt_interval
    real(kind=8), allocatable :: tchk(:)

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
