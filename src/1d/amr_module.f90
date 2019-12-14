module amr_module

    implicit none
       
    save
       
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::   data structure info.
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    integer, parameter :: rsize = 3
    integer, parameter :: nsize = 17

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
    integer, parameter :: storeaux  = 11
    integer, parameter :: storeflags  = 12
    integer, parameter :: numflags  = 13
    integer, parameter :: domflags_base  = 14
    integer, parameter :: domflags2  = 15
    integer, parameter :: bndListSt  = 16
    integer, parameter :: bndListNum = 17

    ! :::::::  real part of node descriptor
    integer, parameter :: cornxlo  = 1
    integer, parameter :: cornxhi  = 2
    integer, parameter :: timemult = 3

    ! :::::::   for linking nodes
    integer, parameter :: nextfree = 2
    integer, parameter :: null = 0
    integer, parameter :: nil  = 0

    integer, parameter :: gridNbor = 1 !use 1st col, 2nd col is nextfree - the link

    ! :::::::  for flagging points   
    real(kind=8), parameter :: UNSET = -1.0
    real(kind=8), parameter :: DONTFLAG = 0.0
    real(kind=8), parameter :: DOFLAG  = 2.0
    real(kind=8), parameter :: badpro = 3.0

    real(kind=8), parameter :: NEEDS_TO_BE_SET = 10.e33
    real(kind=8), parameter :: rinfinity = 10.e32
    integer, parameter :: iinfinity = 999999999
    integer, parameter :: horizontal = 1
    integer, parameter :: maxgr = 15000
    integer, parameter :: maxlv = 10
    integer, parameter :: maxcl = 5000

    ! The max1d parameter controls the number of grid cells in 
    ! any single grid patch before breaking up. 
    ! Might want to set smaller to break up into more patches for OpenMP
    ! max1d should now be set in setrun.py
    integer:: max1d

    integer, parameter :: maxvar = 10
    integer, parameter :: maxaux = 20
    integer, parameter :: maxwave = 10


    ! put linked list of grids into array and save.
    ! order is coarsest level to finest. is redone after regridding
    ! and on restarting.  note use of sentinel in listStart
    integer :: listOfGrids(maxgr),listStart(0:maxlv+1)
    integer,parameter :: bndListSize = 8*maxgr
    integer :: bndList(bndListSize,2)  ! guess size, average # nbors 4? manage as linked list

    real(kind=8) hxposs(maxlv),possk(maxlv),rnode(rsize, maxgr)



    real(kind=8) tol, tolsp
    integer ibuff,  mstart, ndfree, ndfree_bnd, lfine, node(nsize, maxgr), &
            icheck(maxlv),lstart(maxlv),newstl(maxlv), &
            listsp(maxlv),intratx(maxlv), &
            kratio(maxlv), iregsz(maxlv), &
            iregst(maxlv), &
            iregend(maxlv), &
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
    !real(kind=8), pointer, dimension(:) :: alloc   ! old way, changed mjb Sept. 2014
    real(kind=8), allocatable, dimension(:) :: alloc    ! new way, use allocatable, not pointer
    integer memsize
       
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\
    ! :::::   for space management of alloc array
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    integer, parameter :: lfdim=5000
    integer lfree(lfdim,2),lenf

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::  domain description variables
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    logical xperdom
    real(kind=8) :: xupper, xlower
    integer :: nghost, mthbc(2)

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::  collect stats
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    real(kind=8)  rvoll(maxlv),evol,rvol,avenumgrids(maxlv)
    integer ::  iregridcount(maxlv), tvoll(maxlv)
    integer :: timeRegridding, timeUpdating, timeValout
    integer :: timeFlglvl,timeGrdfit2,timeGrdfit3,timeGrdfitAll
    integer :: timeSetaux,timeFilval,timeBound,timeStepgrid,timeFilvalTot
    integer :: timeFlagger, timeBufnst, timeTick, tick_clock_start
    real(kind=8) tvollCPU(maxlv), timeTickCPU
    real(kind=8) timeBoundCPU,timeStepgridCPU,timeSetauxCPU,timeRegriddingCPU
    real(kind=8) timeValoutCPU

    integer lentot,lenmax,lendim

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! :::::  method parameters
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    character(len=10), allocatable :: auxtype(:)
    integer  method(4), mwaves, mcapa
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
    real(kind=8) tmass0

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
