
       module amr_module
       implicit none
       
       save
       
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::   data structure info.
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       integer, parameter :: rsize = 5
       integer, parameter :: nsize = 13

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
       integer, parameter :: storeaux  = 13

! :::::::  real part of node descriptor
       integer, parameter :: cornxlo  = 1
       integer, parameter :: cornylo  = 2
       integer, parameter :: cornxhi  = 3
       integer, parameter :: cornyhi  = 4
       integer, parameter :: timemult = 5

! :::::::   for linking nodes
       integer, parameter :: nextfree = 2
       integer, parameter :: null = 0
       integer, parameter :: nil  = 0

! :::::::  for flagging points
       
       real(kind=8), parameter :: goodpt = 0.0
       real(kind=8), parameter :: badpt  = 2.0
       real(kind=8), parameter :: badpro = 3.0

       real(kind=8), parameter :: rinfinity = 10.e32
       integer, parameter :: iinfinity = 999999999
       integer, parameter :: horizontal = 1
       integer, parameter :: vertical = 2
       integer, parameter :: maxgr = 5000
       integer, parameter :: maxlv = 10
       integer, parameter :: maxcl = 500

!      The max1d parameter should be changed if using OpenMP grid based 
!      looping, usually set to max1d = 60
!#ifdef MAX1D
       !integer, parameter :: max1d = MAX1D
!#else
       integer, parameter :: max1d = 600
!#endif

       integer, parameter :: maxvar = 10
       integer, parameter :: maxaux = 20
       integer, parameter :: maxout = 5000


       real(kind=8) hxposs(maxlv), hyposs(maxlv),possk(maxlv),&
               rnode(rsize, maxgr) 



       real(kind=8) tol, tolsp
       integer ibuff,  mstart, ndfree, lfine, node(nsize, maxgr), &
               icheck(maxlv),lstart(maxlv),newstl(maxlv), &
               listsp(maxlv),intratx(maxlv),intraty(maxlv), &
               kratio(maxlv), iregsz(maxlv),jregsz(maxlv), &
               iregst(maxlv),jregst(maxlv), &
               iregend(maxlv),jregend(maxlv), &
               numgrids(maxlv),numcells(maxlv), &
               iorder,mxnest,kcheck,nghost

       integer ngrids
!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      ::::  for alloc array/memory
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!       Static memory implementation
!        parameter  (memsize = 10000000)
!        common  /calloc/   alloc(memsize)

!      Dynamic memory: 
       real(kind=8), allocatable, target, dimension(:) :: storage
       real(kind=8), pointer, dimension(:) :: alloc
       integer memsize
       
       
!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::   for space management of alloc array
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       integer, parameter :: lfdim=5000

       integer lfree(lfdim,2),lenf

!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::  domain description variables
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       logical xperdom, yperdom, spheredom
       real(kind=8) xupper,yupper,xlower,ylower

!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::  collect stats
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       real(kind=8) rvoll(10),evol,rvol
       integer lentot,lenmax,lendim

!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::  method parameters
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       integer, parameter :: maxwave = 10
       character * 10 auxtype(maxaux)
       integer  method(7), mthlim(maxwave), mwaves, mcapa
       real(kind=8) cfl,cflmax,cflv1,cfl_level


       logical :: fwave
       logical :: flag_richardson,flag_gradient
       integer :: verbosity_regrid


!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      ::::: Parameters and variables related to I/O
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


       logical    printout,matlabout,ncarout

!      variables for conservation checking:
       real(kind=8) tstart,tmass0

!      variables for specifying output format
       integer :: nq_components,naux_components,output_format
       integer, dimension(maxvar) :: output_q_components
       integer, dimension(maxaux) :: output_aux_components
       logical :: output_aux_onlyonce
       
       integer :: matlabu

       integer, parameter :: parmunit = 12
       integer, parameter :: chkunit = 10
       integer, parameter :: inunit  = 5
       integer, parameter :: outunit = 66
       integer, parameter :: pltunit1 = 3
       integer, parameter :: rstunit = 9
       integer, parameter :: dbugunit = 11
       integer, parameter :: matunit = 70

!      ::::  Debugging flags (verbose output)

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


!     Restart file name:
      character*12 :: rstfile

       end module amr_module
