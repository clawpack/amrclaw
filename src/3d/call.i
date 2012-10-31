
!  :::::::   number of spatial dimensions
       parameter(numdim = 3)

!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::   data structure info.
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       integer    cornxlo,cornylo,cornxhi,cornyhi,timemult
       integer    cornzlo,        cornzhi
       integer    store1,store2,storeaux
       integer    tempptr,errptr,ffluxptr,cfluxptr
       integer    rsize
       integer    iplane, jplane, kplane

       parameter (rsize =  7)
       parameter (nsize = 15)

!  :::::::   integer part of node descriptor
       parameter (levelptr  = 1)
       parameter (tempptr   = 2)
       parameter (errptr    = 3)
       parameter (nestlevel = 4)
       parameter (cfluxptr  = 5)
       parameter (ffluxptr  = 6)
       parameter (store1    = 7)
       parameter (store2    = 8)
       parameter (ndilo     = 9)
       parameter (ndihi     = 10)
       parameter (ndjlo     = 11)
       parameter (ndjhi     = 12)
       parameter (ndklo     = 13)
       parameter (ndkhi     = 14)
       parameter (storeaux  = 15)

! :::::::  real part of node descriptor
       parameter (cornxlo  = 1)
       parameter (cornylo  = 2)
       parameter (cornzlo  = 3)
       parameter (cornxhi  = 4)
       parameter (cornyhi  = 5)
       parameter (cornzhi  = 6)
       parameter (timemult = 7)

! :::::::   for linking nodes
       parameter (nextfree = 2)
       parameter (null = 0)
       parameter (nil  = 0)

! :::::::  for flagging points
       parameter (goodpt = 0.0)
       parameter (badpt  = 2.0)
       parameter (badpro = 3.0)

       parameter (rinfinity = 10.e32)
       parameter (iinfinity = 999999)
       parameter (iplane    = 1)
       parameter (jplane    = 2)
       parameter (kplane    = 3)
       parameter  (maxgr = 852)
       parameter  (maxlv = 6)
       parameter  (maxcl = 852)
!      parameter  (max1d = 120)
       parameter  (max1d = 60)
       parameter  (maxvar = 9)
       parameter  (maxaux = 20)
       parameter  (maxout = 50)

       logical    printout, matlabout, ncarout

       common  /nodal/
     1         hxposs(maxlv),hyposs(maxlv),hzposs(maxlv),possk(maxlv),
     2         rnode(rsize, maxgr), node(nsize, maxgr),
     3         icheck(maxlv),lstart(maxlv),newstl(maxlv),listsp(maxlv),
     4         intratx(maxlv), intraty(maxlv),intratz(maxlv),
     5         kratio(maxlv),iregsz(maxlv),jregsz(maxlv),kregsz(maxlv),
     6         iregst(maxlv),jregst(maxlv),kregst(maxlv),
     7         iregend(maxlv),jregend(maxlv),kregend(maxlv),
     8         numgrids(maxlv),numcells(maxlv),
     9         tol, tolsp, ibuff,  mstart, ndfree, lfine,
     1         iorder,mxnest,kcheck,nghost,printout,matlabout,ncarout

       common /cmatlab/ matlabu,ngrids

!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      ::::  for alloc array/memory
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!       Static memory implementation
!        parameter  (memsize = 10000000)
!        common  /calloc/   alloc(memsize)

!      Dynamic memmory
       double precision, pointer, dimension(:) :: alloc
       common /calloc/ alloc, memsize
       

!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::   for space management of alloc array
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       parameter (lfdim=500)

       common /space/
     1               lfree(lfdim,2),lenf

!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::  domain description variables
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       logical xperdom, yperdom, zperdom

       common /cdom/ xupper,yupper,zupper,xlower,ylower,zlower,
     1                 xperdom, yperdom, zperdom


!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::  collect stats
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       common   /stats/  rvoll(maxlv),avenumgrids(maxlv),
     1                   iregridcount(maxlv),evol,rvol,
     2                   lentot,lenmax,lendim

!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::  method parameters
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       parameter (maxwave = 10)
       character * 10 auxtype(maxaux)
       common /cmeth/ method(7), mthlim(maxwave), mwaves, mcapa
       common /auxstuff/ auxtype
       common /ccfl/ cfl,cflmax, cflv1
!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      ::::  for i/o assignments
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!

       integer chkunit,dbugunit,inunit,outunit,pltunit1,rstunit
       integer matunit

       parameter (chkunit = 10)
       parameter (inunit  = 5)
       parameter (outunit = 66)
       parameter (pltunit1 = 3)
       parameter (rstunit = 9)
       parameter (dbugunit = 11)
       parameter (matunit = 70)


!      ::::  common for debugging flags (verbose output)

       logical
     .         dprint,     !  domain flags output
     .         eprint,     !  error estimation output
     .         edebug,     !  even more error estimation output
     .         gprint,     !  verbose grid gen. (clustering,colating)
     .         nprint,     !  nestck reporting
     .         pprint,     !  projec tagged pts.
     .         rprint,     !  regridding -  summary of new grids
     .         sprint,     !  space (memory) output
     .         tprint,     !  tick (time stepping) reporting
     .         uprint      !  updating/upbnding reporting


       common /bugflags/ dprint, eprint, edebug, gprint, nprint, pprint,
     .                   rprint, sprint, tprint, uprint


c
c      ::::  common for conservation check
      common /ctstart/ tstart
      common /comconck/ tmass0

