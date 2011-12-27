
       module amr_module

c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c      :::::   data structure info.
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
       integer    cornxlo,cornylo,cornxhi,cornyhi,timemult
       integer    store1,store2,storeaux
       integer    tempptr,errptr,ffluxptr,cfluxptr
       integer    rsize

       parameter (rsize =  5)
       parameter (nsize = 13)

c  :::::::   integer part of node descriptor
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
       parameter (storeaux  = 13)

c :::::::  real part of node descriptor
       parameter (cornxlo  = 1)
       parameter (cornylo  = 2)
       parameter (cornxhi  = 3)
       parameter (cornyhi  = 4)
       parameter (timemult = 5)

c :::::::   for linking nodes
       parameter (nextfree = 2)
       parameter (null = 0)
       parameter (nil  = 0)

c :::::::  for flagging points
       parameter (goodpt = 0.0)
       parameter (badpt  = 2.0)
       parameter (badpro = 3.0)

       parameter (rinfinity = 10.e32)
       parameter (iinfinity = 999999)
       parameter (horizontal = 1)
       parameter (vertical = 2)
       parameter  (maxgr = 500)
       parameter  (maxlv = 10)
       parameter  (maxcl = 500)
       parameter  (max1d = 500)
       parameter  (maxvar = 10)
       parameter  (maxaux = 20)
       parameter  (maxout = 5000)

       logical    printout,matlabout,ncarout

       double precision hxposs(maxlv), hyposs(maxlv),possk(maxlv),
     &         rnode(rsize, maxgr) 



       double precision tol, tolsp
       integer ibuff,  mstart, ndfree, lfine, node(nsize, maxgr),
     &         icheck(maxlv),lstart(maxlv),newstl(maxlv),
     &         listsp(maxlv),intratx(maxlv),intraty(maxlv),
     &         kratio(maxlv), iregsz(maxlv),jregsz(maxlv),
     &         iregst(maxlv),jregst(maxlv),
     &         iregend(maxlv),jregend(maxlv),
     &         iorder,mxnest,kcheck,nghost

       integer matlabu,ngrids
c
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c      ::::  for alloc array/memory
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c       Static memory implementation
c        parameter  (memsize = 10000000)
c        common  /calloc/   alloc(memsize)

c      Dynamic memory: 
       double precision, save, allocatable, target, dimension(:) ::
     &        storage
       double precision, save, pointer, dimension(:) :: alloc
       integer memsize
       
       
c
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c      :::::   for space management of alloc array
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
       parameter (lfdim=5000)

       integer lfree(lfdim,2),lenf

c
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c      :::::  domain description variables
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
       logical xperdom, yperdom, spheredom
       double precision xupper,yupper,xlower,ylower

c
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c      :::::  collect stats
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
       double precision rvoll(10),evol,rvol
       integer lentot,lenmax,lendim

c
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c      :::::  method parameters
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
       parameter (maxwave = 10)
       character * 10 auxtype(maxaux)
       integer  method(7), mthlim(maxwave), mwaves, mcapa
       double precision cfl,cflmax,cflv1
c
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c      ::::  for i/o assignments
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

       integer chkunit,dbugunit,inunit,outunit,pltunit1,rstunit
       integer matunit,parmunit

       parameter (parmunit = 12)
       parameter (chkunit = 10)
       parameter (inunit  = 5)
       parameter (outunit = 66)
       parameter (pltunit1 = 3)
       parameter (rstunit = 9)
       parameter (dbugunit = 11)
       parameter (matunit = 70)


c      ::::  common for debugging flags (verbose output)

       logical
     &         dprint,     !  domain flags output
     &         eprint,     !  error estimation output
     &         edebug,     !  even more error estimation output
     &         gprint,     !  verbose grid generation (clustering,colating...)
     &         nprint,     !  nestck reporting
     &         pprint,     !  projec tagged pts.
     &         rprint,     !  regridding -  summary of new grids
     &         sprint,     !  space (memory) output
     &         tprint,     !  tick (time stepping) reporting
     &         uprint      !  updating/upbnding reporting


       end module amr_module
