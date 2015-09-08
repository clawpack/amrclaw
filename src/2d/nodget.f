c
c ------------------------------------------------------------
c
      integer function nodget(dummy)
c
      use amr_module
      implicit double precision (a-h,o-z)

c
c ::::::::::::::::: NODGET ::::::::::::::::::::::::::::::::::::;
c nodget =  get first free node of the linked list kept in node
c            array. adjust pointers accordingly.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      if (ndfree .ne. null) go to 10
          write(outunit,100) maxgr
          write(*,100)       maxgr
100       format(' out of nodal space - allowed ',i5,' grids')
          stop
c
 10     nodget         = ndfree
        ndfree         = node(nextfree,ndfree)
c
c  initialize nodal block
c
        do 20 i        = 1, nsize
           node(i,nodget) = 0
 20     continue
c
        do 30 i         = 1, rsize
           rnode(i,nodget) = 0.0d0
 30     continue
c
      return
      end
c
c ------------------------------------------------------------
c
      integer function nodget_bnd(dummy)
c
      use amr_module
      implicit double precision (a-h,o-z)

c
c ::::::::::::::::: NODGET_BND ::::::::::::::::::::::::::::::::::::;
c nodget_bnd =  same as above but for bndry list
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      if (ndfree_bnd .ne. null) go to 10
          write(outunit,100) bndListSize
          write(*,100)       bndListSize
100       format(' out of bndry space - allowed ',i5,' bndry grids')
          stop
c
 10   nodget_bnd      = ndfree_bnd
      ndfree_bnd      = bndList(ndfree_bnd,nextfree)
c
c     ##  initialize to 0
      bndList(nodget_bnd,1) = 0
      bndList(nodget_bnd,2) = 0
c
      return
      end
