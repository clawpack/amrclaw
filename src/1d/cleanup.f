c
c ---------------------------------------------------------
c
      subroutine cleanup(nvar,naux)
c
c :::::::::::::::::::::: CLEANUP ::::::::::::::::::::::::::::::::;
c   this is just a check to make sure all storage was accounted for.
c   routine is called after all the data has been checkpointed.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      use amr_module
      implicit double precision (a-h,o-z)
c
c      ## clean up storage to double check that everything taken care of
c      ## done after the checkpoint so pointers sitll work on restart
       do  120 level = 1, lfine
         call putsp(1,level,nvar,naux)
         mptr =  lstart(level)
 110        nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
            mitot  = nx + 2*nghost
            nwords  = mitot*nvar
            call reclam(node(store1, mptr), nwords)
            if (level .lt. mxnest) 
     .         call reclam(node(store2, mptr), nwords)
            if (naux .gt. 0) 
     .         call reclam(node(storeaux, mptr), mitot*naux)
        mptr = node(levelptr, mptr)
        if (mptr .ne. 0) go to 110
120    continue 

      return
      end
