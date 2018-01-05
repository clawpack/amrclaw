c
c ---------------------------------------------------------
c
      subroutine cleanup(nvar,naux)
c
c :::::::::::::::::::::: CLEANUP ::::::::::::::::::::::::::::::::;
!> this is just a check to make sure all storage was accounted for.
!! routine is called after all the data has been checkpointed.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      use amr_module
#ifdef CUDA
      use memory_module, only: gpu_deallocate, cpu_deallocated_pinned
      use cuda_module, only: device_id
#endif
      implicit double precision (a-h,o-z)
c
c      ## clean up storage to double check that everything taken care of
c      ## done after the checkpoint so pointers sitll work on restart
       do  120 level = 1, lfine
         call putsp(1,level,nvar,naux)
         mptr =  lstart(level)
 110        nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
            ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
            mitot  = nx + 2*nghost
            mjtot  = ny + 2*nghost
            nwords  = mitot*mjtot*nvar
            call reclam(node(store1, mptr), nwords)
#ifdef CUDA
            if (associated(grid_data(mptr)%ptr)) then
                call cpu_deallocated_pinned(grid_data(mptr)%ptr)
            else
                print *, 'grid_data of grid: ',mptr,' is already freed'
                stop
            endif

            if (associated(grid_data_d(mptr)%ptr)) then
                call gpu_deallocate(grid_data_d(mptr)%ptr, device_id)
            else
                print *, "grid_data_d of grid: ",mptr,
     &              " is already freed"
                stop
            endif

            if (associated(fms_d(mptr)%ptr)) then
                call gpu_deallocate(fms_d(mptr)%ptr, device_id)
            else
                print *, "fms_d of grid: ", mptr, " is already freed"
                stop
            endif

            if (associated(fps_d(mptr)%ptr)) then
                call gpu_deallocate(fps_d(mptr)%ptr, device_id)
            else
                print *, "fps_d of grid: ", mptr, " is already freed"
                stop
            endif

            if (associated(gms_d(mptr)%ptr)) then
                call gpu_deallocate(gms_d(mptr)%ptr, device_id)
            else
                print *, "gms_d of grid: ", mptr, " is already freed"
                stop
            endif

            if (associated(gps_d(mptr)%ptr)) then
                call gpu_deallocate(gps_d(mptr)%ptr, device_id)
            else
                print *, "gps_d of grid: ", mptr, " is already freed"
                stop
            endif
            call gpu_deallocate(sx_d(mptr)%ptr, device_id)
            call gpu_deallocate(sy_d(mptr)%ptr, device_id)
            call gpu_deallocate(wave_x_d(mptr)%ptr, device_id)
            call gpu_deallocate(wave_y_d(mptr)%ptr, device_id)
#endif
            if (level .lt. mxnest) 
     .         call reclam(node(store2, mptr), nwords)
            if (naux .gt. 0) 
     .         call reclam(node(storeaux, mptr), mitot*mjtot*naux)
        mptr = node(levelptr, mptr)
        if (mptr .ne. 0) go to 110
120    continue 

      return
      end
