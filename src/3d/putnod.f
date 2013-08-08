c
c -------------------------------------------------------------
c
      subroutine putnod (mptr)
c
      use amr_module
      implicit double precision (a-h,o-z)


c :::::::::::::::::::::::::::::: PUTNOD :::::::::::::::::::::;
c
c  return mptr node to the linked list kept in node array.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      node(nextfree, mptr) = ndfree
      ndfree               = mptr
c
      return
      end
