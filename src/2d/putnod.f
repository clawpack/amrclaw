c
c -------------------------------------------------------------
!> Return **mptr** node to the linked list kept in node array.
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
c
c -------------------------------------------------------------
c
      subroutine putnod_bnd (mcell)
c
      use amr_module
      implicit double precision (a-h,o-z)


c :::::::::::::::::::::::::::::: PUTNOD_BND :::::::::::::::::::::;
c
c  return bndry list node to the linked list 
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      bndList(mcell, nextfree) = ndfree_bnd
      ndfree_bnd               = mcell
c
      return
      end
