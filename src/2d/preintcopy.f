c
c --------------------------------------------------------------
c
      subroutine preintcopy(val,mitot,mjtot,nvar,ilo,ihi,jlo,jhi,
     1                      level,fliparray)
c
      use amr_module
      implicit double precision (a-h,o-z)

      dimension fliparray((mitot+mjtot)*nghost*nvar)
      dimension val(nvar,mitot,mjtot)
      dimension ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)

c  NEW INDEXING ORDER SWITCHED
      iadd(ivar,i,j)  = locflip + ivar-1 + nvar*((j-1)*nr+i-1)
c
c  :::::::::::::: PREINTCOPY :::::::::::::::::::::::::::::::::::::::::::
c     For periodic boundary conditions more work needed to initialize a
c     new grid that sticks out. This routine was
c     called because the patch sticks out of the domain,
c     and has periodic bc.s preprocess the patch before calling
c     intcopy to shift the patch periodically back into the domain.
c
c     Inputs to this routine:
c     ilo,ihi,jlo,jhi = the location in index space of this patch.
c
c     Outputs from this routine:
c     The values of the grid are inserted
c     directly into the enlarged val array for this piece.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
        locflip = 1        

c
c     # will divide patch into 9 possibilities (some empty): 
c       x sticks out left, x interior, x sticks out right
c       same for y. for example, the max. would be
c       i from (ilo,-1), (0,iregsz(level)-1), (iregsz(level),ihi)
        
        if (xperdom) then
           ist(1) = ilo
           ist(2) = 0
           ist(3) = iregsz(level)
           iend(1) = -1
           iend(2) = iregsz(level)-1
           iend(3) = ihi
           ishift(1) = iregsz(level)
           ishift(2) = 0
           ishift(3) = -iregsz(level)
        else    ! if not periodic, set vals to only have nonnull intersection for interior region
           ist(1)    = iregsz(level)
           ist(2)    = ilo
           ist(3)    = iregsz(level)
           iend(1)   = -1
           iend(2)   = ihi
           iend(3)   = -1
           ishift(1) = 0
           ishift(2) = 0
           ishift(3) = 0
        endif


        if (yperdom .or. spheredom) then
           jst(1) = jlo
           jst(2) = 0
           jst(3) = jregsz(level)
           jend(1) = -1
           jend(2) = jregsz(level)-1
           jend(3) = jhi
           jshift(1) = jregsz(level)
           jshift(2) = 0
           jshift(3) = -jregsz(level)
        else
           jst(1)    = jregsz(level)
           jst(2)    = jlo
           jst(3)    = jregsz(level)
           jend(1)   = -1
           jend(2)   = jhi
           jend(3)   = -1
           jshift(1) = 0
           jshift(2) = 0
           jshift(3) = 0
        endif


       do 20 i = 1, 3
          i1 = max(ilo,  ist(i))
          i2 = min(ihi, iend(i))
          if (i1 .gt. i2) go to 20
       do 10 j = 1, 3
          j1 = max(jlo,  jst(j))
          j2 = min(jhi, jend(j))


          if (j1 .le. j2) then ! part of patch in this region
c
c check if special mapping needed for spherical bc. 
c(j=2 is interior,nothing special needed)
            if (.not. spheredom .or. j .eq. 2) then
            iputst = (i1 - ilo) + 1
            jputst = (j1 - jlo) + 1
            call intcopy(val,mitot,mjtot,nvar,
     2                    i1+ishift(i),i2+ishift(i),
     3                    j1+jshift(j),j2+jshift(j),level,
     4                    iputst,jputst)
          else
              nr = i2 - i1 + 1
              nc = j2 - j1 + 1
              ng = 0    ! no ghost cells in this little patch. fill everything.

              jbump = 0
              if (j1 < 0)   jbump = abs(j1)  ! bump up so new bottom is 0
              if (j2 >= jregsz(level)) jbump = -(j2+1-jregsz(level)) ! bump down

c             next 2 lines would take care of periodicity in x
              iwrap1 = i1 + ishift(i)
              iwrap2 = i2 + ishift(i)
c             next 2 lines take care of mapped sphere bcs
              iwrap1 = iregsz(level) - iwrap1 -1
              iwrap2 = iregsz(level) - iwrap2 -1
c             swap so that smaller one is left index, etc since mapping reflects
              tmp = iwrap1
              iwrap1 = iwrap2
              iwrap2 = tmp

              jwrap1 = j1 + jbump
              xlwrap = xlower + iwrap1*hxposs(level)
              ybwrap = ylower + jwrap1*hyposs(level)
              jwrap2 = j2 + jbump

c             write(dbugunit,101) i1,i2,j1,j2
c             write(dbugunit6,102) iwrap1,iwrap2,j1+jbump,j2+jbump
 101          format(" actual patch from i:",2i5," j :",2i5)
 102          format(" intcopy called w i:",2i5," j :",2i5)
              call intcopy(fliparray,nr,nc,nvar,
     1                     iwrap1,iwrap2,jwrap1,jwrap2,level,1,1)

c             copy back using weird mapping for spherical folding
              nrowst = 1   ! start filling up val at (1,1) - no additional offset
              ncolst = 1
              do ii = i1, i2
              do jj = j1, j2
c            write(dbugunit6,100)nrowst+ii-ilo,ncolst+jj-jlo,nr-(ii-i1),
c    1                            nc-jj+j1
 100          format(" filling loc ",2i5," with ",2i5)

              do ivar = 1, nvar
                 iindex = nr-(ii-i1)
                 jindex = nc-(jj-j1)
                 index = iadd(ivar,nr-(ii-i1),nc-(jj-j1))
                 val(ivar,nrowst+(ii-ilo),ncolst+(jj-jlo)) = 
     1                  fliparray(iadd(ivar,nr-(ii-i1),nc-(jj-j1)))
              end do
              end do
              end do
             

            endif

          endif

 10    continue
 20    continue
      
     
      return
      end
