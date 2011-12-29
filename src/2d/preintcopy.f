c
c --------------------------------------------------------------
c
      subroutine preintcopy(val,mitot,mjtot,nvar,ilo,ihi,jlo,jhi,
     1                      level,locflip)
c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension val(nvar,mitot,mjtot)
      dimension ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)

c  OLD INDEXING
c     iadd(i,j,ivar)  = locflip + i - 1 + nr*((ivar-1)*nc+j-1)
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
c     # will divide patch into 9 possibilities (some empty): 
c       x sticks out left, x interior, x sticks out right
c       same for y. for example, the max. would be
c       i from (ilo,-1), (0,iregsz(level)-1), (iregsz(level),ihi)
        
        ist(1) = ilo
        ist(2) = 0
        ist(3) = iregsz(level)
        iend(1) = -1
        iend(2) = iregsz(level)-1
        iend(3) = ihi

        jst(1) = jlo
        jst(2) = 0
        jst(3) = jregsz(level)
        jend(1) = -1
        jend(2) = jregsz(level)-1
        jend(3) = jhi

        ishift(1) = iregsz(level)
        ishift(2) = 0
        ishift(3) = -iregsz(level)
        jshift(1) = jregsz(level)
        jshift(2) = 0
        jshift(3) = -jregsz(level)

       do 20 i = 1, 3
          i1 = max(ilo,  ist(i))
          i2 = min(ihi, iend(i))
       do 10 j = 1, 3
          j1 = max(jlo,  jst(j))
          j2 = min(jhi, jend(j))


          if ((i1 .le. i2) .and. (j1 .le. j2)) then ! part of patch in this region
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

c             ## with dynamic memory alloc had to get scratch space above
c             ## in case alloc moved
c             locflip = igetsp(nr*nc*nvar)

c             write(dbugunit,101) i1,i2,j1,j2
c             write(dbugunit6,102) iwrap1,iwrap2,j1+jbump,j2+jbump
 101          format(" actual patch from i:",2i5," j :",2i5)
 102          format(" intcopy called w i:",2i5," j :",2i5)
              call intcopy(alloc(locflip),nr,nc,nvar,
     1                     iwrap1,iwrap2,jwrap1,jwrap2,level,1,1)

c             copy back using weird mapping for spherical folding
              nrowst = 1   ! start filling up val at (1,1) - no additional offset
              ncolst = 1
              do 15 ii = i1, i2
              do 15 jj = j1, j2
c            write(dbugunit6,100)nrowst+ii-ilo,ncolst+jj-jlo,nr-(ii-i1),
c    1                            nc-jj+j1
 100          format(" filling loc ",2i5," with ",2i5)

              do 15 ivar = 1, nvar
                 val(ivar,nrowst+(ii-ilo),ncolst+(jj-jlo)) = 
     1                  alloc(iadd(ivar,nr-(ii-i1),nc-(jj-j1)))
 15           continue
             
c             call reclam(locflip, (nr*nc*nvar))

            endif

          endif

 10    continue
 20    continue
      
     
      return
      end
