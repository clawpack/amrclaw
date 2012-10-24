c
c ---------------------------------------------------------------
c
        subroutine trimbd(used,set,il ,ir ,jf ,jr ,kb ,kt ,
     &                             ilo,ihi,jlo,jhi,klo,khi)
c
c ::::::::::::::::::::: TRIMBD :::::::::::::::::::::::::::::::::::
c  if used array is completely set (=1.) then return set=true, 
c  otherwise return false, and the dimensions of the smallest rectangular
c  parallelopiped containing all unset points in il,ir,jf,jr,kb,kt.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
        implicit double precision (a-h,o-z)
        dimension  used(ilo:ihi,jlo:jhi,klo:khi)
        logical   set

        utot = 0.
        nrow = ihi - ilo + 1
        ncol = jhi - jlo + 1
        nfil = khi - klo + 1

        do 100 k = klo, khi
        do 100 j = jlo, jhi
        do 100 i = ilo, ihi
100        utot = utot + used(i,j,k)

        if (utot .ge. dble(nrow*ncol*nfil)) then
                set = .true.
                return
        endif
 
c some unset cells were found

        set = .false.
 
        uleft = 1.
        do 200 i = ilo, ihi
           do 220 k = klo, khi
           do 220 j = jlo, jhi
              uleft = dmin1(uleft,used(i,j,k))
220        continue
           il = i -ilo + 1
           if (uleft .eq. 0.) go to 230
200     continue

230     uright = 1.
        do 300 i = ihi, ilo, -1
           do 320 k = klo, khi
           do 320 j = jlo, jhi
              uright = dmin1(uright,used(i,j,k))
320        continue
           ir = i - ilo + 1
           if (uright .eq. 0.) go to 330
300     continue

330     ufront = 1.
        do 400 j = jlo, jhi
           do 420 k = klo, khi
           do 420 i = ilo, ihi
              ufront = dmin1(ufront,used(i,j,k))
420        continue
           jf = j - jlo + 1
           if (ufront .eq. 0.) go to 430
400        continue
 
430     urear = 1.
        do 500 j = jhi, jlo, -1
           do 520 k = klo, khi
           do 520 i = ilo, ihi
              urear = dmin1(urear,used(i,j,k))
520        continue
           jr = j - jlo + 1
           if (urear .eq. 0.) go to 530
500     continue
 
530     ubot = 1.
        do 600 k = klo, khi
           do 620 j = jlo, jhi
           do 620 i = ilo, ihi
              ubot = dmin1(ubot,used(i,j,k))
620        continue
           kb = k - klo + 1
           if (ubot .eq. 0.) go to 630
600     continue
 
630     utop = 1.
        do 700 k = khi, klo, -1
           do 720 j = jlo, jhi
           do 720 i = ilo, ihi
              utop = dmin1(utop,used(i,j,k))
720        continue
           kt = k - klo + 1
           if (utop .eq. 0.) go to 730
700     continue

730     continue

        return
        end
