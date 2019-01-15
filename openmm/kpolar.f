c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kpolar  --  assign polarizability parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kpolar" assigns atomic dipole polarizabilities to the atoms
c     within the structure and processes any new or changed values
c
c     literature references:
c
c     A. C. Simmonett, F. C. Pickard IV, J. W. Ponder and B. R. Brooks,
c     "An Empirical Extrapolation Scheme for Efficient Treatment of
c     Induced Dipoles", Journal of Chemical Physics, 145, 164101 (2016)
c     [OPT method]
c
c     F. Aviat, L. Lagardere and J.-P. Piquemal, "The Truncated
c     Conjugate Gradient (TCG), a Non-Iterative/Fixed-Cost Strategy for
c     Computing Polarization in Molecular Dynamics: Fast Evaluation of
c     Analytical Forces", Journal of Chemical Physics, 147, 161724
c     (2018)  [TCG method]
c
c
      subroutine kpolar
      use atoms
      use chgpen
      use inform
      use iounit
      use keys
      use kpolr
      use mplpot
      use mpole
      use neigh
      use polar
      use polopt
      use polpot
      use polpcg
      use poltcg
      use potent
      implicit none
      integer i,j,k
      integer next,size
      integer nlist,npg
      integer pg(maxval)
      integer, allocatable :: list(:)
      real*8 pol,thl
      real*8 sixth
      logical header
      character*1 digit
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(n))
c
c     set defaults for numbers and lists of polarizable atoms
c
      nlist = 0
      do i = 1, n
         list(i) = 0
      end do
c
c     set defaults for PCG induced dipole parameters
c
      pcgprec = .true.
      pcgguess = .true.
      pcgpeek = 1.0d0
c
c     set defaults for TCG induced dipole parameters
c
      tcgorder = 0
      tcgprec = .true.
      tcgguess = .true.
      tcgpeek = 1.0d0
      if (poltyp .eq. 'TCG   ')  poltyp = 'TCG2  '
      if (poltyp .eq. 'TCG0  ') then
         poltyp = 'DIRECT'
      else if (poltyp .eq. 'TCG1  ') then
         tcgorder = 1
      else if (poltyp .eq. 'TCG2  ') then
         tcgorder = 2
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(copt))  deallocate (copt)
      if (allocated(copm))  deallocate (copm)
      allocate (copt(0:maxopt))
      allocate (copm(0:maxopt))
c
c     set defaults for OPT induced dipole coefficients
c
      do i = 0, maxopt
         copt(i) = 0.0d0
         copm(i) = 0.0d0
      end do
      if (poltyp .eq. 'OPT   ')  poltyp = 'OPT4  '
      if (poltyp .eq. 'OPT1  ') then
         copt(0) = 0.530d0
         copt(1) = 0.604d0
      else if (poltyp .eq. 'OPT2  ') then
         copt(0) = 0.042d0
         copt(1) = 0.635d0
         copt(2) = 0.414d0
      else if (poltyp .eq. 'OPT3  ') then
         copt(0) = -0.132d0
         copt(1) = 0.218d0
         copt(2) = 0.637d0
         copt(3) = 0.293d0
      else if (poltyp .eq. 'OPT4  ') then
         copt(0) = -0.071d0
         copt(1) = -0.096d0
         copt(2) = 0.358d0
         copt(3) = 0.587d0
         copt(4) = 0.216d0
      else if (poltyp .eq. 'OPT5  ') then
         copt(0) = -0.005d0
         copt(1) = -0.129d0
         copt(2) = -0.026d0
         copt(3) = 0.465d0
         copt(4) = 0.528d0
         copt(5) = 0.161d0
      else if (poltyp .eq. 'OPT6  ') then
         copt(0) = 0.014d0
         copt(1) = -0.041d0
         copt(2) = -0.172d0
         copt(3) = 0.073d0
         copt(4) = 0.535d0
         copt(5) = 0.467d0
         copt(6) = 0.122d0
      end if
c
c     get keywords containing polarization-related options
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:12) .eq. 'POLARIZABLE ') then
            read (string,*,err=10,end=10)  (list(i),i=nlist+1,n)
   10       continue
            do while (list(nlist+1) .ne. 0)
               nlist = nlist + 1
            end do
         else if (keyword(1:12) .eq. 'PCG-PRECOND ') then
            pcgprec = .true.
         else if (keyword(1:14) .eq. 'PCG-NOPRECOND ') then
            pcgprec = .false.
         else if (keyword(1:10) .eq. 'PCG-GUESS ') then
            pcgguess = .true.
         else if (keyword(1:12) .eq. 'PCG-NOGUESS ') then
            pcgguess = .false.
         else if (keyword(1:9) .eq. 'PCG-PEEK ') then
            read (string,*,err=20,end=20)  pcgpeek
         else if (keyword(1:12) .eq. 'TCG-PRECOND ') then
            tcgprec = .true.
         else if (keyword(1:14) .eq. 'TCG-NOPRECOND ') then
            tcgprec = .false.
         else if (keyword(1:10) .eq. 'TCG-GUESS ') then
            tcgguess = .true.
         else if (keyword(1:12) .eq. 'TCG-NOGUESS ') then
            tcgguess = .false.
         else if (keyword(1:9) .eq. 'TCG-PEEK ') then
            read (string,*,err=20,end=20)  tcgpeek
         else if (keyword(1:10) .eq. 'OPT-COEFF ') then
            do i = 0, maxopt
               copt(i) = 0.0d0
            end do
            read (string,*,err=20,end=20)  (copt(i),i=0,maxopt)
         end if
   20    continue
      end do
c
c     get maximum coefficient order for OPT induced dipoles
c
      if (poltyp(1:3) .eq. 'OPT') then
         coptmax = 0
         do i = 1, maxopt
            if (copt(i) .ne. 0.0d0)  coptmax = max(i,coptmax)
         end do
         size = 1
         call numeral (coptmax,digit,size)
         poltyp = 'OPT'//digit//'  '
         do i = 0, coptmax
            do j = coptmax, i, -1
               copm(i) = copm(i) + copt(j)
            end do
         end do
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ipolar))  deallocate (ipolar)
      if (allocated(polarity))  deallocate (polarity)
      if (allocated(thole))  deallocate (thole)
      if (allocated(pdamp))  deallocate (pdamp)
      if (allocated(udir))  deallocate (udir)
      if (allocated(udirp))  deallocate (udirp)
      if (allocated(uind))  deallocate (uind)
      if (allocated(uinp))  deallocate (uinp)
      if (allocated(douind))  deallocate (douind)
      allocate (ipolar(n))
      allocate (polarity(n))
      allocate (thole(n))
      allocate (pdamp(n))
      allocate (udir(3,n))
      allocate (udirp(3,n))
      allocate (uind(3,n))
      allocate (uinp(3,n))
      allocate (douind(n))
      if (allocated(uopt))  deallocate (uopt)
      if (allocated(uoptp))  deallocate (uoptp)
      if (allocated(fopt))  deallocate (fopt)
      if (allocated(foptp))  deallocate (foptp)
      if (poltyp(1:3) .eq. 'OPT') then
         allocate (uopt(0:coptmax,3,n))
         allocate (uoptp(0:coptmax,3,n))
         allocate (fopt(0:coptmax,10,n))
         allocate (foptp(0:coptmax,10,n))
      end if
c
c     set the atoms allowed to have nonzero induced dipoles
c
      do i = 1, n
         douind(i) = .true.
      end do
      i = 1
      do while (list(i) .ne. 0)
         if (i .eq. 1) then
            do j = 1, n
               douind(j) = .false.
            end do
         end if
         if (list(i).gt.0 .and. list(i).le.n) then
            j = list(i)
            if (.not. douind(j)) then
               douind(j) = .true.
            end if
         else if (list(i).lt.0 .and. list(i).ge.-n) then
            do j = abs(list(i)), abs(list(i+1))
               if (.not. douind(j)) then
                  douind(j) = .true.
               end if
            end do
            i = i + 1
         end if
         i = i + 1
      end do
c
c     perform dynamic allocation of some local arrays
c
      deallocate (list)
c
c     process keywords containing polarizability parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            thl = -1.0d0
            do j = 1, maxval
               pg(j) = 0
            end do
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=30,end=40)  pol,(pg(j),j=1,maxval)
            goto 40
   30       continue
            read (string,*,err=40,end=40)  pol,thl,(pg(j),j=1,maxval)
   40       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Atomic Dipole',
     &                       ' Polarizability Parameters :')
                  if (thl .ge. 0.0d0) then
                     write (iout,60)
   60                format (/,5x,'Atom Type',11x,'Alpha',7x,
     &                          'Thole',5x,'Group Atom Types',/)
                  else
                     write (iout,70)
   70                format (/,5x,'Atom Type',11x,'Alpha',5x,
     &                          'Group Atom Types',/)
                  end if
               end if
               if (k .le. maxtyp) then
                  polr(k) = pol
                  athl(k) = max(0.0d0,thl)
                  do j = 1, maxval
                     pgrp(j,k) = pg(j)
                     if (pg(j) .eq. 0) then
                        npg = j - 1
                        goto 80
                     end if
                  end do
   80             continue
                  if (.not. silent) then
                     if (thl .ge. 0.0d0) then
                        write (iout,90)  k,pol,thl,(pg(j),j=1,npg)
   90                   format (6x,i6,8x,f10.3,2x,f10.3,7x,20i5)
                     else
                        write (iout,100)  k,pol,(pg(j),j=1,npg)
  100                   format (6x,i6,8x,f10.3,7x,20i5)
                     end if
                  end if
               else
                  write (iout,110)
  110             format (/,' KPOLAR  --  Too many Dipole',
     &                       ' Polarizability Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     find and store the atomic dipole polarizability parameters
c
      do i = 1, n
         polarity(i) = polr(type(i))
         thole(i) = athl(type(i))
      end do
c
c     process keywords containing atom specific polarizabilities
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            thl = -1.0d0
            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:240)
               read (string,*,err=120,end=120)  pol,thl
  120          continue
               if (header) then
                  header = .false.
                  write (iout,130)
  130             format (/,' Additional Dipole Polarizabilities',
     &                       ' for Specific Atoms :')
                  if (thl .ge. 0.0d0) then
                     write (iout,140)
  140                format (/,6x,'Atom',15x,'Alpha',7x,'Thole',/)
                  else
                     write (iout,150)
  150                format (/,6x,'Atom',15x,'Alpha',/)
                  end if
               end if
               if (.not. silent) then
                  if (thl .ge. 0.0d0) then
                     write (iout,160)  k,pol,thl
  160                format (6x,i6,8x,f10.3,2x,f10.3)
                  else
                     write (iout,170)  k,pol
  170                format (6x,i6,8x,f10.3)
                  end if
               end if
               polarity(k) = pol
               thole(k) = max(0.0d0,thl)
            end if
         end if
      end do
c
c     remove zero and undefined polarizable sites from the list
c
      npolar = 0
      if (use_polar .or. use_solv) then
         npole = 0
         do i = 1, n
            if (polarity(i) .eq. 0.0d0)  douind(i) = .false.
            if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0) then
               npole = npole + 1
               ipole(npole) = i
               pollist(i) = npole
               zaxis(npole) = zaxis(i)
               xaxis(npole) = xaxis(i)
               yaxis(npole) = yaxis(i)
               polaxe(npole) = polaxe(i)
               do k = 1, maxpole
                  pole(k,npole) = pole(k,i)
               end do
               if (palpha(i) .ne. 0.0d0)  ncp = ncp + 1
               pcore(npole) = pcore(i)
               pval(npole) = pval(i)
               palpha(npole) = palpha(i)
               if (polarity(i) .ne. 0.0d0) then
                  npolar = npolar + 1
                  ipolar(npolar) = npole
                  douind(i) = .true.
               end if
               polarity(npole) = polarity(i)
               thole(npole) = thole(i)
            end if
         end do
      end if
c
c     set the values used in the scaling of the polarizability
c
      sixth = 1.0d0 / 6.0d0
      do i = 1, npole
         if (thole(i) .eq. 0.0d0) then
            pdamp(i) = 0.0d0
         else
            pdamp(i) = polarity(i)**sixth
         end if
      end do
c
c     assign polarization group connectivity of each atom
c
      call polargrp
c
c     test multipoles at chiral sites and invert if necessary
c
      call chkpole
c
c     turn off polarizable multipole potentials if not used
c
      if (npole .eq. 0)  use_mpole = .false.
      if (ncp .ne. 0)  use_chgpen = .true.
      if (npolar .eq. 0)  use_polar = .false.
      if (use_polar .and. ncp.eq.0)  use_thole = .true.
c
c     perform dynamic allocation of some global arrays
c
      if (use_polar) then
         if (allocated(mindex))  deallocate (mindex)
         if (allocated(minv))  deallocate (minv)
         allocate (mindex(npole))
         allocate (minv(3*maxulst*npole))
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine polargrp  --  polarization group connectivity  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "polargrp" generates members of the polarization group of
c     each atom and separate lists of the 1-2, 1-3 and 1-4 group
c     connectivities
c
c
      subroutine polargrp
      use atoms
      use couple
      use inform
      use iounit
      use kpolr
      use mpole
      use polgrp
      implicit none
      integer i,j,k,m
      integer it,jt
      integer jj,kk
      integer start,stop
      integer nkeep,nlist
      integer maxkeep,maxlist
      integer, allocatable :: keep(:)
      integer, allocatable :: list(:)
      integer, allocatable :: mask(:)
      logical done
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(np11))  deallocate (np11)
      if (allocated(np12))  deallocate (np12)
      if (allocated(np13))  deallocate (np13)
      if (allocated(np14))  deallocate (np14)
      if (allocated(ip11))  deallocate (ip11)
      if (allocated(ip12))  deallocate (ip12)
      if (allocated(ip13))  deallocate (ip13)
      if (allocated(ip14))  deallocate (ip14)
      allocate (np11(n))
      allocate (np12(n))
      allocate (np13(n))
      allocate (np14(n))
      allocate (ip11(maxp11,n))
      allocate (ip12(maxp12,n))
      allocate (ip13(maxp13,n))
      allocate (ip14(maxp14,n))
c
c     find the directly connected group members for each atom
c
      do i = 1, n
         np11(i) = 1
         ip11(1,i) = i
         it = type(i)
         do j = 1, n12(i)
            jj = i12(j,i)
            jt = type(jj)
            do k = 1, maxval
               kk = pgrp(k,it)
               if (kk .eq. 0)  goto 20
               if (pgrp(k,it) .eq. jt) then
                  np11(i) = np11(i) + 1
                  if (np11(i) .le. maxp11) then
                     ip11(np11(i),i) = jj
                  else
                     write (iout,10)
   10                format (/,' POLARGRP  --  Too many Atoms',
     &                          ' in Polarization Group')
                     abort = .true.
                     goto 30
                  end if
               end if
            end do
   20       continue
         end do
      end do
   30 continue
c
c     make sure all connected group members are bidirectional
c
      do i = 1, n
         do j = 1, np11(i)
            k = ip11(j,i)
            do m = 1, np11(k)
               if (ip11(m,k) .eq. i)  goto 50
            end do
            write (iout,40)  min(i,k),max(i,k)
   40       format (/,' POLARGRP  --  Check Polarization Groups for',
     &                 ' Atoms',i9,' and',i9)
            abort = .true.
   50       continue
         end do
      end do
      if (abort)  call fatal
c
c     perform dynamic allocation of some local arrays
c
      maxkeep = 100
      maxlist = 10000
      allocate (keep(maxkeep))
      allocate (list(maxlist))
      allocate (mask(n))
c
c     find any other group members for each atom in turn
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         done = .false.
         start = 1
         stop = np11(i)
         do j = start, stop
            jj = ip11(j,i)
            if (jj .lt. i) then
               done = .true.
               np11(i) = np11(jj)
               do k = 1, np11(i)
                  ip11(k,i) = ip11(k,jj)
               end do
            else
               mask(jj) = i
            end if
         end do
         do while (.not. done)
            done = .true.
            do j = start, stop
               jj = ip11(j,i)
               do k = 1, np11(jj)
                  kk = ip11(k,jj)
                  if (mask(kk) .ne. i) then
                     np11(i) = np11(i) + 1
                     if (np11(i) .le. maxp11) then
                        ip11(np11(i),i) = kk
                     else
                        write (iout,60)
   60                   format (/,' POLARGRP  --  Too many Atoms',
     &                             ' in Polarization Group')
                        abort = .true.
                        goto 70
                     end if
                     mask(kk) = i
                  end if
               end do
            end do
            if (np11(i) .ne. stop) then
               done = .false.
               start = stop + 1
               stop = np11(i)
            end if
         end do
         call sort (np11(i),ip11(1,i))
      end do
   70 continue
c
c     loop over atoms finding all the 1-2 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         nkeep = 0
         do j = 1, np11(i)
            jj = ip11(j,i)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               if (mask(kk) .ne. i) then
                  nkeep = nkeep + 1
                  keep(nkeep) = kk
               end if
            end do
         end do
         nlist = 0
         do j = 1, nkeep
            jj = keep(j)
            do k = 1, np11(jj)
               kk = ip11(k,jj)
               nlist = nlist + 1
               list(nlist) = kk
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp12) then
            np12(i) = nlist
            do j = 1, nlist
               ip12(j,i) = list(j)
            end do
         else
            write (iout,80)
   80       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-2 Polarization Group')
            abort = .true.
            goto 90
         end if
      end do
   90 continue
c
c     loop over atoms finding all the 1-3 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np12(i)
            jj = ip12(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp13) then
            np13(i) = nlist
            do j = 1, nlist
               ip13(j,i) = list(j)
            end do
         else
            write (iout,100)
  100       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-3 Polarization Group')
            abort = .true.
            goto 110
         end if
      end do
  110 continue
c
c     loop over atoms finding all the 1-4 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         do j = 1, np13(i)
            jj = ip13(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np13(i)
            jj = ip13(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp14) then
            np14(i) = nlist
            do j = 1, nlist
               ip14(j,i) = list(j)
            end do
         else
            write (iout,120)
  120       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-4 Polarization Group')
            abort = .true.
            goto 130
         end if
      end do
  130 continue
c
c     perform deallocation of some local arrays
c
      deallocate (keep)
      deallocate (list)
      deallocate (mask)
      return
      end
