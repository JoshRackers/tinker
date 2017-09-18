c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program exchfit  --  fit parameters to structure & energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "exchfit" determines optimized van der Waals and electrostatic
c     parameters by fitting to crystal structures, lattice energies,
c     and dimer structures and interaction energies
c
c
      program exchfit
      use sizes
      use bound
      use boxes
      use files
      use iounit
      use molcul
      use potent
      use vdwpot
      use xtals
      implicit none
      integer i,ixtal
      integer atom1,atom2
      integer nresid,prmtyp
      real*8 grdmin
      real*8 exchmin1
      real*8, allocatable :: xx(:)
      real*8, allocatable :: resid(:)
      real*8, allocatable :: g(:)
      real*8, allocatable :: xlo(:)
      real*8, allocatable :: xhi(:)
      real*8, allocatable :: fjac(:,:)
      logical exist,query
      character*5 vindex
      character*16 label(15)
      character*120 record
      character*120 string
      character*10 optimizer
      external xtalerr,xtalwrt
      external exchmin1,paramwrt1
c
c
c     initialize some variables to be used during fitting
c
      call initial
      nvary = 0
      nresid = 0
c
c     print informational header about available parameters
c
      write (iout,10)
   10 format (/,' The Following Parameters can be Fit for',
     &           ' each Atom Type :',
     &        //,4x,'(1) Exchange-Repulsion A',
     &        /,4x,'(2) Exchange-Repulsion B/r',
     &        /,4x,'(3) Exchange-Repulsion sigma',
     &        /,4x,'(4) Gordon Damping Parameter',
     &        /,4x,'(5) Single Ex-Rep Scalar',
     &        /,4x,'(6) C6-COEFFICIENT',
     &        /,4x,'(7) Exchange-Repulsion C/r2',
     &        /,4x,'(8) Hydrogen Reduction Factor',
     &        /,4x,'(11) vdW radius',
     &        /,4x,'(12) vdW eps',
     &        /,4x,'(13) Number of Valence Electrons',
     &        /,4x,'(14) Exchange Alpha',
     &        /,4x,'(15) Exchange Valence Electrons')
c
c     get types of potential parameters to be optimized
c
      query = .true.
      do while (query)
         prmtyp = -1
         atom1 = 0
         atom2 = 0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=20,end=20)  prmtyp
         call nextarg (string,exist)
         if (exist)  read (string,*,err=20,end=20)  atom1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=20,end=20)  atom2
   20    continue
         if (prmtyp .ne. 0) then
            prmtyp = 0
            write (iout,30)
   30       format (/,' Enter Parameter Type then Atom Class',
     &                 ' or Type(s) :  ',$)
            read (input,40)  record
   40       format (a120)
            read (record,*,err=50,end=50)  prmtyp,atom1,atom2
   50       continue
         end if
         if (prmtyp .eq. 0) then
            query = .false.
         else
            query = .true.
            nvary = nvary + 1
            ivary(nvary) = prmtyp
            vary(1,nvary) = atom1
c            if (prmtyp.eq.5 .or. prmtyp.eq.6) then
c               vary(1,nvary) = min(atom1,atom2)
c               vary(2,nvary) = max(atom1,atom2)
c            end if
         end if
      end do
c
c     determine which optimizer will be used
c
      write (iout,51)
 51   format (/,' Enter type of optimizer to use :  ',$)
      read (input,52)  optimizer
 52   format (A10)
      write (iout,53)
 53   format (/,' Are you fitting an exponential?  ',$)
      read (input,54)  exponential
 54   format (A10)
c
c     get termination criterion as RMS gradient over parameters
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  grdmin
   60 continue
      if (grdmin .le. 0.0d0) then
         write (iout,70)
   70    format (/,' Enter RMS Gradient Termination Criterion',
     &              ' [0.1] :  ',$)
         read (input,80)  grdmin
   80    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 0.1d0
c
c     get the number of structures to use in optimization
c
      nxtal = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  nxtal
   90 continue
      if (nxtal .le. 0) then
         write (iout,100)
  100    format (/,' Enter Number of Structures to be Used [1] :  ',$)
         read (input,110)  nxtal
  110    format (i10)
      end if
c
c     check for too few or too many molecular structures
c
      if (nxtal .eq. 0)  nxtal = 1
      if (nxtal .gt. maxref) then
         write (iout,120)
  120    format (/,' XTALFIT  --  Too many Structures,',
     &              ' Increase the Value of MAXREF')
         call fatal
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(nvary))
c
c     get coordinates and parameters for current structure
c
      do ixtal = 1, nxtal
         call initial
         call getxyz
         call mechanic
c
c     get ideal value for intermolecular or lattice energy
c
         e0_exchange = 0.0d0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=130,end=130)  e0_exchange
            query = .false.
         end if
  130    continue
         if (query) then
            write (iout,140)
  140       format (/,' Enter Target Exchange Value',
     &                 ' [<CR> to omit] :  ',$)
            read (input,150)  e0_exchange
  150       format (f20.0)
         end if
cccccccccccccccccccccccccc
         e0_mpole = 0.0d0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=151,end=151)  e0_mpole
            query = .false.
         end if
 151     continue
         if (query) then
            write (iout,152)
 152        format (/,' Enter Target MPOLE Value',
     &           ' [<CR> to omit] :  ',$)
            read (input,153)  e0_mpole
 153        format (f20.0)
         end if
ccccccccccccccccccccccccccccc
         e0_disp = 0.0d0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=154,end=154)  e0_disp
            query = .false.
         end if
 154     continue
         if (query) then
            write (iout,155)
 155        format (/,' Enter Target Dispersion Value',
     &           ' [<CR> to omit] :  ',$)
            read (input,156)  e0_disp
 156        format (f20.0)
         end if
ccccccccccccccccccccccccc
         e0_vdw = 0.0d0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=157,end=157)  e0_vdw
            query = .false.
         end if
 157     continue
         if (query) then
            write (iout,158)
c     co-opted this to do total
c 158        format (/,' Enter Any Non-zero number for VDW',
c     &           ' [<CR> to omit] :  ',$)
 158        format (/,' Enter Total Interaction Energy',
     &           ' [<CR> to omit] :  ',$)
            read (input,159)  e0_vdw
 159        format (f20.0)
         end if
c
         e0_xfer = 0.0d0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=160,end=160)  e0_xfer
            query = .false.
         end if
 160     continue
         if (query) then
            write (iout,161)
 161        format (/,' Enter Target Charge Transfer Value',
     &           ' [<CR> to omit] :  ',$)
            read (input,162)  e0_xfer
 162        format (f20.0)
         end if
ccccccccccccccccccccccccc
c
c     print molecules per structure, energy and dipole values
c
         write (iout,169)  ixtal,filename(1:35),nmol
 169     format (/,' File Name of Target Structure',i4,' :',8x,a35,
     &        /,' Number of Molecules per Structure :',i13)
         if (e0_exchange .ne. 0.0d0) then
            nresid = nresid + 1
            iresid(nresid) = ixtal
            if (use_bounds) then
               rsdtyp(nresid) = 'Exchange Energy'
            else
               rsdtyp(nresid) = 'E Exchange     '
            end if
            write (iout,170)  e0_exchange
  170       format (' Target E-Exchange or E-Inter Value :  ',f13.2)
         end if
ccccccccccccccccccccccc
         if (e0_mpole .ne. 0.0d0) then
            nresid = nresid + 1
            iresid(nresid) = ixtal
            if (use_bounds) then
               rsdtyp(nresid) = 'MPOLE     Energy'
            else
               rsdtyp(nresid) = 'E MPOLE         '
            end if
            write (iout,171)  e0_mpole
 171        format (' Target E-MPOLE or E-Inter Value :  ',f13.2)
         end if
cccccccccccccccccccccccc
         if (e0_disp .ne. 0.0d0) then
            nresid = nresid + 1
            iresid(nresid) = ixtal
            if (use_bounds) then
               rsdtyp(nresid) = 'DISP     Energy'
            else
               rsdtyp(nresid) = 'E DISP         '
            end if
            write (iout,172)  e0_disp
 172        format (' Target E-DISP or E-Inter Value :  ',f13.2)
         end if
cccccccccccccccccccccccccc
         if (e0_vdw .ne. 0.0d0) then
            nresid = nresid + 1
            iresid(nresid) = ixtal
            if (use_bounds) then
               rsdtyp(nresid) = 'VDW      Energy'
            else
               rsdtyp(nresid) = 'E VDW          '
            end if
            write (iout,173)  e0_vdw
 173        format (' Target E-VDW or E-Inter Value :  ',f13.2)
         end if
ccccccccccccccccccccccc
         if (e0_xfer .ne. 0.0d0) then
            nresid = nresid + 1
            iresid(nresid) = ixtal
            if (use_bounds) then
               rsdtyp(nresid) = 'Charge Transfer Energy'
            else
               rsdtyp(nresid) = 'E ChgTransfer  '
            end if
            write (iout,174)  e0_xfer
 174        format (' Target E-XFER or E-Inter Value :  ',f13.2)
         end if
ccccccccccccccccccccccccc
c
c     set the initial values of the parameters
c
         call xtalprm ('STORE',ixtal,xx)
      end do
c
c     turn off all local interactions
c
      call potoff
      use_vdw = .true.
      use_charge = .true.
      use_chgdpl = .true.
      use_dipole = .true.
      use_mpole = .true.
      use_polar = .true.
      use_extra = .true.
c
c     types of variables for use in optimization
c
      label(1) = 'EX-OVERLAP    '
      label(2) = 'EX-BOVERLAP   '
      label(3) = 'EX-OVERLAPR   '
      label(4) = 'ALPHA-PERM    '
      label(5) = 'Single Scalar '
      label(6) = 'C6 Coefficient'
      label(7) = 'Overlap C/r2  '
      label(8) = 'Hydrogen Red  '
      label(11) = 'vdw radius    '
      label(12) = 'vdw eps       '
      label(13) = 'VALENCE-ELE   '
      label(14) = 'ALPHA-EXCH    '
      label(15) = 'VALENCE-ELE-EXCH '
      do i = 1, nvary
         vartyp(i) = label(ivary(i))
      end do
      vindex = 'Class'
      if (vdwindex .eq. 'TYPE ')  vindex = 'Type '
c
c     print the initial parameter values
c
      write (iout,180)
  180 format (/,' Initial Values of the Parameters :',/)
      do i = 1, nvary
         write (iout,190)  i,vartyp(i),vindex,vary(1,i),xx(i)
 190     format (3x,'(',i2,')',2x,a16,4x,'Atom ',a5,i5,4x,f12.4)
c         if (ivary(i) .le. 3) then
c            write (iout,190)  i,vartyp(i),vindex,vary(1,i),xx(i)
c  190       format (3x,'(',i2,')',2x,a16,4x,'Atom ',a5,i5,4x,f12.4)
c         else if (ivary(i) .ne. 6) then
c            write (iout,200)  i,vartyp(i),vary(1,i),xx(i)
c  200       format (3x,'(',i2,')',2x,a16,4x,'Atom Type ',i5,4x,f12.4)
c         else
c            write (iout,210)  i,vartyp(i),vary(1,i),vary(2,i),xx(i)
c  210       format (3x,'(',i2,')',2x,a16,4x,'Bond Type ',2i5,f12.4)
c         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (resid(nresid))
      allocate (g(nvary))
      allocate (xlo(nvary))
      allocate (xhi(nvary))
      allocate (fjac(nresid,nvary))
c
c     set upper and lower bounds based on the parameter type
c
      do i = 1, nvary
         if (xx(i).ge.0.0d0) then
            xlo(i) = 0.25d0 * xx(i)
            xhi(i) = 4.0d0 * xx(i)
         else
            xlo(i) = 4.0d0 * xx(i)
            xhi(i) = 0.25d0 * xx(i)
         end if
         if (ivary(i).eq.8) then
            xlo(i) = 0.8d0
            xhi(i) = 1.0d0
         end if
         if (ivary(i).eq.13) then
            xlo(i) = 0.0d0
            print *,"WARNING: Using vdwclasses"
            if (vary(1,i).ge.1) xhi(i) = 1.0d0
            if (vary(1,i).ge.7) xhi(i) = 6.0d0
            if (vary(1,i).ge.14) xhi(i) = 7.0d0
            if (vary(1,i).ge.17) xhi(i) = 8.0d0
            if (vary(1,i).ge.23) xhi(i) = 15.0d0
            if (vary(1,i).ge.24) xhi(i) = 16.0d0
            if (vary(1,i).ge.26) xhi(i) = 9.0d0
            if (vary(1,i).ge.27) xhi(i) = 17.0d0
            if (vary(1,i).ge.28) xhi(i) = 35.0d0
         end if
      end do
c
c     use nonlinear least squares to refine the parameters
c
      if (optimizer .eq. "lbfgs") then
         call lbfgs (nvary,xx,resid,grdmin,exchmin1,paramwrt1)
      else
         call square (nvary,nresid,xlo,xhi,xx,resid,g,fjac,
     &        grdmin,xtalerr,xtalwrt)
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (xlo)
      deallocate (xhi)
      deallocate (fjac)
c
c     print final values of parameters and scaled derivatives
c
      write (iout,220)
  220 format (/,' Final Values of Parameters and Scaled',
     &           ' Derivatives :',/)
      do i = 1, nvary
         write (iout,230)  i,vartyp(i),vindex,vary(1,i),xx(i),g(i)
 230     format (3x,'(',i2,')',2x,a16,4x,'Atom ',a5,i5,2x,2f14.4)
c         if (ivary(i) .le. 3) then
c            write (iout,230)  i,vartyp(i),vindex,vary(1,i),xx(i),g(i)
c  230       format (3x,'(',i2,')',2x,a16,4x,'Atom ',a5,i5,2x,2f14.4)
c         else if (ivary(i) .ne. 6) then
c            write (iout,240)  i,vartyp(i),vary(1,i),xx(i),g(i)
c  240       format (3x,'(',i2,')',2x,a16,4x,'Atom Type ',i5,2x,2f14.4)
c         else
c            write (iout,250)  i,vartyp(i),vary(1,i),vary(2,i),xx(i),g(i)
c  250       format (3x,'(',i2,')',2x,a16,4x,'Bond Type ',2i5,
c     &                 f11.4,f14.4)
c         end if
      end do
c
c     print final values of the individual residual functions
c
      write (iout,260)
  260 format (/,' Final Residual Error Function Values :',/)
      do i = 1, nresid
         if (i .lt. 100) then
            write (iout,270)  i,rsdtyp(i),iresid(i),resid(i)
  270       format (3x,'(',i2,')',2x,a16,6x,'Structure',i4,4x,f12.4)
         else
            write (iout,280)  i,rsdtyp(i),iresid(i),resid(i)
  280       format (2x,'(',i3,')',2x,a16,6x,'Structure',i4,4x,f12.4)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      deallocate (resid)
      deallocate (g)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine xtalprm  --  energy/optimization conversion  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "xtalprm" stores or retrieves a molecular structure; used to
c     make a previously stored structure the active structure, or to
c     store a structure for later use
c
c     the current version only provides for intermolecular potential
c     energy terms
c
c
      subroutine xtalprm (mode,ixtal,xx)
      use sizes
      use atoms
      use atomid
      use bound
c      use boxes
      use charge
      use dipole
      use files
      use fracs
      use inform
      use kvdws
      use molcul
      use mpole
      use polar
      use potderivs
      use vdw
      use vdwpot
      use xtals
      use chgpen
      implicit none
      integer i,j,k
      integer ixtal,prmtyp
      integer jclass
      real*8 rd,ep,sixth
      real*8 xmid,ymid,zmid
      real*8 e0_exchanges(maxref)
      real*8 e0_mpoles(maxref)
      real*8 e0_disps(maxref)
      real*8 e0_vdws(maxref)
      real*8 e0_xfers(maxref)
      real*8 xx(*)
      logical first
      character*5 mode
      save e0_exchanges
      save e0_mpoles
      save e0_disps
      save e0_vdws
      save e0_xfers
      save first
      data first  / .true. /
c
c
c     save or restore the key values for the current crystal
c
      if (mode .eq. 'STORE') then
         call makeref (ixtal)
      else if (mode .eq. 'RESET') then
      if (allocated(pot)) deallocate (pot)
      if (allocated(field)) deallocate (field)
      if (allocated(gradfield)) deallocate (gradfield)
c
c     permanent - permanent exclusion rule
c
      if (allocated(potm)) deallocate (potm)
      if (allocated(fieldm)) deallocate (fieldm)
      if (allocated(gradfieldm)) deallocate (gradfieldm)
c
c     ewald damping
c
      if (allocated(pot_ewald)) deallocate (pot_ewald)
      if (allocated(field_ewald)) deallocate (field_ewald)
      if (allocated(gradfield_ewald)) 
     &     deallocate (gradfield_ewald)
c
c     thole damping with p and d exclusion rules needed for
c     computing induced dipoles
c
      if (allocated(fieldd_thole)) deallocate (fieldd_thole)
      if (allocated(fieldp_thole)) deallocate (fieldp_thole)
c
      if (allocated(fieldd_func)) deallocate (fieldd_func)
      if (allocated(fieldp_func)) deallocate (fieldp_func)
c
c     gordon charge penetration damping
c
      if (allocated(potm_gordon)) deallocate (potm_gordon)
      if (allocated(fieldm_gordon)) 
     &     deallocate (fieldm_gordon)
      if (allocated(gradfieldm_gordon))
     &     deallocate (gradfieldm_gordon)
c
c     gordon overlap damping
c
      if (allocated(potm_goverlap)) deallocate (potm_goverlap)
      if (allocated(fieldm_goverlap))
     &     deallocate (fieldm_goverlap)
      if (allocated(gradfieldm_goverlap))
     &     deallocate (gradfieldm_goverlap)
c
c     gordon charge penetration damping - nuclei
c
      if (allocated(nucpotm_gordon))
     &     deallocate (nucpotm_gordon)
      if (allocated(nucfieldm_gordon))
     &     deallocate (nucfieldm_gordon)
c
c     gordon damping for polarization
c
      if (allocated(fieldd_gordon))
     &     deallocate (fieldd_gordon)
      if (allocated(fieldp_gordon))
     &     deallocate (fieldp_gordon)
c
c     gordon damping with nuclear regularization
c
      if (allocated(fieldd_gordonreg))
     &     deallocate (fieldd_gordonreg)
      if (allocated(fieldp_gordonreg))
     &     deallocate (fieldp_gordonreg)
c
         call getref (ixtal)
         call basefile (filename)
         silent = .true.
         call mechanic
         silent = .false.
         if (use_bounds)  call bounds
      end if
c
c     perform dynamic allocation of some global arrays
c
c      if (mode .eq. 'RESET') then
c         if (first) then
c            first = .false.
c            allocate (xfrac(nmol))
c            allocate (yfrac(nmol))
c            allocate (zfrac(nmol))
c         end if
c      end if
c
c     coordinates of molecular centers of mass
c
c      if (mode .eq. 'RESET') then
c         do i = 1, nmol
c            xmid = 0.0d0
c            ymid = 0.0d0
c            zmid = 0.0d0
c            do j = imol(1,i), imol(2,i)
c               k = kmol(j)
c               xmid = xmid + x(k)*mass(k)
c               ymid = ymid + y(k)*mass(k)
c               zmid = zmid + z(k)*mass(k)
c            end do
c            zmid = zmid / gamma_term
c            ymid = (ymid - zmid*beta_term) / gamma_sin
c            xmid = xmid - ymid*gamma_cos - zmid*beta_cos
c            xfrac(i) = xmid / (xbox * molmass(i))
c            yfrac(i) = ymid / (ybox * molmass(i))
c            zfrac(i) = zmid / (zbox * molmass(i))
c         end do
c      end if
c
c     values of ideal intermolecular or lattice energy
c
      if (mode .eq. 'STORE') then
         e0_exchanges(ixtal) = e0_exchange
         e0_mpoles(ixtal) = e0_mpole
         e0_disps(ixtal) = e0_disp
         e0_vdws(ixtal) = e0_vdw
         e0_xfers(ixtal) = e0_xfer
      else if (mode .eq. 'RESET') then
         e0_exchange = e0_exchanges(ixtal)
         e0_mpole = e0_mpoles(ixtal)
         e0_disp = e0_disps(ixtal)
         e0_vdw = e0_vdws(ixtal)
         e0_xfer = e0_xfers(ixtal)
      end if
c
c     store or reset values of the optimization variables
c
      do j = 1, nvary
         prmtyp = ivary(j)
         jclass = vary(1,j)
ccccccccccccccccccccccccccccccccc
         if (prmtyp .eq. 1) then
            if (mode .eq. 'STORE') then
c               xx(j) = overlap(jclass)
               xx(j) = overlap(jclass)
            else if (mode .eq. 'RESET') then
c               overlap(jclass) = xx(j)
               overlap(jclass) = xx(j)
            end if
         else if (prmtyp .eq. 2) then
            if (mode .eq. 'STORE') then
               xx(j) = boverlap(jclass)
            else if (mode .eq. 'RESET') then
               boverlap(jclass) = xx(j)
            end if
         else if (prmtyp .eq. 3) then
            if (mode .eq. 'STORE') then
               xx(j) = overlapr(jclass)
            else if (mode .eq. 'RESET') then
               overlapr(jclass) = xx(j)
            end if
         else if (prmtyp .eq. 4) then
            if (mode .eq. 'STORE') then
c               xx(j) = alpha(jclass)
               xx(j) = alpha(jclass)
            else if (mode .eq. 'RESET') then
c               alpha(jclass) = xx(j)
               alpha(jclass) = xx(j)
            end if
         else if (prmtyp .eq. 5) then
            if (mode .eq. 'STORE') then
               xx(j) = soverlap
            else if (mode .eq. 'RESET') then
               soverlap = xx(j)
            end if
         else if (prmtyp .eq. 6) then
            if (mode .eq. 'STORE') then
               xx(j) = csix(jclass)
            else if (mode .eq. 'RESET') then
               csix(jclass) = xx(j)
            end if
         else if (prmtyp .eq. 7) then
            if (mode .eq. 'STORE') then
               xx(j) = coverlap(jclass)
            else if (mode .eq. 'RESET') then
               coverlap(jclass) = xx(j)
            end if
         else if (prmtyp .eq. 8) then
            if (mode .eq. 'STORE') then
               do i = 1, n
                  if (class(i) .eq. jclass) then
                     xx(j) = kred(i)
                     goto 10
                  end if
               end do
            else if (mode .eq. 'RESET') then
               do i = 1, n
                  if (class(i) .eq. jclass)  kred(i) = xx(j)
               end do
            end if
c     vdw radius and eps
         else if (prmtyp .eq. 11) then
            if (mode .eq. 'STORE') then
               xx(j) = rad(jclass)
            else if (mode .eq. 'RESET') then
               rad(jclass) = xx(j)
               do i = 1, maxclass
                  if (rad(i).eq.0.0d0 .and. rad(jclass).eq.0.0d0) then
                     rd = 0.0d0
                  else if (radrule(1:10) .eq. 'ARITHMETIC') then
                     rd = rad(i) + rad(jclass)
                  else if (radrule(1:9) .eq. 'GEOMETRIC') then
                     rd = 2.0d0 * sqrt(rad(i) * rad(jclass))
                  else if (radrule(1:10) .eq. 'CUBIC-MEAN') then
                     rd = 2.0d0 * (rad(i)**3+rad(jclass)**3)
     &                       / (rad(i)**2+rad(jclass)**2)
                  else
                     rd = rad(i) + rad(jclass)
                  end if
                  radmin(i,jclass) = rd
                  radmin(jclass,i) = rd
               end do
            end if
         else if (prmtyp .eq. 12) then
            if (mode .eq. 'STORE') then
               xx(j) = eps(jclass)
            else if (mode .eq. 'RESET') then
               eps(jclass) = abs(xx(j))
               do i = 1, maxclass
                  if (eps(i).eq.0.0d0 .and. eps(jclass).eq.0.0d0) then
                     ep = 0.0d0
                  else if (epsrule(1:10) .eq. 'ARITHMETIC') then
                     ep = 0.5d0 * (eps(i) + eps(jclass))
                  else if (epsrule(1:9) .eq. 'GEOMETRIC') then
                     ep = sqrt(eps(i) * eps(jclass))
                  else if (epsrule(1:8) .eq. 'HARMONIC') then
                     ep = 2.0d0 * (eps(i)*eps(jclass))
     &                       / (eps(i)+eps(jclass))
                  else if (epsrule(1:3) .eq. 'HHG') then
                     ep = 4.0d0 * (eps(i)*eps(jclass))
     &                      / (sqrt(eps(i))+sqrt(eps(jclass)))**2
                  else if (epsrule(1:3) .eq. 'W-H') then
                     ep = 2.0d0 * sqrt(eps(i)) * sqrt(eps(jclass)) * 
     &                    (rad(i)*rad(jclass))**3 /
     &                    (rad(i)**6 + rad(jclass)**6)
                  else
                     ep = sqrt(eps(i) * eps(jclass))
                  end if
                  epsilon(i,jclass) = ep
                  epsilon(jclass,i) = ep
               end do
            end if
cccccccccccccccccccccccccccccccccc
         else if (prmtyp .eq. 13) then
            if (mode .eq. 'STORE') then
               xx(j) = val_ele(jclass)
            else if (mode .eq. 'RESET') then
               val_ele(jclass) = xx(j)
            end if
         else if (prmtyp .eq. 14) then
            if (mode .eq. 'STORE') then
               xx(j) = exch_alpha(jclass)
            else if (mode .eq. 'RESET') then
               exch_alpha(jclass) = xx(j)
            end if
         else if (prmtyp .eq. 15) then
            if (mode .eq. 'STORE') then
               xx(j) = exch_val_ele(jclass)
            else if (mode .eq. 'RESET') then
               exch_val_ele(jclass) = xx(j)
            end if
         end if
   10    continue
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine xtalerr  --  error function for xtalfit  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "xtalerr" computes an error function value derived from
c     lattice energies, dimer intermolecular energies and the
c     gradient with respect to structural parameters
c
c
      subroutine xtalerr (nvaried,nresid,xx,resid)
      use sizes
      use atoms
      use boxes
      use bound
      use charge
      use dipole
      use energi
      use limits
      use math
      use molcul
      use mpole
      use polar
      use vdw
      use xtals
c
      use chgpen
      use inter
c
      implicit none
      integer i,k,ixtal
      integer nresid,nvaried
      real*8 energy
      real*8 eps,temp
      real*8 e,e0
      real*8 e_monomer
      real*8 e_exchange
      real*8 e_mpole
      real*8 e_disp
      real*8 e_vdw
      real*8 e_cp
      real*8 e_xfer
      real*8 e0_exchange_log
      real*8 e0_xfer_log
      real*8 e0_mpole_log
      real*8 e0_disp_log
      real*8 e0_vdw_log
      real*8 dmol,big
      real*8 e1,e2,e3
      real*8 e4,e5,e6
      real*8 g1,g2,g3
      real*8 g4,g5,g6
      real*8 xx(*)
      real*8 resid(*)
c
c
c     zero out number of residuals and set numerical step size
c
      nresid = 0
      eps = 1.0d-4
c
c     set force field parameter values and find the base energy
c
      do ixtal = 1, nxtal
         call xtalprm ('RESET',ixtal,xx)
         call analysis (energy)
         if (nanflag) then
            print *,"dimer number: ",ixtal
            print *,"parameters"
            do i = 1, nvary
               print *,i,xx(i)
            end do
            call fatal
         end if
c         print *,"ixtal",ixtal,ex
         e_exchange = ex
c     this only works if intramolecular multipoles are turned off
         e_mpole = em
         e_disp = edis
         e_cp = ecp
ccccccccccccccccccccccccccc
c     co-opted for total
c         e_vdw = e_exchange + e_disp
         e_vdw = einter
ccccccccccccccccccccccccccc
         e_xfer = exfr
c     get exchange from ehal3
         if (exmodel .eq. "14-ONLY") e_exchange = ev
c
         if (exponential .eq. "YES") then
            e_exchange = log(ex)
            if (exmodel .eq. "14-ONLY") e_exchange = log(ev)
            e0_exchange_log = log(e0_exchange)
c
            e_xfer = log(-exfr)
            if (e0_xfer .le. 0.0d0) then
               e0_xfer_log = log(-e0_xfer)
            else
               e0_xfer_log = e_xfer
            end if
         end if
c
c     compute residual due to intermolecular or lattice energy;
c     weight energies more heavily, since there are fewer of them
c
         if (e0_exchange .ne. 0.0d0) then
            nresid = nresid + 1
            if (exponential .eq. "YES") then
               resid(nresid) = e_exchange - e0_exchange_log
            else
               resid(nresid) = e_exchange - e0_exchange
            end if
         end if
         if (e0_mpole .ne. 0.0d0) then
            nresid = nresid + 1
            resid(nresid) = e_mpole - e0_mpole
         end if
         if (e0_disp .ne. 0.0d0) then
            nresid = nresid + 1
            resid(nresid) = e_disp - e0_disp
         end if
         if (e0_vdw .ne. 0.0d0) then
            nresid = nresid + 1
            resid(nresid) = e_vdw - e0_vdw
         end if
         if (e0_xfer .ne. 0.0d0) then
            nresid = nresid + 1
            if (exponential .eq. "YES") then
               resid(nresid) = e_xfer - e0_xfer_log
            else
               resid(nresid) = e_xfer - e0_xfer
            end if
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine xtalmove  --  translation of rigid molecules  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "xtalmove" converts fractional to Cartesian coordinates for
c     rigid molecules during optimization of force field parameters
c
c
      subroutine xtalmove
      use sizes
      use atoms
      use atomid
      use boxes
      use fracs
      use molcul
      implicit none
      integer i,j,k
      integer init,stop
      real*8 weigh
      real*8 xmid,ymid,zmid
      real*8, allocatable :: xoff(:)
      real*8, allocatable :: yoff(:)
      real*8, allocatable :: zoff(:)
c
c
c     get values for fractional coordinate interconversion
c
      call lattice
c
c     perform dynamic allocation of some local arrays
c
      allocate (xoff(n))
      allocate (yoff(n))
      allocate (zoff(n))
c
c     locate the center of mass of each molecule
c
      do i = 1, nmol
         init = imol(1,i)
         stop = imol(2,i)
         xmid = 0.0d0
         ymid = 0.0d0
         zmid = 0.0d0
         do j = init, stop
            k = kmol(j)
            weigh = mass(k)
            xmid = xmid + x(k)*weigh
            ymid = ymid + y(k)*weigh
            zmid = zmid + z(k)*weigh
         end do
         weigh = molmass(i)
         xmid = xmid / weigh
         ymid = ymid / weigh
         zmid = zmid / weigh
c
c     save atomic coordinates relative to center of mass
c
         do j = init, stop
            k = kmol(j)
            xoff(k) = x(k) - xmid
            yoff(k) = y(k) - ymid
            zoff(k) = z(k) - zmid
         end do
c
c     convert fractional center of mass to Cartesian coordinates
c
         xmid = xfrac(i)*xbox + yfrac(i)*ybox*gamma_cos
     &             + zfrac(i)*zbox*beta_cos
         ymid = yfrac(i)*ybox*gamma_sin + zfrac(i)*zbox*beta_term
         zmid = zfrac(i)*zbox*gamma_term
c
c     translate coordinates via offset from center of mass
c
         do j = init, stop
            k = kmol(j)
            x(k) = xoff(k) + xmid
            y(k) = yoff(k) + ymid
            z(k) = zoff(k) + zmid
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xoff)
      deallocate (yoff)
      deallocate (zoff)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine xtalwrt  --  output optimization parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "xtalwrt" prints intermediate results during fitting of
c     force field parameters to structures and energies
c
c
      subroutine xtalwrt (niter,nresid,xx,gs,resid)
      use iounit
      use vdwpot
      use xtals
      implicit none
      integer i,niter
      integer nresid
      real*8 xx(*)
      real*8 gs(*)
      real*8 resid(*)
      character*5 vindex
c
c
c     print the values of parameters and scaled derivatives
c
      vindex = 'Class'
      if (vdwindex .eq. 'TYPE ')  vindex = 'Type '
      write (iout,10)  niter
   10 format (/,' Parameters and Scaled Derivatives at',
     &          ' Iteration',i4,' :',/)
      do i = 1, nvary
         write (iout,20)  i,vartyp(i),vindex,vary(1,i),xx(i),gs(i)
 20      format (3x,'(',i2,')',2x,a16,4x,'Atom ',a5,i5,2x,2f14.4)
c         if (ivary(i) .le. 3) then
c            write (iout,20)  i,vartyp(i),vindex,vary(1,i),xx(i),gs(i)
c   20       format (3x,'(',i2,')',2x,a16,4x,'Atom ',a5,i5,2x,2f14.4)
c         else if (ivary(i) .ne. 6) then
c            write (iout,30)  i,vartyp(i),vary(1,i),xx(i),gs(i)
c   30       format (3x,'(',i2,')',2x,a16,4x,'Atom Type ',i5,2x,2f14.4)
c         else
c            write (iout,40)  i,vartyp(i),vary(1,i),vary(2,i),xx(i),gs(i)
c   40       format (3x,'(',i2,')',2x,a16,4x,'Bond Type ',2i5,
c     &                 f11.4,f14.4)
c         end if
      end do
c
c     print the values of the individual residual functions
c
      write (iout,50)  niter
   50 format (/,' Residual Error Function Values at Iteration',
     &           i4,' :',/)
      do i = 1, nresid
         if (i .lt. 100) then
            write (iout,60)  i,rsdtyp(i),iresid(i),resid(i)
   60       format (3x,'(',i2,')',2x,a16,6x,'Structure',i4,4x,f12.4)
         else
            write (iout,70)  i,rsdtyp(i),iresid(i),resid(i)
   70       format (2x,'(',i3,')',2x,a16,6x,'Structure',i4,4x,f12.4)
         end if
      end do
      write (iout,80)
   80 format ()
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  function exchmin1  --  error function for fitpen     ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "exchmin1" computes an error function value derived from
c     derivatives on each atom as well as the intermolecular
c     energy.  
c
c
      function exchmin1 (xx,g)
      use sizes
      use atoms
      use boxes
      use bound
      use charge
      use dipole
      use energi
      use limits
      use math
      use molcul
      use mpole
      use polar
      use vdw
      use xtals
      use chgpen
c
c      use deriv
c      use energi
c      use fitprm
c      use inter
c      use molcul
c      use mplpot
c      use polar
c      use chgpen
      implicit none
      integer i,j,ixtal
      real*8 exchmin1
      real*8 energy,gi,eps,xx0
      real*8 ediff
      real*8 value,value0
      real*8 e,e0
      real*8 e_monomer
      real*8 e_exchange,e0_exchange_log
      real*8 e_mpole
      real*8 e_disp
      real*8 e_vdw
      real*8 resid_ex
      real*8 resid_mpole
      real*8 resid_disp
      real*8 resid_vdw
      real*8 resid_tot
      real*8 xx(*)
      real*8 g(*)
c
c     zero out error function
c
      exchmin1 = 0.0d0
      do i = 1, nvary
         g(i) = 0.0d0
c         print *,"xx",xx(i)
      end do
      do ixtal = 1, nxtal
         call xtalprm ('RESET',ixtal,xx)
         call analysis (energy)
c         print *,"check",energy
         e_exchange = ex
         e_mpole = em
         e_disp = edis
         e_vdw = e_exchange + e_disp
c
         if (exponential .eq. "YES") then
            e_exchange = log(ex)
            if (exmodel .eq. "14-ONLY") e_exchange = log(ev)
            e0_exchange_log = log(e0_exchange)
         end if
c
         if (e0_exchange .ne. 0.0d0) then
            if (exponential .eq. "YES") then
               resid_ex = (e_exchange - e0_exchange_log)**2
               exchmin1 = exchmin1 + resid_ex
            else
               resid_ex = (e_exchange - e0_exchange)**2
               exchmin1 = exchmin1 + resid_ex
            end if
         end if
         if (e0_mpole .ne. 0.0d0) then
            resid_mpole = (e_mpole - e0_mpole)**2
            exchmin1 = exchmin1 + resid_mpole
         end if
         if (e0_disp .ne. 0.0d0) then
            resid_disp = (e_disp - e0_disp)**2
            exchmin1 = exchmin1 + resid_disp
         end if
         if (e0_vdw .ne. 0.0d0) then
            resid_vdw = (e_vdw - e0_vdw)**2
            exchmin1 = exchmin1 + resid_vdw
         end if
c         print *,"residual",ixtal,exchmin1
         do i = 1, nvary
            eps = 0.000001d0
            xx0 = xx(i)
            xx(i) = xx(i) - 0.5d0*eps
            call xtalprm ('RESET',ixtal,xx)
            call analysis (energy)
            e_exchange = ex
            e_mpole = em
            e_disp = edis
            e_vdw = e_exchange + e_disp
c            print *,"check 1",energy
c
            value0 = 0.0d0
            if (e0_exchange .ne. 0.0d0) then
               if (exponential .eq. "YES") then
                  resid_ex = (e_exchange - e0_exchange_log)**2
                  value0 = value0 + resid_ex
               else
                  resid_ex = (e_exchange - e0_exchange)**2
                  value0 = value0 + resid_ex
               end if
            end if
            if (e0_mpole .ne. 0.0d0) then
               resid_mpole = (e_mpole - e0_mpole)**2
               value0 = value0 + resid_mpole
            end if
            if (e0_disp .ne. 0.0d0) then
               resid_disp = (e_disp - e0_disp)**2
               value0 = value0 + resid_disp
            end if
            if (e0_vdw .ne. 0.0d0) then
               resid_vdw = (e_vdw - e0_vdw)**2
               value0 = value0 + resid_vdw
            end if
c
            xx(i) = xx(i) + eps
            call xtalprm ('RESET',ixtal,xx)
            call analysis (energy)
            e_exchange = ex
            e_mpole = em
            e_disp = edis
            e_vdw = e_exchange + e_disp
c            print *,"check 2",energy
c     
            value = 0.0d0
            if (e0_exchange .ne. 0.0d0) then
               if (exponential .eq. "YES") then
                  resid_ex = (e_exchange - e0_exchange_log)**2
                  value = value + resid_ex
               else
                  resid_ex = (e_exchange - e0_exchange)**2
                  value = value + resid_ex
               end if
            end if
            if (e0_mpole .ne. 0.0d0) then
               resid_mpole = (e_mpole - e0_mpole)**2
               value = value + resid_mpole
            end if
            if (e0_disp .ne. 0.0d0) then
               resid_disp = (e_disp - e0_disp)**2
               value = value + resid_disp
            end if
            if (e0_vdw .ne. 0.0d0) then
               resid_vdw = (e_vdw - e0_vdw)**2
               value = value + resid_vdw
            end if
c
            gi = (value - value0) / eps
            g(i) = g(i) + gi
            xx(i) = xx0
         end do
      end do
      return
      end
c
c     ###############################################################
c     ##                                                           ## 
c     ##  subroutine paramwrt1  --  output optimization parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c     "paramwrt1" is a utility that prints intermediate results     
c     during fitting of force field parameters
c
      subroutine paramwrt1 (niter,f,xx)
      use xtals
      use iounit
      implicit none
      integer i,niter
      real*8 f
      real*8 xx(*)
c 
c
c     write the values of parameters and scaled derivatives
c
      write (iout,10)  niter
 10       format (/,' Parameter Values at Iteration',i6,' :',/)
      do i = 1, nvary
c         write (iout,20)  i,vartyp(i),vary(1,i),abs(xx(i))
c 20          format (3x,'(',i2,')',2x,a16,7x,'Atom Class',i5,4x,f12.4)
         write (iout,20)  vartyp(i),vary(1,i),abs(xx(i))
 20            format (3x,a16,7x,i5,4x,f12.4)
      end do
      write (iout,30)
 30       format ()
      return
      end
