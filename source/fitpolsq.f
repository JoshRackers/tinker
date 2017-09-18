c
c
c     ####################################################################
c     ##  COPYRIGHT (C)  2014  by  Jay William Ponder and Josh Rackers  ##
c     ##                      All Rights Reserved                       ##
c     ####################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program fitpolsq  --  fit parameters to dimer structures  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fitpol" computes an optimized set of potential energy
c     parameters for user specified polarization energy
c     interactions by fitting to dimer structures and energies.
c
      program fitpolsq
      use sizes
      use atoms
      use files
      use fitprm
      use iounit
      implicit none
      integer i,j,idim,imon
      integer atom1,atom2
      integer nresid,prmtyp
      real*8 grdmin
      real*8, allocatable :: xx(:)
      real*8, allocatable :: f(:)
      real*8, allocatable :: g(:)
      real*8, allocatable :: xlo(:)
      real*8, allocatable :: xhi(:)
      real*8, allocatable :: jacobian(:,:)
      real*8 e,derivs
      real*8 minimiz1
      real*8 e_mon
      logical exist,query
      character*16 blank
      character*16 label(10)
      character*120 record
      character*120 string
      external dimerr,dimwrt
      external minimiz1
      external optsave
c
c
c     initialize some variables to be used during fitting
c
      call initial
      nvary = 0
      nresid = 0
      blank = '                '
      do i = 1, maxlsq
         vartyp(i) = blank
      end do
      do i = 1, maxrsd
         rsdtyp(i) = blank
      end do
c
c     print informational header about available parameters
c
      write (iout,10)
   10 format (/,' The Following Parameters can be Fit for',
     &           ' each Atom Type :',
     &        //,4x,'(1) Atomic Polarizability',
     &        /,4x,'(2) Thole Damping Factor',
     &        /,4x,'(3) Gordon Induce Dipole Damping',
     &        /,4x,'(4) Gordon Perm Elec Damping',
     &        /,4x,'(6) Regularization')
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
            read (record,*,err=41,end=41)  prmtyp,atom1,atom2
   41         continue
         end if
         if (prmtyp .eq. 0) then
            query = .false.
         else
            query = .true.
            nvary = nvary + 1
            ivary(nvary) = prmtyp
            vary(1,nvary) = atom1
         end if
      end do
c
c     choose scale factors
c
      factor = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=47,end=47)  factor
 47      continue
      if (factor .le. 0.0d0) then
         write (iout,48)
 48            format (/,' Enter Scale Factor for Total Intermolecular',
     &        ' Energy :  ',$)
         read (input,49)  factor
 49            format (f20.0)
      end if
      if (factor .le. 0) factor = 0.0d0
      factor_force = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=51,end=51)  factor_force
 51         continue
      if (factor_force .le. 0.0d0) then
         write (iout,52)
 52         format (/,' Enter Scale Factor for Forces on Atoms',
     &              ' :  ',$)
         read (input,53)  factor_force
 53         format (f20.0)
      end if
      if (factor_force .le. 0) factor_force = 0.0d0
      factor_pol = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=54,end=54)  factor_pol
 54         continue
      if (factor_pol .le. 0.0d0) then
         write (iout,55)
 55                  format (/,' Enter Scale Factor for Intermolecular',
     &              ' Polarization Energy :  ',$)
         read (input,56)  factor_pol
 56      format (f20.0)
      end if
      if (factor_pol .le. 0) factor_pol = 0.0d0
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
c     get number of dimer structures to use in optimization
c
      ndim = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  ndim
   90 continue
      if (ndim .le. 0) then
         write (iout,100)
  100    format (/,' Enter Number of Structures to be Used [1] :  ',$)
         read (input,110)  ndim
  110    format (i10)
      end if
      if (ndim .eq. 0)  ndim = 1
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(nvary))
c
c     get the structural data for each dimer in turn
c
      do idim = 1, ndim
         call getxyz
c
c     get an ideal value for the intermolecular energy
c
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=120,end=120)  e0_dimer
            query = .false.
         end if
  120    continue
         if (query) then
            write (iout,130)
  130       format (/,' Enter Intermolecular Energy Value'
     &           '[<CR> to omit] :  ',$)
            read (input,140)  e0_dimer
  140       format (f20.0)
         end if
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=150,end=150)  e0_pol
            query = .false.
         end if
 150     continue
         if (query) then
            write (iout,160)
 160        format (/,' Enter Intermolecular Polarization Energy Value'
     &           ' [<CR> to omit] :  ',$)
            read (input, 170) e0_pol
 170        format (f20.0)
         end if
c
c     print number of structures and intermolecular energy
c
         write (iout,180)  idim,filename(1:35)
  180    format (/,' File Name of Dimer Structure',i6,' :',7x,a)   
c
c     set the types of residuals for use in optimization
c
c         do i = 1, n
c            nresid = nresid + 1
c            iresid(nresid) = idim
c            rsdtyp(nresid) = 'Force on atom i'
c         end do
c         if (e0_dimer .ne. 0.0d0) then
c            nresid = nresid + 1
c            iresid(nresid) = idim
c            rsdtyp(nresid) = 'Einter '
c            write (iout,190)  e0_dimer
c  190       format (' Value of Intermolecular Energy :  ',f13.2)
c         end if
         if (e0_pol .ne. 0.0d0) then
            nresid = nresid + 1
            iresid(nresid) = idim
            rsdtyp(nresid) = 'Einter (pol) '
            write (iout,191)  e0_pol
 191        format (' Value of Intermolecular Polarization Energy'
     &           ':  ',f13.2)
         end if
c
c     set the initial values of the parameters
c
         call setstruc ('STORE',idim)
      end do
c
c     types of variables for use in optimization
c
      call setparam ('STORE',xx)
      label(1) = 'Polarizability     '
      label(2) = 'Thole Damping      '
      label(3) = 'Gordon Damping     '
      label(4) = 'Gordon Perm Damping'
      label(5) = 'Piquemal Factor '
      label(6) = 'Regularization     '
      label(7) = 'Truhlar Zeta    '
      do i = 1, nvary
         vartyp(i) = label(ivary(i))
      end do
c
c     print the initial parameter values
c
      write (iout,210)
  210 format (/,' Initial Values of the Parameters :',/)
      do i = 1, nvary
         write (iout,220)  i,vartyp(i),vary(1,i),xx(i)
  220    format (3x,'(',i4,')',2x,a16,4x,'Atom Type',i5,4x,f12.4)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (f(nresid))
      allocate (g(nvary))
      allocate (xlo(nvary))
      allocate (xhi(nvary))
      allocate (jacobian(nresid,nvary))
c
c     set upper and lower bounds based on the parameter type
c
      do i = 1, nvary
         if (ivary(i) .eq. 1) then
            xlo(i) = 0.25d0 * xx(i)
            xhi(i) = 4.0d0 * xx(i)
         else if (ivary(i) .eq. 2) then
            xlo(i) = 0.25d0 * xx(i)
            xhi(i) = 4.0d0 * xx(i)
         else if (ivary(i) .eq. 3) then
            xlo(i) = 0.25d0 * xx(i)
            xhi(i) = 4.0d0 * xx(i)
         else if (ivary(i) .eq. 4) then
            xlo(i) = 0.25d0 * xx(i)
            xhi(i) = 4.0d0 * xx(i)
         else if (ivary(i) .eq. 5) then
            xlo(i) = 0.25d0 * xx(i)
            xhi(i) = 4.0d0 * xx(i)
         else if (ivary(i) .eq. 6) then
            xlo(i) = 0.25d0 * xx(i)
            xhi(i) = 4.0d0 * xx(i)
         else if (ivary(i) .eq. 7) then
            xlo(i) = 0.25d0 * xx(i)
            xhi(i) = 4.0d0 * xx(i)
         end if
      end do
c
c     use nonlinear least squares to refine the parameters
c
      call square (nvary,nresid,xlo,xhi,xx,f,g,jacobian,
     &                  grdmin,dimerr,dimwrt)
c      call square (nresid,nvary,xlo,xhi,xx,f,g,jacobian,
c     &                maxrsd,grdmin,dimerr,dimwrt)
c
c     perform deallocation of some local arrays
c
      deallocate (xlo)
      deallocate (xhi)
      deallocate (jacobian)
c
c     print the final parameter values
c
      write (iout,250)
  250 format (/,' Final Values of the Parameters and Scaled',
     &            ' Derivatives :',/)
      do i = 1, nvary
         write (iout,260)  i,vartyp(i),vary(1,i),xx(i),g(i)
  260    format (3x,'(',i3,')',2x,a16,4x,'Atom Type',i5,4x,2f12.4)
      end do
c
c     write out residual function values
c
      write (iout,270)  
 270  format (/,' Final Error Function Values ')
      do i = 1, nresid
         write (iout,280)  i,rsdtyp(i),iresid(i),f(i)
 280     format (3x,'(',i2,')',2x,a16,4x,2x,'Dimer',i4,4x,f12.4)
      end do
      write (iout,290)
 290  format ()
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      deallocate (f)
      deallocate (g)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine setstruc  --  manage optimization structures  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "setstruc" stores or retrieves a reference structure and
c     its associated energy values
c
c
      subroutine setstruc (mode,idim)
      use fitprm
      use chgpen
      implicit none
      integer maxstruc
      parameter (maxstruc=1000)
      integer idim
      real*8 e0_dimers(maxstruc)
      real*8 e0_pols(maxstruc)
      character*5 mode
      save e0_dimers
      save e0_pols
c
c
c     store or reset the atomic coordinates and other info
c
      if (mode .eq. 'STORE') then
         call makeref (idim)
      else if (mode .eq. 'RESET') then
         call getref (idim)
      end if
c
c     store or reset the ideal intermolecular energy value
c
      if (mode .eq. 'STORE') then
         e0_dimers(idim) = e0_dimer
         e0_pols(idim) = e0_pol
      else if (mode .eq. 'RESET') then
         e0_dimer = e0_dimers(idim)
         e0_pol = e0_pols(idim)
      end if
c
c     setup mechanics for structures with different constitutions
c
      call mechanic
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine setparam  --  manage optimization parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "setparam" stores or resets the force field parameters
c     being optimized
c
c
      subroutine setparam (mode,xx)
      use sizes
      use atomid
      use atoms
      use fitprm
      use kvdws
      use vdw
      use kpolr
      use polar
      use chgpen
      implicit none
      integer i,j
      integer prmtyp
      integer jclass
      real*8 sepsi,sepsj
      real*8 xx(*)
      character*5 mode
c
c
c     store or reset the force field parameters to be optimized
c
c     note jclass refers to either class or type depending on parameter
c
      do j = 1, nvary
         prmtyp = ivary(j)
         jclass = vary(1,j)
c     atomic polarizability
         if (prmtyp .eq. 1) then
            if (mode .eq. 'STORE') then
               xx(j) = xpolr(jclass)
            else if (mode .eq. 'RESET') then
               xpolr(jclass) = abs(xx(j))
c               do i = 1, n
c                  if (cpclass(i).eq.jclass) then
c                     polarity(i) = xpolr(jclass)
c                     alphap(jclass) = polfactor(jclass) *
c     &                 polarity(i)**(1.0d0/3.0d0)
c                  end if
c               end do
            end if
c     thole damping factor
         else if (prmtyp .eq. 2) then
            if (mode .eq. 'STORE') then
               xx(j) = xathl(jclass)
            else if (mode .eq. 'RESET') then
               xathl(jclass) = abs(xx(j))
c               do i = 1, n
c                  if (cpclass(i) .eq. jclass)  thole(i) = xathl(jclass)
c               end do
            end if
c     gordon damping factor - polarizability
         else if (prmtyp .eq. 3) then
            if (mode .eq. 'STORE') then
               xx(j) = alphap(jclass)
            else if (mode .eq. 'RESET') then
               alphap(jclass) = abs(xx(j))
            end if
c     gordon damping factor - permanent electrostatics
         else if (prmtyp .eq. 4) then
            if (mode .eq. 'STORE') then
               xx(j) = alpha(jclass)
            else if (mode .eq. 'RESET') then
               alpha(jclass) = abs(xx(j))
            end if
c     gordon damping proportional to polarizability
         else if (prmtyp .eq. 5) then
            if (mode .eq. 'STORE') then
               xx(j) = polfactor(jclass)
            else if (mode .eq. 'RESET') then
               polfactor(jclass) = abs(xx(j))
c               do i = 1, n
c                  if (cpclass(i) .eq. jclass) then
c                     alphap(jclass) = polfactor(jclass) * 
c     &                 polarity(i)**(1.0d0/3.0d0)
c                  end if
c               end do
            end if
         else if (prmtyp .eq. 6) then
            if (mode .eq. 'STORE') then
               xx(j) = regular(jclass)
            else if (mode .eq. 'RESET') then
               regular(jclass) = abs(xx(j))
            end if
         end if      
      end do
      do j = 1, nvary
         prmtyp = ivary(j)
         jclass = vary(1,j)
c     atomic polarizability                                                
         if (prmtyp .eq. 1) then
            do i = 1, n
               if (cpclass(i).eq.jclass) then
                  polarity(i) = xpolr(jclass)
c     comment out if you want alphap fixed
c                  alphap(jclass) = polfactor(jclass) *
c     &                 polarity(i)**(1.0d0/3.0d0)
               end if
            end do
c     thole damping
         else if (prmtyp .eq. 2) then
            do i = 1, n
               if (cpclass(i) .eq. jclass)  thole(i) = xathl(jclass)
            end do
c     gordon damping proportional to polarizability
         else if (prmtyp .eq. 5) then
            do i = 1, n
               if (cpclass(i).eq.jclass) then
                  alphap(jclass) = polfactor(jclass) *
     &                 polarity(i)**(1.0d0/3.0d0)
               end if
            end do
         end if   
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine dimerr  --  error function for fitpen    ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "dimerr" computes an error function value derived from
c     derivatives on each atom as well as the intermolecular
c     energy.  
c
c
      subroutine dimerr (nvaried,nresid,xx,resid)
      use sizes
      use atoms
      use deriv
      use energi
      use fitprm
      use inter
      use molcul
      use mplpot
      use polar
      use chgpen
      implicit none
      integer i,j,i1,i2,idim
      integer nresid,nvaried
      real*8 minimiz1
      real*8 energy,e,e1
      real*8 xforce,yforce,zforce
      real*8 force2,forcesum
      real*8 e_dimer,e_inter
      real*8 xx(*)
      real*8 resid(*)
      real*8, allocatable :: derivs(:,:)
c
c     zero out the number of residual functions
c
      nresid = 0
c      forcesum = 0.0d0
c
c     set the values of the potential energy parameters
c     and get the energy of the base structure
c
      do idim = 1, ndim
         call setstruc ('RESET',idim)
         call setparam ('RESET',xx)
         call analysis (energy)
c
c     perform dynamic allocation of some local arrays
c
c         allocate (derivs(3,n))
c
c     get energy derivatives for all atoms
c
c         call gradient (e,derivs)
c         do j = 1, n
c            nresid = nresid + 1
c            xforce = -derivs(1,j)
c            yforce = -derivs(2,j)
c            zforce = -derivs(3,j)
c            force2 = xforce*xforce + yforce*yforce + 
c     &               zforce*zforce
c            resid(nresid) = factor_force * sqrt(force2) / dble(n)
c         end do
c
c     compute the intermolecular energy
c     note: getting the monomer energy seperately is possible when
c           different size xyz files can be swapped in and out
c
c         e_inter = einter
c
c     compute the intermolecular electrostatic energy
c     note: getting the monomer electrostatic energy interactively
c           is also possible when different size xyz files
c           can be swapped in and out
c
c         e_elec = em + ex - (7.0282*2)
c         e_elec = em + ex
c
c
c     compute the residual vector as a collection
c     of the forces on all the atoms as well as
c     the intermolecular energy
c
c         if (e0_dimer .ne. 0.0d0) then
c            nresid = nresid + 1
c            resid(nresid) = factor * (e_inter - e0_dimer)
c         end if
         if (e0_pol .ne. 0.0d0) then
            nresid = nresid + 1
            resid(nresid) = factor_pol * (ep - e0_pol)
         end if
c         print *,idim,ep,e0_pol,resid(nresid)
c
c     perform deallocation of some local arrays
c
c         deallocate (derivs)
      end do
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine dimwrt  --  write current parameters   ##
c     ##                                                    ##
c     ########################################################
c
c
c     "dimwrt" is a utility that prints intermediate results
c     during fitting of force field parameters
c
c
      subroutine dimwrt (niter,nresid,xx,gs,f)
      use iounit
      use fitprm
      implicit none
      integer i,niter
      integer nresid
      real*8 xx(*)
      real*8 gs(*)
      real*8 f(*)
c
c     write the values of parameters and scaled derivatives
c
      write (iout,10)  niter
   10 format (/,' Parameters and Scaled Derivatives at',
     &          ' Iteration',i4,' :',/)
      do i = 1, nvary
         write (iout,20)  i,vartyp(i),vary(1,i),xx(i),gs(i)
   20    format (3x,'(',i2,')',2x,a16,4x,'Atom Class',i5,4x,2f12.4)
      end do
c
c     write the values of the residual functions
c
      write (iout,50)  niter
   50 format (/,' Residual Error Function Values at Iteration',
     &           i4,' :',/)
      do i = 1, nresid
         write (iout,60)  i,rsdtyp(i),iresid(i),f(i)
   60    format (3x,'(',i3,')',2x,a16,4x,2x,'Dimer',i4,4x,f12.4)
      end do
      write (iout,70)
   70 format ()
      return
      end
