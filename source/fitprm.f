c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2014  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  fitprm.i  --  force field parameters used in fitting   ##
c     ##                                                         ##
c     #############################################################
c
c
c     maxlsq       maximum number of least squares variables
c     maxrsd       maximum number of least squares residuals
c
c     ndim         number of dimer structures to be stored
c     nvary        number of potential parameters to optimize
c     ivary        index for the types of potential parameters
c     vary         atom numbers involved in potential parameters
c     iresid       dimer structure to which each residual refers
c     e0_dimer     ideal intermolecular energy for the current dimer
c     e0_monomer   sum of the energies of the two monomers
c     e0_elec      ideal intermolecular electrostaic energy for current dimer
c     factor       weight of total intermolecular energy in optimization function
c     factor_force weight of forces in optimazation function
c     factor_ele   weight of intermolecular electrostatic energy in optimization
c     rsdtyp       experimental variable for each of the residuals
c     vartyp       type of potential parameter to be optimized
cccccccc
c     e0_pol       ideal polarization energy
c     e0_polab     ideal a <- b polarization energy
c     e0_polba     ideal b <- a polarization energy
c     poltotal     total interactive molecular polarizability
c     polx         x polarizability
c     poly         y polarizability 
c     polz         z polarizability
c     factor_pol
c     factor_total
c     factor_xyz
c
c
      module fitprm
      implicit none
      integer maxlsq
      integer maxrsd
      parameter (maxlsq=1000)
      parameter (maxrsd=100000)
      integer ndim,nvary
      integer ivary(maxlsq)
      integer vary(2,maxlsq)
      integer iresid(maxrsd)
      real*8 e0_dimer
      real*8 e0_monomer
      real*8 e0_elec
      real*8 factor
      real*8 factor_force
      real*8 factor_ele
ccccccc
      real*8 e0_pol
      real*8 e0_polab
      real*8 e0_polba
      real*8 poltotal
      real*8 polx,poly,polz
      real*8 factor_pol
      real*8 factor_total
      real*8 factor_xyz
ccccccc
      character*16 rsdtyp(maxrsd)
      character*16 vartyp(maxlsq)
      save
      end
