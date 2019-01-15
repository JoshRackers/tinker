c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module mpole  --  atomic multipoles in current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxpole   max components (monopole=1,dipole=4,quadrupole=13)
c
c     npole     total number of multipole sites in the system
c     ipole     number of the atom for each multipole site
c     polsiz    number of multipole components at each atom
c     pollist   multipole site for each atom (0=no multipole)
c     zaxis     number of the z-axis defining atom for each site
c     xaxis     number of the x-axis defining atom for each site
c     yaxis     number of the y-axis defining atom for each site
c     pole      traceless Cartesian multipoles in the local frame
c     rpole     traceless Cartesian multipoles in the global frame
c     spole     spherical harmonic multipoles in the local frame
c     srpole    spherical harmonic multipoles in the global frame
c     polaxe    local axis type for each multipole site
c
c
      module mpole
      implicit none
      integer maxpole
      parameter (maxpole=13)
      integer npole
      integer, allocatable :: ipole(:)
      integer, allocatable :: polsiz(:)
      integer, allocatable :: pollist(:)
      integer, allocatable :: zaxis(:)
      integer, allocatable :: xaxis(:)
      integer, allocatable :: yaxis(:)
      real*8, allocatable :: pole(:,:)
      real*8, allocatable :: rpole(:,:)
      real*8, allocatable :: spole(:,:)
      real*8, allocatable :: srpole(:,:)
      character*8, allocatable :: polaxe(:)
      save
      end
