c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine epolar  --  polarization energy              ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "epolar" calculates the electrostatic energy due to
c     dipole polarizability interactions
c
c
      subroutine epolar
      use sizes
      use analyz
      use energi
      use limits
      use mpole
      use potent
      implicit none
      integer i,ii
c
c
c     choose the method for summing over multipole interactions
c
      if (use_ewald) then
         call epolar0b
      else
         call epolar0a
      end if
c
c
c     zero out energy terms and analysis which are not in use
c
      if (.not. use_polar) then
         ep = 0.0d0
         do i = 1, npole
            ii = ipole(i)
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine epolar0a  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "epolar0a" calculates the total dipole           
c     polarizability interaction energy and partitions the 
c     energy among the atoms
c
c
      subroutine epolar0a
      use sizes
      use action
      use analyz
c      use atomid
      use atoms
      use bound
      use boxes
      use cell
      use chgpen
      use chgpot
c      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use limits
c      use math
      use molcul
      use mplpot      
      use mpole
      use polar
c      use polgrp
c      use polpot
      use potderivs
      use potent
c      use shunt
      use usage
      implicit none
      integer i
      integer ii
      integer ix,iy,iz
      real*8 e,ei,fgrp
      real*8 f
      real*8 uix,uiy,uiz
      real*8, allocatable :: fieldd_damp(:,:)
      real*8, allocatable :: fieldp_damp(:,:)
      logical header
c
c
c     zero out the polarization energy and partitioning
c
      ep = 0.0d0
      header = .true.
      if (npole .eq. 0)  return
c
c     compute the induced dipoles at each polarizable atom
c
      damp_thole = .false.
      damp_gordon = .false.
      damp_piquemal = .false.
      call induce
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldd_damp(3,npole))
      allocate (fieldp_damp(3,npole))
c
c     check what kind of damping was used
c
      if (directdamp .eq. "GORDON") then
         fieldp_damp = fieldp_gordon
      else if (directdamp .eq. "PIQUEMAL") then
         fieldp_damp = fieldp_piquemal
      else if (directdamp .eq. "THOLE") then
         fieldp_damp = fieldp_thole
      end if
c
c     set conversion factor
c
      f = electric / dielec
c
c     calculate the total multipole interaction energy
c
      do i = 1, npole
         ii = ipole(i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
c
c     contract induced dipoles with damp damped permanent field
c
         ei = 0.5*(uix*fieldp_damp(1,ii) + uiy*fieldp_damp(2,ii) + 
     &        uiz*fieldp_damp(3,ii))
c
c     apply f constant
c
         ei = f * ei
c
c     increment the overall polarization energy
c
         ep = ep + ei
      end do
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar0b  --  ewald summation multipole analysis  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar0b" calculates the atomic multipole and dipole 
c     polarizability interaction energy using a particle mesh
c     Ewald summation
c
c
      subroutine epolar0b
      use sizes
      use action
      use analyz
      use atoms
      use boxes
      use chgpen
      use chgpot
      use energi
      use ewald
      use inter
      use limits
      use math
      use mpole
      use polar
      use potderivs
      use potent
      implicit none
      integer i,ii
      real*8 e,ei,eintra
      real*8 ereal,eself,efix,erecip
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8, allocatable :: fieldd_damp(:,:)
      real*8, allocatable :: fieldp_damp(:,:)
c
c
c     zero out the polarization energy and partitioning
c
      ep = 0.0d0
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     terms needed for self energy
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
c
c     compute the induced dipoles at each polarizable atom
c
      damp_thole = .false.
      damp_gordon = .false.
      damp_piquemal = .false.
      call induce
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldd_damp(3,npole))
      allocate (fieldp_damp(3,npole))
c
c     check what kind of damping was used
c
      if (directdamp .eq. "GORDON") then
         fieldp_damp = fieldp_gordon
      else if (directdamp .eq. "PIQUEMAL") then
         fieldp_damp = fieldp_piquemal
      else if (directdamp .eq. "THOLE") then
         fieldp_damp = fieldp_thole
      end if
c
c     compute the real space, reciprocal space and self-energy 
c     parts of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
c
c     terms needed for self energy
c
         uii = dix*uix + diy*uiy + diz*uiz
c
c     real space
c
         ereal = 0.5d0*(
     &        uix*field_ewald(1,ii) + uiy*field_ewald(2,ii) + 
     &        uiz*field_ewald(3,ii))
c
c     full real space energies needed for scaled interactions
c
         efix = 0.5d0*(
     &        uix*field(1,ii) + uiy*field(2,ii) +
     &        uiz*field(3,ii))
c
         efix = efix - 0.5d0*(
     &        uix*fieldp_damp(1,ii) + uiy*fieldp_damp(2,ii) +
     &        uiz*fieldp_damp(3,ii))
c
c     self-energy
c
         eself = fterm * term * uii / 3.0d0
c
c     reciprocal space
c
         erecip = 0.5d0*(uix*field_recip(1,ii) +
     &        uiy*field_recip(2,ii)+uiz*field_recip(3,ii))
c
c     apply f constant
c
         ereal = f * ereal
         efix = f * efix
         erecip = f * erecip
c
c     increment the overall polarization energy
c
         ep = ep + ereal + eself - efix + erecip
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            uix = uind(1,i)
            uiy = uind(2,i)
            uiz = uind(3,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
            xu = xu + uix
            yu = yu + uiy
            zu = zu + uiz
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
      end if
      return
      end
