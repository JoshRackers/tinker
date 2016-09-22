c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole3  --  mpole/polar energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole3" calculates the electrostatic energy due to
c     atomic multipole and dipole polarizability interactions,
c     and partitions the energy among the atoms
c
c
      subroutine empole3
      use sizes
      use analyz
      use chgpen
      use energi
      use limits
      use mpole
      use potent
      use potderivs
      implicit none
      integer i,ii
c
c
c     choose polarization damping if necessary
c     NEED TO RESET IN EPOLAR ROUTINES IF USING DIFFERENT TYPES OF DAMPING
c
      if (use_polar) then
         if (directdamp .eq. "GORDON") then
            damp_gordon = .true.
            if (regularize .eq. "YES") damp_gordonreg = .true.
         else if (directdamp .eq. "PIQUEMAL") then
            damp_piquemal = .true.
         else if (directdamp .eq. "THOLE") then
            damp_thole = .true.
         end if
      end if
      if (penetration .eq. "GORDON") damp_gordon = .true.
      if (penetration .eq. "PIQUEMAL") damp_piquemal = .true.
c
c     choose the method for summing over multipole interactions 
c
      if (use_ewald) then
         damp_none = .true.
         damp_ewald = .true.
         call empole3b
      else
         damp_none = .true.
         damp_ewald = .false.
         call empole3a
      end if
c
c
c     zero out energy terms and analysis which are not in use
c
      if (.not. use_mpole) then
         em = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            aem(ii) = 0.0d0
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole3a  --  nonewald multipole analysis     ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole3a" calculates the atomic multipole and dipole           
c     polarizability interaction energy and partitions the 
c     energy among the atoms
c
c
      subroutine empole3a
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use cell
      use chgpot
      use chgpen
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use limits
      use math
      use molcul
      use mplpot      
      use mpole
      use polar
      use polgrp
      use polpot
      use potderivs
      use potent
      use shunt
      use usage
      implicit none
      integer i
      integer ii
      integer ix,iy,iz
      real*8 e,fgrp
      real*8 eold
      real*8 f,fm,fp
      real*8 ci,zi,qi
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      logical proceed
      logical header,huge
      logical usei,usek
      logical muse
      character*6 mode
c
c
c     zero out multipole and polarization energy and partitioning
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
      header = .true.
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the permanent electric potential,
c     field and field gradient at each multipole site
c
      call permfield2
c
c     set conversion factor
c
      f = electric / dielec
c
c     calculate the total multipole interaction energy
c
      do i = 1, npole
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         qixy = 2.0d0*qixy
         qixz = 2.0d0*qixz
         qiyz = 2.0d0*qiyz
c
c     split nuclear and electronic charge if applying charge penetration
c
         if (penetration .eq. "GORDON") then
            zi = atomic(i)
            if (num_ele .eq. "VALENCE") then
               if (atomic(i) .gt. 2)  zi = zi - 2.0d0
               if (atomic(i) .gt. 10)  zi = zi - 8.0d0
               if (atomic(i) .gt. 18)  zi = zi - 8.0d0
               if (atomic(i) .gt. 20)  zi = zi - 10.0d0
            end if
            qi = ci - zi
         end if
c     
c     contract multipoles with permanent potential, field and field gradient
c
c     charge
         eold = 0.5d0*ci*potm(ii)
c     permanent dipole
         eold = eold + 0.5d0*
     &        (dix*fieldm(1,ii)+diy*fieldm(2,ii)+diz*fieldm(3,ii))
c     quadrupole
         eold = eold + 0.5d0*
     &        (qixx*gradfieldm(1,1,ii) + qixy*gradfieldm(2,1,ii) +
     &        qixz*gradfieldm(3,1,ii) + qiyy*gradfieldm(2,2,ii) +
     &        qiyz*gradfieldm(3,2,ii) + qizz*gradfieldm(3,3,ii))
c
c     apply f constant
c
         eold = f * eold
c
c     compute charge penetration corrected permanent multipole energy
c
c     charge - nucleus
         e = 0.5d0*zi*nucpotm_gordon(ii)
c     charge - electrons
         e = e + 0.5d0*qi*potm_gordon(ii)
c     permanent dipole
         e = e + 0.5d0*
     &        (dix*fieldm_gordon(1,ii)+diy*fieldm_gordon(2,ii)+
     &        diz*fieldm_gordon(3,ii))
c     quadrupole
         e = e + 0.5d0*
     &        (qixx*gradfieldm_gordon(1,1,ii) + 
     &        qixy*gradfieldm_gordon(2,1,ii) +
     &        qixz*gradfieldm_gordon(3,1,ii) + 
     &        qiyy*gradfieldm_gordon(2,2,ii) +
     &        qiyz*gradfieldm_gordon(3,2,ii) + 
     &        qizz*gradfieldm_gordon(3,3,ii))
c
c     apply f constant
c
         e = f * e
c
c     increment the overall multipole and polarization energies
c
         em = em + e
         nem = nem + 1
         aem(ii) = aem(ii) + e
      end do
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine empole3b  --  ewald summation multipole analysis  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "empole3b" calculates the atomic multipole and dipole 
c     polarizability interaction energy using a particle mesh
c     Ewald summation
c
c
      subroutine empole3b
      use sizes
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use inter
      use limits
      use math
      use mpole
      use polar
      use potderivs
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
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the permanent electric potential,
c     field and field gradient at each multipole site
c
      call permfield2
c
c     get reciprocal space potential, field and field gradient
c
      call permrecip1
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
c
c     compute the real space, reciprocal space and self-energy 
c     parts of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
c
c     terms needed for self energy
c
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
c
c     terms for real space energy
c 
         qixy = 2.0d0*qixy
         qixz = 2.0d0*qixz
         qiyz = 2.0d0*qiyz
c
c     real space
c
c     charge
         ereal = 0.5d0*ci*pot_ewald(ii)
c     permanent dipole
         ereal = ereal + 0.5d0*(dix*field_ewald(1,ii) + 
     &        diy*field_ewald(2,ii) + diz*field_ewald(3,ii))
c     quadrupole
         ereal = ereal + 0.5d0*(
     &    qixx*gradfield_ewald(1,1,ii)+ qixy*gradfield_ewald(2,1,ii) +
     &    qixz*gradfield_ewald(3,1,ii)+ qiyy*gradfield_ewald(2,2,ii) +
     &    qiyz*gradfield_ewald(3,2,ii)+ qizz*gradfield_ewald(3,3,ii))
c
c     self-energy
c
         eself = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
c
c     full real space energies needed for scaled interactions
c
         efix = 0.5d0*ci*pot(ii)
         efix = efix + 0.5d0*(dix*field(1,ii) +
     &        diy*field(2,ii) + diz*field(3,ii))
         efix = efix + 0.5d0*(
     &    qixx*gradfield(1,1,ii)+ qixy*gradfield(2,1,ii) +
     &    qixz*gradfield(3,1,ii)+ qiyy*gradfield(2,2,ii) +
     &    qiyz*gradfield(3,2,ii)+ qizz*gradfield(3,3,ii))
c
c     scaled interactions
c
         efix = efix - 0.5d0*ci*potm(ii)
         efix = efix - 0.5d0*(dix*fieldm(1,ii) +
     &        diy*fieldm(2,ii) + diz*fieldm(3,ii))
         efix = efix - 0.5d0*(
     &    qixx*gradfieldm(1,1,ii)+ qixy*gradfieldm(2,1,ii) +
     &    qixz*gradfieldm(3,1,ii)+ qiyy*gradfieldm(2,2,ii) +
     &    qiyz*gradfieldm(3,2,ii)+ qizz*gradfieldm(3,3,ii))
c
c     reciprocal space
c
c     charge
         erecip = 0.5d0*ci*pot_recip(i) 
c     permanent dipole
         erecip = erecip + 0.5d0*(dix*field_recip(1,ii) + 
     &        diy*field_recip(2,ii)+diz*field_recip(3,ii))
c     quadrupole
         erecip = erecip + 0.5d0*(
     &     qixx*gradfield_recip(1,1,ii) + qixy*gradfield_recip(2,1,ii) +
     &     qixz*gradfield_recip(3,1,ii) + qiyy*gradfield_recip(2,2,ii) +
     &     qiyz*gradfield_recip(3,2,ii) + qizz*gradfield_recip(3,3,ii))
c
c     apply f constant
c
         ereal = f * ereal
         efix = f * efix
         erecip = f * erecip
c
c     increment the overall multipole energy
c
         em = em + ereal + eself - efix + erecip
         nem = nem + 1
         aem(i) = aem(i) + e
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
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
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
         nem = nem + 1
         em = em + term*(xd*xd+yd*yd+zd*zd)
      end if
      return
      end

