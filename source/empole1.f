c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole1  --  mpole/polar energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1" calculates the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
      subroutine empole1
      use sizes
      use chgpen
      use deriv
      use energi
      use limits
      use mpole
      use potent
      use potderivs
      implicit none
      integer i,j,ii
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
         call empole1b
      else
         damp_none = .true.
         damp_ewald = .false.
         call empole1a
      end if
c
c     zero out energy and derivative terms which are not in use
c
      if (.not. use_mpole) then
         em = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            do j = 1, 3
               dem(j,ii) = 0.0d0
            end do
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole1a  --  nonewald multipole derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole1a" calculates the multipole
c     energy and derivatives with respect to Cartesian coordinates
c
c
      subroutine empole1a
      use sizes
      use atomid
      use atoms
      use bound
      use boxes
      use cell
      use chgpot
      use chgpen
      use couple
      use deriv
      use energi
      use group
      use inter
      use limits
      use molcul
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use potderivs
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,ii
      real*8 e,ei,f
      real*8 fx,fy,fz
      real*8 ci,zi,qi
      real*8 dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8, allocatable :: frc(:,:)
      real*8, allocatable :: trq(:,:)
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (frc(3,n))
      allocate (trq(3,npole))
c
c     zero out multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
      end do
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
c     field, field gradient and field hessian 
c     at each multipole site
c
      call permfield3
c
c     set conversion factor
c
      f = electric / dielec
c
c     calculate the total multipole energy and derivatives
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
c     increment the overall multipole energy
c
         em = em + e
c
c     compute charge penetration corrected multipole gradient
c
c     charge - nucleus
         fx = zi*nucfieldm_gordon(1,ii)
         fy = zi*nucfieldm_gordon(2,ii)
         fz = zi*nucfieldm_gordon(3,ii)
c     charge - electrons
         fx = fx + qi*fieldm_gordon(1,ii)
         fy = fy + qi*fieldm_gordon(2,ii)
         fz = fz + qi*fieldm_gordon(3,ii)
c     dipole
         fx = fx + dix*gradfieldm_gordon(1,1,ii) + 
     &        diy*gradfieldm_gordon(2,1,ii) +
     &        diz*gradfieldm_gordon(3,1,ii)
         fy = fy + dix*gradfieldm_gordon(2,1,ii) + 
     &        diy*gradfieldm_gordon(2,2,ii) +
     &        diz*gradfieldm_gordon(3,2,ii)
         fz = fz + dix*gradfieldm_gordon(3,1,ii) + 
     &        diy*gradfieldm_gordon(3,2,ii) +
     &        diz*gradfieldm_gordon(3,3,ii)
c     quadrupole
         fx = fx + qixx*hessfieldm_gordon(1,1,1,ii) +
     &        qixy*hessfieldm_gordon(2,1,1,ii) +
     &        qixz*hessfieldm_gordon(3,1,1,ii) + 
     &        qiyy*hessfieldm_gordon(2,2,1,ii) +
     &        qiyz*hessfieldm_gordon(3,2,1,ii) + 
     &        qizz*hessfieldm_gordon(3,3,1,ii)
         fy = fy + qixx*hessfieldm_gordon(2,1,1,ii) + 
     &        qixy*hessfieldm_gordon(2,2,1,ii) +
     &        qixz*hessfieldm_gordon(3,2,1,ii) + 
     &        qiyy*hessfieldm_gordon(2,2,2,ii) +
     &        qiyz*hessfieldm_gordon(3,2,2,ii) + 
     &        qizz*hessfieldm_gordon(3,3,2,ii)
         fz = fz + qixx*hessfieldm_gordon(3,1,1,ii) +
     &        qixy*hessfieldm_gordon(3,2,1,ii) +
     &        qixz*hessfieldm_gordon(3,3,1,ii) + 
     &        qiyy*hessfieldm_gordon(3,2,2,ii) +
     &        qiyz*hessfieldm_gordon(3,3,2,ii) + 
     &        qizz*hessfieldm_gordon(3,3,3,ii)
c
c     increment the permanent multipole gradient
c
         dem(1,ii) = f * fx
         dem(2,ii) = f * fy
         dem(3,ii) = f * fz
c
c     calculate permanent multipole torques (no torques on charges)
c
c     dipole
         trq(1,ii) = diz*fieldm_gordon(2,ii) - diy*fieldm_gordon(3,ii)
         trq(2,ii) = dix*fieldm_gordon(3,ii) - diz*fieldm_gordon(1,ii)
         trq(3,ii) = diy*fieldm_gordon(1,ii) - dix*fieldm_gordon(2,ii)
c     quadrupole
         trq(1,ii) = trq(1,ii) + 
     &        2.0d0*(qizz - qiyy)*gradfieldm_gordon(3,2,ii)
     &        + qixz*gradfieldm_gordon(2,1,ii) + 
     &        qiyz*gradfieldm_gordon(2,2,ii)
     &        - qixy*gradfieldm_gordon(3,1,ii) - 
     &        qiyz*gradfieldm_gordon(3,3,ii)
         trq(2,ii) = trq(2,ii) + 
     &        2.0d0*(qixx - qizz)*gradfieldm_gordon(3,1,ii)
     &        + qixy*gradfieldm_gordon(3,2,ii) + 
     &        qixz*gradfieldm_gordon(3,3,ii)
     &        - qixz*gradfieldm_gordon(1,1,ii) - 
     &        qiyz*gradfieldm_gordon(2,1,ii)
         trq(3,ii) = trq(3,ii) + 
     &        2.0d0*(qiyy - qixx)*gradfieldm_gordon(2,1,ii)
     &        + qixy*gradfieldm_gordon(1,1,ii) + 
     &        qiyz*gradfieldm_gordon(3,1,ii)
     &        - qixy*gradfieldm_gordon(2,2,ii) - 
     &        qixz*gradfieldm_gordon(3,2,ii)
         trq(1,ii) = f * trq(1,ii)
         trq(2,ii) = f * trq(2,ii)
         trq(3,ii) = f * trq(3,ii)
      end do
c
c     distribute torques into permanent multipole gradient
c
      do i = 1, n
         frc(1,i) = 0.0d0
         frc(2,i) = 0.0d0
         frc(3,i) = 0.0d0
      end do
      call torque2 (trq,frc)
      do i = 1, n
         dem(1,i) = dem(1,i) + frc(1,i)
         dem(2,i) = dem(2,i) + frc(2,i)
         dem(3,i) = dem(3,i) + frc(3,i)
      end do
      return
      end
c
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine empole1b  --  ewald permanent multipole energy and  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "empole1b" calculates the atomic multipole energy and 
c     gradient using a particle mesh Ewald summation
c
c
      subroutine empole1b
      use sizes
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use inter
      use limits
      use math
      use mpole
      use potderivs
      implicit none
      integer i,j,ii
      real*8 e,ei,eintra
      real*8 ereal,eself,efix,erecip
      real*8 fxreal,fxfix,fxrecip
      real*8 fyreal,fyfix,fyrecip
      real*8 fzreal,fzfix,fzrecip
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield
      real*8 ydfield
      real*8 zdfield
      real*8, allocatable :: frc(:,:)
      real*8, allocatable :: trq(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (frc(3,n))
      allocate (trq(3,npole))
c
c     zero out multipole and polarization energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
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
c     field, field gradient and field hessian at each multipole site
c
      call permfield3
c
c     get reciprocal space potential, field and field gradient
c     and field hessian
c
      call permrecip1
c
c     terms for self energy
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
     &        + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
c
c     terms for real space energy
c 
         qixy = 2.0d0*qixy
         qixz = 2.0d0*qixz
         qiyz = 2.0d0*qiyz
c
c     real space energy
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
c
c     calculate the permanent multipole forces
c
c
c     real space forces
c
c     charge
         fxreal = ci*field_ewald(1,ii)
         fyreal = ci*field_ewald(2,ii)
         fzreal = ci*field_ewald(3,ii)
c     dipole
         fxreal = fxreal + dix*gradfield_ewald(1,1,ii) + 
     &        diy*gradfield_ewald(2,1,ii) + diz*gradfield_ewald(3,1,ii)
         fyreal = fyreal + dix*gradfield_ewald(2,1,ii) +
     &        diy*gradfield_ewald(2,2,ii) + diz*gradfield_ewald(3,2,ii)
         fzreal = fzreal + dix*gradfield_ewald(3,1,ii) +
     &        diy*gradfield_ewald(3,2,ii) + diz*gradfield_ewald(3,3,ii)
c     quadrupole
         fxreal = fxreal + qixx*hessfield_ewald(1,1,1,ii) +
     &        qixy*hessfield_ewald(2,1,1,ii) +
     &        qixz*hessfield_ewald(3,1,1,ii) + 
     &        qiyy*hessfield_ewald(2,2,1,ii) +
     &        qiyz*hessfield_ewald(3,2,1,ii) + 
     &        qizz*hessfield_ewald(3,3,1,ii)
         fyreal = fyreal + qixx*hessfield_ewald(2,1,1,ii) +
     &        qixy*hessfield_ewald(2,2,1,ii) +
     &        qixz*hessfield_ewald(3,2,1,ii) + 
     &        qiyy*hessfield_ewald(2,2,2,ii) +
     &        qiyz*hessfield_ewald(3,2,2,ii) + 
     &        qizz*hessfield_ewald(3,3,2,ii)
         fzreal = fzreal + qixx*hessfield_ewald(3,1,1,ii) +
     &        qixy*hessfield_ewald(3,2,1,ii) +
     &        qixz*hessfield_ewald(3,3,1,ii) + 
     &        qiyy*hessfield_ewald(3,2,2,ii) +
     &        qiyz*hessfield_ewald(3,3,2,ii) + 
     &        qizz*hessfield_ewald(3,3,3,ii)
c
c     full real space forces needed for scaled interactions
c
c     charge
         fxfix = ci*field(1,ii)
         fyfix = ci*field(2,ii)
         fzfix = ci*field(3,ii)
c     dipole
         fxfix = fxfix + dix*gradfield(1,1,ii) +
     &        diy*gradfield(2,1,ii) + diz*gradfield(3,1,ii)
         fyfix = fyfix + dix*gradfield(2,1,ii) +
     &        diy*gradfield(2,2,ii) + diz*gradfield(3,2,ii)
         fzfix = fzfix + dix*gradfield(3,1,ii) +
     &        diy*gradfield(3,2,ii) + diz*gradfield(3,3,ii)
c     quadrupole
         fxfix = fxfix + qixx*hessfield(1,1,1,ii) + 
     &        qixy*hessfield(2,1,1,ii) +
     &        qixz*hessfield(3,1,1,ii) + qiyy*hessfield(2,2,1,ii) +
     &        qiyz*hessfield(3,2,1,ii) + qizz*hessfield(3,3,1,ii)
         fyfix = fyfix + qixx*hessfield(2,1,1,ii) + 
     &        qixy*hessfield(2,2,1,ii) +
     &        qixz*hessfield(3,2,1,ii) + qiyy*hessfield(2,2,2,ii) +
     &        qiyz*hessfield(3,2,2,ii) + qizz*hessfield(3,3,2,ii)
         fzfix = fzfix + qixx*hessfield(3,1,1,ii) + 
     &        qixy*hessfield(3,2,1,ii) +
     &        qixz*hessfield(3,3,1,ii) + qiyy*hessfield(3,2,2,ii) +
     &        qiyz*hessfield(3,3,2,ii) + qizz*hessfield(3,3,3,ii)
c
c     scaled interactions
c
c     charge
         fxfix = fxfix - ci*fieldm(1,ii)
         fyfix = fyfix - ci*fieldm(2,ii)
         fzfix = fzfix - ci*fieldm(3,ii)
c     dipole
         fxfix = fxfix - dix*gradfieldm(1,1,ii) -
     &        diy*gradfieldm(2,1,ii) - diz*gradfieldm(3,1,ii)
         fyfix = fyfix - dix*gradfieldm(2,1,ii) -
     &        diy*gradfieldm(2,2,ii) - diz*gradfieldm(3,2,ii)
         fzfix = fzfix - dix*gradfieldm(3,1,ii) -
     &        diy*gradfieldm(3,2,ii) - diz*gradfieldm(3,3,ii)
c     quadrupole
         fxfix = fxfix - qixx*hessfieldm(1,1,1,ii) -
     &        qixy*hessfieldm(2,1,1,ii) -
     &        qixz*hessfieldm(3,1,1,ii) - qiyy*hessfieldm(2,2,1,ii) -
     &        qiyz*hessfieldm(3,2,1,ii) - qizz*hessfieldm(3,3,1,ii)
         fyfix = fyfix - qixx*hessfieldm(2,1,1,ii) -
     &        qixy*hessfieldm(2,2,1,ii) -
     &        qixz*hessfieldm(3,2,1,ii) - qiyy*hessfieldm(2,2,2,ii) -
     &        qiyz*hessfieldm(3,2,2,ii) - qizz*hessfieldm(3,3,2,ii)
         fzfix = fzfix - qixx*hessfieldm(3,1,1,ii) -
     &        qixy*hessfieldm(3,2,1,ii) -
     &        qixz*hessfieldm(3,3,1,ii) - qiyy*hessfieldm(3,2,2,ii) -
     &        qiyz*hessfieldm(3,3,2,ii) - qizz*hessfieldm(3,3,3,ii)
c
c     reciprocal space forces
c
c     charge
         fxrecip = ci*field_recip(1,ii)
         fyrecip = ci*field_recip(2,ii)
         fzrecip = ci*field_recip(3,ii)
c     dipole
         fxrecip = fxrecip + dix*gradfield_recip(1,1,ii) +
     &        diy*gradfield_recip(2,1,ii) + diz*gradfield_recip(3,1,ii)
         fyrecip = fyrecip + dix*gradfield_recip(2,1,ii) +
     &        diy*gradfield_recip(2,2,ii) + diz*gradfield_recip(3,2,ii)
         fzrecip = fzrecip + dix*gradfield_recip(3,1,ii) +
     &        diy*gradfield_recip(3,2,ii) + diz*gradfield_recip(3,3,ii)
c     quadrupole
         fxrecip = fxrecip + qixx*hessfield_recip(1,1,1,ii) + 
     &        qixy*hessfield_recip(2,1,1,ii) +
     &        qixz*hessfield_recip(3,1,1,ii) + 
     &        qiyy*hessfield_recip(2,2,1,ii) +
     &        qiyz*hessfield_recip(3,2,1,ii) + 
     &        qizz*hessfield_recip(3,3,1,ii)
         fyrecip = fyrecip + qixx*hessfield_recip(2,1,1,ii) +
     &        qixy*hessfield_recip(2,2,1,ii) +
     &        qixz*hessfield_recip(3,2,1,ii) + 
     &        qiyy*hessfield_recip(2,2,2,ii) +
     &        qiyz*hessfield_recip(3,2,2,ii) + 
     &        qizz*hessfield_recip(3,3,2,ii)
         fzrecip = fzrecip + qixx*hessfield_recip(3,1,1,ii) +
     &        qixy*hessfield_recip(3,2,1,ii) +
     &        qixz*hessfield_recip(3,3,1,ii) + 
     &        qiyy*hessfield_recip(3,2,2,ii) +
     &        qiyz*hessfield_recip(3,3,2,ii) + 
     &        qizz*hessfield_recip(3,3,3,ii)
c
c     apply f constant
c     
         fxreal = f * fxreal
         fyreal = f * fyreal
         fzreal = f * fzreal 
         fxfix = f * fxfix
         fyfix = f * fyfix
         fzfix = f * fzfix
         fxrecip = f * fxrecip
         fyrecip = f * fyrecip
         fzrecip = f * fzrecip
c
c     accumulate permanent multipole forces
c
         dem(1,ii) = fxreal - fxfix + fxrecip
         dem(2,ii) = fyreal - fyfix + fyrecip
         dem(3,ii) = fzreal - fzfix + fzrecip
c
c     calculate permanent multipole torques
c
c     dipole
         trq(1,ii) = diz*field_ewald(2,ii) - diy*field_ewald(3,ii)
         trq(2,ii) = dix*field_ewald(3,ii) - diz*field_ewald(1,ii)
         trq(3,ii) = diy*field_ewald(1,ii) - dix*field_ewald(2,ii)
c     quadrupole
         trq(1,ii) = trq(1,ii) + 
     &        2.0d0*(qizz - qiyy)*gradfield_ewald(3,2,ii)
     &        + qixz*gradfield_ewald(2,1,ii) + 
     &        qiyz*gradfield_ewald(2,2,ii)
     &        - qixy*gradfield_ewald(3,1,ii) - 
     &        qiyz*gradfield_ewald(3,3,ii)
         trq(2,ii) = trq(2,ii) + 
     &        2.0d0*(qixx - qizz)*gradfield_ewald(3,1,ii)
     &        + qixy*gradfield_ewald(3,2,ii) + 
     &        qixz*gradfield_ewald(3,3,ii)
     &        - qixz*gradfield_ewald(1,1,ii) - 
     &        qiyz*gradfield_ewald(2,1,ii)
         trq(3,ii) = trq(3,ii) + 
     &        2.0d0*(qiyy - qixx)*gradfield_ewald(2,1,ii)
     &        + qixy*gradfield_ewald(1,1,ii) + 
     &        qiyz*gradfield_ewald(3,1,ii)
     &        - qixy*gradfield_ewald(2,2,ii) - 
     &        qixz*gradfield_ewald(3,2,ii)
c     dipole
         trq(1,ii) = trq(1,ii) - (diz*field(2,ii) - diy*field(3,ii))
         trq(2,ii) = trq(2,ii) - (dix*field(3,ii) - diz*field(1,ii))
         trq(3,ii) = trq(3,ii) - (diy*field(1,ii) - dix*field(2,ii))
c     quadrupole
         trq(1,ii) = trq(1,ii) - (2.0d0*(qizz - qiyy)*gradfield(3,2,ii)
     &        + qixz*gradfield(2,1,ii) + qiyz*gradfield(2,2,ii)
     &        - qixy*gradfield(3,1,ii) - qiyz*gradfield(3,3,ii))
         trq(2,ii) = trq(2,ii) - (2.0d0*(qixx - qizz)*gradfield(3,1,ii)
     &        + qixy*gradfield(3,2,ii) + qixz*gradfield(3,3,ii)
     &        - qixz*gradfield(1,1,ii) - qiyz*gradfield(2,1,ii))
         trq(3,ii) = trq(3,ii) - (2.0d0*(qiyy - qixx)*gradfield(2,1,ii)
     &        + qixy*gradfield(1,1,ii) + qiyz*gradfield(3,1,ii)
     &        - qixy*gradfield(2,2,ii) - qixz*gradfield(3,2,ii))
c     dipole
         trq(1,ii) = trq(1,ii) + diz*fieldm(2,ii) - diy*fieldm(3,ii)
         trq(2,ii) = trq(2,ii) + dix*fieldm(3,ii) - diz*fieldm(1,ii)
         trq(3,ii) = trq(3,ii) + diy*fieldm(1,ii) - dix*fieldm(2,ii)
c     quadrupole
         trq(1,ii) = trq(1,ii) + 2.0d0*(qizz - qiyy)*gradfieldm(3,2,ii)
     &        + qixz*gradfieldm(2,1,ii) + qiyz*gradfieldm(2,2,ii)
     &        - qixy*gradfieldm(3,1,ii) - qiyz*gradfieldm(3,3,ii)
         trq(2,ii) = trq(2,ii) + 2.0d0*(qixx - qizz)*gradfieldm(3,1,ii)
     &        + qixy*gradfieldm(3,2,ii) + qixz*gradfieldm(3,3,ii)
     &        - qixz*gradfieldm(1,1,ii) - qiyz*gradfieldm(2,1,ii)
         trq(3,ii) = trq(3,ii) + 2.0d0*(qiyy - qixx)*gradfieldm(2,1,ii)
     &        + qixy*gradfieldm(1,1,ii) + qiyz*gradfieldm(3,1,ii)
     &        - qixy*gradfieldm(2,2,ii) - qixz*gradfieldm(3,2,ii)
c     dipole
         trq(1,ii) = trq(1,ii) + diz*field_recip(2,ii) - 
     &        diy*field_recip(3,ii)
         trq(2,ii) = trq(2,ii) + dix*field_recip(3,ii) - 
     &        diz*field_recip(1,ii)
         trq(3,ii) = trq(3,ii) + diy*field_recip(1,ii) - 
     &        dix*field_recip(2,ii)
c     quadrupole
         trq(1,ii) = trq(1,ii) + 
     &        2.0d0*(qizz - qiyy)*gradfield_recip(3,2,ii)
     &        + qixz*gradfield_recip(2,1,ii) + 
     &        qiyz*gradfield_recip(2,2,ii)
     &        - qixy*gradfield_recip(3,1,ii) - 
     &        qiyz*gradfield_recip(3,3,ii)
         trq(2,ii) = trq(2,ii) + 
     &        2.0d0*(qixx - qizz)*gradfield_recip(3,1,ii)
     &        + qixy*gradfield_recip(3,2,ii) + 
     &        qixz*gradfield_recip(3,3,ii)
     &        - qixz*gradfield_recip(1,1,ii) - 
     &        qiyz*gradfield_recip(2,1,ii)
         trq(3,ii) = trq(3,ii) + 
     &        2.0d0*(qiyy - qixx)*gradfield_recip(2,1,ii)
     &        + qixy*gradfield_recip(1,1,ii) + 
     &        qiyz*gradfield_recip(3,1,ii)
     &        - qixy*gradfield_recip(2,2,ii) - 
     &        qixz*gradfield_recip(3,2,ii)
c
         trq(1,ii) = f * trq(1,ii)
         trq(2,ii) = f * trq(2,ii)
         trq(3,ii) = f * trq(3,ii)
      end do
c
c     distribute torques into permanent multipole gradient
c
      do i = 1, n
         frc(1,i) = 0.0d0
         frc(2,i) = 0.0d0
         frc(3,i) = 0.0d0
      end do
      call torque2 (trq,frc)
      do i = 1, n
         dem(1,i) = dem(1,i) + frc(1,i)
         dem(2,i) = dem(2,i) + frc(2,i)
         dem(3,i) = dem(3,i) + frc(3,i)
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         do i = 1, npole
            ii = ipole(i)
            dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
            dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
            dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         do i = 1, npole
            trq(1,i) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
            trq(2,i) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
            trq(3,i) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
         end do
         do i = 1, n
            frc(1,i) = 0.0d0
            frc(2,i) = 0.0d0
            frc(3,i) = 0.0d0
         end do
         call torque2 (trq,frc)
         do i = 1, n
            dem(1,i) = dem(1,i) + frc(1,i)
            dem(2,i) = dem(2,i) + frc(2,i)
            dem(3,i) = dem(3,i) + frc(3,i)
         end do
      end if
      return
      end
