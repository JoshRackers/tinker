c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar1  --  polarization energy and derivatives ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar1" calculates the electrostatic energy and forces due to
c     dipole polarizability interactions
c
c
      subroutine epolar1
      use sizes
      use deriv
      use energi
      use limits
      use mpole
      use potent
      implicit none
      integer i,j,ii
c
c
c     choose the method for summing over multipole interactions
c
      if (use_ewald) then
         call epolar1b
      else
         call epolar1a
      end if
c
c
c     zero out energy terms and analysis which are not in use
c
      if (.not. use_polar) then
         ep = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            do j = 1, 3
               dep(j,ii) = 0.0d0
            end do
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine epolar1a  --  polarization energy and         ##
c     ##                           derivatives via a double loop   ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "epolar1a" calculates the total dipole           
c     polarizability interaction energy and forces
c
c
      subroutine epolar1a
      use sizes
      use atoms
      use atomid
      use bound
      use boxes
      use cell
      use chgpen
      use chgpot
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
      integer i,j,k
      integer ii
      integer ix,iy,iz
      real*8 e,ei,fgrp
      real*8 f
      real*8 fx,fy,fz
      real*8 ci,zi,qi
      real*8 dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8, allocatable :: frc(:,:)
      real*8, allocatable :: trq(:,:)
      real*8, allocatable :: fieldd_damp(:,:)
      real*8, allocatable :: fieldp_damp(:,:)
      real*8, allocatable :: gradfieldd_damp(:,:,:)
      real*8, allocatable :: gradfieldp_damp(:,:,:)
      real*8, allocatable :: udfield_damp(:,:)
      real*8, allocatable :: upfield_damp(:,:)
      real*8, allocatable :: udgradfield_damp(:,:,:)
      real*8, allocatable :: upgradfield_damp(:,:,:)
      real*8, allocatable :: udhessfield_damp(:,:,:,:)
      real*8, allocatable :: uphessfield_damp(:,:,:,:)
      real*8, allocatable :: udnucfieldd_damp(:,:)
      real*8, allocatable :: upnucfieldp_damp(:,:)
      real*8, allocatable :: udfieldd_damp(:,:)
      real*8, allocatable :: upfieldp_damp(:,:)
      real*8, allocatable :: udgradfieldd_damp(:,:,:)
      real*8, allocatable :: upgradfieldp_damp(:,:,:)
      real*8, allocatable :: udhessfieldd_damp(:,:,:,:)
      real*8, allocatable :: uphessfieldp_damp(:,:,:,:)
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (frc(3,n))
      allocate (trq(3,npole))
c
c     zero out the polarization energy and partitioning
c
      ep = 0.0d0
      do i = 1, n
         do j = 1, 3
            dep(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     compute the induced dipoles at each polarizable atom
c
      damp_thole = .false.
      damp_gordon = .false.
      damp_piquemal = .false.
      call induce
c
c     get induced field gradient hessian for forces
c
      damp_thole = .false.
      damp_gordon = .false.
      damp_piquemal = .false.
      if (mutualdamp .eq. "GORDON") damp_gordon = .true.
      if (mutualdamp .eq. "PIQUEMAL") damp_piquemal = .true.
      if (mutualdamp .eq. "THOLE") damp_thole = .true.
      call mutualfield3
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldd_damp(3,npole))
      allocate (fieldp_damp(3,npole))
      allocate (gradfieldd_damp(3,3,npole))
      allocate (gradfieldp_damp(3,3,npole))
c
      allocate (udfield_damp(3,npole))
      allocate (upfield_damp(3,npole))
      allocate (udgradfield_damp(3,3,npole))
      allocate (upgradfield_damp(3,3,npole))
      allocate (udhessfield_damp(3,3,3,npole))
      allocate (uphessfield_damp(3,3,3,npole))
c
      allocate (udnucfieldd_damp(3,npole))
      allocate (upnucfieldp_damp(3,npole))
      allocate (udfieldd_damp(3,npole))
      allocate (upfieldp_damp(3,npole))
      allocate (udgradfieldd_damp(3,3,npole))
      allocate (upgradfieldp_damp(3,3,npole))
      allocate (udhessfieldd_damp(3,3,3,npole))
      allocate (uphessfieldp_damp(3,3,3,npole))
c
c     check what kind of direct damping was used
c
      if (directdamp .eq. "GORDON") then
         fieldd_damp = fieldd_gordon
         fieldp_damp = fieldp_gordon
         gradfieldd_damp = gradfieldd_gordon
         gradfieldp_damp = gradfieldp_gordon
      else if (directdamp .eq. "PIQUEMAL") then
         fieldd_damp = fieldd_piquemal
         fieldp_damp = fieldp_piquemal
         gradfieldd_damp = gradfieldd_piquemal
         gradfieldp_damp = gradfieldp_piquemal
      else if (directdamp .eq. "THOLE") then
         fieldd_damp = fieldd_thole
         fieldp_damp = fieldp_thole
         gradfieldd_damp = gradfieldd_thole
         gradfieldp_damp = gradfieldp_thole
      end if
c
c     check what kind of mutual damping was used
c
      if (mutualdamp .eq. "GORDON") then
         udfield_damp = udfield_gordon
         upfield_damp = upfield_gordon
         udgradfield_damp = udgradfield_gordon
         upgradfield_damp = upgradfield_gordon
         udhessfield_damp = udhessfield_gordon
         uphessfield_damp = uphessfield_gordon
c
         udnucfieldd_damp = udnucfieldd_gordon
         upnucfieldp_damp = upnucfieldp_gordon
         udfieldd_damp = udfieldd_gordon
         upfieldp_damp = upfieldp_gordon
         udgradfieldd_damp = udgradfieldd_gordon
         upgradfieldp_damp = upgradfieldp_gordon
         udhessfieldd_damp = udhessfieldd_gordon
         uphessfieldp_damp = uphessfieldp_gordon
      else if (mutualdamp .eq. "PIQUEMAL") then
         udfield_damp = udfield_piquemal
         upfield_damp = upfield_piquemal
         udgradfield_damp = udgradfield_piquemal
         upgradfield_damp = upgradfield_piquemal
         udhessfield_damp = udhessfield_piquemal
         uphessfield_damp = uphessfield_piquemal
c
         udfieldd_damp = udfieldd_piquemal
         upfieldp_damp = upfieldp_piquemal
         udgradfieldd_damp = udgradfieldd_piquemal
         upgradfieldp_damp = upgradfieldp_piquemal
         udhessfieldd_damp = udhessfieldd_piquemal
         uphessfieldp_damp = uphessfieldp_piquemal
      else if (mutualdamp .eq. "THOLE") then
         udfield_damp = udfield_thole
         upfield_damp = upfield_thole
         udgradfield_damp = udgradfield_thole
         upgradfield_damp = upgradfield_thole
         udhessfield_damp = udhessfield_thole
         uphessfield_damp = uphessfield_thole
c
         udfieldd_damp = udfieldd_thole
         upfieldp_damp = upfieldp_thole
         udgradfieldd_damp = udgradfieldd_thole
         upgradfieldp_damp = upgradfieldp_thole
         udhessfieldd_damp = udhessfieldd_thole
         uphessfieldp_damp = uphessfieldp_thole
      end if
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
c      mode = 'MPOLE'
c      call switch (mode)
c
c     calculate the total multipole interaction energy
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
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
c
c     split nuclear and electronic charge
c     
         if (mutualdamp .eq. "GORDON") then
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
c
c     calculate forces involving induced dipoles
c
c
c     CHECK ABOUT D AND P INDUCED DIPOLES
c
c
c     force on induced dipole from permanent field gradient
c
         fx = uix*gradfieldp_damp(1,1,ii) + 
     &        uiy*gradfieldp_damp(2,1,ii) +
     &        uiz*gradfieldp_damp(3,1,ii)
         fy = uix*gradfieldp_damp(2,1,ii) + 
     &        uiy*gradfieldp_damp(2,2,ii) +
     &        uiz*gradfieldp_damp(3,2,ii)
         fz = uix*gradfieldp_damp(3,1,ii) + 
     &        uiy*gradfieldp_damp(3,2,ii) +
     &        uiz*gradfieldp_damp(3,3,ii)
c
c     force on induced dipole from mutual field gradient
c
         if (poltyp .eq. 'MUTUAL') then
c            fx = fx + uix*udgradfield_damp(1,1,ii) +
c     &           uiy*udgradfield_damp(2,1,ii) +
c     &           uiz*udgradfield_damp(3,1,ii)
c            fy = fy + uix*udgradfield_damp(2,1,ii) +
c     &           uiy*udgradfield_damp(2,2,ii) +
c     &           uiz*udgradfield_damp(3,2,ii)
c            fz = fz + uix*udgradfield_damp(3,1,ii) +
c     &           uiy*udgradfield_damp(3,2,ii) +
c     &           uiz*udgradfield_damp(3,3,ii)
c
            fx = fx + uix*upgradfield_damp(1,1,ii) +
     &           uiy*upgradfield_damp(2,1,ii) +
     &           uiz*upgradfield_damp(3,1,ii)
            fy = fy + uix*upgradfield_damp(2,1,ii) +
     &           uiy*upgradfield_damp(2,2,ii) +
     &           uiz*upgradfield_damp(3,2,ii)
            fz = fz + uix*upgradfield_damp(3,1,ii) +
     &           uiy*upgradfield_damp(3,2,ii) +
     &           uiz*upgradfield_damp(3,3,ii)
         end if
c
c     force on permanent moments from induced dipoles
c 
c
c     should this be p induced dipoles with p exclusion rules??????
c
         if (mutualdamp .eq. "GORDON") then
c     charges - nuclear
            fx = fx + zi*upnucfieldp_damp(1,ii)
            fy = fy + zi*upnucfieldp_damp(2,ii)
            fz = fz + zi*upnucfieldp_damp(3,ii)
c     charges - electrons
            fx = fx + qi*upfieldp_damp(1,ii)
            fy = fy + qi*upfieldp_damp(2,ii)
            fz = fz + qi*upfieldp_damp(3,ii)
         else if (mutualdamp .eq. "THOLE") then
            fx = fx + ci*upfieldp_damp(1,ii)
            fy = fy + ci*upfieldp_damp(2,ii)
            fz = fz + ci*upfieldp_damp(3,ii)
         end if
c     dipoles
         fx = fx + dix*upgradfieldp_damp(1,1,ii) +
     &        diy*upgradfieldp_damp(2,1,ii) +
     &        diz*upgradfieldp_damp(3,1,ii)
         fy = fy + dix*upgradfieldp_damp(2,1,ii) +
     &        diy*upgradfieldp_damp(2,2,ii) +
     &        diz*upgradfieldp_damp(3,2,ii)
         fz = fz + dix*upgradfieldp_damp(3,1,ii) +
     &        diy*upgradfieldp_damp(3,2,ii) +
     &        diz*upgradfieldp_damp(3,3,ii)
c     quadrupoles
         fx = fx + qixx*uphessfieldp_damp(1,1,1,ii) +
     &        qixy*uphessfieldp_damp(2,1,1,ii) +
     &        qixz*uphessfieldp_damp(3,1,1,ii) + 
     &        qiyy*uphessfieldp_damp(2,2,1,ii) + 
     &        qiyz*uphessfieldp_damp(3,2,1,ii) + 
     &        qizz*uphessfieldp_damp(3,3,1,ii)
         fy = fy + qixx*uphessfieldp_damp(2,1,1,ii) + 
     &        qixy*uphessfieldp_damp(2,2,1,ii) +
     &        qixz*uphessfieldp_damp(3,2,1,ii) + 
     &        qiyy*uphessfieldp_damp(2,2,2,ii) +
     &        qiyz*uphessfieldp_damp(3,2,2,ii) + 
     &        qizz*uphessfieldp_damp(3,3,2,ii)
         fz = fz + qixx*uphessfieldp_damp(3,1,1,ii) + 
     &        qixy*uphessfieldp_damp(3,2,1,ii) +
     &        qixz*uphessfieldp_damp(3,3,1,ii) + 
     &        qiyy*uphessfieldp_damp(3,2,2,ii) +
     &        qiyz*uphessfieldp_damp(3,3,2,ii) + 
     &        qizz*uphessfieldp_damp(3,3,3,ii)
c
c     increment the induced dipole gradient
c
         dep(1,ii) = f * fx
         dep(2,ii) = f * fy
         dep(3,ii) = f * fz
c
c     calculate the induced dipole torques
c
         trq(1,ii) = uiz*fieldp_damp(2,ii) - uiy*fieldp_damp(3,ii)
         trq(2,ii) = uix*fieldp_damp(3,ii) - uiz*fieldp_damp(1,ii)
         trq(3,ii) = uiy*fieldp_damp(1,ii) - uix*fieldp_damp(2,ii)
c
c     induced dipole torques from mutual field
c
         if (poltyp .eq. 'MUTUAL') then
            trq(1,ii) = trq(1,ii) + uiz*upfield_damp(2,ii) - 
     &           uiy*upfield_damp(3,ii)
            trq(2,ii) = trq(2,ii) + uix*upfield_damp(3,ii) - 
     &           uiz*upfield_damp(1,ii)
            trq(3,ii) = trq(3,ii) + uiy*upfield_damp(1,ii) - 
     &           uix*upfield_damp(2,ii)
         end if
c
c     calculate torques on permanent moments from mutual field
c
c     dipole
         trq(1,ii) = trq(1,ii) + diz*upfieldp_damp(2,ii) -
     &        diy*upfieldp_damp(3,ii)
         trq(2,ii) = trq(2,ii) + dix*upfieldp_damp(3,ii) -
     &        diz*upfieldp_damp(1,ii)
         trq(3,ii) = trq(3,ii) + diy*upfieldp_damp(1,ii) -
     &        dix*upfieldp_damp(2,ii)
c     quadrupole
         trq(1,ii) = trq(1,ii) + 2.0d0*(qizz - qiyy)*
     &        upgradfieldp_damp(3,2,ii)
     &        + qixz*upgradfieldp_damp(2,1,ii) + 
     &        qiyz*upgradfieldp_damp(2,2,ii)
     &        - qixy*upgradfieldp_damp(3,1,ii) - 
     &        qiyz*upgradfieldp_damp(3,3,ii)
         trq(2,ii) = trq(2,ii) + 2.0d0*(qixx - qizz)*
     &        upgradfieldp_damp(3,1,ii)
     &        + qixy*upgradfieldp_damp(3,2,ii) + 
     &        qixz*upgradfieldp_damp(3,3,ii)
     &        - qixz*upgradfieldp_damp(1,1,ii) - 
     &        qiyz*upgradfieldp_damp(2,1,ii)
         trq(3,ii) = trq(3,ii) + 2.0d0*(qiyy - qixx)*
     &        upgradfieldp_damp(2,1,ii)
     &        + qixy*upgradfieldp_damp(1,1,ii) + 
     &        qiyz*upgradfieldp_damp(3,1,ii)
     &        - qixy*upgradfieldp_damp(2,2,ii) - 
     &        qixz*upgradfieldp_damp(3,2,ii)
c
c     apply f constant to torques
c
         trq(1,ii) = f * trq(1,ii)
         trq(2,ii) = f * trq(2,ii)
         trq(3,ii) = f * trq(3,ii)
      end do
c
c     distribute torques into polarization gradient
c
      do i = 1, n
         frc(1,i) = 0.0d0
         frc(2,i) = 0.0d0
         frc(3,i) = 0.0d0
      end do
      call torque2 (trq,frc)
      do i = 1, n
         dep(1,i) = dep(1,i) + frc(1,i)
         dep(2,i) = dep(2,i) + frc(2,i)
         dep(3,i) = dep(3,i) + frc(3,i)
      end do
      deallocate (frc)
      deallocate (trq)
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar1b  --  ewald summation polarization energy ##
c     ##                           and derivatives                     ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar1b" calculates the dipole polarizability interaction 
c     energy and forces using a particle mesh Ewald summation
c
c
      subroutine epolar1b
      use sizes
      use action
      use analyz
      use atoms
      use boxes
      use chgpen
      use chgpot
      use deriv
      use energi
      use ewald
      use inter
      use limits
      use math
      use mpole
      use polar
      use polpot
      use potderivs
      use potent
      implicit none
      integer i,ii,j
      real*8 e,ei,eintra
      real*8 ereal,eself,efix,erecip
      real*8 fx,fy,fz
      real*8 f,term,term2,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 ci,qi,zi
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 field_tot(3),upfield_tot(3),upnucfield_tot(3)
      real*8 gradfield_tot(3,3),upgradfield_tot(3,3)
      real*8 uphessfield_tot(3,3,3)
      real*8, allocatable :: frc(:,:)
      real*8, allocatable :: trq(:,:)
      real*8, allocatable :: fieldd_damp(:,:)
      real*8, allocatable :: fieldp_damp(:,:)
      real*8, allocatable :: gradfieldd_damp(:,:,:)
      real*8, allocatable :: gradfieldp_damp(:,:,:)
      real*8, allocatable :: udfield_damp(:,:)
      real*8, allocatable :: upfield_damp(:,:)
      real*8, allocatable :: udgradfield_damp(:,:,:)
      real*8, allocatable :: upgradfield_damp(:,:,:)
      real*8, allocatable :: udhessfield_damp(:,:,:,:)
      real*8, allocatable :: uphessfield_damp(:,:,:,:)
      real*8, allocatable :: udnucfieldd_damp(:,:)
      real*8, allocatable :: upnucfieldp_damp(:,:)
      real*8, allocatable :: udfieldd_damp(:,:)
      real*8, allocatable :: upfieldp_damp(:,:)
      real*8, allocatable :: udgradfieldd_damp(:,:,:)
      real*8, allocatable :: upgradfieldp_damp(:,:,:)
      real*8, allocatable :: udhessfieldd_damp(:,:,:,:)
      real*8, allocatable :: uphessfieldp_damp(:,:,:,:)
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (frc(3,n))
      allocate (trq(3,npole))
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
      term2 = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      fterm = -f * aewald / sqrtpi
c
c     compute the induced dipoles at each polarizable atom
c
      damp_thole = .false.
      damp_gordon = .false.
      damp_piquemal = .false.
      call induce
c
c     get mutual field, field gradient and field hessian
c
      call mutualfield3
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldd_damp(3,npole))
      allocate (fieldp_damp(3,npole))
      allocate (gradfieldd_damp(3,3,npole))
      allocate (gradfieldp_damp(3,3,npole))
c
      allocate (udfield_damp(3,npole))
      allocate (upfield_damp(3,npole))
      allocate (udgradfield_damp(3,3,npole))
      allocate (upgradfield_damp(3,3,npole))
      allocate (udhessfield_damp(3,3,3,npole))
      allocate (uphessfield_damp(3,3,3,npole))
c
      allocate (udnucfieldd_damp(3,npole))
      allocate (upnucfieldp_damp(3,npole))
      allocate (udfieldd_damp(3,npole))
      allocate (upfieldp_damp(3,npole))
      allocate (udgradfieldd_damp(3,3,npole))
      allocate (upgradfieldp_damp(3,3,npole))
      allocate (udhessfieldd_damp(3,3,3,npole))
      allocate (uphessfieldp_damp(3,3,3,npole))
c
c     check what kind of direct damping was used
c
      if (directdamp .eq. "GORDON") then
         fieldd_damp = fieldd_gordon
         fieldp_damp = fieldp_gordon
         gradfieldd_damp = gradfieldd_gordon
         gradfieldp_damp = gradfieldp_gordon
      else if (directdamp .eq. "PIQUEMAL") then
         fieldd_damp = fieldd_piquemal
         fieldp_damp = fieldp_piquemal
         gradfieldd_damp = gradfieldd_piquemal
         gradfieldp_damp = gradfieldp_piquemal
      else if (directdamp .eq. "THOLE") then
         fieldd_damp = fieldd_thole
         fieldp_damp = fieldp_thole
         gradfieldd_damp = gradfieldd_thole
         gradfieldp_damp = gradfieldp_thole
      end if
c
c     check what kind of mutual damping was used
c
      if (mutualdamp .eq. "GORDON") then
         udfield_damp = udfield_gordon
         upfield_damp = upfield_gordon
         udgradfield_damp = udgradfield_gordon
         upgradfield_damp = upgradfield_gordon
         udhessfield_damp = udhessfield_gordon
         uphessfield_damp = uphessfield_gordon
c
         udnucfieldd_damp = udnucfieldd_gordon
         upnucfieldp_damp = upnucfieldp_gordon
         udfieldd_damp = udfieldd_gordon
         upfieldp_damp = upfieldp_gordon
         udgradfieldd_damp = udgradfieldd_gordon
         upgradfieldp_damp = upgradfieldp_gordon
         udhessfieldd_damp = udhessfieldd_gordon
         uphessfieldp_damp = uphessfieldp_gordon
      else if (mutualdamp .eq. "PIQUEMAL") then
         udfield_damp = udfield_piquemal
         upfield_damp = upfield_piquemal
         udgradfield_damp = udgradfield_piquemal
         upgradfield_damp = upgradfield_piquemal
         udhessfield_damp = udhessfield_piquemal
         uphessfield_damp = uphessfield_piquemal
c
         udfieldd_damp = udfieldd_piquemal
         upfieldp_damp = upfieldp_piquemal
         udgradfieldd_damp = udgradfieldd_piquemal
         upgradfieldp_damp = upgradfieldp_piquemal
         udhessfieldd_damp = udhessfieldd_piquemal
         uphessfieldp_damp = uphessfieldp_piquemal
      else if (mutualdamp .eq. "THOLE") then
         udfield_damp = udfield_thole
         upfield_damp = upfield_thole
         udgradfield_damp = udgradfield_thole
         upgradfield_damp = upgradfield_thole
         udhessfield_damp = udhessfield_thole
         uphessfield_damp = uphessfield_thole
c
         udfieldd_damp = udfieldd_thole
         upfieldp_damp = upfieldp_thole
         udgradfieldd_damp = udgradfieldd_thole
         upgradfieldp_damp = upgradfieldp_thole
         udhessfieldd_damp = udhessfieldd_thole
         uphessfieldp_damp = uphessfieldp_thole
      end if
c
c     get reciprocal space field, field gradient and field hessian
c
      call mutualrecip3
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
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         qixy = 2.0d0*qixy
         qixz = 2.0d0*qixz
         qiyz = 2.0d0*qiyz
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
c
c     calculate induce dipole forces and torques
c
c
c     induced dipole with permanent field gradient
c
         field_tot(:) = field_ewald(:,ii) - field(:,ii) +
     &        fieldp_damp(:,ii) + field_recip(:,ii)
         gradfield_tot(:,:) = gradfield_ewald(:,:,ii) - 
     &        gradfield(:,:,ii) + gradfieldp_damp(:,:,ii) + 
     &        gradfield_recip(:,:,ii)
         fx = uix*gradfield_tot(1,1) + 
     &        uiy*gradfield_tot(2,1) +
     &        uiz*gradfield_tot(3,1)
         fy = uix*gradfield_tot(2,1) + 
     &        uiy*gradfield_tot(2,2) +
     &        uiz*gradfield_tot(3,2)
         fz = uix*gradfield_tot(3,1) + 
     &        uiy*gradfield_tot(3,2) +
     &        uiz*gradfield_tot(3,3)
c     torques
         trq(1,ii) = uiz*field_tot(2) - uiy*field_tot(3)
         trq(2,ii) = uix*field_tot(3) - uiz*field_tot(1)
         trq(3,ii) = uiy*field_tot(1) - uix*field_tot(2)
c
c     induced dipole with induced dipole field gradient
c
         if (poltyp .eq. 'MUTUAL') then
            upfield_tot(:) = upfield_ewald(:,ii) - upfield(:,ii) +
     &           upfield_damp(:,ii) + upfield_recip(:,ii)
            upgradfield_tot(:,:) = upgradfield_ewald(:,:,ii) - 
     &           upgradfield(:,:,ii) + upgradfield_damp(:,:,ii) + 
     &           upgradfield_recip(:,:,ii)
            fx = fx + uix*upgradfield_tot(1,1) +
     &           uiy*upgradfield_tot(2,1) +
     &           uiz*upgradfield_tot(3,1)
            fy = fy + uix*upgradfield_tot(2,1) +
     &           uiy*upgradfield_tot(2,2) +
     &           uiz*upgradfield_tot(3,2)
            fz = fz + uix*upgradfield_tot(3,1) +
     &           uiy*upgradfield_tot(3,2) +
     &           uiz*upgradfield_tot(3,3)
c     torques
            trq(1,ii) = trq(1,ii) + uiz*upfield_tot(2) - 
     &           uiy*upfield_tot(3)
            trq(2,ii) = trq(2,ii) + uix*upfield_tot(3) - 
     &           uiz*upfield_tot(1)
            trq(3,ii) = trq(3,ii) + uiy*upfield_tot(1) - 
     &           uix*upfield_tot(2)
         end if
c
c     permanent moments with induced dipole field, 
c     field gradient, field hessian
c
         upnucfield_tot(:) = upfield_ewald(:,ii) - upfield(:,ii) +
     &        upnucfieldp_damp(:,ii) + upfield_recip(:,ii)
         upfield_tot(:) = upfield_ewald(:,ii) - upfield(:,ii) + 
     &        upfieldp_damp(:,ii) + upfield_recip(:,ii)
         upgradfield_tot(:,:) = upgradfield_ewald(:,:,ii) - 
     &        upgradfield(:,:,ii) + upgradfieldp_damp(:,:,ii) + 
     &        upgradfield_recip(:,:,ii)
         uphessfield_tot(:,:,:) = uphessfield_ewald(:,:,:,ii) - 
     &        uphessfield(:,:,:,ii) + uphessfieldp_damp(:,:,:,ii) + 
     &        uphessfield_recip(:,:,:,ii)
         if (mutualdamp .eq. "GORDON") then
c     charges - nuclear
            fx = fx + zi*upnucfield_tot(1)
            fy = fy + zi*upnucfield_tot(2)
            fz = fz + zi*upnucfield_tot(3)
c     charges - electrons
            fx = fx + qi*upfield_tot(1)
            fy = fy + qi*upfield_tot(2)
            fz = fz + qi*upfield_tot(3)
         else if (mutualdamp .eq. "THOLE") then
            fx = fx + ci*upfield_tot(1)
            fy = fy + ci*upfield_tot(2)
            fz = fz + ci*upfield_tot(3)
         end if
c     dipoles
         fx = fx + dix*upgradfield_tot(1,1) +
     &        diy*upgradfield_tot(2,1) +
     &        diz*upgradfield_tot(3,1)
         fy = fy + dix*upgradfield_tot(2,1) +
     &        diy*upgradfield_tot(2,2) +
     &        diz*upgradfield_tot(3,2)
         fz = fz + dix*upgradfield_tot(3,1) +
     &        diy*upgradfield_tot(3,2) +
     &        diz*upgradfield_tot(3,3)
c     quadrupoles
         fx = fx + qixx*uphessfield_tot(1,1,1) +
     &        qixy*uphessfield_tot(2,1,1) +
     &        qixz*uphessfield_tot(3,1,1) + 
     &        qiyy*uphessfield_tot(2,2,1) + 
     &        qiyz*uphessfield_tot(3,2,1) + 
     &        qizz*uphessfield_tot(3,3,1)
         fy = fy + qixx*uphessfield_tot(2,1,1) + 
     &        qixy*uphessfield_tot(2,2,1) +
     &        qixz*uphessfield_tot(3,2,1) + 
     &        qiyy*uphessfield_tot(2,2,2) +
     &        qiyz*uphessfield_tot(3,2,2) + 
     &        qizz*uphessfield_tot(3,3,2)
         fz = fz + qixx*uphessfield_tot(3,1,1) + 
     &        qixy*uphessfield_tot(3,2,1) +
     &        qixz*uphessfield_tot(3,3,1) + 
     &        qiyy*uphessfield_tot(3,2,2) +
     &        qiyz*uphessfield_tot(3,3,2) + 
     &        qizz*uphessfield_tot(3,3,3)
c     torques
c     dipole
         trq(1,ii) = trq(1,ii) + diz*upfield_tot(2) -
     &        diy*upfield_tot(3)
         trq(2,ii) = trq(2,ii) + dix*upfield_tot(3) -
     &        diz*upfield_tot(1)
         trq(3,ii) = trq(3,ii) + diy*upfield_tot(1) -
     &        dix*upfield_tot(2)
c     quadrupole
         trq(1,ii) = trq(1,ii) + 2.0d0*(qizz - qiyy)*
     &        upgradfield_tot(3,2)
     &        + qixz*upgradfield_tot(2,1) + 
     &        qiyz*upgradfield_tot(2,2)
     &        - qixy*upgradfield_tot(3,1) - 
     &        qiyz*upgradfield_tot(3,3)
         trq(2,ii) = trq(2,ii) + 2.0d0*(qixx - qizz)*
     &        upgradfield_tot(3,1)
     &        + qixy*upgradfield_tot(3,2) + 
     &        qixz*upgradfield_tot(3,3)
     &        - qixz*upgradfield_tot(1,1) - 
     &        qiyz*upgradfield_tot(2,1)
         trq(3,ii) = trq(3,ii) + 2.0d0*(qiyy - qixx)*
     &        upgradfield_tot(2,1)
     &        + qixy*upgradfield_tot(1,1) + 
     &        qiyz*upgradfield_tot(3,1)
     &        - qixy*upgradfield_tot(2,2) - 
     &        qixz*upgradfield_tot(3,2)
c
c     self torque due to induced dipole
c
c     WHY DO I GET A GOOD ANSWER WHEN THIS IS TURNED OFF?????
c         trq(1,ii) = trq(1,ii) + term2 * (diy*uiz-diz*uiy)
c         trq(2,ii) = trq(2,ii) + term2 * (diz*uix-dix*uiz)
c         trq(3,ii) = trq(3,ii) + term2 * (dix*uiy-diy*uix)
c
c     apply f constant to forces and torques
c
         dep(1,ii) = f * fx
         dep(2,ii) = f * fy
         dep(3,ii) = f * fz
         trq(1,ii) = f * trq(1,ii)
         trq(2,ii) = f * trq(2,ii)
         trq(3,ii) = f * trq(3,ii)
      end do
c
c     distribute torques into polarization gradient
c
      do i = 1, n
         frc(1,i) = 0.0d0
         frc(2,i) = 0.0d0
         frc(3,i) = 0.0d0
      end do
      call torque2 (trq,frc)
      do i = 1, n
         dep(1,i) = dep(1,i) + frc(1,i)
         dep(2,i) = dep(2,i) + frc(2,i)
         dep(3,i) = dep(3,i) + frc(3,i)
      end do
      deallocate (frc)
      deallocate (trq)
c
c     calculate cell boundary correction term
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
         nep = nep + 1
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
      end if
      return
      end
