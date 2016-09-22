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
      use bound
      use boxes
      use cell
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
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
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
      call induce
c
c     get induced field gradient hessian for forces
c
      call mutualfield3
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
c     contract induced dipoles with thole damped permanent field
c
         ei = 0.5*(uix*fieldp_thole(1,ii) + uiy*fieldp_thole(2,ii) + 
     &        uiz*fieldp_thole(3,ii))
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
         fx = uix*gradfieldp_thole(1,1,ii) + 
     &        uiy*gradfieldp_thole(2,1,ii) +
     &        uiz*gradfieldp_thole(3,1,ii)
         fy = uix*gradfieldp_thole(2,1,ii) + 
     &        uiy*gradfieldp_thole(2,2,ii) +
     &        uiz*gradfieldp_thole(3,2,ii)
         fz = uix*gradfieldp_thole(3,1,ii) + 
     &        uiy*gradfieldp_thole(3,2,ii) +
     &        uiz*gradfieldp_thole(3,3,ii)
c
c     force on induced dipole from mutual field gradient
c
         if (poltyp .eq. 'MUTUAL') then
c            fx = fx + uix*udgradfield_thole(1,1,ii) +
c     &           uiy*udgradfield_thole(2,1,ii) +
c     &           uiz*udgradfield_thole(3,1,ii)
c            fy = fy + uix*udgradfield_thole(2,1,ii) +
c     &           uiy*udgradfield_thole(2,2,ii) +
c     &           uiz*udgradfield_thole(3,2,ii)
c            fz = fz + uix*udgradfield_thole(3,1,ii) +
c     &           uiy*udgradfield_thole(3,2,ii) +
c     &           uiz*udgradfield_thole(3,3,ii)
c
            fx = fx + uix*upgradfield_thole(1,1,ii) +
     &           uiy*upgradfield_thole(2,1,ii) +
     &           uiz*upgradfield_thole(3,1,ii)
            fy = fy + uix*upgradfield_thole(2,1,ii) +
     &           uiy*upgradfield_thole(2,2,ii) +
     &           uiz*upgradfield_thole(3,2,ii)
            fz = fz + uix*upgradfield_thole(3,1,ii) +
     &           uiy*upgradfield_thole(3,2,ii) +
     &           uiz*upgradfield_thole(3,3,ii)
         end if
c
c     force on permanent moments from induced dipoles
c 
c
c     should this be p induced dipoles with p exclusion rules??????
c
c     charges
         fx = fx + ci*upfieldp_thole(1,ii)
         fy = fy + ci*upfieldp_thole(2,ii)
         fz = fz + ci*upfieldp_thole(3,ii)
c     dipoles
         fx = fx + dix*upgradfieldp_thole(1,1,ii) +
     &        diy*upgradfieldp_thole(2,1,ii) +
     &        diz*upgradfieldp_thole(3,1,ii)
         fy = fy + dix*upgradfieldp_thole(2,1,ii) +
     &        diy*upgradfieldp_thole(2,2,ii) +
     &        diz*upgradfieldp_thole(3,2,ii)
         fz = fz + dix*upgradfieldp_thole(3,1,ii) +
     &        diy*upgradfieldp_thole(3,2,ii) +
     &        diz*upgradfieldp_thole(3,3,ii)
c     quadrupoles
         fx = fx + qixx*uphessfieldp_thole(1,1,1,ii) +
     &        qixy*uphessfieldp_thole(2,1,1,ii) +
     &        qixz*uphessfieldp_thole(3,1,1,ii) + 
     &        qiyy*uphessfieldp_thole(2,2,1,ii) + 
     &        qiyz*uphessfieldp_thole(3,2,1,ii) + 
     &        qizz*uphessfieldp_thole(3,3,1,ii)
         fy = fy + qixx*uphessfieldp_thole(2,1,1,ii) + 
     &        qixy*uphessfieldp_thole(2,2,1,ii) +
     &        qixz*uphessfieldp_thole(3,2,1,ii) + 
     &        qiyy*uphessfieldp_thole(2,2,2,ii) +
     &        qiyz*uphessfieldp_thole(3,2,2,ii) + 
     &        qizz*uphessfieldp_thole(3,3,2,ii)
         fz = fz + qixx*uphessfieldp_thole(3,1,1,ii) + 
     &        qixy*uphessfieldp_thole(3,2,1,ii) +
     &        qixz*uphessfieldp_thole(3,3,1,ii) + 
     &        qiyy*uphessfieldp_thole(3,2,2,ii) +
     &        qiyz*uphessfieldp_thole(3,3,2,ii) + 
     &        qizz*uphessfieldp_thole(3,3,3,ii)
c
c     increment the induced dipole gradient
c
         dep(1,ii) = f * fx
         dep(2,ii) = f * fy
         dep(3,ii) = f * fz
c
c     calculate the induced dipole torques
c
         trq(1,ii) = uiz*fieldp_thole(2,ii) - uiy*fieldp_thole(3,ii)
         trq(2,ii) = uix*fieldp_thole(3,ii) - uiz*fieldp_thole(1,ii)
         trq(3,ii) = uiy*fieldp_thole(1,ii) - uix*fieldp_thole(2,ii)
c
c     induced dipole torques from mutual field
c
         if (poltyp .eq. 'MUTUAL') then
            trq(1,ii) = trq(1,ii) + uiz*upfield_thole(2,ii) - 
     &           uiy*upfield_thole(3,ii)
            trq(2,ii) = trq(2,ii) + uix*upfield_thole(3,ii) - 
     &           uiz*upfield_thole(1,ii)
            trq(3,ii) = trq(3,ii) + uiy*upfield_thole(1,ii) - 
     &           uix*upfield_thole(2,ii)
         end if
c
c     calculate torques on permanent moments from mutual field
c
c     dipole
         trq(1,ii) = trq(1,ii) + diz*upfieldp_thole(2,ii) -
     &        diy*upfieldp_thole(3,ii)
         trq(2,ii) = trq(2,ii) + dix*upfieldp_thole(3,ii) -
     &        diz*upfieldp_thole(1,ii)
         trq(3,ii) = trq(3,ii) + diy*upfieldp_thole(1,ii) -
     &        dix*upfieldp_thole(2,ii)
c     quadrupole
         trq(1,ii) = trq(1,ii) + 2.0d0*(qizz - qiyy)*
     &        upgradfieldp_thole(3,2,ii)
     &        + qixz*upgradfieldp_thole(2,1,ii) + 
     &        qiyz*upgradfieldp_thole(2,2,ii)
     &        - qixy*upgradfieldp_thole(3,1,ii) - 
     &        qiyz*upgradfieldp_thole(3,3,ii)
         trq(2,ii) = trq(2,ii) + 2.0d0*(qixx - qizz)*
     &        upgradfieldp_thole(3,1,ii)
     &        + qixy*upgradfieldp_thole(3,2,ii) + 
     &        qixz*upgradfieldp_thole(3,3,ii)
     &        - qixz*upgradfieldp_thole(1,1,ii) - 
     &        qiyz*upgradfieldp_thole(2,1,ii)
         trq(3,ii) = trq(3,ii) + 2.0d0*(qiyy - qixx)*
     &        upgradfieldp_thole(2,1,ii)
     &        + qixy*upgradfieldp_thole(1,1,ii) + 
     &        qiyz*upgradfieldp_thole(3,1,ii)
     &        - qixy*upgradfieldp_thole(2,2,ii) - 
     &        qixz*upgradfieldp_thole(3,2,ii)
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
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 field_tot(3),upfield_tot(3)
      real*8 gradfield_tot(3,3),upgradfield_tot(3,3)
      real*8 uphessfield_tot(3,3,3)
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
      call induce
c
c     get mutual field, field gradient and field hessian
c
      call mutualfield3
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
     &        uix*fieldp_thole(1,ii) + uiy*fieldp_thole(2,ii) +
     &        uiz*fieldp_thole(3,ii))
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
     &        fieldp_thole(:,ii) + field_recip(:,ii)
         gradfield_tot(:,:) = gradfield_ewald(:,:,ii) - 
     &        gradfield(:,:,ii) + gradfieldp_thole(:,:,ii) + 
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
     &           upfield_thole(:,ii) + upfield_recip(:,ii)
            upgradfield_tot(:,:) = upgradfield_ewald(:,:,ii) - 
     &           upgradfield(:,:,ii) + upgradfield_thole(:,:,ii) + 
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
         upfield_tot(:) = upfield_ewald(:,ii) - upfield(:,ii) + 
     &        upfieldp_thole(:,ii) + upfield_recip(:,ii)
         upgradfield_tot(:,:) = upgradfield_ewald(:,:,ii) - 
     &        upgradfield(:,:,ii) + upgradfieldp_thole(:,:,ii) + 
     &        upgradfield_recip(:,:,ii)
         uphessfield_tot(:,:,:) = uphessfield_ewald(:,:,:,ii) - 
     &        uphessfield(:,:,:,ii) + uphessfieldp_thole(:,:,:,ii) + 
     &        uphessfield_recip(:,:,:,ii)
c     charges
         fx = fx + ci*upfield_tot(1)
         fy = fy + ci*upfield_tot(2)
         fz = fz + ci*upfield_tot(3)
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
