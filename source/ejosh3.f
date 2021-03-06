c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ejosh3  --  atomic multipole energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ejosh3" calculates the electrostatic energy due to atomic
c     multipole interactions, and partitions the energy among atoms
c
c
      subroutine ejosh3
      use limits
      use mpole
      implicit none
c
c     save permanent electric field
c     
      savefield = .true.
c
c
c     choose the method for summing over multipole interactions
c
      if (use_ewald) then
         if (use_mlist) then
            call ejosh3d
         else
            call ejosh3c
         end if
      else
         if (use_mlist) then
            call ejosh3b
         else
            call ejosh3a
         end if
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ejosh3a  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ejosh3a" calculates the atomic multipole interaction energy
c     using a double loop, and partitions the energy among atoms
c
c
      subroutine ejosh3a
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use polgrp
      use polpot
      use potent
      use shunt
      use usage
      use chgpen
      use disp
      use pauli
      implicit none
      integer i,j,k
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      integer rorder
      real*8 e,f,fgrp
      real*8 ecc,ecv,evc,evv
      real*8 e_ele,e_disp,e_pauli
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 alphai,alphak
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 corei,vali
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 corek,valk
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 term1i,term2i,term3i
      real*8 term1k,term2k,term3k
      real*8 rr6
      real*8 c6i,c6k,c6ik
      real*8 displam
      real*8 pvali,pvalk
      real*8 overlapi,overlapk,oik
      real*8 apauli,apaulk
      real*8 fid(3),fkd(3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: lambdai(:)
      real*8, allocatable :: lambdak(:)
      real*8, allocatable :: lambdaik(:)
      logical proceed
      logical header,huge
      logical usei,usek
      character*6 mode
c
c
c     zero out total atomic multipole energy and partitioning
c
      nem = 0
      nedis = 0
      nepr = 0
      em = 0.0d0
      edis = 0.0d0
      epr = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aedis(i) = 0.0d0
         aepr(i) = 0.0d0
         do j = 1, 3
            permfield(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     set maximum power of 1/r for damping
c     (9 for quadrupole-quadrupole energy)
c
      rorder = 9
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (dscale(n))
      allocate (lambdai(rorder))
      allocate (lambdak(rorder))
      allocate (lambdaik(rorder))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     initialize charge penetration scale factors
c
      do i = 1, rorder
         lambdai(i) = 1.0d0
         lambdak(i) = 1.0d0
         lambdaik(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     calculate the multipole interaction energy term
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         corei = monopole(1,i)
         vali = monopole(2,i)
         alphai = alphaele(i)
         c6i = csix(i)
         overlapi = overpauli(i)
         apauli = alphapauli(i)
         pvali = monopauli(i)
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
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do k = i+1, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (proceed) then
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  corek = monopole(1,k)
                  valk = monopole(2,k)
                  alphak = alphaele(k)
                  c6k = csix(k)
                  overlapk = overpauli(k)
                  apaulk = alphapauli(k)
                  pvalk = monopauli(k)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
c                  rr1 = f * mscale(kk) / r
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
c
c     intermediates involving moments and distance separation
c
                  dri = dix*xr + diy*yr + diz*zr
                  drk = dkx*xr + dky*yr + dkz*zr
                  dik = dix*dkx + diy*dky + diz*dkz
                  qrix = qixx*xr + qixy*yr + qixz*zr
                  qriy = qixy*xr + qiyy*yr + qiyz*zr
                  qriz = qixz*xr + qiyz*yr + qizz*zr
                  qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qrky = qkxy*xr + qkyy*yr + qkyz*zr
                  qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qrri = qrix*xr + qriy*yr + qriz*zr
                  qrrk = qrkx*xr + qrky*yr + qrkz*zr
                  qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                  qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                  diqrk = dix*qrkx + diy*qrky + diz*qrkz
                  dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
c
c     calculate valence - valence interaction intermediate terms
c
c                  term1 = ci*ck
                  term1ik = vali*valk
                  term2ik = valk*dri - vali*drk + dik
                  term3ik = vali*qrrk + valk*qrri - dri*drk
     &                       + 2.0d0*(dkqri-diqrk+qik)
                  term4ik = dri*qrrk - drk*qrri - 4.0d0*qrrik
                  term5ik = qrri*qrrk
c
c     calculate core - valence interaction intermediate terms
c
                  term1i = corek*vali
                  term2i = corek*dri
                  term3i = corek*qrri
c
c     calculate valence - core interaction intermediate terms
c
                  term1k = corei*valk
                  term2k = -corei*drk
                  term3k = corei*qrrk
c
c     calculate core - core interaction intermediate terms
c
                  term1 = corei*corek
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     electrostatic interaction energy!
c
                  call damphlike(r,rorder,alphai,alphak,
     &                 lambdai,lambdak,lambdaik)
c
c     compute the valence - valence energy contribution for this interaction
c
                  evv = term1ik*rr1*lambdaik(1) + 
     &                  term2ik*rr3*lambdaik(3) + 
     &                  term3ik*rr5*lambdaik(5) +
     &                  term4ik*rr7*lambdaik(7) + 
     &                  term5ik*rr9*lambdaik(9)
c
c     compute the core - valence energy contribution for this interaction
c
                  ecv = term1i*rr1*lambdai(1) +
     &                  term2i*rr3*lambdai(3) +
     &                  term3i*rr5*lambdai(5)
c
                  evc = term1k*rr1*lambdak(1) +
     &                  term2k*rr3*lambdak(3) +
     &                  term3k*rr5*lambdak(5)
c
c     compute the core - core energy contribution for this interaction
c
                  ecc = term1*rr1
c
c     add together energy components for total damped energy
c
                  e_ele = f * mscale(kk) * (evv + ecv + evc + ecc)
c
                  if (use_group)  e_ele = e_ele * fgrp
c
c     save permanent electric field for induced dipole calculation
c
               fid(1) = -xr*(rr3*corek + rr3*lambdak(3)*valk -
     &              rr5*lambdak(5)*drk + rr7*lambdak(7)*qrrk)
     &              - rr3*lambdak(3)*dkx + 2.0d0*rr5*lambdak(5)*qrkx
               fid(2) = -yr*(rr3*corek + rr3*lambdak(3)*valk -
     &              rr5*lambdak(5)*drk+rr7*lambdak(7)*qrrk)
     &              - rr3*lambdak(3)*dky + 2.0d0*rr5*lambdak(5)*qrky
               fid(3) = -zr*(rr3*corek + rr3*lambdak(3)*valk -
     &              rr5*lambdak(5)*drk+rr7*lambdak(7)*qrrk)
     &              - rr3*lambdak(3)*dkz + 2.0d0*rr5*lambdak(5)*qrkz
               fkd(1) = xr*(rr3*corei + rr3*lambdai(3)*vali +
     &              rr5*lambdai(5)*dri + rr7*lambdai(7)*qrri)
     &              - rr3*lambdai(3)*dix - 2.0d0*rr5*lambdai(5)*qrix
               fkd(2) = yr*(rr3*corei + rr3*lambdai(3)*vali +
     &              rr5*lambdai(5)*dri + rr7*lambdai(7)*qrri)
     &              - rr3*lambdai(3)*diy - 2.0d0*rr5*lambdai(5)*qriy
               fkd(3) = zr*(rr3*corei + rr3*lambdai(3)*vali +
     &              rr5*lambdai(5)*dri + rr7*lambdai(7)*qrri)
     &              - rr3*lambdai(3)*diz - 2.0d0*rr5*lambdai(5)*qriz
c
c     increment electric field on both sites
c
               do j = 1, 3
                  permfield(j,i) = permfield(j,i) + fid(j)*dscale(kk)
                  permfield(j,k) = permfield(j,k) + fkd(j)*dscale(kk)
               end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     dispersion energy!
c
                  rr6 = rr3**2
c
c     c6 multiplicative combining rule
c
                  c6ik = c6i*c6k
c
c     dispersion damping factor
c
                  displam = 0.5d0*(3.0d0*lambdaik(5) - lambdaik(3))
c
c     compute damped 1/r^6 energy term
c
                  e_disp = -c6ik * rr6 * mscale(kk) * displam**2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     pauli repulsion energy!
c
                  call damppauli(r,r2,rr1,rr3,rr5,rr7,rr9,rr11,rorder,
     &                 apauli,apaulk,lambdaik)                  
c
c     recompute terms with number of pauli valence electrons
c
                  term1ik = pvali*pvalk
                  term2ik = pvalk*dri - pvali*drk + dik
                  term3ik = pvali*qrrk + pvalk*qrri - dri*drk
     &                       + 2.0d0*(dkqri-diqrk+qik)
c
c     compute valence - valence energy contribution for this interaction
c     (pauli repulsion has no terms involving the core)
c
                  evv = term1ik*lambdaik(1) +
     &                  term2ik*lambdaik(3) +
     &                  term3ik*lambdaik(5) +
     &                  term4ik*lambdaik(7) +
     &                  term5ik*lambdaik(9)
c
c     combining rule for pauli repulsion prefactor
c
                  oik = overlapi*overlapk
c
c     total pauli repulsion energy
c
                  e_pauli = oik * mscale(kk) * evv * rr1
c
c     print components
c
c                  if (debug .and. (e_pauli.ne.0.0d0)) then
c                     print *,"Repulsion Components for:",ii,name(ii),
c     &                    kk,name(kk),r,e_pauli
c                     print *,"Monopole-Monople",
c     &                    oik*mscale(kk)*rr1*term1ik*lambdaik(1)
c                     print *,"Order3",
c     &                    oik*mscale(kk)*rr1*term2ik*lambdaik(3)
c                     print *,"Order5",
c     &                    oik*mscale(kk)*rr1*term3ik*lambdaik(5)
c                     print *,"Charge-Quadrupole",
c     &                    oik*mscale(kk)*rr1*(pvali*qrrk + pvalk*qrri)*
c     &                    lambdaik(5)
cc                     if ((ii .eq.2).and.(kk .eq. 4)) then
c                        print *,"quaddollarsign",oik*mscale(kk)*rr1*
c     &                       (pvali*qrrk)*lambdaik(5)
cc                     end if
c                     print *,"Order7",
c     &                    oik*mscale(kk)*rr1*term4ik*lambdaik(7)
c                     print *,"Quadrupole-Quadrupole",
c     &                    oik*mscale(kk)*rr1*term5ik*lambdaik(9)
c                  end if
c
                  if (use_group)  e_pauli = e_pauli * fgrp
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     increment the overall multipole energy components
c
                  if (e_ele .ne. 0.0d0) then
                     nem = nem + 1
                     em = em + e_ele
                     aem(ii) = aem(ii) + 0.5d0*e_ele
                     aem(kk) = aem(kk) + 0.5d0*e_ele
                     if (molcule(ii) .ne. molcule(kk))
     &                  einter = einter + e_ele
                  end if
c
c     increment the overall dispersion energy components
c
                  if (e_disp .ne. 0.0d0) then
                     nedis = nedis + 1
                     edis = edis + e_disp
                     aedis(ii) = aedis(ii) + 0.5d0*e_disp
                     aedis(kk) = aedis(kk) + 0.5d0*e_disp
                     if (molcule(ii) .ne. molcule(kk))
     &                  einter = einter + e_disp
                  end if
c
c     increment the overall pauli repulsion energy components
c
                  if (e_pauli .ne. 0.0d0) then
                     nepr = nepr + 1
                     epr = epr + e_pauli
                     aepr(ii) = aepr(ii) + 0.5d0*e_pauli
                     aepr(kk) = aepr(kk) + 0.5d0*e_pauli
                     if (molcule(ii) .ne. molcule(kk))
     &                  einter = einter + e_pauli
                  end if
c
c     print message if the energy of multipole interaction is large
c
                  huge = (abs(e_ele) .gt. 100.0d0)
                  if ((debug.and.e_ele.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Atomic Multipole',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  ii,name(ii),kk,name(kk),r,e_ele
   30                format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                          f10.4,2x,f12.4)
                  end if
c
c     print message if the energy of the dispersion interaction is large
c     
                  huge = (abs(e_disp) .gt. 100.0d0)
                  if ((debug.and.e_disp.ne.0.0d0)
     &                 .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,21)
 21                     format (/,' Individual Dispersion',
     &                       ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                       15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,31)  ii,name(ii),kk,name(kk),r,e_disp
 31                  format (' Disp  ',4x,2(i7,'-',a3),9x,
     &                    f10.4,2x,f12.4)
                  end if
c
c     print message if the energy of the repulsion interaction is large 
c
                  huge = (abs(e_pauli) .gt. 100.0d0)
                  if ((debug.and.e_pauli.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,22)
 22                     format (/,' Individual Pauli Repulsion',
     &                       ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                       15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,32)  ii,name(ii),kk,name(kk),r,e_pauli
 32                  format (' Repulsion',1x,2(i7,'-',a3),9x,
     &                    f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction energy with other unit cells
c
         do i = 1, npole
            ii = ipole(i)
            iz = zaxis(i)
            ix = xaxis(i)
            iy = yaxis(k)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
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
            usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = m2scale
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = m3scale
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = m4scale
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = m5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do k = i, npole
               kk = ipole(k)
               kz = zaxis(k)
               kx = xaxis(k)
               ky = yaxis(k)
               usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
               if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
               proceed = .true.
               if (proceed)  proceed = (usei .or. usek)
               if (proceed) then
                  do j = 1, ncell
                     xr = x(kk) - xi
                     yr = y(kk) - yi
                     zr = z(kk) - zi
                     call imager (xr,yr,zr,j)
                     r2 = xr*xr + yr* yr + zr*zr
                     if (.not. (use_polymer .and. r2.le.polycut2))
     &                  mscale(kk) = 1.0d0
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        ck = rpole(1,k)
                        dkx = rpole(2,k)
                        dky = rpole(3,k)
                        dkz = rpole(4,k)
                        qkxx = rpole(5,k)
                        qkxy = rpole(6,k)
                        qkxz = rpole(7,k)
                        qkyy = rpole(9,k)
                        qkyz = rpole(10,k)
                        qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
                        rr1 = f / r
                        rr3 = rr1 / r2
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr9 = 7.0d0 * rr7 / r2
c
c     intermediates involving moments and distance separation
c
                        dri = dix*xr + diy*yr + diz*zr
                        drk = dkx*xr + dky*yr + dkz*zr
                        dik = dix*dkx + diy*dky + diz*dkz
                        qrix = qixx*xr + qixy*yr + qixz*zr
                        qriy = qixy*xr + qiyy*yr + qiyz*zr
                        qriz = qixz*xr + qiyz*yr + qizz*zr
                        qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                        qrky = qkxy*xr + qkyy*yr + qkyz*zr
                        qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                        qrri = qrix*xr + qriy*yr + qriz*zr
                        qrrk = qrkx*xr + qrky*yr + qrkz*zr
                        qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                        qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                           + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                        diqrk = dix*qrkx + diy*qrky + diz*qrkz
                        dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
                        term1 = ci*ck
                        term2 = ck*dri - ci*drk + dik
                        term3 = ci*qrrk + ck*qrri - dri*drk
     &                             + 2.0d0*(dkqri-diqrk+qik)
                        term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
                        term5 = qrri*qrrk
c
c     compute the energy contribution for this interaction
c
                        e = term1*rr1 + term2*rr3 + term3*rr5
     &                         + term4*rr7 + term5*rr9
                        e = e * mscale(kk)
                        if (use_group)  e = e * fgrp
                        if (ii .eq. kk)  e = 0.5d0 * e
c
c     increment the overall multipole energy components
c
                        if (e .ne. 0.0d0) then
                           nem = nem + 1
                           em = em + e
                           aem(ii) = aem(ii) + 0.5d0*e
                           aem(kk) = aem(kk) + 0.5d0*e
                           einter = einter + e
                        end if
c
c     print message if the energy of this interaction is large
c
                        huge = (abs(e) .gt. 100.0d0)
                        if ((debug.and.e.ne.0.0d0)
     &                        .or. (verbose.and.huge)) then
                           if (header) then
                              header = .false.
                              write (iout,40)
   40                         format (/,' Individual Atomic Multipole',
     &                                   ' Interactions :',
     &                                //,' Type',14x,'Atom Names',
     &                                   15x,'Distance',8x,'Energy',/)
                           end if
                           write (iout,50)  ii,name(ii),kk,name(kk),r,e
   50                      format (' M-Pole',4x,2(i7,'-',a3),1x,
     &                                '(XTAL)',2x,f10.4,2x,f12.4)
                        end if
                     end if
                  end do
               end if
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ejosh3b  --  neighbor list multipole analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ejosh3b" calculates the atomic multipole interaction energy
c     using a neighbor list, and partitions the energy among the atoms
c
c
      subroutine ejosh3b
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use potent
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8, allocatable :: mscale(:)
      logical proceed
      logical header,huge
      logical usei,usek
      character*6 mode
c
c
c     zero out total atomic multipole energy and partitioning
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,xaxis,yaxis,zaxis,rpole,use,
!$OMP& n12,i12,n13,i13,n14,i14,n15,i15,m2scale,m3scale,m4scale,
!$OMP& m5scale,nelst,elst,use_group,use_intra,use_bounds,off2,
!$OMP& f,molcule,name,verbose,debug,header,iout)
!$OMP& firstprivate(mscale) shared (em,einter,nem,aem)
!$OMP DO reduction(+:em,einter,nem,aem) schedule(guided)
c
c     calculate the multipole interaction energy term
c
      do i = 1, npole
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
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
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (proceed) then
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
                  rr1 = f * mscale(kk) / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
c
c     intermediates involving moments and distance separation
c
                  dri = dix*xr + diy*yr + diz*zr
                  drk = dkx*xr + dky*yr + dkz*zr
                  dik = dix*dkx + diy*dky + diz*dkz
                  qrix = qixx*xr + qixy*yr + qixz*zr
                  qriy = qixy*xr + qiyy*yr + qiyz*zr
                  qriz = qixz*xr + qiyz*yr + qizz*zr
                  qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qrky = qkxy*xr + qkyy*yr + qkyz*zr
                  qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qrri = qrix*xr + qriy*yr + qriz*zr
                  qrrk = qrkx*xr + qrky*yr + qrkz*zr
                  qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                  qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                  diqrk = dix*qrkx + diy*qrky + diz*qrkz
                  dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
                  term1 = ci*ck
                  term2 = ck*dri - ci*drk + dik
                  term3 = ci*qrrk + ck*qrri - dri*drk
     &                       + 2.0d0*(dkqri-diqrk+qik)
                  term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
                  term5 = qrri*qrrk
c
c     compute the energy contribution for this interaction
c
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
                  if (use_group)  e = e * fgrp
c
c     increment the overall multipole energy components
c
                  if (e .ne. 0.0d0) then
                     nem = nem + 1
                     em = em + e
                     aem(ii) = aem(ii) + 0.5d0*e
                     aem(kk) = aem(kk) + 0.5d0*e
                     if (molcule(ii) .ne. molcule(kk))
     &                  einter = einter + e
                  end if
c
c     print message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Atomic Multipole',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  ii,name(ii),kk,name(kk),r,e
   30                format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                          f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ejosh3c  --  Ewald multipole analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ejosh3c" calculates the atomic multipole interaction energy
c     using a particle mesh Ewald summation and double loop, and
c     partitions the energy among the atoms
c
c
      subroutine ejosh3c
      use sizes
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      implicit none
      integer i,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      em = 0.0d0
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
c     compute the real space part of the Ewald summation
c
      call emreal3c
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip
c
c     compute the self-energy part of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
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
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
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
         do i = 1, npole
            ii = ipole(i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         e = term * (xd*xd+yd*yd+zd*zd)
         em = em + e
         nem = nem + 1
         do i = 1, npole
            ii = ipole(i)
            aem(ii) = aem(ii) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine ejoshreal3c  --  real space mpole analysis via loop  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "ejoshreal3c" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions and partitions
c     the energy among the atoms
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine ejoshreal3c
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use energi
      use ewald
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 e,efull,f
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 bn(0:4)
      real*8, allocatable :: mscale(:)
      logical header,huge
      character*6 mode
      external erfc
c
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
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
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
               qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                  + qixx*qkxx + qiyy*qkyy + qizz*qkzz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0d0*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
               term5 = qrri*qrrk
c
c     compute the full undamped energy for this interaction
c
               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                    + term4*rr7 + term5*rr9
               efull = mscale(kk) * efull
               if (efull .ne. 0.0d0) then
                  nem = nem + 1
                  aem(ii) = aem(ii) + 0.5d0*efull
                  aem(kk) = aem(kk) + 0.5d0*efull
                  if (molcule(ii) .ne. molcule(kk))
     &               einter = einter + efull
               end if
c
c     modify error function terms to account for scaling
c
               scalekk = 1.0d0 - mscale(kk)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               em = em + e
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 100.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Atomic Multipole',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  ii,name(ii),kk,name(kk),r,efull
   30             format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction energy with other unit cells
c
         do i = 1, npole
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
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
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = m2scale
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = m3scale
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = m4scale
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = m5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do k = i, npole
               kk = ipole(k)
               do jcell = 1, ncell
                  xr = x(kk) - xi
                  yr = y(kk) - yi
                  zr = z(kk) - zi
                  call imager (xr,yr,zr,jcell)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (.not. (use_polymer .and. r2.le.polycut2))
     &               mscale(kk) = 1.0d0
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     ck = rpole(1,k)
                     dkx = rpole(2,k)
                     dky = rpole(3,k)
                     dkz = rpole(4,k)
                     qkxx = rpole(5,k)
                     qkxy = rpole(6,k)
                     qkxz = rpole(7,k)
                     qkyy = rpole(9,k)
                     qkyz = rpole(10,k)
                     qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
                     rr1 = f / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
                     rr9 = 7.0d0 * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) / r
                     alsq2 = 2.0d0 * aewald**2
                     alsq2n = 0.0d0
                     if (aewald .gt. 0.0d0)
     &                  alsq2n = 1.0d0 / (sqrtpi*aewald)
                     exp2a = exp(-ralpha**2)
                     do j = 1, 4
                        bfac = dble(j+j-1)
                        alsq2n = alsq2 * alsq2n
                        bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                     end do
                     do j = 0, 4
                        bn(j) = f * bn(j)
                     end do
c
c     intermediates involving moments and distance separation
c
                     dri = dix*xr + diy*yr + diz*zr
                     drk = dkx*xr + dky*yr + dkz*zr
                     dik = dix*dkx + diy*dky + diz*dkz
                     qrix = qixx*xr + qixy*yr + qixz*zr
                     qriy = qixy*xr + qiyy*yr + qiyz*zr
                     qriz = qixz*xr + qiyz*yr + qizz*zr
                     qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qrky = qkxy*xr + qkyy*yr + qkyz*zr
                     qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qrri = qrix*xr + qriy*yr + qriz*zr
                     qrrk = qrkx*xr + qrky*yr + qrkz*zr
                     qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                     qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                     diqrk = dix*qrkx + diy*qrky + diz*qrkz
                     dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
                     term1 = ci*ck
                     term2 = ck*dri - ci*drk + dik
                     term3 = ci*qrrk + ck*qrri - dri*drk
     &                          + 2.0d0*(dkqri-diqrk+qik)
                     term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
                     term5 = qrri*qrrk
c
c     compute the full undamped energy for this interaction
c
                     efull = term1*rr1 + term2*rr3 + term3*rr5
     &                          + term4*rr7 + term5*rr9
                     efull = mscale(kk) * efull
                     if (ii .eq. kk)  efull = 0.5d0 * efull
                     if (efull .ne. 0.0d0) then
                        nem = nem + 1
                        aem(ii) = aem(ii) + 0.5d0*efull
                        aem(kk) = aem(kk) + 0.5d0*efull
                        einter = einter + efull
                     end if
c
c     modify distances to account for Ewald and exclusions
c
                     scalekk = 1.0d0 - mscale(kk)
                     rr1 = bn(0) - scalekk*rr1
                     rr3 = bn(1) - scalekk*rr3
                     rr5 = bn(2) - scalekk*rr5
                     rr7 = bn(3) - scalekk*rr7
                     rr9 = bn(4) - scalekk*rr9
c
c     compute the energy contribution for this interaction
c
                     e = term1*rr1 + term2*rr3 + term3*rr5
     &                      + term4*rr7 + term5*rr9
                     if (ii .eq. kk)  e = 0.5d0 * e
                     em = em + e
c
c     print message if the energy of this interaction is large
c
                     huge = (abs(efull) .gt. 100.0d0)
                     if ((debug.and.efull.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Atomic Multipole',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',8x,'Energy',/)
                        end if
                        write (iout,50)  ii,name(ii),kk,name(kk),r,efull
   50                   format (' M-Pole',4x,2(i7,'-',a3),1x,
     &                             '(XTAL)',2x,f10.4,2x,f12.4)
                     end if
                  end if
               end do
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ejosh3d  --  Ewald multipole analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ejosh3d" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a neighbor list, and
c     partitions the energy among the atoms
c
c
      subroutine ejosh3d
      use sizes
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use disp
      implicit none
      integer i,ii,j
      real*8 e,f
      real*8 pre
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      nedis = 0
      nepr = 0
      em = 0.0d0
      edis = 0.0d0
      epr = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aedis(i) = 0.0d0
         aepr(i) = 0.0d0
         do j = 1, 3
            permfield(j,i) = 0.0d0
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
c     compute the real space part of the Ewald summation
c
      call ejoshreal3d
c
c     compute the reciprocal space part of multipole Ewald summation
c
      call emrecip
c
c     compute the reciprocal space part of dispersion Ewald summation 
c
      call edisprecip
c
c     compute the self-energy part of the multipole Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
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
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
         nem = nem + 1
         aem(i) = aem(i) + e
      end do
c
c     compute the self term for the dispersion Ewald summation
c
      do i = 1, npole
         pre = adewald**6/12.0d0
         edis = edis + pre*csix(i)*csix(i)
      end do
c
c     compute self energy portion of electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
            permfield(j,i) = permfield(j,i) + term*rpole(j+1,i)
         end do
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
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         e = term * (xd*xd+yd*yd+zd*zd)
         em = em + e
         nem = nem + 1
         do i = 1, npole
            ii = ipole(i)
            aem(ii) = aem(ii) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine ejoshreal3d  --  real space mpole analysis via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "ejoshreal3d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions, and partitions
c     the energy among the atoms using a pairwise neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine ejoshreal3d
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpot
      use couple
      use energi
      use ewald
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use polgrp
      use polpot
      use potent
      use shunt
      use chgpen
      use disp
      use pauli
      use tarray
      use openmp
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      integer nlocal,nchunk
      integer tid,maxlocal
!$    integer omp_get_thread_num
      integer, allocatable :: toffset(:)
      integer, allocatable :: ilocal(:,:)
      real*8 e,f
      real*8 ecc,ecv,evc,evv
      real*8 e_ele,e_disp,e_pauli
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 alphai,alphak
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 rr1core,rr3core
      real*8 rr1i,rr3i,rr5i
      real*8 rr7i,rr9i
      real*8 rr1k,rr3k,rr5k
      real*8 rr7k,rr9k
      real*8 rr1ik,rr3ik,rr5ik
      real*8 rr7ik,rr9ik
      real*8 corei,vali
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 corek,valk
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 term1i,term2i,term3i
      real*8 term1k,term2k,term3k
      real*8 rr6
      real*8 c6i,c6k,c6ik
      real*8 displam
      real*8 damp,term,expterm
      real*8 ralpha2
      real*8 pvali,pvalk
      real*8 overlapi,overlapk,oik
      real*8 apauli,apaulk
      real*8 fid(3),fkd(3)
      real*8 bn(0:4)
      real*8 lambdai(9),lambdak(9),lambdaik(9)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: muscale(:)
      real*8, allocatable :: dlocal(:,:)
      logical header,huge
      character*6 mode
      external erfc
c
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c
c     values for storage of mutual polarization intermediates
c
      nchunk = int(0.5d0*dble(npole)/dble(nthread)) + 1
      maxlocal = int(dble(npole)*dble(maxelst)/dble(nthread))
      nlocal = 0
      ntpair = 0
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (dscale(n))
      allocate (muscale(n))
      allocate (toffset(0:nthread-1))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
         dscale(i) = 1.0d0
         muscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,monopole,alphaele,csix,overpauli,
!$OMP& alphapauli,monopauli,rpole,n12,i12,n13,i13,n14,i14,n15,
!$OMP& np11,np12,np13,np14,ip11,ip12,ip13,ip14, 
!$OMP& i15,m2scale,m3scale,m4scale,m5scale,
!$OMP& d1scale,d2scale,d3scale,d4scale,
!$OMP& mu2scale,mu3scale,mu4scale,mu5scale,atomic,
!$OMP& nelst,elst,use_bounds,
!$OMP& f,off2,aewald,molcule,name,verbose,debug,header,iout,
!$OMP& adewald,ntpair,tindex,                         
!$OMP& tdipdip,toffset,maxlocal,maxelst,                     
!$OMP& nthread,nchunk) 
!$OMP& firstprivate(mscale,dscale,muscale,nlocal) shared (em,edis,epr,
!$OMP& einter,nem,aem,permfield)
c
c     perform dynamic allocation of some local arrays
c
      allocate (ilocal(2,maxlocal))
      allocate (dlocal(6,maxlocal))
c
!$OMP DO reduction(+:em,edis,epr,einter,nem,aem,
!$OMP& permfield) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         corei = monopole(1,i)
         vali = monopole(2,i)
         alphai = alphaele(i)
         c6i = csix(i)
         overlapi = overpauli(i)
         apauli = alphapauli(i)
         pvali = monopauli(i)
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
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            muscale(i12(j,ii)) = mu2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            muscale(i13(j,ii)) = mu3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            muscale(i14(j,ii)) = mu4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            muscale(i15(j,ii)) = mu5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               corek = monopole(1,k)
               valk = monopole(2,k)
               alphak = alphaele(k)
               c6k = csix(k)
               overlapk = overpauli(k)
               apaulk = alphapauli(k)
               pvalk = monopauli(k)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c               do j = 0, 4
c                  bn(j) = f * bn(j)
c               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
               qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                  + qixx*qkxx + qiyy*qkyy + qizz*qkzz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c
c     calculate intermediate terms for multipole interaction
c
c
c     calculate valence - valence interaction intermediate terms
c
c                  term1 = ci*ck
               term1ik = vali*valk
               term2ik = valk*dri - vali*drk + dik
               term3ik = vali*qrrk + valk*qrri - dri*drk
     &              + 2.0d0*(dkqri-diqrk+qik)
               term4ik = dri*qrrk - drk*qrri - 4.0d0*qrrik
               term5ik = qrri*qrrk
c     
c     calculate core - valence interaction intermediate terms
c
               term1i = corek*vali
               term2i = corek*dri
               term3i = corek*qrri
c
c     calculate valence - core interaction intermediate terms
c
               term1k = corei*valk
               term2k = -corei*drk
               term3k = corei*qrrk
c
c     calculate core - core interaction intermediate terms
c
               term1 = corei*corek
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     electrostatic interaction energy!
c
               call damphlike(r,9,alphai,alphak,
     &              lambdai,lambdak,lambdaik)
c
c     compute the valence - valence energy contribution for this interaction
c
               evv = term1ik*rr1*lambdaik(1) + 
     &              term2ik*rr3*lambdaik(3) + 
     &              term3ik*rr5*lambdaik(5) +
     &              term4ik*rr7*lambdaik(7) + 
     &              term5ik*rr9*lambdaik(9)
c
c     compute the core - valence energy contribution for this interaction
c
               ecv = term1i*rr1*lambdai(1) +
     &              term2i*rr3*lambdai(3) +
     &              term3i*rr5*lambdai(5)
c     
               evc = term1k*rr1*lambdak(1) +
     &              term2k*rr3*lambdak(3) +
     &              term3k*rr5*lambdak(5)
c
c     compute the core - core energy contribution for this interaction
c
               ecc = term1*rr1
c
c     add together energy components for total damped energy
c
               e_ele = f * mscale(kk) * (evv + ecv + evc + ecc)
c
c     compute the full undamped energy for this interaction
c
c               efull = term1*rr1 + term2*rr3 + term3*rr5
c     &                    + term4*rr7 + term5*rr9
c               efull = mscale(kk) * efull
c               print *,"e_ele",i,k,e_ele
               if (e_ele .ne. 0.0d0) then
                  nem = nem + 1
                  aem(ii) = aem(ii) + 0.5d0*e_ele
                  aem(kk) = aem(kk) + 0.5d0*e_ele
                  if (molcule(ii) .ne. molcule(kk))
     &               einter = einter + e_ele
               end if
c
c     modify error function terms to account for scaling
c
c     ewald - naked + damped
c
c               scalekk = 1.0d0 - mscale(kk)
               rr1core = bn(0) - (1.0d0 - mscale(kk))*rr1
               rr3core = bn(1) - (1.0d0 - mscale(kk))*rr3
c
               rr1i = bn(0) - (1.0d0 - mscale(kk)*lambdai(1))*rr1
               rr3i = bn(1) - (1.0d0 - mscale(kk)*lambdai(3))*rr3
               rr5i = bn(2) - (1.0d0 - mscale(kk)*lambdai(5))*rr5
               rr7i = bn(3) - (1.0d0 - mscale(kk)*lambdai(7))*rr7
               rr9i = bn(4) - (1.0d0 - mscale(kk)*lambdai(9))*rr9
c     
               rr1k = bn(0) - (1.0d0 - mscale(kk)*lambdak(1))*rr1
               rr3k = bn(1) - (1.0d0 - mscale(kk)*lambdak(3))*rr3
               rr5k = bn(2) - (1.0d0 - mscale(kk)*lambdak(5))*rr5
               rr7k = bn(3) - (1.0d0 - mscale(kk)*lambdak(7))*rr7
               rr9k = bn(4) - (1.0d0 - mscale(kk)*lambdak(9))*rr9
c
               rr1ik = bn(0) - (1.0d0 - mscale(kk)*lambdaik(1))*rr1
               rr3ik = bn(1) - (1.0d0 - mscale(kk)*lambdaik(3))*rr3
               rr5ik = bn(2) - (1.0d0 - mscale(kk)*lambdaik(5))*rr5
               rr7ik = bn(3) - (1.0d0 - mscale(kk)*lambdaik(7))*rr7
               rr9ik = bn(4) - (1.0d0 - mscale(kk)*lambdaik(9))*rr9
c
c     compute the valence - valence energy contribution for this interaction
c
               evv = term1ik*rr1ik + 
     &              term2ik*rr3ik + 
     &              term3ik*rr5ik +
     &              term4ik*rr7ik + 
     &              term5ik*rr9ik
c     
c     compute the core - valence energy contribution for this interaction
c     
               ecv = term1i*rr1i +
     &              term2i*rr3i +
     &              term3i*rr5i
c     
               evc = term1k*rr1k +
     &              term2k*rr3k +
     &              term3k*rr5k
c     
c     compute the core - core energy contribution for this interaction
c     
               ecc = term1*rr1core
c
c     compute the energy contribution for this interaction
c
               e = evv + ecv + evc + ecc
               em = em + f*e
c
c     save permanent electric field for induced dipole calculation
c
               rr1core = bn(0) - (1.0d0 - dscale(kk))*rr1
               rr3core = bn(1) - (1.0d0 - dscale(kk))*rr3
c
               rr1i = bn(0) - (1.0d0 - dscale(kk)*lambdai(1))*rr1
               rr3i = bn(1) - (1.0d0 - dscale(kk)*lambdai(3))*rr3
               rr5i = bn(2) - (1.0d0 - dscale(kk)*lambdai(5))*rr5
               rr7i = bn(3) - (1.0d0 - dscale(kk)*lambdai(7))*rr7
c
               rr1k = bn(0) - (1.0d0 - dscale(kk)*lambdak(1))*rr1
               rr3k = bn(1) - (1.0d0 - dscale(kk)*lambdak(3))*rr3
               rr5k = bn(2) - (1.0d0 - dscale(kk)*lambdak(5))*rr5
               rr7k = bn(3) - (1.0d0 - dscale(kk)*lambdak(7))*rr7
c
               fid(1) = -xr*(rr3core*corek + rr3k*valk -
     &              rr5k*drk + rr7k*qrrk)
     &              - rr3k*dkx + 2.0d0*rr5k*qrkx
               fid(2) = -yr*(rr3core*corek + rr3k*valk -
     &              rr5k*drk+rr7k*qrrk)
     &              - rr3k*dky + 2.0d0*rr5k*qrky
               fid(3) = -zr*(rr3core*corek + rr3k*valk -
     &              rr5k*drk+rr7k*qrrk)
     &              - rr3k*dkz + 2.0d0*rr5k*qrkz
               fkd(1) = xr*(rr3core*corei + rr3i*vali +
     &              rr5i*dri + rr7i*qrri)
     &              - rr3i*dix - 2.0d0*rr5i*qrix
               fkd(2) = yr*(rr3core*corei + rr3i*vali +
     &              rr5i*dri + rr7i*qrri)
     &              - rr3i*diy - 2.0d0*rr5i*qriy
               fkd(3) = zr*(rr3core*corei + rr3i*vali +
     &              rr5i*dri + rr7i*qrri)
     &              - rr3i*diz - 2.0d0*rr5i*qriz
c
c     increment electric field on both sites
c
               do j = 1, 3
                  permfield(j,i) = permfield(j,i) + fid(j)
                  permfield(j,k) = permfield(j,k) + fkd(j)
               end do
c
c     save dipole - dipole t matrix for mutual induction
c
c     INSERT MUTUAL EXCLUSION RULES HERE!!!!
c
               rr3ik = bn(1) - (1.0d0 - muscale(kk)*lambdaik(3))*rr3
               rr5ik = bn(2) - (1.0d0 - muscale(kk)*lambdaik(5))*rr5
               nlocal = nlocal + 1
               ilocal(1,nlocal) = i
               ilocal(2,nlocal) = k
               dlocal(1,nlocal) = -rr3ik + rr5ik*xr*xr
               dlocal(2,nlocal) = rr5ik*xr*yr
               dlocal(3,nlocal) = rr5ik*xr*zr
               dlocal(4,nlocal) = -rr3ik + rr5ik*yr*yr
               dlocal(5,nlocal) = rr5ik*yr*zr
               dlocal(6,nlocal) = -rr3ik + rr5ik*zr*zr
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     dispersion energy!
c
               rr6 = rr3**2
c
c     c6 multiplicative combining rule
c
               c6ik = c6i*c6k
c
c     dispersion damping factor
c
               displam = 0.5d0*(3.0d0*lambdaik(5) - lambdaik(3))
c
c     dispersion ewald damping
c
               ralpha2 = r2 * adewald**2
               damp = 1.0d0
               if (ralpha2 .lt. 50.0d0) then
                  expterm = exp(-ralpha2)
                  term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                  damp = term*expterm
               end if
c               print *,i,k,damp
c
c     compute damped 1/r^6 energy term
c
               e_disp = -c6ik*rr6*(damp +
     &              (mscale(kk)*displam**2 - 1.0d0))
c
c
c     accumulate dispersion energy
c
               edis = edis + e_disp
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     pauli repulsion energy!
c
                  call damppauli(r,r2,rr1,rr3,rr5,rr7,rr9,rr11,9,
     &                 apauli,apaulk,lambdaik)                  
c
c     recompute terms with number of pauli valence electrons
c
                  term1ik = pvali*pvalk
                  term2ik = pvalk*dri - pvali*drk + dik
                  term3ik = pvali*qrrk + pvalk*qrri - dri*drk
     &                       + 2.0d0*(dkqri-diqrk+qik)
c
c     compute valence - valence energy contribution for this interaction
c     (pauli repulsion has no terms involving the core)
c
                  evv = term1ik*lambdaik(1) +
     &                  term2ik*lambdaik(3) +
     &                  term3ik*lambdaik(5) +
     &                  term4ik*lambdaik(7) +
     &                  term5ik*lambdaik(9)
c
c     combining rule for pauli repulsion prefactor
c
                  oik = overlapi*overlapk
c
c     total pauli repulsion energy
c
                  e_pauli = oik * mscale(kk) * evv * rr1
c
c     accumulate pauli repulsion energy
c
                  epr = epr + e_pauli







c
c
c     print a message if the energy of this interaction is large
c
c               huge = (abs(e_ele) .gt. 100.0d0)
c               if ((debug.and.e_ele.ne.0.0d0)
c     &               .or. (verbose.and.huge)) then
c                  if (header) then
c                     header = .false.
c                     write (iout,20)
c   20                format (/,' Individual Atomic Multipole',
c     &                          ' Interactions :',
c     &                       //,' Type',14x,'Atom Names',
c     &                          15x,'Distance',8x,'Energy',/)
c                  end if
c                  write (iout,30)  ii,name(ii),kk,name(kk),r,e_ele
c   30             format (' M-Pole',4x,2(i7,'-',a3),9x,
c     &                       f10.4,2x,f12.4)
c               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            muscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            muscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            muscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            muscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
c
c     find offset into global arrays for the current thread
c
!$OMP CRITICAL
      tid = 0
!$    tid = omp_get_thread_num ()
      toffset(tid) = ntpair
      ntpair = ntpair + nlocal
!$OMP END CRITICAL
c
c     store terms used later for mutual polarization
c
      k = toffset(tid)
      do i = 1, nlocal
         m = k + i
         tindex(1,m) = ilocal(1,i)
         tindex(2,m) = ilocal(2,i)
         do j = 1, 6
            tdipdip(j,m) = dlocal(j,i)
         end do
      end do
      deallocate (ilocal)
      deallocate (dlocal)
c
c
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (muscale)
      return
      end
