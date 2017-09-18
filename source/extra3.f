c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine extra3  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "extra3" calculates any additional user defined potential
c     contribution and also partitions the energy among the atoms
c
c
      subroutine extra3
      use sizes
      use action
      use analyz
      use atoms
      use energi
c
      use atomid
      use inform
      use math
      use molcul
      use mpole
      use couple
      use chgpen
      use vdw
      use potent
      use chgpot
      use inter
c      use vdwpot
      implicit none
      integer i,k,j,m
      integer ii,kk
      integer rorder
      real*8 r,r2
      real*8 xr,yr,zr
      real*8 ci,qi,zi
      real*8 ck,qk,zk
      real*8 alphai,alphak,alphamix
      real*8 alphai2,alphak2
      real*8 alpham,dampm,expdampm
      real*8 oldalpha3,oldalpha4
      real*8 invalphai,invalphak
      real*8 dampi,dampk,damp
      real*8 dampi2,dampi3,dampi4
      real*8 dampk2,dampk3,dampk4
      real*8 expdampi,expdampk,expdamp
      real*8 expdampi2,expdampk2
      real*8 overlapi,overlapk
      real*8 boverlapi,boverlapk
      real*8 coverlapi,coverlapk
      real*8 overlapri,overlaprk
      real*8 sigma
      real*8 oik,boik,coik,roik,oik2,term,pre
      real*8 e,e2,es2r
ccccccccccccccccc
      real*8 ee,en
      real*8 eundamp,edamp
      real*8 eundamp_ele,edamp_ele
      real*8 eover,eoverexp,eover2
      real*8 f,fm
      real*8 termi,termk
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 scale1i,scale1k,scale1ik
      real*8 scale3i,scale3k,scale3ik
      real*8 scale5i,scale5k,scale5ik
      real*8 scale7ik
      real*8 scale9ik
      real*8 oscale1ik,oscale3ik,oscale5ik
      real*8 oscale7ik,oscale9ik
      real*8 expscale1ik,expscale3ik,expscale5ik
      real*8 expscale7ik,expscale9ik
      real*8 dip_mag,quad_mag
c
      real*8 s,ds,dds,ddds,dddds
      real*8 sc(10)
      real*8 gl(0:4),gln(0:4)
      real*8 iv,rdn
      real*8 diff
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      character*20 list(6)
      logical proceed
c
c
c     zero out energy and partitioning due to extra potential terms
c
c      nex = 0
      ex = 0.0d0
c      do i = 1, n
c         aex(i) = 0.0d0
c      end do
      f = electric / dielec
      if (.not. use_mpole) then
         call kmpole
         call chkpole
         call rotpole
      end if
c      if (.not. use_vdw) then 
c         call kvdw
c      end if
c
c     i do not want any hydrogen reductions
c
c      allocate (xred(n))
c      allocate (yred(n))
c      allocate (zred(n))
c      do k = 1, nvdw
c         i = ivdw(k)
c         iv = ired(i)
c         rdn = kred(i)
c         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
c         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
c         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
c      end do
c
c     add any user-defined extra potentials and partitioning
c
      do i = 1, npole-1
         ii = ivdw(i)
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
         zi = dble(atomic(i))
         if (exch_num_ele .eq. "VALENCE") then
            if (atomic(i) .gt. 2)  zi = zi - 2.0d0
            if (atomic(i) .gt. 10)  zi = zi - 8.0d0
            if (atomic(i) .gt. 18)  zi = zi - 8.0d0
            if (atomic(i) .gt. 20)  zi = zi - 10.0d0
         else if (exch_num_ele .eq. "TWO") then
            if (atomic(i) .ne. 1) zi = ci + 2.0d0
         else if (exch_num_ele .eq. "TWOPLUS") then
            if (atomic(i) .ne. 1) zi = 2.0d0
         else if (exch_num_ele .eq. "TWO-VARIABLE") then
            if (atomic(i) .ne. 1) zi = ci + 2.0d0
         else if (exch_num_ele .eq. "TWO-VARIABLE-H") then
            if (atomic(i) .ne. 1) zi = ci + 2.0d0
         else if (exch_num_ele .eq. "ZERO") then
            zi = ci
         else if (exch_num_ele .eq. "VARIABLE") then
            if (use_vdwclass) then
               zi = ci + exch_val_ele(vdwclass(i))
            else
               zi = ci + exch_val_ele(cpclass(i))
            end if
         else if (exch_num_ele .eq. "VARIABLE-H") then
            if (use_vdwclass) then
               if (atomic(i) .ne. 1) zi = ci + exch_val_ele(vdwclass(i))
            else
               if (atomic(i) .ne. 1) zi = ci + exch_val_ele(cpclass(i))
            end if
         end if
         qi = ci - zi
c
         do k = i+1, npole
            kk = ivdw(k)
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
            proceed = .true.
            if (molcule(kk) .eq. molcule(ii))  proceed = .false.
            if (proceed) then
c
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
c     
               r2 = xr*xr + yr*yr + zr*zr
               r = sqrt(r2)
               if (use_vdwclass) then
                  overlapi = overlap(vdwclass(ii))
                  alphai = exch_alpha(vdwclass(ii))
                  overlapk = overlap(vdwclass(kk))
                  alphak = exch_alpha(vdwclass(kk))
               else
                  overlapi = overlap(cpclass(ii))
                  alphai = exch_alpha(cpclass(ii))
                  overlapk = overlap(cpclass(kk))
                  alphak = exch_alpha(cpclass(kk))
               end if
               oik = overlapi * overlapk
               if (exmix .eq. 'ARITHMETIC') then
                  oik = 0.5d0 * (overlapi + overlapk)
               else if (exmix .eq. 'GEOMETRIC') then
                  oik = sqrt(overlapi) * sqrt(overlapk)
               else if (exmix .eq. 'HARMONIC') then
                  oik = 2.0d0 * (overlapi*overlapk)/(overlapi+overlapk)
               else if (exmix .eq. 'HHG') then
                  oik = 4.0d0 * (overlapi*overlapk) / 
     &                 (sqrt(overlapi) + sqrt(overlapk))**2
               else if (exmix .eq. "W-H") then
                  oik2 = overlapi * overlapk
                  oik = sqrt(oik2)
               else if (exmix .eq."W-H1") then
                  termi = 2.0d0*(alphai**3)*(alphak**3)
                  termk = (alphai**6) + (alphak**6)
                  oik = sqrt(overlapi*overlapk)*termi/termk
c     bohm ahlrichs (reference?)
               else if (exmix .eq. "BA") then
                  alphamix = 2.0d0 * alphai * alphak / (alphai + alphak)
                  termi = overlapi**(1.0d0/alphai)
                  termk = overlapk**(1.0d0/alphak)
                  oik = (termi*termk)**(alphamix/2.0d0)
               else if (exmix .eq. "EPS") then
                  oik = epsilon(k,i)
                  print *,"eps oik",oik
               end if
c     
               zk = dble(atomic(k))
               if (exch_num_ele .eq. "VALENCE") then
                  if (atomic(k) .gt. 2)  zk = zk - 2.0d0
                  if (atomic(k) .gt. 10)  zk = zk - 8.0d0
                  if (atomic(k) .gt. 18)  zk = zk - 8.0d0
                  if (atomic(k) .gt. 20)  zk = zk - 10.0d0
               else if (exch_num_ele .eq. "TWO") then
                  if (atomic(k) .ne. 1) zk = ck + 2.0d0
               else if (exch_num_ele .eq. "TWOPLUS") then
                  if (atomic(k) .ne. 1) zk = 2.0d0
               else if (exch_num_ele .eq. "TWO-VARIABLE") then
                  if (atomic(k) .ne. 1) zk = ck + 2.0d0
               else if (exch_num_ele .eq. "TWO-VARIABLE-H") then
                  if (atomic(k) .ne. 1) zk = ck + 2.0d0
               else if (exch_num_ele .eq. "ZERO") then
                  zk = ck
               else if (exch_num_ele .eq. "VARIABLE") then
                  if (use_vdwclass) then
                     zk = ck + exch_val_ele(vdwclass(k))
                  else
                     zk = ck + exch_val_ele(cpclass(k))
                  end if
               else if (exch_num_ele .eq. "VARIABLE-H") then
                  if (use_vdwclass) then
                     if (atomic(k) .ne. 1) zk = ck +
     &                    exch_val_ele(vdwclass(k))
                  else
                     if (atomic(k) .ne. 1) zk = ck + 
     &                    exch_val_ele(cpclass(k))
                  end if
               end if
               qk = ck - zk
c
c     
c     construct some intermediate quadrupole values
c
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
c     
c     calculate the scalar products for permanent multipoles
c     
               sc(2) = dix*dkx + diy*dky + diz*dkz
               sc(3) = dix*xr + diy*yr + diz*zr
               sc(4) = dkx*xr + dky*yr + dkz*zr
               sc(5) = qix*xr + qiy*yr + qiz*zr
               sc(6) = qkx*xr + qky*yr + qkz*zr
               sc(7) = qix*dkx + qiy*dky + qiz*dkz
               sc(8) = qkx*dix + qky*diy + qkz*diz
               sc(9) = qix*qkx + qiy*qky + qiz*qkz
               sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &              + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c     
c     calculate the gl functions for permanent multipoles - electrons
c     
               gl(0) = qi*qk
               gl(1) = qk*sc(3) - qi*sc(4) + sc(2)
               gl(2) = qi*sc(6) + qk*sc(5) - sc(3)*sc(4)
     &              + 2.0d0*(sc(7)-sc(8)+sc(10))
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
               gl(4) = sc(5)*sc(6)
c     
c     calculate the gl functions for permanent multipoles - nuclei
c     
               gln(0) = zi*zk + zi*qk + zk*qi
               gln(1) = zk*sc(3) - zi*sc(4)
               gln(2) = zi*sc(6) + zk*sc(5)
c     
c     compute the energy contributions for this interaction
c     
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     get exchange parameters
c     
               if (exmodel .eq. "ORBITAL2ONLYR") then
                  alphai2 = 0.5d0*alphai
                  alphak2 = 0.5d0*alphak
                  dampi = alphai2*r
                  dampk = alphak2*r
                  expdampi = exp(-dampi)
                  expdampk = exp(-dampk)
                  diff = abs(alphai - alphak)
                  if (diff .ge. 0.001d0) then
c                  if (alphai2 .ne. alphak2) then
                     term = alphai2**2 - alphak2**2
                     pre = 64.0d0*(alphai**3)*(alphak**3)/(term**4)
c
c    new
c
      s = (dampi - 0.4D1 / term * alphak2 * alphai2) * expdampk + (dampk
     & + 0.4D1 / term * alphak2 * alphai2) * expdampi
c
      ds = (alphai2 * alphak2 * r2 - 0.4D1 * alphai2 * alphak2 ** 2 / te
     &rm * r - 0.4D1 / term * alphak2 * alphai2) * expdampk + (alphai2 *
     & alphak2 * r2 + 0.4D1 / term * alphak2 * alphai2 + 0.4D1 * alphai2
     & ** 2 * alphak2 / term * r) * expdampi
c
      dds = (alphai2 * alphak2 * r2 / 0.3D1 - 0.4D1 * alphai2 * alphak2 
     &** 2 / term * r + alphai2 * alphak2 ** 2 * r ** 3 / 0.3D1 - 0.4D1 
     &/ 0.3D1 * alphai2 * alphak2 ** 3 / term * r2 - 0.4D1 / term * alph
     &ak2 * alphai2) * expdampk + (alphai2 * alphak2 * r2 / 0.3D1 + 0.4D
     &1 / term * alphak2 * alphai2 + alphai2 ** 2 * alphak2 * r ** 3 / 0
     &.3D1 + 0.4D1 * alphai2 ** 2 * alphak2 / term * r + 0.4D1 / 0.3D1 *
     & alphai2 ** 3 * alphak2 / term * r2) * expdampi

c
      ddds = (-0.4D1 / term * alphak2 * alphai2 - 0.4D1 * alphai2 * alph
     &ak2 ** 2 / term * r - 0.8D1 / 0.5D1 * alphai2 * alphak2 ** 3 / ter
     &m * r2 + alphai2 * alphak2 * r2 / 0.5D1 + alphai2 * alphak2 ** 2 *
     & r ** 3 / 0.5D1 + alphai2 * alphak2 ** 3 * r ** 4 / 0.15D2 - 0.4D1
     & / 0.15D2 * alphai2 * alphak2 ** 4 / term * r ** 3) * expdampk + (
     &alphai2 ** 3 * alphak2 * r ** 4 / 0.15D2 + 0.4D1 / 0.15D2 * alphai
     &2 ** 4 * alphak2 / term * r ** 3 + 0.4D1 / term * alphak2 * alphai
     &2 + 0.4D1 * alphai2 ** 2 * alphak2 / term * r + 0.8D1 / 0.5D1 * al
     &phai2 ** 3 * alphak2 / term * r2 + alphai2 * alphak2 * r2 / 0.5D1 
     &+ alphai2 ** 2 * alphak2 * r ** 3 / 0.5D1) * expdampi
c
      dddds = (-0.12D2 / 0.7D1 * alphai2 * alphak2 ** 3 / term * r2 - 0.
     &8D1 / 0.21D2 * alphai2 * alphak2 ** 4 / term * r ** 3 + alphai2 * 
     &alphak2 * r2 / 0.7D1 + alphai2 * alphak2 ** 2 * r ** 3 / 0.7D1 - 0
     &.4D1 / term * alphak2 * alphai2 + 0.2D1 / 0.35D2 * alphai2 * alpha
     &k2 ** 3 * r ** 4 - 0.4D1 / 0.105D3 * alphai2 * alphak2 ** 5 / term
     & * r ** 4 + alphai2 * alphak2 ** 4 * r ** 5 / 0.105D3 - 0.4D1 * al
     &phai2 * alphak2 ** 2 / term * r) * expdampk + (alphai2 ** 2 * alph
     &ak2 * r ** 3 / 0.7D1 + 0.4D1 / term * alphak2 * alphai2 + 0.2D1 / 
     &0.35D2 * alphai2 ** 3 * alphak2 * r ** 4 + 0.4D1 / 0.105D3 * alpha
     &i2 ** 5 * alphak2 / term * r ** 4 + alphai2 ** 4 * alphak2 * r ** 
     &5 / 0.105D3 + 0.4D1 * alphai2 ** 2 * alphak2 / term * r + 0.12D2 /
     & 0.7D1 * alphai2 ** 3 * alphak2 / term * r2 + 0.8D1 / 0.21D2 * alp
     &hai2 ** 4 * alphak2 / term * r ** 3 + alphai2 * alphak2 * r2 / 0.7
     &D1) * expdampi
c
      s = s*rr1
      ds = ds*rr3
      dds = dds*rr5
      ddds = ddds*rr7
      dddds = dddds*rr9
c
      oscale1ik = pre*s**2
      oscale3ik = 2.0d0*pre*s*ds
      oscale5ik = 2.0d0*pre*(ds*ds + s*dds)
      oscale7ik = 2.0d0*pre*(dds*ds + ds*dds + ds*dds + s*ddds)
      oscale9ik = 2.0d0*pre*(ddds*ds + dds*dds + dds*dds + ds*ddds +
     &     dds*dds + ds*ddds + ds*ddds + s*dddds)               
c
                  else
                     alphai = min(alphai,alphak)
                     alphai2 = min(alphai2,alphak2)
                     dampi = alphai2*r
                     expdampi = exp(-dampi)
c                     dampi = min(dampi,dampk)
c                     expdampi = max(expdampi,expdampk)
                     pre = (alphai**6)/(alphai2**6)
c
      s = (r + r2 * alphai2 + r ** 3 * alphai2 ** 2 / 0.3D1) * expdampi
      ds = (r ** 3 * alphai2 ** 2 / 0.3D1 + r ** 4 * alphai2 ** 3 / 0.3D
     &1) * expdampi
      dds = alphai2 ** 4 * expdampi * r ** 5 / 0.9D1
      ddds = alphai2 ** 5 * expdampi * r ** 6 / 0.45D2
      dddds = (alphai2 ** 5 * r ** 6 / 0.315D3 + alphai2 ** 6 * r ** 7 /
     & 0.315D3) * expdampi
c
      s = s*rr1
      ds = ds*rr3
      dds = dds*rr5
      ddds = ddds*rr7
      dddds = dddds*rr9
c
      oscale1ik = pre*s**2
      oscale3ik = 2.0d0*pre*s*ds
      oscale5ik = 2.0d0*pre*(ds*ds + s*dds)
      oscale7ik = 2.0d0*pre*(dds*ds + ds*dds + ds*dds + s*ddds)
      oscale9ik = 2.0d0*pre*(ddds*ds + dds*dds + dds*dds + ds*ddds +
     &     dds*dds + ds*ddds + ds*ddds + s*dddds)
c                  
                  end if
c     
c     undamped energy
c
                  ee = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                 + gl(3)*rr7 + gl(4)*rr9
                  en = gln(0)*rr1 + gln(1)*rr3 + gln(2)*rr5
                  eundamp = ee + en
c
c     exchange energy
c                  
                  es2r = gl(0)*oscale1ik
     &                 + gl(1)*oscale3ik
     &                 + gl(2)*oscale5ik
     &                 + gl(3)*oscale7ik
     &                 + gl(4)*oscale9ik
c
c                  print *,"energy contributions"
c                  print *,atomic(i),atomic(k),i,k,r
c                  print *,"charge-charge",gl(0)*oscale1ik
c                  print *,"charge-dipole",gl(1)*oscale3ik
c                  print *,"d-d and q-Q  ",gl(2)*oscale5ik
c                  print *,"Q-d          ",gl(3)*oscale7ik
c                  print *,"Q-Q          ",gl(4)*oscale9ik
c                  print *,"total        ",es2r
c
c     total exchange energy for this pair
c
                  if (cptype .eq. "CHARGE-CHARGE") then
                     es2r = gl(0)*oscale1ik
                  end if
                  e = oik*zi*zk*es2r/r
c
c     increment the total exchange energy
c
                  ex = ex + e
c
c     increment the total intermolecular energy
c
                  if (molcule(ii) .ne. molcule(kk)) then
                     einter = einter + e
                  end if
c
c     use combining rule instead of two exponentials
c
               else if (exmodel .eq. "S2-COMBINE") then
                  alphai2 = 0.5d0*alphai
                  alphak2 = 0.5d0*alphak
                  pre = (alphai**6)/(alphai2**6)
c
c     alpha combining rule
c
                  alpham = (alphai2 + alphak2)/2.0d0
                  dampm = alpham*r
                  expdampm = exp(-dampm)
                  pre = ((2.0d0*alpham)**6)/(alpham**6)
c
c     compute scaling factors for combined exchange function
c
      s = (r + r2 * alpham + r ** 3 * alpham ** 2 / 0.3D1) * expdampm
      ds = (r ** 3 * alpham ** 2 / 0.3D1 + r ** 4 * alpham ** 3 / 0.3D
     &1) * expdampm
      dds = alpham ** 4 * expdampm * r ** 5 / 0.9D1
      ddds = alpham ** 5 * expdampm * r ** 6 / 0.45D2
      dddds = (alpham ** 5 * r ** 6 / 0.315D3 + alpham ** 6 * r ** 7 /
     & 0.315D3) * expdampm
c
      s = s*rr1
      ds = ds*rr3
      dds = dds*rr5
      ddds = ddds*rr7
      dddds = dddds*rr9
c
      oscale1ik = pre*s**2
      oscale3ik = 2.0d0*pre*s*ds
      oscale5ik = 2.0d0*pre*(ds*ds + s*dds)
      oscale7ik = 2.0d0*pre*(dds*ds + ds*dds + ds*dds + s*ddds)
      oscale9ik = 2.0d0*pre*(ddds*ds + dds*dds + dds*dds + ds*ddds +
     &     dds*dds + ds*ddds + ds*ddds + s*dddds)
c
c     exchange energy
c                  
                  es2r = gl(0)*oscale1ik
     &                 + gl(1)*oscale3ik
     &                 + gl(2)*oscale5ik
     &                 + gl(3)*oscale7ik
     &                 + gl(4)*oscale9ik
c
c     compute total exchange energy
c
                  e = oik*zi*zk*es2r/r
c
c     increment
c
                  ex = ex + e
c
c     increment the total intermolecular energy
c
                  if (molcule(ii) .ne. molcule(kk)) then
                     einter = einter + e
                  end if
               else
c                  print *,"ERROR: SET EXCH-MODEL    ORBITAL2ONLYR"
               end if
            end if
         end do
      end do
      return
      end
