c
c
c     #############################################################
c     ## COPYRIGHT (C) 2016 by Josh Rackers & Jay William Ponder ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  fieldstuff  --  contains routines to compute potential,    ##
c     ##  field, field gradient and field double gradient            ##
c     ##                                                             ##
c     #################################################################
c
c
c     "fieldstuff" has routines to compute the electric potential, field
c     field gradient or field double gradient
c
c
c     #################################################
c     ##                                             ##
c     ##  subroutine potik  --  electric potential   ##
c     ##                                             ##
c     #################################################
c
c
c     "potik" computes the electric potential given two multipoles
c     and damped or undamped tmatrix
c
      subroutine potik(i,k,t0,t1,t2,poti,potk)
      use mpole
      implicit none
      integer i,k
      real*8 ci
      real*8 dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck
      real*8 dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 t0,t1(3),t2(3,3)
      real*8 poti,potk
c
c     read in multipole values
c
      ci = rpole(1,i)
      dix = rpole(2,i)
      diy = rpole(3,i)
      diz = rpole(4,i)
      qixx = rpole(5,i)
      qixy = rpole(6,i)*2.0d0
      qixz = rpole(7,i)*2.0d0
      qiyy = rpole(9,i)
      qiyz = rpole(10,i)*2.0d0
      qizz = rpole(13,i)
      ck = rpole(1,k)
      dkx = rpole(2,k)
      dky = rpole(3,k)
      dkz = rpole(4,k)
      qkxx = rpole(5,k)
      qkxy = rpole(6,k)*2.0d0
      qkxz = rpole(7,k)*2.0d0
      qkyy = rpole(9,k)
      qkyz = rpole(10,k)*2.0d0
      qkzz = rpole(13,k)
c
c     calculate potential
c
c     from charges
      poti = t0 * ck
      potk = t0 * ci
c     from dipoles
      poti = poti + t1(1)*dkx + t1(2)*dky + t1(3)*dkz
      potk = potk - t1(1)*dix - t1(2)*diy - t1(3)*diz
c     from quadrupoles
      poti = poti + t2(1,1)*qkxx + t2(2,1)*qkxy + 
     &     t2(3,1)*qkxz + t2(2,2)*qkyy + 
     &     t2(3,2)*qkyz + t2(3,3)*qkzz
      potk = potk + t2(1,1)*qixx + t2(2,1)*qixy +
     &     t2(3,1)*qixz + t2(2,2)*qiyy + 
     &     t2(3,2)*qiyz + t2(3,3)*qizz
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine cp_potik  --  electric potential   ##
c     ##                                                ##
c     ####################################################
c
c
c     "cp_potik" computes the electric potential given two multipoles
c     and damped or undamped tmatrix
c
      subroutine cp_potik(i,k,
     &     t0,t1,t2,t0i,t1i,t2i,t0k,t1k,t2k,t0ik,t1ik,t2ik,
     &     nucpoti,nucpotk,elepoti,elepotk)
      use atomid
      use chgpen
      use mpole
      implicit none
      integer i,k
      real*8 ci,zi,qi
      real*8 dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,zk,qk
      real*8 dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 t0,t1(3),t2(3,3)
      real*8 t0i,t1i(3),t2i(3,3)
      real*8 t0k,t1k(3),t2k(3,3)
      real*8 t0ik,t1ik(3),t2ik(3,3)
      real*8 nucpoti,nucpotk
      real*8 elepoti,elepotk
c
c     read in multipole values
c
      ci = rpole(1,i)
      dix = rpole(2,i)
      diy = rpole(3,i)
      diz = rpole(4,i)
      qixx = rpole(5,i)
      qixy = rpole(6,i)*2.0d0
      qixz = rpole(7,i)*2.0d0
      qiyy = rpole(9,i)
      qiyz = rpole(10,i)*2.0d0
      qizz = rpole(13,i)
      ck = rpole(1,k)
      dkx = rpole(2,k)
      dky = rpole(3,k)
      dkz = rpole(4,k)
      qkxx = rpole(5,k)
      qkxy = rpole(6,k)*2.0d0
      qkxz = rpole(7,k)*2.0d0
      qkyy = rpole(9,k)
      qkyz = rpole(10,k)*2.0d0
      qkzz = rpole(13,k)
c
c     split nuclear and electronic charges
c
      zi = atomic(i)
      zk = atomic(k)
      if (num_ele .eq. "VALENCE") then
         if (atomic(i) .gt. 2)  zi = zi - 2.0d0
         if (atomic(i) .gt. 10)  zi = zi - 8.0d0
         if (atomic(i) .gt. 18)  zi = zi - 8.0d0
         if (atomic(i) .gt. 20)  zi = zi - 10.0d0
         if (atomic(k) .gt. 2)  zk = zk - 2.0d0
         if (atomic(k) .gt. 10)  zk = zk - 8.0d0
         if (atomic(k) .gt. 18)  zk = zk - 8.0d0
         if (atomic(k) .gt. 20)  zk = zk - 10.0d0
      end if
      qi = ci - zi
      qk = ck - zk
c
c     calculate potential at nuclei
c
c     from nuclei
      nucpoti = t0 * zk
      nucpotk = t0 * zi
c     from electrons
      nucpoti = nucpoti + t0k * qk
      nucpotk = nucpotk + t0i * qi
c     from dipoles
      nucpoti = nucpoti + t1k(1)*dkx + t1k(2)*dky + t1k(3)*dkz
      nucpotk = nucpotk - t1i(1)*dix - t1i(2)*diy - t1i(3)*diz
c     from quadrupoles
      nucpoti = nucpoti + t2k(1,1)*qkxx + t2k(2,1)*qkxy +
     &     t2k(3,1)*qkxz + t2k(2,2)*qkyy +
     &     t2k(3,2)*qkyz + t2k(3,3)*qkzz
      nucpotk = nucpotk + t2i(1,1)*qixx + t2i(2,1)*qixy +
     &     t2i(3,1)*qixz + t2i(2,2)*qiyy +
     &     t2i(3,2)*qiyz + t2i(3,3)*qizz
c
c     calculate potential at electrons
c
c     from nuclei
      elepoti = t0i * zk
      elepotk = t0k * zi
c     from electrons
      elepoti = elepoti + t0ik * qk
      elepotk = elepotk + t0ik * qi
c     from dipoles
      elepoti = elepoti + t1ik(1)*dkx + t1ik(2)*dky + t1ik(3)*dkz
      elepotk = elepotk - t1ik(1)*dix - t1ik(2)*diy - t1ik(3)*diz
c     from quadrupoles
      elepoti = elepoti + t2ik(1,1)*qkxx + t2ik(2,1)*qkxy + 
     &     t2ik(3,1)*qkxz + t2ik(2,2)*qkyy + 
     &     t2ik(3,2)*qkyz + t2ik(3,3)*qkzz
      elepotk = elepotk + t2ik(1,1)*qixx + t2ik(2,1)*qixy +
     &     t2ik(3,1)*qixz + t2ik(2,2)*qiyy + 
     &     t2ik(3,2)*qiyz + t2ik(3,3)*qizz
      return
      end
c
c
c     #################################################
c     ##                                             ##
c     ##  subroutine upotik  --  electric potential   ##
c     ##                                             ##
c     #################################################
c
c
c     "upotik" computes the electric potential given two induced dipoles
c     and damped or undamped tmatrix
c
      subroutine upotik(i,k,t1,potid,potkd,potip,potkp)
      use polar
      implicit none
      integer i,k
      real*8 duix,duiy,duiz
      real*8 dukx,duky,dukz
      real*8 puix,puiy,puiz
      real*8 pukx,puky,pukz
      real*8 t1(3)
      real*8 potid,potkd
      real*8 potip,potkp
c
c     read in induced dipole values
c
      duix = uind(1,i)
      duiy = uind(2,i)
      duiz = uind(3,i)
      dukx = uind(1,k)
      duky = uind(2,k)
      dukz = uind(3,k)
      puix = uinp(1,i)
      puiy = uinp(2,i)
      puiz = uinp(3,i)
      pukx = uinp(1,k)
      puky = uinp(2,k)
      pukz = uinp(3,k)
c
c     calculate potential
c
c     from dipoles
      potid = t1(1)*dukx + t1(2)*duky + t1(3)*dukz
      potkd = -t1(1)*duix - t1(2)*duiy - t1(3)*duiz
      potip = t1(1)*pukx + t1(2)*puky + t1(3)*pukz
      potkp = -t1(1)*puix - t1(2)*puiy - t1(3)*puiz
      return
      end
c
c
c     #################################################
c     ##                                             ##
c     ##  subroutine fieldik  --  electric field     ##
c     ##                                             ##
c     #################################################
c
c
c     "fieldik" computes the electric potential given two multipoles
c     and damped or undamped tmatrix
c
      subroutine fieldik(i,k,t1,t2,t3,fieldi,fieldk)
      use mpole
      implicit none
      integer i,k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 t1(3),t2(3,3),t3(3,3,3)
      real*8 fieldi(3),fieldk(3)
c
c     read in multipole values
c
      ci = rpole(1,i)
      dix = rpole(2,i)
      diy = rpole(3,i)
      diz = rpole(4,i)
      qixx = rpole(5,i)
      qixy = rpole(6,i)*2.0d0
      qixz = rpole(7,i)*2.0d0
      qiyy = rpole(9,i)
      qiyz = rpole(10,i)*2.0d0
      qizz = rpole(13,i)
      ck = rpole(1,k)
      dkx = rpole(2,k)
      dky = rpole(3,k)
      dkz = rpole(4,k)
      qkxx = rpole(5,k)
      qkxy = rpole(6,k)*2.0d0
      qkxz = rpole(7,k)*2.0d0
      qkyy = rpole(9,k)
      qkyz = rpole(10,k)*2.0d0
      qkzz = rpole(13,k)
c
c     calculate the electric field
c     
c     from charges
      fieldi(1) = -t1(1)*ck
      fieldi(2) = -t1(2)*ck
      fieldi(3) = -t1(3)*ck
c     
      fieldk(1) = t1(1)*ci
      fieldk(2) = t1(2)*ci
      fieldk(3) = t1(3)*ci
c     from dipoles
      fieldi(1) = fieldi(1) - t2(1,1)*dkx - t2(2,1)*dky - 
     &     t2(3,1)*dkz
      fieldi(2) = fieldi(2) - t2(2,1)*dkx - t2(2,2)*dky - 
     &     t2(3,2)*dkz
      fieldi(3) = fieldi(3) - t2(3,1)*dkx - t2(3,2)*dky -
     &     t2(3,3)*dkz
c     
      fieldk(1) = fieldk(1) - t2(1,1)*dix - t2(2,1)*diy -
     &     t2(3,1)*diz
      fieldk(2) = fieldk(2) - t2(2,1)*dix - t2(2,2)*diy -
     &     t2(3,2)*diz
      fieldk(3) = fieldk(3) - t2(3,1)*dix - t2(3,2)*diy -
     &     t2(3,3)*diz
c     from quadrupoles
      fieldi(1) = fieldi(1)- t3(1,1,1)*qkxx - t3(2,1,1)*qkxy
     &     -t3(3,1,1)*qkxz - t3(2,2,1)*qkyy - t3(3,2,1)*qkyz 
     &     -t3(3,3,1)*qkzz
      fieldi(2) = fieldi(2)- t3(2,1,1)*qkxx - t3(2,2,1)*qkxy
     &     -t3(3,2,1)*qkxz - t3(2,2,2)*qkyy - t3(3,2,2)*qkyz
     &     -t3(3,3,2)*qkzz
      fieldi(3) = fieldi(3)- t3(3,1,1)*qkxx - t3(3,2,1)*qkxy
     &     -t3(3,3,1)*qkxz - t3(3,2,2)*qkyy - t3(3,3,2)*qkyz
     &     -t3(3,3,3)*qkzz
c     
      fieldk(1) = fieldk(1)+ t3(1,1,1)*qixx + t3(2,1,1)*qixy
     &     +t3(3,1,1)*qixz + t3(2,2,1)*qiyy + t3(3,2,1)*qiyz
     &     +t3(3,3,1)*qizz
      fieldk(2) = fieldk(2)+ t3(2,1,1)*qixx + t3(2,2,1)*qixy
     &     +t3(3,2,1)*qixz + t3(2,2,2)*qiyy + t3(3,2,2)*qiyz
     &     +t3(3,3,2)*qizz
      fieldk(3) = fieldk(3)+ t3(3,1,1)*qixx + t3(3,2,1)*qixy
     &     +t3(3,3,1)*qixz + t3(3,2,2)*qiyy + t3(3,3,2)*qiyz
     &     +t3(3,3,3)*qizz
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine cp_fieldik  --  electric field     ##
c     ##                                                ##
c     ####################################################
c
c
c     "cp_fieldik" computes the electric potential given two multipoles
c     and damped or undamped tmatrix
c
      subroutine cp_fieldik(i,k,
     &     t1,t2,t3,t1i,t2i,t3i,t1k,t2k,t3k,t1ik,t2ik,t3ik,
     &     nucfieldi,nucfieldk,elefieldi,elefieldk)
      use atomid
      use chgpen
      use mpole
      implicit none
      integer i,k
      real*8 ci,zi,qi
      real*8 dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,zk,qk
      real*8 dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 t1(3),t2(3,3),t3(3,3,3)
      real*8 t1i(3),t2i(3,3),t3i(3,3,3)
      real*8 t1k(3),t2k(3,3),t3k(3,3,3)
      real*8 t1ik(3),t2ik(3,3),t3ik(3,3,3)
      real*8 nucfieldi(3),nucfieldk(3)
      real*8 elefieldi(3),elefieldk(3)
c
c     read in multipole values
c
      ci = rpole(1,i)
      dix = rpole(2,i)
      diy = rpole(3,i)
      diz = rpole(4,i)
      qixx = rpole(5,i)
      qixy = rpole(6,i)*2.0d0
      qixz = rpole(7,i)*2.0d0
      qiyy = rpole(9,i)
      qiyz = rpole(10,i)*2.0d0
      qizz = rpole(13,i)
      ck = rpole(1,k)
      dkx = rpole(2,k)
      dky = rpole(3,k)
      dkz = rpole(4,k)
      qkxx = rpole(5,k)
      qkxy = rpole(6,k)*2.0d0
      qkxz = rpole(7,k)*2.0d0
      qkyy = rpole(9,k)
      qkyz = rpole(10,k)*2.0d0
      qkzz = rpole(13,k)
c
c     split nuclear and electronic charges
c
      zi = atomic(i)
      zk = atomic(k)
      if (num_ele .eq. "VALENCE") then
         if (atomic(i) .gt. 2)  zi = zi - 2.0d0
         if (atomic(i) .gt. 10)  zi = zi - 8.0d0
         if (atomic(i) .gt. 18)  zi = zi - 8.0d0
         if (atomic(i) .gt. 20)  zi = zi - 10.0d0
         if (atomic(k) .gt. 2)  zk = zk - 2.0d0
         if (atomic(k) .gt. 10)  zk = zk - 8.0d0
         if (atomic(k) .gt. 18)  zk = zk - 8.0d0
         if (atomic(k) .gt. 20)  zk = zk - 10.0d0
      end if
      qi = ci - zi
      qk = ck - zk
c
c     calculate the electric field at nuclei
c     
c     from nuclei
      nucfieldi(1) = -t1(1)*zk
      nucfieldi(2) = -t1(2)*zk
      nucfieldi(3) = -t1(3)*zk
c     
      nucfieldk(1) = t1(1)*zi
      nucfieldk(2) = t1(2)*zi
      nucfieldk(3) = t1(3)*zi
c     from electrons
      nucfieldi(1) = nucfieldi(1) - t1k(1)*qk
      nucfieldi(2) = nucfieldi(2) - t1k(2)*qk
      nucfieldi(3) = nucfieldi(3) - t1k(3)*qk
c
      nucfieldk(1) = nucfieldk(1) + t1i(1)*qi
      nucfieldk(2) = nucfieldk(2) + t1i(2)*qi
      nucfieldk(3) = nucfieldk(3) + t1i(3)*qi
c     from dipoles
      nucfieldi(1) = nucfieldi(1) - t2k(1,1)*dkx - t2k(2,1)*dky - 
     &     t2k(3,1)*dkz
      nucfieldi(2) = nucfieldi(2) - t2k(2,1)*dkx - t2k(2,2)*dky - 
     &     t2k(3,2)*dkz
      nucfieldi(3) = nucfieldi(3) - t2k(3,1)*dkx - t2k(3,2)*dky -
     &     t2k(3,3)*dkz
c     
      nucfieldk(1) = nucfieldk(1) - t2i(1,1)*dix - t2i(2,1)*diy -
     &     t2i(3,1)*diz
      nucfieldk(2) = nucfieldk(2) - t2i(2,1)*dix - t2i(2,2)*diy -
     &     t2i(3,2)*diz
      nucfieldk(3) = nucfieldk(3) - t2i(3,1)*dix - t2i(3,2)*diy -
     &     t2i(3,3)*diz
c     from quadrupoles
      nucfieldi(1) = nucfieldi(1)- t3k(1,1,1)*qkxx - t3k(2,1,1)*qkxy
     &     -t3k(3,1,1)*qkxz - t3k(2,2,1)*qkyy - t3k(3,2,1)*qkyz 
     &     -t3k(3,3,1)*qkzz
      nucfieldi(2) = nucfieldi(2)- t3k(2,1,1)*qkxx - t3k(2,2,1)*qkxy
     &     -t3k(3,2,1)*qkxz - t3k(2,2,2)*qkyy - t3k(3,2,2)*qkyz
     &     -t3k(3,3,2)*qkzz
      nucfieldi(3) = nucfieldi(3)- t3k(3,1,1)*qkxx - t3k(3,2,1)*qkxy
     &     -t3k(3,3,1)*qkxz - t3k(3,2,2)*qkyy - t3k(3,3,2)*qkyz
     &     -t3k(3,3,3)*qkzz
c     
      nucfieldk(1) = nucfieldk(1)+ t3i(1,1,1)*qixx + t3i(2,1,1)*qixy
     &     +t3i(3,1,1)*qixz + t3i(2,2,1)*qiyy + t3i(3,2,1)*qiyz
     &     +t3i(3,3,1)*qizz
      nucfieldk(2) = nucfieldk(2)+ t3i(2,1,1)*qixx + t3i(2,2,1)*qixy
     &     +t3i(3,2,1)*qixz + t3i(2,2,2)*qiyy + t3i(3,2,2)*qiyz
     &     +t3i(3,3,2)*qizz
      nucfieldk(3) = nucfieldk(3)+ t3i(3,1,1)*qixx + t3i(3,2,1)*qixy
     &     +t3i(3,3,1)*qixz + t3i(3,2,2)*qiyy + t3i(3,3,2)*qiyz
     &     +t3i(3,3,3)*qizz
c
c     calculate the electric field at electrons
c     
c     from nuclei
      elefieldi(1) = -t1i(1)*zk
      elefieldi(2) = -t1i(2)*zk
      elefieldi(3) = -t1i(3)*zk
c     
      elefieldk(1) = t1k(1)*zi
      elefieldk(2) = t1k(2)*zi
      elefieldk(3) = t1k(3)*zi
c     from electrons
      elefieldi(1) = elefieldi(1) - t1ik(1)*qk
      elefieldi(2) = elefieldi(2) - t1ik(2)*qk
      elefieldi(3) = elefieldi(3) - t1ik(3)*qk
c
      elefieldk(1) = elefieldk(1) + t1ik(1)*qi
      elefieldk(2) = elefieldk(2) + t1ik(2)*qi
      elefieldk(3) = elefieldk(3) + t1ik(3)*qi
c     from dipoles
      elefieldi(1) = elefieldi(1) - t2ik(1,1)*dkx - t2ik(2,1)*dky - 
     &     t2ik(3,1)*dkz
      elefieldi(2) = elefieldi(2) - t2ik(2,1)*dkx - t2ik(2,2)*dky - 
     &     t2ik(3,2)*dkz
      elefieldi(3) = elefieldi(3) - t2ik(3,1)*dkx - t2ik(3,2)*dky -
     &     t2ik(3,3)*dkz
c     
      elefieldk(1) = elefieldk(1) - t2ik(1,1)*dix - t2ik(2,1)*diy -
     &     t2ik(3,1)*diz
      elefieldk(2) = elefieldk(2) - t2ik(2,1)*dix - t2ik(2,2)*diy -
     &     t2ik(3,2)*diz
      elefieldk(3) = elefieldk(3) - t2ik(3,1)*dix - t2ik(3,2)*diy -
     &     t2ik(3,3)*diz
c     from quadrupoles
      elefieldi(1) = elefieldi(1)- t3ik(1,1,1)*qkxx - t3ik(2,1,1)*qkxy
     &     -t3ik(3,1,1)*qkxz - t3ik(2,2,1)*qkyy - t3ik(3,2,1)*qkyz 
     &     -t3ik(3,3,1)*qkzz
      elefieldi(2) = elefieldi(2)- t3ik(2,1,1)*qkxx - t3ik(2,2,1)*qkxy
     &     -t3ik(3,2,1)*qkxz - t3ik(2,2,2)*qkyy - t3ik(3,2,2)*qkyz
     &     -t3ik(3,3,2)*qkzz
      elefieldi(3) = elefieldi(3)- t3ik(3,1,1)*qkxx - t3ik(3,2,1)*qkxy
     &     -t3ik(3,3,1)*qkxz - t3ik(3,2,2)*qkyy - t3ik(3,3,2)*qkyz
     &     -t3ik(3,3,3)*qkzz
c     
      elefieldk(1) = elefieldk(1)+ t3ik(1,1,1)*qixx + t3ik(2,1,1)*qixy
     &     +t3ik(3,1,1)*qixz + t3ik(2,2,1)*qiyy + t3ik(3,2,1)*qiyz
     &     +t3ik(3,3,1)*qizz
      elefieldk(2) = elefieldk(2)+ t3ik(2,1,1)*qixx + t3ik(2,2,1)*qixy
     &     +t3ik(3,2,1)*qixz + t3ik(2,2,2)*qiyy + t3ik(3,2,2)*qiyz
     &     +t3ik(3,3,2)*qizz
      elefieldk(3) = elefieldk(3)+ t3ik(3,1,1)*qixx + t3ik(3,2,1)*qixy
     &     +t3ik(3,3,1)*qixz + t3ik(3,2,2)*qiyy + t3ik(3,3,2)*qiyz
     &     +t3ik(3,3,3)*qizz
      return
      end
c
c
c     #################################################
c     ##                                             ##
c     ##  subroutine ufieldik  --  electric field    ##
c     ##                                             ##
c     #################################################
c
c
c     "ufieldik" computes the electric field given two induced dipoles
c     and damped or undamped tmatrix
c
      subroutine ufieldik(i,k,t2,fieldid,fieldkd,fieldip,fieldkp)
      use polar
      implicit none
      integer i,k
      real*8 duix,duiy,duiz
      real*8 dukx,duky,dukz
      real*8 puix,puiy,puiz
      real*8 pukx,puky,pukz
      real*8 t2(3,3)
      real*8 fieldid(3),fieldkd(3)
      real*8 fieldip(3),fieldkp(3)
c
c     read in induced dipole values
c
      duix = uind(1,i)
      duiy = uind(2,i)
      duiz = uind(3,i)
      dukx = uind(1,k)
      duky = uind(2,k)
      dukz = uind(3,k)
      puix = uinp(1,i)
      puiy = uinp(2,i)
      puiz = uinp(3,i)
      pukx = uinp(1,k)
      puky = uinp(2,k)
      pukz = uinp(3,k)
c
c     calculate the electric field
c     
      fieldid(1) = -t2(1,1)*dukx - t2(2,1)*duky - 
     &     t2(3,1)*dukz
      fieldid(2) = -t2(2,1)*dukx - t2(2,2)*duky - 
     &     t2(3,2)*dukz
      fieldid(3) = -t2(3,1)*dukx - t2(3,2)*duky -
     &     t2(3,3)*dukz
c     
      fieldkd(1) = -t2(1,1)*duix - t2(2,1)*duiy -
     &     t2(3,1)*duiz
      fieldkd(2) = -t2(2,1)*duix - t2(2,2)*duiy -
     &     t2(3,2)*duiz
      fieldkd(3) = -t2(3,1)*duix - t2(3,2)*duiy -
     &     t2(3,3)*duiz
c
      fieldip(1) = -t2(1,1)*pukx - t2(2,1)*puky -
     &     t2(3,1)*pukz
      fieldip(2) = -t2(2,1)*pukx - t2(2,2)*puky -
     &     t2(3,2)*pukz
      fieldip(3) = -t2(3,1)*pukx - t2(3,2)*puky -
     &     t2(3,3)*pukz
c
      fieldkp(1) = -t2(1,1)*puix - t2(2,1)*puiy -
     &     t2(3,1)*puiz
      fieldkp(2) = -t2(2,1)*puix - t2(2,2)*puiy -
     &     t2(3,2)*puiz
      fieldkp(3) = -t2(3,1)*puix - t2(3,2)*puiy -
     &     t2(3,3)*puiz
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine gradfieldik  --  electric field gradient ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "gradfieldik" computes the electric field gradient given two multipoles
c     and damped or undamped tmatrix
c
      subroutine gradfieldik(i,k,t2,t3,t4,gradfieldi,gradfieldk)
      use mpole
      implicit none
      integer i,k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 t2(3,3),t3(3,3,3),t4(3,3,3,3)
      real*8 gradfieldi(3,3),gradfieldk(3,3)
c
c     read in multipole values
c
      ci = rpole(1,i)
      dix = rpole(2,i)
      diy = rpole(3,i)
      diz = rpole(4,i)
      qixx = rpole(5,i)
      qixy = rpole(6,i)*2.0d0
      qixz = rpole(7,i)*2.0d0
      qiyy = rpole(9,i)
      qiyz = rpole(10,i)*2.0d0
      qizz = rpole(13,i)
      ck = rpole(1,k)
      dkx = rpole(2,k)
      dky = rpole(3,k)
      dkz = rpole(4,k)
      qkxx = rpole(5,k)
      qkxy = rpole(6,k)*2.0d0
      qkxz = rpole(7,k)*2.0d0
      qkyy = rpole(9,k)
      qkyz = rpole(10,k)*2.0d0
      qkzz = rpole(13,k)
c
c     calculate electric field gradient
c
c     from charges
      gradfieldi(1,1) = t2(1,1)*ck
      gradfieldi(2,1) = t2(2,1)*ck
      gradfieldi(3,1) = t2(3,1)*ck
c      gradfieldi(2,1) = t2(2,1)*ck
      gradfieldi(2,2) = t2(2,2)*ck
      gradfieldi(3,2) = t2(3,2)*ck
c      gradfieldi(3,1) = t2(3,1)*ck
c      gradfieldi(3,2) = t2(3,2)*ck
      gradfieldi(3,3) = t2(3,3)*ck
c     
      gradfieldk(1,1) = t2(1,1)*ci
      gradfieldk(2,1) = t2(2,1)*ci
      gradfieldk(3,1) = t2(3,1)*ci
c      gradfieldk(2,1) = t2(2,1)*ci
      gradfieldk(2,2) = t2(2,2)*ci
      gradfieldk(3,2) = t2(3,2)*ci
c      gradfieldk(3,1) = t2(3,1)*ci
c      gradfieldk(3,2) = t2(3,2)*ci
      gradfieldk(3,3) = t2(3,3)*ci
c     from dipoles
      gradfieldi(1,1) = gradfieldi(1,1) + t3(1,1,1)*dkx + 
     &     t3(2,1,1)*dky + t3(3,1,1)*dkz
      gradfieldi(2,1) = gradfieldi(2,1) + t3(2,1,1)*dkx +
     &     t3(2,2,1)*dky + t3(3,2,1)*dkz
      gradfieldi(3,1) = gradfieldi(3,1) + t3(3,1,1)*dkx +
     &     t3(3,2,1)*dky + t3(3,3,1)*dkz
c      gradfieldi(2,1) = gradfieldi(2,1) + t3(2,1,1)*dkx +
c     &     t3(2,2,1)*dky + t3(3,2,1)*dkz
      gradfieldi(2,2) = gradfieldi(2,2) + t3(2,2,1)*dkx +
     &     t3(2,2,2)*dky + t3(3,2,2)*dkz
      gradfieldi(3,2) = gradfieldi(3,2) + t3(3,2,1)*dkx +
     &     t3(3,2,2)*dky + t3(3,3,2)*dkz
c      gradfieldi(3,1) = gradfieldi(3,1) + t3(3,1,1)*dkx +
c     &     t3(3,2,1)*dky + t3(3,3,1)*dkz
c      gradfieldi(3,2) = gradfieldi(3,2) + t3(3,2,1)*dkx +
c     &     t3(3,2,2)*dky + t3(3,3,2)*dkz
      gradfieldi(3,3) = gradfieldi(3,3) + t3(3,3,1)*dkx +
     &     t3(3,3,2)*dky + t3(3,3,3)*dkz
c     
      gradfieldk(1,1) = gradfieldk(1,1) - t3(1,1,1)*dix -
     &     t3(2,1,1)*diy - t3(3,1,1)*diz
      gradfieldk(2,1) = gradfieldk(2,1) - t3(2,1,1)*dix -
     &     t3(2,2,1)*diy - t3(3,2,1)*diz
      gradfieldk(3,1) = gradfieldk(3,1) - t3(3,1,1)*dix -
     &     t3(3,2,1)*diy - t3(3,3,1)*diz
c      gradfieldk(2,1) = gradfieldk(2,1) - t3(2,1,1)*dix -
c     &                 t3(2,2,1)*diy - t3(3,2,1)*diz
      gradfieldk(2,2) = gradfieldk(2,2) - t3(2,2,1)*dix -
     &     t3(2,2,2)*diy - t3(3,2,2)*diz
      gradfieldk(3,2) = gradfieldk(3,2) - t3(3,2,1)*dix -
     &     t3(3,2,2)*diy - t3(3,3,2)*diz
c      gradfieldk(3,1) = gradfieldk(3,1) - t3(3,1,1)*dix -
c     &     t3(3,2,1)*diy - t3(3,3,1)*diz
c      gradfieldk(3,2) = gradfieldk(3,2) - t3(3,2,1)*dix -
c     &     t3(3,2,2)*diy - t3(3,3,2)*diz
      gradfieldk(3,3) = gradfieldk(3,3) - t3(3,3,1)*dix -
     &     t3(3,3,2)*diy - t3(3,3,3)*diz
c     from quadrupoles
      gradfieldi(1,1) = gradfieldi(1,1) + t4(1,1,1,1)*qkxx +
     &     t4(2,1,1,1)*qkxy + t4(3,1,1,1)*qkxz + 
     &     t4(2,2,1,1)*qkyy + t4(3,2,1,1)*qkyz + 
     &     t4(3,3,1,1)*qkzz
      gradfieldi(2,1) = gradfieldi(2,1) + t4(2,1,1,1)*qkxx +
     &     t4(2,2,1,1)*qkxy + t4(3,2,1,1)*qkxz +
     &     t4(2,2,2,1)*qkyy + t4(3,2,2,1)*qkyz +
     &     t4(3,3,2,1)*qkzz
      gradfieldi(3,1) = gradfieldi(3,1) + t4(3,1,1,1)*qkxx +
     &     t4(3,2,1,1)*qkxy + t4(3,3,1,1)*qkxz +
     &     t4(3,2,2,1)*qkyy + t4(3,3,2,1)*qkyz +
     &     t4(3,3,3,1)*qkzz
c      gradfieldi(2,1) = gradfieldi(2,1) + t4(2,1,1,1)*qkxx +
c     &     t4(2,2,1,1)*qkxy + t4(3,2,1,1)*qkxz +
c     &     t4(2,2,2,1)*qkyy + t4(3,2,2,1)*qkyz +
c     &     t4(3,3,2,1)*qkzz
      gradfieldi(2,2) = gradfieldi(2,2) + t4(2,2,1,1)*qkxx +
     &     t4(2,2,2,1)*qkxy + t4(3,2,2,1)*qkxz +
     &     t4(2,2,2,2)*qkyy + t4(3,2,2,2)*qkyz +
     &     t4(3,3,2,2)*qkzz
      gradfieldi(3,2) = gradfieldi(3,2) + t4(3,2,1,1)*qkxx +
     &     t4(3,2,2,1)*qkxy + t4(3,3,2,1)*qkxz +
     &     t4(3,2,2,2)*qkyy + t4(3,3,2,2)*qkyz +
     &     t4(3,3,3,2)*qkzz
c      gradfieldi(3,1) = gradfieldi(3,1) + t4(3,1,1,1)*qkxx +
c     &     t4(3,2,1,1)*qkxy + t4(3,3,1,1)*qkxz +
c     &     t4(3,2,2,1)*qkyy + t4(3,3,2,1)*qkyz +
c     &     t4(3,3,3,1)*qkzz
c      gradfieldi(3,2) = gradfieldi(3,2) + t4(3,2,1,1)*qkxx +
c     &     t4(3,2,2,1)*qkxy + t4(3,3,2,1)*qkxz +
c     &     t4(3,2,2,2)*qkyy + t4(3,3,2,2)*qkyz +
c     &     t4(3,3,3,2)*qkzz
      gradfieldi(3,3) = gradfieldi(3,3) + t4(3,3,1,1)*qkxx +
     &     t4(3,3,2,1)*qkxy + t4(3,3,3,1)*qkxz +
     &     t4(3,3,2,2)*qkyy + t4(3,3,3,2)*qkyz +
     &     t4(3,3,3,3)*qkzz
c     
      gradfieldk(1,1) = gradfieldk(1,1) + t4(1,1,1,1)*qixx +
     &     t4(2,1,1,1)*qixy + t4(3,1,1,1)*qixz +
     &     t4(2,2,1,1)*qiyy + t4(3,2,1,1)*qiyz +
     &     t4(3,3,1,1)*qizz
      gradfieldk(2,1) = gradfieldk(2,1) + t4(2,1,1,1)*qixx +
     &     t4(2,2,1,1)*qixy + t4(3,2,1,1)*qixz +
     &     t4(2,2,2,1)*qiyy + t4(3,2,2,1)*qiyz +
     &     t4(3,3,2,1)*qizz
      gradfieldk(3,1) = gradfieldk(3,1) + t4(3,1,1,1)*qixx +
     &     t4(3,2,1,1)*qixy + t4(3,3,1,1)*qixz +
     &     t4(3,2,2,1)*qiyy + t4(3,3,2,1)*qiyz +
     &     t4(3,3,3,1)*qizz
c      gradfieldk(2,1) = gradfieldk(2,1) + t4(2,1,1,1)*qixx +
c     &     t4(2,2,1,1)*qixy + t4(3,2,1,1)*qixz +
c     &     t4(2,2,2,1)*qiyy + t4(3,2,2,1)*qiyz +
c     &     t4(3,3,2,1)*qizz
      gradfieldk(2,2) = gradfieldk(2,2) + t4(2,2,1,1)*qixx +
     &     t4(2,2,2,1)*qixy + t4(3,2,2,1)*qixz +
     &     t4(2,2,2,2)*qiyy + t4(3,2,2,2)*qiyz +
     &     t4(3,3,2,2)*qizz
      gradfieldk(3,2) = gradfieldk(3,2) + t4(3,2,1,1)*qixx +
     &     t4(3,2,2,1)*qixy + t4(3,3,2,1)*qixz +
     &     t4(3,2,2,2)*qiyy + t4(3,3,2,2)*qiyz +
     &     t4(3,3,3,2)*qizz
c      gradfieldk(3,1) = gradfieldk(3,1) + t4(3,1,1,1)*qixx +
c     &     t4(3,2,1,1)*qixy + t4(3,3,1,1)*qixz +
c     &     t4(3,2,2,1)*qiyy + t4(3,3,2,1)*qiyz +
c     &     t4(3,3,3,1)*qizz
c      gradfieldk(3,2) = gradfieldk(3,2) + t4(3,2,1,1)*qixx +
c     &     t4(3,2,2,1)*qixy + t4(3,3,2,1)*qixz +
c     &     t4(3,2,2,2)*qiyy + t4(3,3,2,2)*qiyz +
c     &     t4(3,3,3,2)*qizz
      gradfieldk(3,3) = gradfieldk(3,3) + t4(3,3,1,1)*qixx +
     &     t4(3,3,2,1)*qixy + t4(3,3,3,1)*qixz +
     &     t4(3,3,2,2)*qiyy + t4(3,3,3,2)*qiyz +
     &     t4(3,3,3,3)*qizz
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine cp_gradfieldik  --  electric field gradient ##
c     ##                                                         ##
c     #############################################################
c
c
c     "cp_gradfieldik" computes the electric field gradient given two multipoles
c     and damped or undamped tmatrix
c
      subroutine cp_gradfieldik(i,k,
     &     t2,t3,t4,t2i,t3i,t4i,t2k,t3k,t4k,t2ik,t3ik,t4ik,
     &     elegradfieldi,elegradfieldk)
      use atomid
      use chgpen
      use mpole
      implicit none
      integer i,k
      real*8 ci,zi,qi
      real*8 dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,zk,qk
      real*8 dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 t2(3,3),t3(3,3,3),t4(3,3,3,3)
      real*8 t2i(3,3),t3i(3,3,3),t4i(3,3,3,3)
      real*8 t2k(3,3),t3k(3,3,3),t4k(3,3,3,3)
      real*8 t2ik(3,3),t3ik(3,3,3),t4ik(3,3,3,3)
      real*8 elegradfieldi(3,3),elegradfieldk(3,3)
c
c     read in multipole values
c
      ci = rpole(1,i)
      dix = rpole(2,i)
      diy = rpole(3,i)
      diz = rpole(4,i)
      qixx = rpole(5,i)
      qixy = rpole(6,i)*2.0d0
      qixz = rpole(7,i)*2.0d0
      qiyy = rpole(9,i)
      qiyz = rpole(10,i)*2.0d0
      qizz = rpole(13,i)
      ck = rpole(1,k)
      dkx = rpole(2,k)
      dky = rpole(3,k)
      dkz = rpole(4,k)
      qkxx = rpole(5,k)
      qkxy = rpole(6,k)*2.0d0
      qkxz = rpole(7,k)*2.0d0
      qkyy = rpole(9,k)
      qkyz = rpole(10,k)*2.0d0
      qkzz = rpole(13,k)
c
c     split nuclear and electronic charges
c
      zi = atomic(i)
      zk = atomic(k)
      if (num_ele .eq. "VALENCE") then
         if (atomic(i) .gt. 2)  zi = zi - 2.0d0
         if (atomic(i) .gt. 10)  zi = zi - 8.0d0
         if (atomic(i) .gt. 18)  zi = zi - 8.0d0
         if (atomic(i) .gt. 20)  zi = zi - 10.0d0
         if (atomic(k) .gt. 2)  zk = zk - 2.0d0
         if (atomic(k) .gt. 10)  zk = zk - 8.0d0
         if (atomic(k) .gt. 18)  zk = zk - 8.0d0
         if (atomic(k) .gt. 20)  zk = zk - 10.0d0
      end if
      qi = ci - zi
      qk = ck - zk
c
c     calculate electric field gradient at electrons
c
c     from nuclei
      elegradfieldi(1,1) = t2i(1,1)*zk
      elegradfieldi(2,1) = t2i(2,1)*zk
      elegradfieldi(3,1) = t2i(3,1)*zk
c      elegradfieldi(2,1) = t2i(2,1)*zk
      elegradfieldi(2,2) = t2i(2,2)*zk
      elegradfieldi(3,2) = t2i(3,2)*zk
c      elegradfieldi(3,1) = t2i(3,1)*zk
c      elegradfieldi(3,2) = t2i(3,2)*zk
      elegradfieldi(3,3) = t2i(3,3)*zk
c     
      elegradfieldk(1,1) = t2k(1,1)*zi
      elegradfieldk(2,1) = t2k(2,1)*zi
      elegradfieldk(3,1) = t2k(3,1)*zi
c      elegradfieldk(2,1) = t2k(2,1)*zi
      elegradfieldk(2,2) = t2k(2,2)*zi
      elegradfieldk(3,2) = t2k(3,2)*zi
c      elegradfieldk(3,1) = t2k(3,1)*zi
c      elegradfieldk(3,2) = t2k(3,2)*zi
      elegradfieldk(3,3) = t2k(3,3)*zi
c     from electrons
      elegradfieldi(1,1) = elegradfieldi(1,1) + t2ik(1,1)*qk
      elegradfieldi(2,1) = elegradfieldi(2,1) + t2ik(2,1)*qk
      elegradfieldi(3,1) = elegradfieldi(3,1) + t2ik(3,1)*qk
c      elegradfieldi(2,1) = elegradfieldi(2,1) + t2ik(2,1)*qk
      elegradfieldi(2,2) = elegradfieldi(2,2) + t2ik(2,2)*qk
      elegradfieldi(3,2) = elegradfieldi(3,2) + t2ik(3,2)*qk
c      elegradfieldi(3,1) = elegradfieldi(3,1) + t2ik(3,1)*qk
c      elegradfieldi(3,2) = elegradfieldi(3,2) + t2ik(3,2)*qk
      elegradfieldi(3,3) = elegradfieldi(3,3) + t2ik(3,3)*qk
c
      elegradfieldk(1,1) = elegradfieldk(1,1) + t2ik(1,1)*qi
      elegradfieldk(2,1) = elegradfieldk(2,1) + t2ik(2,1)*qi
      elegradfieldk(3,1) = elegradfieldk(3,1) + t2ik(3,1)*qi
c      elegradfieldk(2,1) = elegradfieldk(2,1) + t2ik(2,1)*qi
      elegradfieldk(2,2) = elegradfieldk(2,2) + t2ik(2,2)*qi
      elegradfieldk(3,2) = elegradfieldk(3,2) + t2ik(3,2)*qi
c      elegradfieldk(3,1) = elegradfieldk(3,1) + t2ik(3,1)*qi
c      elegradfieldk(3,2) = elegradfieldk(3,2) + t2ik(3,2)*qi
      elegradfieldk(3,3) = elegradfieldk(3,3) + t2ik(3,3)*qi
c     from dipoles
      elegradfieldi(1,1) = elegradfieldi(1,1) + t3ik(1,1,1)*dkx + 
     &     t3ik(2,1,1)*dky + t3ik(3,1,1)*dkz
      elegradfieldi(2,1) = elegradfieldi(2,1) + t3ik(2,1,1)*dkx +
     &     t3ik(2,2,1)*dky + t3ik(3,2,1)*dkz
      elegradfieldi(3,1) = elegradfieldi(3,1) + t3ik(3,1,1)*dkx +
     &     t3ik(3,2,1)*dky + t3ik(3,3,1)*dkz
c      elegradfieldi(2,1) = elegradfieldi(2,1) + t3ik(2,1,1)*dkx +
c     &     t3ik(2,2,1)*dky + t3ik(3,2,1)*dkz
      elegradfieldi(2,2) = elegradfieldi(2,2) + t3ik(2,2,1)*dkx +
     &     t3ik(2,2,2)*dky + t3ik(3,2,2)*dkz
      elegradfieldi(3,2) = elegradfieldi(3,2) + t3ik(3,2,1)*dkx +
     &     t3ik(3,2,2)*dky + t3ik(3,3,2)*dkz
c      elegradfieldi(3,1) = elegradfieldi(3,1) + t3ik(3,1,1)*dkx +
c     &     t3ik(3,2,1)*dky + t3ik(3,3,1)*dkz
c      elegradfieldi(3,2) = elegradfieldi(3,2) + t3ik(3,2,1)*dkx +
c     &     t3ik(3,2,2)*dky + t3ik(3,3,2)*dkz
      elegradfieldi(3,3) = elegradfieldi(3,3) + t3ik(3,3,1)*dkx +
     &     t3ik(3,3,2)*dky + t3ik(3,3,3)*dkz
c     
      elegradfieldk(1,1) = elegradfieldk(1,1) - t3ik(1,1,1)*dix -
     &     t3ik(2,1,1)*diy - t3ik(3,1,1)*diz
      elegradfieldk(2,1) = elegradfieldk(2,1) - t3ik(2,1,1)*dix -
     &     t3ik(2,2,1)*diy - t3ik(3,2,1)*diz
      elegradfieldk(3,1) = elegradfieldk(3,1) - t3ik(3,1,1)*dix -
     &     t3ik(3,2,1)*diy - t3ik(3,3,1)*diz
c      elegradfieldk(2,1) = elegradfieldk(2,1) - t3ik(2,1,1)*dix -
c     &                 t3ik(2,2,1)*diy - t3ik(3,2,1)*diz
      elegradfieldk(2,2) = elegradfieldk(2,2) - t3ik(2,2,1)*dix -
     &     t3ik(2,2,2)*diy - t3ik(3,2,2)*diz
      elegradfieldk(3,2) = elegradfieldk(3,2) - t3ik(3,2,1)*dix -
     &     t3ik(3,2,2)*diy - t3ik(3,3,2)*diz
c      elegradfieldk(3,1) = elegradfieldk(3,1) - t3ik(3,1,1)*dix -
c     &     t3ik(3,2,1)*diy - t3ik(3,3,1)*diz
c      elegradfieldk(3,2) = elegradfieldk(3,2) - t3ik(3,2,1)*dix -
c     &     t3ik(3,2,2)*diy - t3ik(3,3,2)*diz
      elegradfieldk(3,3) = elegradfieldk(3,3) - t3ik(3,3,1)*dix -
     &                        t3ik(3,3,2)*diy - t3ik(3,3,3)*diz
c     from quadrupoles
      elegradfieldi(1,1) = elegradfieldi(1,1) + t4ik(1,1,1,1)*qkxx +
     &     t4ik(2,1,1,1)*qkxy + t4ik(3,1,1,1)*qkxz + 
     &     t4ik(2,2,1,1)*qkyy + t4ik(3,2,1,1)*qkyz + 
     &     t4ik(3,3,1,1)*qkzz
      elegradfieldi(2,1) = elegradfieldi(2,1) + t4ik(2,1,1,1)*qkxx +
     &     t4ik(2,2,1,1)*qkxy + t4ik(3,2,1,1)*qkxz +
     &     t4ik(2,2,2,1)*qkyy + t4ik(3,2,2,1)*qkyz +
     &     t4ik(3,3,2,1)*qkzz
      elegradfieldi(3,1) = elegradfieldi(3,1) + t4ik(3,1,1,1)*qkxx +
     &     t4ik(3,2,1,1)*qkxy + t4ik(3,3,1,1)*qkxz +
     &     t4ik(3,2,2,1)*qkyy + t4ik(3,3,2,1)*qkyz +
     &     t4ik(3,3,3,1)*qkzz
c      elegradfieldi(2,1) = elegradfieldi(2,1) + t4ik(2,1,1,1)*qkxx +
c     &     t4ik(2,2,1,1)*qkxy + t4ik(3,2,1,1)*qkxz +
c     &     t4ik(2,2,2,1)*qkyy + t4ik(3,2,2,1)*qkyz +
c     &     t4ik(3,3,2,1)*qkzz
      elegradfieldi(2,2) = elegradfieldi(2,2) + t4ik(2,2,1,1)*qkxx +
     &     t4ik(2,2,2,1)*qkxy + t4ik(3,2,2,1)*qkxz +
     &     t4ik(2,2,2,2)*qkyy + t4ik(3,2,2,2)*qkyz +
     &     t4ik(3,3,2,2)*qkzz
      elegradfieldi(3,2) = elegradfieldi(3,2) + t4ik(3,2,1,1)*qkxx +
     &     t4ik(3,2,2,1)*qkxy + t4ik(3,3,2,1)*qkxz +
     &     t4ik(3,2,2,2)*qkyy + t4ik(3,3,2,2)*qkyz +
     &     t4ik(3,3,3,2)*qkzz
c      elegradfieldi(3,1) = elegradfieldi(3,1) + t4ik(3,1,1,1)*qkxx +
c     &     t4ik(3,2,1,1)*qkxy + t4ik(3,3,1,1)*qkxz +
c     &     t4ik(3,2,2,1)*qkyy + t4ik(3,3,2,1)*qkyz +
c     &     t4ik(3,3,3,1)*qkzz
c      elegradfieldi(3,2) = elegradfieldi(3,2) + t4ik(3,2,1,1)*qkxx +
c     &     t4ik(3,2,2,1)*qkxy + t4ik(3,3,2,1)*qkxz +
c     &     t4ik(3,2,2,2)*qkyy + t4ik(3,3,2,2)*qkyz +
c     &     t4ik(3,3,3,2)*qkzz
      elegradfieldi(3,3) = elegradfieldi(3,3) + t4ik(3,3,1,1)*qkxx +
     &     t4ik(3,3,2,1)*qkxy + t4ik(3,3,3,1)*qkxz +
     &     t4ik(3,3,2,2)*qkyy + t4ik(3,3,3,2)*qkyz +
     &     t4ik(3,3,3,3)*qkzz
c     
      elegradfieldk(1,1) = elegradfieldk(1,1) + t4ik(1,1,1,1)*qixx +
     &     t4ik(2,1,1,1)*qixy + t4ik(3,1,1,1)*qixz +
     &     t4ik(2,2,1,1)*qiyy + t4ik(3,2,1,1)*qiyz +
     &     t4ik(3,3,1,1)*qizz
      elegradfieldk(2,1) = elegradfieldk(2,1) + t4ik(2,1,1,1)*qixx +
     &     t4ik(2,2,1,1)*qixy + t4ik(3,2,1,1)*qixz +
     &     t4ik(2,2,2,1)*qiyy + t4ik(3,2,2,1)*qiyz +
     &     t4ik(3,3,2,1)*qizz
      elegradfieldk(3,1) = elegradfieldk(3,1) + t4ik(3,1,1,1)*qixx +
     &     t4ik(3,2,1,1)*qixy + t4ik(3,3,1,1)*qixz +
     &     t4ik(3,2,2,1)*qiyy + t4ik(3,3,2,1)*qiyz +
     &     t4ik(3,3,3,1)*qizz
c      elegradfieldk(2,1) = elegradfieldk(2,1) + t4ik(2,1,1,1)*qixx +
c     &     t4ik(2,2,1,1)*qixy + t4ik(3,2,1,1)*qixz +
c     &     t4ik(2,2,2,1)*qiyy + t4ik(3,2,2,1)*qiyz +
c     &     t4ik(3,3,2,1)*qizz
      elegradfieldk(2,2) = elegradfieldk(2,2) + t4ik(2,2,1,1)*qixx +
     &     t4ik(2,2,2,1)*qixy + t4ik(3,2,2,1)*qixz +
     &     t4ik(2,2,2,2)*qiyy + t4ik(3,2,2,2)*qiyz +
     &     t4ik(3,3,2,2)*qizz
      elegradfieldk(3,2) = elegradfieldk(3,2) + t4ik(3,2,1,1)*qixx +
     &     t4ik(3,2,2,1)*qixy + t4ik(3,3,2,1)*qixz +
     &     t4ik(3,2,2,2)*qiyy + t4ik(3,3,2,2)*qiyz +
     &     t4ik(3,3,3,2)*qizz
c      elegradfieldk(3,1) = elegradfieldk(3,1) + t4ik(3,1,1,1)*qixx +
c     &     t4ik(3,2,1,1)*qixy + t4ik(3,3,1,1)*qixz +
c     &     t4ik(3,2,2,1)*qiyy + t4ik(3,3,2,1)*qiyz +
c     &     t4ik(3,3,3,1)*qizz
c      elegradfieldk(3,2) = elegradfieldk(3,2) + t4ik(3,2,1,1)*qixx +
c     &     t4ik(3,2,2,1)*qixy + t4ik(3,3,2,1)*qixz +
c     &     t4ik(3,2,2,2)*qiyy + t4ik(3,3,2,2)*qiyz +
c     &     t4ik(3,3,3,2)*qizz
      elegradfieldk(3,3) = elegradfieldk(3,3) + t4ik(3,3,1,1)*qixx +
     &     t4ik(3,3,2,1)*qixy + t4ik(3,3,3,1)*qixz +
     &     t4ik(3,3,2,2)*qiyy + t4ik(3,3,3,2)*qiyz +
     &     t4ik(3,3,3,3)*qizz
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine ugradfieldik  --  electric field gradient ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "ugradfieldik" computes the electric field gradient given
c     two induced dipoles and damped or undamped tmatrix
c
      subroutine ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &     gradfieldip,gradfieldkp)
      use polar
      implicit none
      integer i,k
      real*8 duix,duiy,duiz
      real*8 dukx,duky,dukz
      real*8 puix,puiy,puiz
      real*8 pukx,puky,pukz
      real*8 t3(3,3,3)
      real*8 gradfieldid(3,3),gradfieldkd(3,3)
      real*8 gradfieldip(3,3),gradfieldkp(3,3)
c
c     read in induced dipole values
c
      duix = uind(1,i)
      duiy = uind(2,i)
      duiz = uind(3,i)
      dukx = uind(1,k)
      duky = uind(2,k)
      dukz = uind(3,k)
      puix = uinp(1,i)
      puiy = uinp(2,i)
      puiz = uinp(3,i)
      pukx = uinp(1,k)
      puky = uinp(2,k)
      pukz = uinp(3,k)
c
c     calculate electric field gradient
c
      gradfieldid(1,1) = t3(1,1,1)*dukx + 
     &     t3(2,1,1)*duky + t3(3,1,1)*dukz
      gradfieldid(2,1) = t3(2,1,1)*dukx +
     &     t3(2,2,1)*duky + t3(3,2,1)*dukz
      gradfieldid(3,1) = t3(3,1,1)*dukx +
     &     t3(3,2,1)*duky + t3(3,3,1)*dukz
c      gradfieldid(2,1) = t3(2,1,1)*dukx +
c     &     t3(2,2,1)*duky + t3(3,2,1)*dukz
      gradfieldid(2,2) = t3(2,2,1)*dukx +
     &     t3(2,2,2)*duky + t3(3,2,2)*dukz
      gradfieldid(3,2) = t3(3,2,1)*dukx +
     &     t3(3,2,2)*duky + t3(3,3,2)*dukz
c      gradfieldid(3,1) = t3(3,1,1)*dukx +
c     &     t3(3,2,1)*duky + t3(3,3,1)*dukz
c      gradfieldid(3,2) = t3(3,2,1)*dukx +
c     &     t3(3,2,2)*duky + t3(3,3,2)*dukz
      gradfieldid(3,3) = t3(3,3,1)*dukx +
     &     t3(3,3,2)*duky + t3(3,3,3)*dukz
c     
      gradfieldkd(1,1) = -t3(1,1,1)*duix -
     &     t3(2,1,1)*duiy - t3(3,1,1)*duiz
      gradfieldkd(2,1) = -t3(2,1,1)*duix -
     &     t3(2,2,1)*duiy - t3(3,2,1)*duiz
      gradfieldkd(3,1) = -t3(3,1,1)*duix -
     &     t3(3,2,1)*duiy - t3(3,3,1)*duiz
c      gradfieldkd(2,1) = -t3(2,1,1)*duix -
c     &     t3(2,2,1)*duiy - t3(3,2,1)*duiz
      gradfieldkd(2,2) = -t3(2,2,1)*duix -
     &     t3(2,2,2)*duiy - t3(3,2,2)*duiz
      gradfieldkd(3,2) = -t3(3,2,1)*duix -
     &     t3(3,2,2)*duiy - t3(3,3,2)*duiz
c      gradfieldkd(3,1) = -t3(3,1,1)*duix -
c     &     t3(3,2,1)*duiy - t3(3,3,1)*duiz
c      gradfieldkd(3,2) = -t3(3,2,1)*duix -
c     &     t3(3,2,2)*duiy - t3(3,3,2)*duiz
      gradfieldkd(3,3) = -t3(3,3,1)*duix -
     &     t3(3,3,2)*duiy - t3(3,3,3)*duiz
c
c     p induced dipoles
c
      gradfieldip(1,1) = t3(1,1,1)*pukx +
     &     t3(2,1,1)*puky + t3(3,1,1)*pukz
      gradfieldip(2,1) = t3(2,1,1)*pukx +
     &     t3(2,2,1)*puky + t3(3,2,1)*pukz
      gradfieldip(3,1) = t3(3,1,1)*pukx +
     &     t3(3,2,1)*puky + t3(3,3,1)*pukz
c      gradfieldip(2,1) = t3(2,1,1)*pukx +
c     &     t3(2,2,1)*puky + t3(3,2,1)*pukz
      gradfieldip(2,2) = t3(2,2,1)*pukx +
     &     t3(2,2,2)*puky + t3(3,2,2)*pukz
      gradfieldip(3,2) = t3(3,2,1)*pukx +
     &     t3(3,2,2)*puky + t3(3,3,2)*pukz
c      gradfieldip(3,1) = t3(3,1,1)*pukx +
c     &     t3(3,2,1)*puky + t3(3,3,1)*pukz
c      gradfieldip(3,2) = t3(3,2,1)*pukx +
c     &     t3(3,2,2)*puky + t3(3,3,2)*pukz
      gradfieldip(3,3) = t3(3,3,1)*pukx +
     &     t3(3,3,2)*puky + t3(3,3,3)*pukz
c
      gradfieldkp(1,1) = -t3(1,1,1)*puix -
     &     t3(2,1,1)*puiy - t3(3,1,1)*puiz
      gradfieldkp(2,1) = -t3(2,1,1)*puix -
     &     t3(2,2,1)*puiy - t3(3,2,1)*puiz
      gradfieldkp(3,1) = -t3(3,1,1)*puix -
     &     t3(3,2,1)*puiy - t3(3,3,1)*puiz
c      gradfieldkp(2,1) = -t3(2,1,1)*puix -
c     &     t3(2,2,1)*puiy - t3(3,2,1)*puiz
      gradfieldkp(2,2) = -t3(2,2,1)*puix -
     &     t3(2,2,2)*puiy - t3(3,2,2)*puiz
      gradfieldkp(3,2) = -t3(3,2,1)*puix -
     &     t3(3,2,2)*puiy - t3(3,3,2)*puiz
c      gradfieldkp(3,1) = -t3(3,1,1)*puix -
c     &     t3(3,2,1)*puiy - t3(3,3,1)*puiz
c      gradfieldkp(3,2) = -t3(3,2,1)*puix -
c     &     t3(3,2,2)*puiy - t3(3,3,2)*puiz
      gradfieldkp(3,3) = -t3(3,3,1)*puix -
     &     t3(3,3,2)*puiy - t3(3,3,3)*puiz
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine hessfieldik  --  electric field hessian  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "hessfieldik" computes the electric field hessian given two multipoles
c     and damped or undamped tmatrix
c
      subroutine hessfieldik(i,k,t3,t4,t5,hessfieldi,hessfieldk)
      use mpole
      implicit none
      integer i,k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 t3(3,3,3),t4(3,3,3,3),t5(3,3,3,3,3)
      real*8 hessfieldi(3,3,3),hessfieldk(3,3,3)
c
c     read in multipole values
c
      ci = rpole(1,i)
      dix = rpole(2,i)
      diy = rpole(3,i)
      diz = rpole(4,i)
      qixx = rpole(5,i)
      qixy = rpole(6,i)*2.0d0
      qixz = rpole(7,i)*2.0d0
      qiyy = rpole(9,i)
      qiyz = rpole(10,i)*2.0d0
      qizz = rpole(13,i)
      ck = rpole(1,k)
      dkx = rpole(2,k)
      dky = rpole(3,k)
      dkz = rpole(4,k)
      qkxx = rpole(5,k)
      qkxy = rpole(6,k)*2.0d0
      qkxz = rpole(7,k)*2.0d0
      qkyy = rpole(9,k)
      qkyz = rpole(10,k)*2.0d0
      qkzz = rpole(13,k)
c
c     calculate electric field hessian
c
c     from charges
      hessfieldi(1,1,1) = -t3(1,1,1)*ck
      hessfieldi(2,1,1) = -t3(2,1,1)*ck
      hessfieldi(3,1,1) = -t3(3,1,1)*ck      
c      hessfieldi(2,1,1) = -t3(2,1,1)*ck
      hessfieldi(2,2,1) = -t3(2,2,1)*ck
      hessfieldi(3,2,1) = -t3(3,2,1)*ck
c      hessfieldi(3,1,1) = -t3(3,1,1)*ck
c      hessfieldi(3,2,1) = -t3(3,2,1)*ck
      hessfieldi(3,3,1) = -t3(3,3,1)*ck
c      hessfieldi(2,1,1) = -t3(2,1,1)*ck
c      hessfieldi(2,2,1) = -t3(2,2,1)*ck
c      hessfieldi(3,2,1) = -t3(3,2,1)*ck
c      hessfieldi(2,2,1) = -t3(2,2,1)*ck
      hessfieldi(2,2,2) = -t3(2,2,2)*ck
      hessfieldi(3,2,2) = -t3(3,2,2)*ck
c      hessfieldi(3,2,1) = -t3(3,2,1)*ck
c      hessfieldi(3,2,2) = -t3(3,2,2)*ck
      hessfieldi(3,3,2) = -t3(3,3,2)*ck
c      hessfieldi(3,1,1) = -t3(3,1,1)*ck
c      hessfieldi(3,2,1) = -t3(3,2,1)*ck
c      hessfieldi(3,3,1) = -t3(3,3,1)*ck
c      hessfieldi(3,2,1) = -t3(3,2,1)*ck
c      hessfieldi(3,2,2) = -t3(3,2,2)*ck
c      hessfieldi(3,3,2) = -t3(3,3,2)*ck
c      hessfieldi(3,3,1) = -t3(3,3,1)*ck
c      hessfieldi(3,3,2) = -t3(3,3,2)*ck
      hessfieldi(3,3,3) = -t3(3,3,3)*ck
c
      hessfieldk(1,1,1) = t3(1,1,1)*ci
      hessfieldk(2,1,1) = t3(2,1,1)*ci
      hessfieldk(3,1,1) = t3(3,1,1)*ci
c      hessfieldk(2,1,1) = t3(2,1,1)*ci
      hessfieldk(2,2,1) = t3(2,2,1)*ci
      hessfieldk(3,2,1) = t3(3,2,1)*ci
c      hessfieldk(3,1,1) = t3(3,1,1)*ci
c      hessfieldk(3,2,1) = t3(3,2,1)*ci
      hessfieldk(3,3,1) = t3(3,3,1)*ci
c      hessfieldk(2,1,1) = t3(2,1,1)*ci
c      hessfieldk(2,2,1) = t3(2,2,1)*ci
c      hessfieldk(3,2,1) = t3(3,2,1)*ci
c      hessfieldk(2,2,1) = t3(2,2,1)*ci
      hessfieldk(2,2,2) = t3(2,2,2)*ci
      hessfieldk(3,2,2) = t3(3,2,2)*ci
c      hessfieldk(3,2,1) = t3(3,2,1)*ci
c      hessfieldk(3,2,2) = t3(3,2,2)*ci
      hessfieldk(3,3,2) = t3(3,3,2)*ci
c      hessfieldk(3,1,1) = t3(3,1,1)*ci
c      hessfieldk(3,2,1) = t3(3,2,1)*ci
c      hessfieldk(3,3,1) = t3(3,3,1)*ci
c      hessfieldk(3,2,1) = t3(3,2,1)*ci
c      hessfieldk(3,2,2) = t3(3,2,2)*ci
c      hessfieldk(3,3,2) = t3(3,3,2)*ci
c      hessfieldk(3,3,1) = t3(3,3,1)*ci
c      hessfieldk(3,3,2) = t3(3,3,2)*ci
      hessfieldk(3,3,3) = t3(3,3,3)*ci
c
c     from dipoles
c
      hessfieldi(1,1,1) = hessfieldi(1,1,1) - t4(1,1,1,1)*dkx - 
     &     t4(2,1,1,1)*dky - t4(3,1,1,1)*dkz
      hessfieldi(2,1,1) = hessfieldi(2,1,1) - t4(2,1,1,1)*dkx -
     &     t4(2,2,1,1)*dky - t4(3,2,1,1)*dkz
      hessfieldi(3,1,1) = hessfieldi(3,1,1) - t4(3,1,1,1)*dkx -
     &     t4(3,2,1,1)*dky - t4(3,3,1,1)*dkz
c      hessfieldi(2,1,1) = hessfieldi(2,1,1) - t4(2,1,1,1)*dkx -
c     &     t4(2,2,1,1)*dky - t4(3,2,1,1)*dkz
      hessfieldi(2,2,1) = hessfieldi(2,2,1) - t4(2,2,1,1)*dkx -
     &     t4(2,2,2,1)*dky - t4(3,2,2,1)*dkz
      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t4(3,2,1,1)*dkx -
     &     t4(3,2,2,1)*dky - t4(3,3,2,1)*dkz
c      hessfieldi(3,1,1) = hessfieldi(3,1,1) - t4(3,1,1,1)*dkx -
c     &     t4(3,2,1,1)*dky - t4(3,3,1,1)*dkz
c      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t4(3,2,1,1)*dkx -
c     &     t4(3,2,2,1)*dky - t4(3,3,2,1)*dkz
      hessfieldi(3,3,1) = hessfieldi(3,3,1) - t4(3,3,1,1)*dkx -
     &     t4(3,3,2,1)*dky - t4(3,3,3,1)*dkz
c      hessfieldi(2,1,1) = hessfieldi(2,1,1) - t4(2,1,1,1)*dkx -
c     &     t4(2,2,1,1)*dky - t4(3,2,1,1)*dkz
c      hessfieldi(2,2,1) = hessfieldi(2,2,1) - t4(2,2,1,1)*dkx -
c     &     t4(2,2,2,1)*dky - t4(3,2,2,1)*dkz
c      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t4(3,2,1,1)*dkx -
c     &     t4(3,2,2,1)*dky - t4(3,3,2,1)*dkz
c      hessfieldi(2,2,1) = hessfieldi(2,2,1) - t4(2,2,1,1)*dkx -
c     &     t4(2,2,2,1)*dky - t4(3,2,2,1)*dkz
      hessfieldi(2,2,2) = hessfieldi(2,2,2) - t4(2,2,2,1)*dkx -
     &     t4(2,2,2,2)*dky - t4(3,2,2,2)*dkz
      hessfieldi(3,2,2) = hessfieldi(3,2,2) - t4(3,2,2,1)*dkx -
     &     t4(3,2,2,2)*dky - t4(3,3,2,2)*dkz
c      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t4(3,2,1,1)*dkx -
c     &     t4(3,2,2,1)*dky - t4(3,3,2,1)*dkz
c      hessfieldi(3,2,2) = hessfieldi(3,2,2) - t4(3,2,2,1)*dkx -
c     &     t4(3,2,2,2)*dky - t4(3,3,2,2)*dkz
      hessfieldi(3,3,2) = hessfieldi(3,3,2) - t4(3,3,2,1)*dkx -
     &     t4(3,3,2,2)*dky - t4(3,3,3,2)*dkz
c      hessfieldi(3,1,1) = hessfieldi(3,1,1) - t4(3,1,1,1)*dkx -
c     &     t4(3,2,1,1)*dky - t4(3,3,1,1)*dkz
c      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t4(3,2,1,1)*dkx -
c     &     t4(3,2,2,1)*dky - t4(3,3,2,1)*dkz
c      hessfieldi(3,3,1) = hessfieldi(3,3,1) - t4(3,3,1,1)*dkx -
c     &     t4(3,3,2,1)*dky - t4(3,3,3,1)*dkz
c      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t4(3,2,1,1)*dkx -
c     &     t4(3,2,2,1)*dky - t4(3,3,2,1)*dkz
c      hessfieldi(3,2,2) = hessfieldi(3,2,2) - t4(3,2,2,1)*dkx -
c     &     t4(3,2,2,2)*dky - t4(3,3,2,2)*dkz
c      hessfieldi(3,3,2) = hessfieldi(3,3,2) - t4(3,3,2,1)*dkx -
c     &     t4(3,3,2,2)*dky - t4(3,3,3,2)*dkz
c      hessfieldi(3,3,1) = hessfieldi(3,3,1) - t4(3,3,1,1)*dkx -
c     &     t4(3,3,2,1)*dky - t4(3,3,3,1)*dkz
c      hessfieldi(3,3,2) = hessfieldi(3,3,2) - t4(3,3,2,1)*dkx -
c     &     t4(3,3,2,2)*dky - t4(3,3,3,2)*dkz
      hessfieldi(3,3,3) = hessfieldi(3,3,3) - t4(3,3,3,1)*dkx -
     &     t4(3,3,3,2)*dky - t4(3,3,3,3)*dkz
c
      hessfieldk(1,1,1) = hessfieldk(1,1,1) - t4(1,1,1,1)*dix - 
     &     t4(2,1,1,1)*diy - t4(3,1,1,1)*diz
      hessfieldk(2,1,1) = hessfieldk(2,1,1) - t4(2,1,1,1)*dix -
     &     t4(2,2,1,1)*diy - t4(3,2,1,1)*diz
      hessfieldk(3,1,1) = hessfieldk(3,1,1) - t4(3,1,1,1)*dix -
     &     t4(3,2,1,1)*diy - t4(3,3,1,1)*diz
c      hessfieldk(2,1,1) = hessfieldk(2,1,1) - t4(2,1,1,1)*dix -
c     &     t4(2,2,1,1)*diy - t4(3,2,1,1)*diz
      hessfieldk(2,2,1) = hessfieldk(2,2,1) - t4(2,2,1,1)*dix -
     &     t4(2,2,2,1)*diy - t4(3,2,2,1)*diz
      hessfieldk(3,2,1) = hessfieldk(3,2,1) - t4(3,2,1,1)*dix -
     &     t4(3,2,2,1)*diy - t4(3,3,2,1)*diz
c      hessfieldk(3,1,1) = hessfieldk(3,1,1) - t4(3,1,1,1)*dix -
c     &     t4(3,2,1,1)*diy - t4(3,3,1,1)*diz
c      hessfieldk(3,2,1) = hessfieldk(3,2,1) - t4(3,2,1,1)*dix -
c     &     t4(3,2,2,1)*diy - t4(3,3,2,1)*diz
      hessfieldk(3,3,1) = hessfieldk(3,3,1) - t4(3,3,1,1)*dix -
     &     t4(3,3,2,1)*diy - t4(3,3,3,1)*diz
c      hessfieldk(2,1,1) = hessfieldk(2,1,1) - t4(2,1,1,1)*dix -
c     &     t4(2,2,1,1)*diy - t4(3,2,1,1)*diz
c      hessfieldk(2,2,1) = hessfieldk(2,2,1) - t4(2,2,1,1)*dix -
c     &     t4(2,2,2,1)*diy - t4(3,2,2,1)*diz
c      hessfieldk(3,2,1) = hessfieldk(3,2,1) - t4(3,2,1,1)*dix -
c     &     t4(3,2,2,1)*diy - t4(3,3,2,1)*diz
c      hessfieldk(2,2,1) = hessfieldk(2,2,1) - t4(2,2,1,1)*dix -
c     &     t4(2,2,2,1)*diy - t4(3,2,2,1)*diz
      hessfieldk(2,2,2) = hessfieldk(2,2,2) - t4(2,2,2,1)*dix -
     &     t4(2,2,2,2)*diy - t4(3,2,2,2)*diz
      hessfieldk(3,2,2) = hessfieldk(3,2,2) - t4(3,2,2,1)*dix -
     &     t4(3,2,2,2)*diy - t4(3,3,2,2)*diz
c      hessfieldk(3,2,1) = hessfieldk(3,2,1) - t4(3,2,1,1)*dix -
c     &     t4(3,2,2,1)*diy - t4(3,3,2,1)*diz
c      hessfieldk(3,2,2) = hessfieldk(3,2,2) - t4(3,2,2,1)*dix -
c     &     t4(3,2,2,2)*diy - t4(3,3,2,2)*diz
      hessfieldk(3,3,2) = hessfieldk(3,3,2) - t4(3,3,2,1)*dix -
     &     t4(3,3,2,2)*diy - t4(3,3,3,2)*diz
c      hessfieldk(3,1,1) = hessfieldk(3,1,1) - t4(3,1,1,1)*dix -
c     &     t4(3,2,1,1)*diy - t4(3,3,1,1)*diz
c      hessfieldk(3,2,1) = hessfieldk(3,2,1) - t4(3,2,1,1)*dix -
c     &     t4(3,2,2,1)*diy - t4(3,3,2,1)*diz
c      hessfieldk(3,3,1) = hessfieldk(3,3,1) - t4(3,3,1,1)*dix -
c     &     t4(3,3,2,1)*diy - t4(3,3,3,1)*diz
c      hessfieldk(3,2,1) = hessfieldk(3,2,1) - t4(3,2,1,1)*dix -
c     &     t4(3,2,2,1)*diy - t4(3,3,2,1)*diz
c      hessfieldk(3,2,2) = hessfieldk(3,2,2) - t4(3,2,2,1)*dix -
c     &     t4(3,2,2,2)*diy - t4(3,3,2,2)*diz
c      hessfieldk(3,3,2) = hessfieldk(3,3,2) - t4(3,3,2,1)*dix -
c     &     t4(3,3,2,2)*diy - t4(3,3,3,2)*diz
c      hessfieldk(3,3,1) = hessfieldk(3,3,1) - t4(3,3,1,1)*dix -
c     &     t4(3,3,2,1)*diy - t4(3,3,3,1)*diz
c      hessfieldk(3,3,2) = hessfieldk(3,3,2) - t4(3,3,2,1)*dix -
c     &     t4(3,3,2,2)*diy - t4(3,3,3,2)*diz
      hessfieldk(3,3,3) = hessfieldk(3,3,3) - t4(3,3,3,1)*dix -
     &     t4(3,3,3,2)*diy - t4(3,3,3,3)*diz
c
c     quadrupoles
c
      hessfieldi(1,1,1) = hessfieldi(1,1,1) - t5(1,1,1,1,1)*qkxx - 
     &     t5(2,1,1,1,1)*qkxy - t5(3,1,1,1,1)*qkxz - t5(2,2,1,1,1)*qkyy-
     &     t5(3,2,1,1,1)*qkyz - t5(3,3,1,1,1)*qkzz
      hessfieldi(2,1,1) = hessfieldi(2,1,1) - t5(2,1,1,1,1)*qkxx -
     &     t5(2,2,1,1,1)*qkxy - t5(3,2,1,1,1)*qkxz - t5(2,2,2,1,1)*qkyy-
     &     t5(3,2,2,1,1)*qkyz - t5(3,3,2,1,1)*qkzz
      hessfieldi(3,1,1) = hessfieldi(3,1,1) - t5(3,1,1,1,1)*qkxx -
     &     t5(3,2,1,1,1)*qkxy - t5(3,3,1,1,1)*qkxz - t5(3,2,2,1,1)*qkyy-
     &     t5(3,3,2,1,1)*qkyz - t5(3,3,3,1,1)*qkzz
c      hessfieldi(2,1,1) = hessfieldi(2,1,1) - t5(2,1,1,1,1)*qkxx -
c     &     t5(2,2,1,1,1)*qkxy - t5(3,2,1,1,1)*qkxz - t5(2,2,2,1,1)*qkyy-
c     &     t5(3,2,2,1,1)*qkyz - t5(3,3,2,1,1)*qkzz
      hessfieldi(2,2,1) = hessfieldi(2,2,1) - t5(2,2,1,1,1)*qkxx -
     &     t5(2,2,2,1,1)*qkxy - t5(3,2,2,1,1)*qkxz - t5(2,2,2,2,1)*qkyy-
     &     t5(3,2,2,2,1)*qkyz - t5(3,3,2,2,1)*qkzz
      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t5(3,2,1,1,1)*qkxx -
     &     t5(3,2,2,1,1)*qkxy - t5(3,3,2,1,1)*qkxz - t5(3,2,2,2,1)*qkyy-
     &     t5(3,3,2,2,1)*qkyz - t5(3,3,3,2,1)*qkzz
c      hessfieldi(3,1,1) = hessfieldi(3,1,1) - t5(3,1,1,1,1)*qkxx -
c     &     t5(3,2,1,1,1)*qkxy - t5(3,3,1,1,1)*qkxz - t5(3,2,2,1,1)*qkyy-
c     &     t5(3,3,2,1,1)*qkyz - t5(3,3,3,1,1)*qkzz
c      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t5(3,2,1,1,1)*qkxx -
c     &     t5(3,2,2,1,1)*qkxy - t5(3,3,2,1,1)*qkxz - t5(3,2,2,2,1)*qkyy-
c     &     t5(3,3,2,2,1)*qkyz - t5(3,3,3,2,1)*qkzz
      hessfieldi(3,3,1) = hessfieldi(3,3,1) - t5(3,3,1,1,1)*qkxx -
     &     t5(3,3,2,1,1)*qkxy - t5(3,3,3,1,1)*qkxz - t5(3,3,2,2,1)*qkyy-
     &     t5(3,3,3,2,1)*qkyz - t5(3,3,3,3,1)*qkzz
c      hessfieldi(2,1,1) = hessfieldi(2,1,1) - t5(2,1,1,1,1)*qkxx -
c     &     t5(2,2,1,1,1)*qkxy - t5(3,2,1,1,1)*qkxz - t5(2,2,2,1,1)*qkyy-
c     &     t5(3,2,2,1,1)*qkyz - t5(3,3,2,1,1)*qkzz
c      hessfieldi(2,2,1) = hessfieldi(2,2,1) - t5(2,2,1,1,1)*qkxx -
c     &     t5(2,2,2,1,1)*qkxy - t5(3,2,2,1,1)*qkxz - t5(2,2,2,2,1)*qkyy-
c     &     t5(3,2,2,2,1)*qkyz - t5(3,3,2,2,1)*qkzz
c      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t5(3,2,1,1,1)*qkxx -
c     &     t5(3,2,2,1,1)*qkxy - t5(3,3,2,1,1)*qkxz - t5(3,2,2,2,1)*qkyy-
c     &     t5(3,3,2,2,1)*qkyz - t5(3,3,3,2,1)*qkzz
c      hessfieldi(2,2,1) = hessfieldi(2,2,1) - t5(2,2,1,1,1)*qkxx -
c     &     t5(2,2,2,1,1)*qkxy - t5(3,2,2,1,1)*qkxz - t5(2,2,2,2,1)*qkyy-
c     &     t5(3,2,2,2,1)*qkyz - t5(3,3,2,2,1)*qkzz
      hessfieldi(2,2,2) = hessfieldi(2,2,2) - t5(2,2,2,1,1)*qkxx -
     &     t5(2,2,2,2,1)*qkxy - t5(3,2,2,2,1)*qkxz - t5(2,2,2,2,2)*qkyy-
     &     t5(3,2,2,2,2)*qkyz - t5(3,3,2,2,2)*qkzz
      hessfieldi(3,2,2) = hessfieldi(3,2,2) - t5(3,2,2,1,1)*qkxx -
     &     t5(3,2,2,2,1)*qkxy - t5(3,3,2,2,1)*qkxz - t5(3,2,2,2,2)*qkyy-
     &     t5(3,3,2,2,2)*qkyz - t5(3,3,3,2,2)*qkzz
c      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t5(3,2,1,1,1)*qkxx -
c     &     t5(3,2,2,1,1)*qkxy - t5(3,3,2,1,1)*qkxz - t5(3,2,2,2,1)*qkyy-
c     &     t5(3,3,2,2,1)*qkyz - t5(3,3,3,2,1)*qkzz
c      hessfieldi(3,2,2) = hessfieldi(3,2,2) - t5(3,2,2,1,1)*qkxx -
c     &     t5(3,2,2,2,1)*qkxy - t5(3,3,2,2,1)*qkxz - t5(3,2,2,2,2)*qkyy-
c     &     t5(3,3,2,2,2)*qkyz - t5(3,3,3,2,2)*qkzz
      hessfieldi(3,3,2) = hessfieldi(3,3,2) - t5(3,3,2,1,1)*qkxx -
     &     t5(3,3,2,2,1)*qkxy - t5(3,3,3,2,1)*qkxz - t5(3,3,2,2,2)*qkyy-
     &     t5(3,3,3,2,2)*qkyz - t5(3,3,3,3,2)*qkzz
c      hessfieldi(3,1,1) = hessfieldi(3,1,1) - t5(3,1,1,1,1)*qkxx -
c     &     t5(3,2,1,1,1)*qkxy - t5(3,3,1,1,1)*qkxz - t5(3,2,2,1,1)*qkyy-
c     &     t5(3,3,2,1,1)*qkyz - t5(3,3,3,1,1)*qkzz
c      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t5(3,2,1,1,1)*qkxx -
c     &     t5(3,2,2,1,1)*qkxy - t5(3,3,2,1,1)*qkxz - t5(3,2,2,2,1)*qkyy-
c     &     t5(3,3,2,2,1)*qkyz - t5(3,3,3,2,1)*qkzz
c      hessfieldi(3,3,1) = hessfieldi(3,3,1) - t5(3,3,1,1,1)*qkxx -
c     &     t5(3,3,2,1,1)*qkxy - t5(3,3,3,1,1)*qkxz - t5(3,3,2,2,1)*qkyy-
c     &     t5(3,3,3,2,1)*qkyz - t5(3,3,3,3,1)*qkzz
c      hessfieldi(3,2,1) = hessfieldi(3,2,1) - t5(3,2,1,1,1)*qkxx -
c     &     t5(3,2,2,1,1)*qkxy - t5(3,3,2,1,1)*qkxz - t5(3,2,2,2,1)*qkyy-
c     &     t5(3,3,2,2,1)*qkyz - t5(3,3,3,2,1)*qkzz
c      hessfieldi(3,2,2) = hessfieldi(3,2,2) - t5(3,2,2,1,1)*qkxx -
c     &     t5(3,2,2,2,1)*qkxy - t5(3,3,2,2,1)*qkxz - t5(3,2,2,2,2)*qkyy-
c     &     t5(3,3,2,2,2)*qkyz - t5(3,3,3,2,2)*qkzz
c      hessfieldi(3,3,2) = hessfieldi(3,3,2) - t5(3,3,2,1,1)*qkxx -
c     &     t5(3,3,2,2,1)*qkxy - t5(3,3,3,2,1)*qkxz - t5(3,3,2,2,2)*qkyy-
c     &     t5(3,3,3,2,2)*qkyz - t5(3,3,3,3,2)*qkzz
c      hessfieldi(3,3,1) = hessfieldi(3,3,1) - t5(3,3,1,1,1)*qkxx -
c     &     t5(3,3,2,1,1)*qkxy - t5(3,3,3,1,1)*qkxz - t5(3,3,2,2,1)*qkyy-
c     &     t5(3,3,3,2,1)*qkyz - t5(3,3,3,3,1)*qkzz
c      hessfieldi(3,3,2) = hessfieldi(3,3,2) - t5(3,3,2,1,1)*qkxx -
c     &     t5(3,3,2,2,1)*qkxy - t5(3,3,3,2,1)*qkxz - t5(3,3,2,2,2)*qkyy-
c     &     t5(3,3,3,2,2)*qkyz - t5(3,3,3,3,2)*qkzz
      hessfieldi(3,3,3) = hessfieldi(3,3,3) -  t5(3,3,3,1,1)*qkxx -
     &     t5(3,3,3,2,1)*qkxy - t5(3,3,3,3,1)*qkxz - t5(3,3,3,2,2)*qkyy-
     &     t5(3,3,3,3,2)*qkyz - t5(3,3,3,3,3)*qkzz
c
      hessfieldk(1,1,1) = hessfieldk(1,1,1) + t5(1,1,1,1,1)*qixx + 
     &     t5(2,1,1,1,1)*qixy + t5(3,1,1,1,1)*qixz + t5(2,2,1,1,1)*qiyy+
     &     t5(3,2,1,1,1)*qiyz + t5(3,3,1,1,1)*qizz
      hessfieldk(2,1,1) = hessfieldk(2,1,1) + t5(2,1,1,1,1)*qixx +
     &     t5(2,2,1,1,1)*qixy + t5(3,2,1,1,1)*qixz + t5(2,2,2,1,1)*qiyy+
     &     t5(3,2,2,1,1)*qiyz + t5(3,3,2,1,1)*qizz
      hessfieldk(3,1,1) = hessfieldk(3,1,1) + t5(3,1,1,1,1)*qixx +
     &     t5(3,2,1,1,1)*qixy + t5(3,3,1,1,1)*qixz + t5(3,2,2,1,1)*qiyy+
     &     t5(3,3,2,1,1)*qiyz + t5(3,3,3,1,1)*qizz
c      hessfieldk(2,1,1) = hessfieldk(2,1,1) + t5(2,1,1,1,1)*qixx +
c     &     t5(2,2,1,1,1)*qixy + t5(3,2,1,1,1)*qixz + t5(2,2,2,1,1)*qiyy+
c     &     t5(3,2,2,1,1)*qiyz + t5(3,3,2,1,1)*qizz
      hessfieldk(2,2,1) = hessfieldk(2,2,1) + t5(2,2,1,1,1)*qixx +
     &     t5(2,2,2,1,1)*qixy + t5(3,2,2,1,1)*qixz + t5(2,2,2,2,1)*qiyy+
     &     t5(3,2,2,2,1)*qiyz + t5(3,3,2,2,1)*qizz
      hessfieldk(3,2,1) = hessfieldk(3,2,1) + t5(3,2,1,1,1)*qixx +
     &     t5(3,2,2,1,1)*qixy + t5(3,3,2,1,1)*qixz + t5(3,2,2,2,1)*qiyy+
     &     t5(3,3,2,2,1)*qiyz + t5(3,3,3,2,1)*qizz
c      hessfieldk(3,1,1) = hessfieldk(3,1,1) + t5(3,1,1,1,1)*qixx +
c     &     t5(3,2,1,1,1)*qixy + t5(3,3,1,1,1)*qixz + t5(3,2,2,1,1)*qiyy+
c     &     t5(3,3,2,1,1)*qiyz + t5(3,3,3,1,1)*qizz
c      hessfieldk(3,2,1) = hessfieldk(3,2,1) + t5(3,2,1,1,1)*qixx +
c     &     t5(3,2,2,1,1)*qixy + t5(3,3,2,1,1)*qixz + t5(3,2,2,2,1)*qiyy+
c     &     t5(3,3,2,2,1)*qiyz + t5(3,3,3,2,1)*qizz
      hessfieldk(3,3,1) = hessfieldk(3,3,1) + t5(3,3,1,1,1)*qixx +
     &     t5(3,3,2,1,1)*qixy + t5(3,3,3,1,1)*qixz + t5(3,3,2,2,1)*qiyy+
     &     t5(3,3,3,2,1)*qiyz + t5(3,3,3,3,1)*qizz
c      hessfieldk(2,1,1) = hessfieldk(2,1,1) + t5(2,1,1,1,1)*qixx +
c     &     t5(2,2,1,1,1)*qixy + t5(3,2,1,1,1)*qixz + t5(2,2,2,1,1)*qiyy+
c     &     t5(3,2,2,1,1)*qiyz + t5(3,3,2,1,1)*qizz
c      hessfieldk(2,2,1) = hessfieldk(2,2,1) + t5(2,2,1,1,1)*qixx +
c     &     t5(2,2,2,1,1)*qixy + t5(3,2,2,1,1)*qixz + t5(2,2,2,2,1)*qiyy+
c     &     t5(3,2,2,2,1)*qiyz + t5(3,3,2,2,1)*qizz
c      hessfieldk(3,2,1) = hessfieldk(3,2,1) + t5(3,2,1,1,1)*qixx +
c     &     t5(3,2,2,1,1)*qixy + t5(3,3,2,1,1)*qixz + t5(3,2,2,2,1)*qiyy+
c     &     t5(3,3,2,2,1)*qiyz + t5(3,3,3,2,1)*qizz
c      hessfieldk(2,2,1) = hessfieldk(2,2,1) + t5(2,2,1,1,1)*qixx +
c     &     t5(2,2,2,1,1)*qixy + t5(3,2,2,1,1)*qixz + t5(2,2,2,2,1)*qiyy+
c     &     t5(3,2,2,2,1)*qiyz + t5(3,3,2,2,1)*qizz
      hessfieldk(2,2,2) = hessfieldk(2,2,2) + t5(2,2,2,1,1)*qixx +
     &     t5(2,2,2,2,1)*qixy + t5(3,2,2,2,1)*qixz + t5(2,2,2,2,2)*qiyy+
     &     t5(3,2,2,2,2)*qiyz + t5(3,3,2,2,2)*qizz
      hessfieldk(3,2,2) = hessfieldk(3,2,2) + t5(3,2,2,1,1)*qixx +
     &     t5(3,2,2,2,1)*qixy + t5(3,3,2,2,1)*qixz + t5(3,2,2,2,2)*qiyy+
     &     t5(3,3,2,2,2)*qiyz + t5(3,3,3,2,2)*qizz
c      hessfieldk(3,2,1) = hessfieldk(3,2,1) + t5(3,2,1,1,1)*qixx +
c     &     t5(3,2,2,1,1)*qixy + t5(3,3,2,1,1)*qixz + t5(3,2,2,2,1)*qiyy+
c     &     t5(3,3,2,2,1)*qiyz + t5(3,3,3,2,1)*qizz
c      hessfieldk(3,2,2) = hessfieldk(3,2,2) + t5(3,2,2,1,1)*qixx +
c     &     t5(3,2,2,2,1)*qixy + t5(3,3,2,2,1)*qixz + t5(3,2,2,2,2)*qiyy+
c     &     t5(3,3,2,2,2)*qiyz + t5(3,3,3,2,2)*qizz
      hessfieldk(3,3,2) = hessfieldk(3,3,2) + t5(3,3,2,1,1)*qixx +
     &     t5(3,3,2,2,1)*qixy + t5(3,3,3,2,1)*qixz + t5(3,3,2,2,2)*qiyy+
     &     t5(3,3,3,2,2)*qiyz + t5(3,3,3,3,2)*qizz
c      hessfieldk(3,1,1) = hessfieldk(3,1,1) + t5(3,1,1,1,1)*qixx +
c     &     t5(3,2,1,1,1)*qixy + t5(3,3,1,1,1)*qixz + t5(3,2,2,1,1)*qiyy+
c     &     t5(3,3,2,1,1)*qiyz + t5(3,3,3,1,1)*qizz
c      hessfieldk(3,2,1) = hessfieldk(3,2,1) + t5(3,2,1,1,1)*qixx +
c     &     t5(3,2,2,1,1)*qixy + t5(3,3,2,1,1)*qixz + t5(3,2,2,2,1)*qiyy+
c     &     t5(3,3,2,2,1)*qiyz + t5(3,3,3,2,1)*qizz
c      hessfieldk(3,3,1) = hessfieldk(3,3,1) + t5(3,3,1,1,1)*qixx +
c     &     t5(3,3,2,1,1)*qixy + t5(3,3,3,1,1)*qixz + t5(3,3,2,2,1)*qiyy+
c     &     t5(3,3,3,2,1)*qiyz + t5(3,3,3,3,1)*qizz
c      hessfieldk(3,2,1) = hessfieldk(3,2,1) + t5(3,2,1,1,1)*qixx +
c     &     t5(3,2,2,1,1)*qixy + t5(3,3,2,1,1)*qixz + t5(3,2,2,2,1)*qiyy+
c     &     t5(3,3,2,2,1)*qiyz + t5(3,3,3,2,1)*qizz
c      hessfieldk(3,2,2) = hessfieldk(3,2,2) + t5(3,2,2,1,1)*qixx +
c     &     t5(3,2,2,2,1)*qixy + t5(3,3,2,2,1)*qixz + t5(3,2,2,2,2)*qiyy+
c     &     t5(3,3,2,2,2)*qiyz + t5(3,3,3,2,2)*qizz
c      hessfieldk(3,3,2) = hessfieldk(3,3,2) + t5(3,3,2,1,1)*qixx +
c     &     t5(3,3,2,2,1)*qixy + t5(3,3,3,2,1)*qixz + t5(3,3,2,2,2)*qiyy+
c     &     t5(3,3,3,2,2)*qiyz + t5(3,3,3,3,2)*qizz
c      hessfieldk(3,3,1) = hessfieldk(3,3,1) + t5(3,3,1,1,1)*qixx +
c     &     t5(3,3,2,1,1)*qixy + t5(3,3,3,1,1)*qixz + t5(3,3,2,2,1)*qiyy+
c     &     t5(3,3,3,2,1)*qiyz + t5(3,3,3,3,1)*qizz
c      hessfieldk(3,3,2) = hessfieldk(3,3,2) + t5(3,3,2,1,1)*qixx +
c     &     t5(3,3,2,2,1)*qixy + t5(3,3,3,2,1)*qixz + t5(3,3,2,2,2)*qiyy+
c     &     t5(3,3,3,2,2)*qiyz + t5(3,3,3,3,2)*qizz
      hessfieldk(3,3,3) = hessfieldk(3,3,3) +  t5(3,3,3,1,1)*qixx +
     &     t5(3,3,3,2,1)*qixy + t5(3,3,3,3,1)*qixz + t5(3,3,3,2,2)*qiyy+
     &     t5(3,3,3,3,2)*qiyz + t5(3,3,3,3,3)*qizz
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine cp_hessfieldik  --  electric field hessian  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "cp_hessfieldik" computes the electric field hessian given two multipoles
c     and damped or undamped tmatrix
c
      subroutine cp_hessfieldik(i,k,
     &     t3,t4,t5,t3i,t4i,t5i,t3k,t4k,t5k,t3ik,t4ik,t5ik,
     &     elehessfieldi,elehessfieldk)
      use atomid
      use chgpen
      use mpole
      implicit none
      integer i,k
      real*8 ci,zi,qi
      real*8 dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,zk,qk
      real*8 dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 t3(3,3,3),t4(3,3,3,3),t5(3,3,3,3,3)
      real*8 t3i(3,3,3),t4i(3,3,3,3),t5i(3,3,3,3,3)
      real*8 t3k(3,3,3),t4k(3,3,3,3),t5k(3,3,3,3,3)
      real*8 t3ik(3,3,3),t4ik(3,3,3,3),t5ik(3,3,3,3,3)
      real*8 elehessfieldi(3,3,3),elehessfieldk(3,3,3)
c
c     read in multipole values
c
      ci = rpole(1,i)
      dix = rpole(2,i)
      diy = rpole(3,i)
      diz = rpole(4,i)
      qixx = rpole(5,i)
      qixy = rpole(6,i)*2.0d0
      qixz = rpole(7,i)*2.0d0
      qiyy = rpole(9,i)
      qiyz = rpole(10,i)*2.0d0
      qizz = rpole(13,i)
      ck = rpole(1,k)
      dkx = rpole(2,k)
      dky = rpole(3,k)
      dkz = rpole(4,k)
      qkxx = rpole(5,k)
      qkxy = rpole(6,k)*2.0d0
      qkxz = rpole(7,k)*2.0d0
      qkyy = rpole(9,k)
      qkyz = rpole(10,k)*2.0d0
      qkzz = rpole(13,k)
c
c     split nuclear and electronic charges
c
      zi = atomic(i)
      zk = atomic(k)
      if (num_ele .eq. "VALENCE") then
         if (atomic(i) .gt. 2)  zi = zi - 2.0d0
         if (atomic(i) .gt. 10)  zi = zi - 8.0d0
         if (atomic(i) .gt. 18)  zi = zi - 8.0d0
         if (atomic(i) .gt. 20)  zi = zi - 10.0d0
         if (atomic(k) .gt. 2)  zk = zk - 2.0d0
         if (atomic(k) .gt. 10)  zk = zk - 8.0d0
         if (atomic(k) .gt. 18)  zk = zk - 8.0d0
         if (atomic(k) .gt. 20)  zk = zk - 10.0d0
      end if
      qi = ci - zi
      qk = ck - zk
c
c     calculate electric field hessian at electrons
c
c     from nuclei
      elehessfieldi(1,1,1) = -t3i(1,1,1)*zk
      elehessfieldi(2,1,1) = -t3i(2,1,1)*zk
      elehessfieldi(3,1,1) = -t3i(3,1,1)*zk      
c      elehessfieldi(2,1,1) = -t3i(2,1,1)*zk
      elehessfieldi(2,2,1) = -t3i(2,2,1)*zk
      elehessfieldi(3,2,1) = -t3i(3,2,1)*zk
c      elehessfieldi(3,1,1) = -t3i(3,1,1)*zk
c      elehessfieldi(3,2,1) = -t3i(3,2,1)*zk
      elehessfieldi(3,3,1) = -t3i(3,3,1)*zk
c      elehessfieldi(2,1,1) = -t3i(2,1,1)*zk
c      elehessfieldi(2,2,1) = -t3i(2,2,1)*zk
c      elehessfieldi(3,2,1) = -t3i(3,2,1)*zk
c      elehessfieldi(2,2,1) = -t3i(2,2,1)*zk
      elehessfieldi(2,2,2) = -t3i(2,2,2)*zk
      elehessfieldi(3,2,2) = -t3i(3,2,2)*zk
c      elehessfieldi(3,2,1) = -t3i(3,2,1)*zk
c      elehessfieldi(3,2,2) = -t3i(3,2,2)*zk
      elehessfieldi(3,3,2) = -t3i(3,3,2)*zk
c      elehessfieldi(3,1,1) = -t3i(3,1,1)*zk
c      elehessfieldi(3,2,1) = -t3i(3,2,1)*zk
c      elehessfieldi(3,3,1) = -t3i(3,3,1)*zk
c      elehessfieldi(3,2,1) = -t3i(3,2,1)*zk
c      elehessfieldi(3,2,2) = -t3i(3,2,2)*zk
c      elehessfieldi(3,3,2) = -t3i(3,3,2)*zk
c      elehessfieldi(3,3,1) = -t3i(3,3,1)*zk
c      elehessfieldi(3,3,2) = -t3i(3,3,2)*zk
      elehessfieldi(3,3,3) = -t3i(3,3,3)*zk
c
      elehessfieldk(1,1,1) = t3k(1,1,1)*zi
      elehessfieldk(2,1,1) = t3k(2,1,1)*zi
      elehessfieldk(3,1,1) = t3k(3,1,1)*zi
c      elehessfieldk(2,1,1) = t3k(2,1,1)*zi
      elehessfieldk(2,2,1) = t3k(2,2,1)*zi
      elehessfieldk(3,2,1) = t3k(3,2,1)*zi
c      elehessfieldk(3,1,1) = t3k(3,1,1)*zi
c      elehessfieldk(3,2,1) = t3k(3,2,1)*zi
      elehessfieldk(3,3,1) = t3k(3,3,1)*zi
c      elehessfieldk(2,1,1) = t3k(2,1,1)*zi
c      elehessfieldk(2,2,1) = t3k(2,2,1)*zi
c      elehessfieldk(3,2,1) = t3k(3,2,1)*zi
c      elehessfieldk(2,2,1) = t3k(2,2,1)*zi
      elehessfieldk(2,2,2) = t3k(2,2,2)*zi
      elehessfieldk(3,2,2) = t3k(3,2,2)*zi
c      elehessfieldk(3,2,1) = t3k(3,2,1)*zi
c      elehessfieldk(3,2,2) = t3k(3,2,2)*zi
      elehessfieldk(3,3,2) = t3k(3,3,2)*zi
c      elehessfieldk(3,1,1) = t3k(3,1,1)*zi
c      elehessfieldk(3,2,1) = t3k(3,2,1)*zi
c      elehessfieldk(3,3,1) = t3k(3,3,1)*zi
c      elehessfieldk(3,2,1) = t3k(3,2,1)*zi
c      elehessfieldk(3,2,2) = t3k(3,2,2)*zi
c      elehessfieldk(3,3,2) = t3k(3,3,2)*zi
c      elehessfieldk(3,3,1) = t3k(3,3,1)*zi
c      elehessfieldk(3,3,2) = t3k(3,3,2)*zi
      elehessfieldk(3,3,3) = t3k(3,3,3)*zi
c
c     from electrons
c
      elehessfieldi(1,1,1) = elehessfieldi(1,1,1) - t3ik(1,1,1)*qk
      elehessfieldi(2,1,1) = elehessfieldi(2,1,1) - t3ik(2,1,1)*qk
      elehessfieldi(3,1,1) = elehessfieldi(3,1,1) - t3ik(3,1,1)*qk      
c      elehessfieldi(2,1,1) = elehessfieldi(2,1,1) - t3ik(2,1,1)*qk
      elehessfieldi(2,2,1) = elehessfieldi(2,2,1) - t3ik(2,2,1)*qk
      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t3ik(3,2,1)*qk
c      elehessfieldi(3,1,1) = elehessfieldi(3,1,1) - t3ik(3,1,1)*qk
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t3ik(3,2,1)*qk
      elehessfieldi(3,3,1) = elehessfieldi(3,3,1) - t3ik(3,3,1)*qk
c      elehessfieldi(2,1,1) = elehessfieldi(2,1,1) - t3ik(2,1,1)*qk
c      elehessfieldi(2,2,1) = elehessfieldi(2,2,1) - t3ik(2,2,1)*qk
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t3ik(3,2,1)*qk
c      elehessfieldi(2,2,1) = elehessfieldi(2,2,1) - t3ik(2,2,1)*qk
      elehessfieldi(2,2,2) = elehessfieldi(2,2,2) - t3ik(2,2,2)*qk
      elehessfieldi(3,2,2) = elehessfieldi(3,2,2) - t3ik(3,2,2)*qk
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t3ik(3,2,1)*qk
c      elehessfieldi(3,2,2) = elehessfieldi(3,2,2) - t3ik(3,2,2)*qk
      elehessfieldi(3,3,2) = elehessfieldi(3,3,2) - t3ik(3,3,2)*qk
c      elehessfieldi(3,1,1) = elehessfieldi(3,1,1) - t3ik(3,1,1)*qk
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t3ik(3,2,1)*qk
c      elehessfieldi(3,3,1) = elehessfieldi(3,3,1) - t3ik(3,3,1)*qk
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t3ik(3,2,1)*qk
c      elehessfieldi(3,2,2) = elehessfieldi(3,2,2) - t3ik(3,2,2)*qk
c      elehessfieldi(3,3,2) = elehessfieldi(3,3,2) - t3ik(3,3,2)*qk
c      elehessfieldi(3,3,1) = elehessfieldi(3,3,1) - t3ik(3,3,1)*qk
c      elehessfieldi(3,3,2) = elehessfieldi(3,3,2) - t3ik(3,3,2)*qk
      elehessfieldi(3,3,3) = elehessfieldi(3,3,3) - t3ik(3,3,3)*qk
c
      elehessfieldk(1,1,1) = elehessfieldk(1,1,1) + t3ik(1,1,1)*qi
      elehessfieldk(2,1,1) = elehessfieldk(2,1,1) + t3ik(2,1,1)*qi
      elehessfieldk(3,1,1) = elehessfieldk(3,1,1) + t3ik(3,1,1)*qi
c      elehessfieldk(2,1,1) = elehessfieldk(2,1,1) + t3ik(2,1,1)*qi
      elehessfieldk(2,2,1) = elehessfieldk(2,2,1) + t3ik(2,2,1)*qi
      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) + t3ik(3,2,1)*qi
c      elehessfieldk(3,1,1) = elehessfieldk(3,1,1) + t3ik(3,1,1)*qi
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) + t3ik(3,2,1)*qi
      elehessfieldk(3,3,1) = elehessfieldk(3,3,1) + t3ik(3,3,1)*qi
c      elehessfieldk(2,1,1) = elehessfieldk(2,1,1) + t3ik(2,1,1)*qi
c      elehessfieldk(2,2,1) = elehessfieldk(2,2,1) + t3ik(2,2,1)*qi
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) + t3ik(3,2,1)*qi
c      elehessfieldk(2,2,1) = elehessfieldk(2,2,1) + t3ik(2,2,1)*qi
      elehessfieldk(2,2,2) = elehessfieldk(2,2,2) + t3ik(2,2,2)*qi
      elehessfieldk(3,2,2) = elehessfieldk(3,2,2) + t3ik(3,2,2)*qi
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) + t3ik(3,2,1)*qi
c      elehessfieldk(3,2,2) = elehessfieldk(3,2,2) + t3ik(3,2,2)*qi
      elehessfieldk(3,3,2) = elehessfieldk(3,3,2) + t3ik(3,3,2)*qi
c      elehessfieldk(3,1,1) = elehessfieldk(3,1,1) + t3ik(3,1,1)*qi
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) + t3ik(3,2,1)*qi
c      elehessfieldk(3,3,1) = elehessfieldk(3,3,1) + t3ik(3,3,1)*qi
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) + t3ik(3,2,1)*qi
c      elehessfieldk(3,2,2) = elehessfieldk(3,2,2) + t3ik(3,2,2)*qi
c      elehessfieldk(3,3,2) = elehessfieldk(3,3,2) + t3ik(3,3,2)*qi
c      elehessfieldk(3,3,1) = elehessfieldk(3,3,1) + t3ik(3,3,1)*qi
c      elehessfieldk(3,3,2) = elehessfieldk(3,3,2) + t3ik(3,3,2)*qi
      elehessfieldk(3,3,3) = elehessfieldk(3,3,3) + t3ik(3,3,3)*qi
c
c     from dipoles
c
      elehessfieldi(1,1,1) = elehessfieldi(1,1,1) - t4ik(1,1,1,1)*dkx - 
     &     t4ik(2,1,1,1)*dky - t4ik(3,1,1,1)*dkz
      elehessfieldi(2,1,1) = elehessfieldi(2,1,1) - t4ik(2,1,1,1)*dkx -
     &     t4ik(2,2,1,1)*dky - t4ik(3,2,1,1)*dkz
      elehessfieldi(3,1,1) = elehessfieldi(3,1,1) - t4ik(3,1,1,1)*dkx -
     &     t4ik(3,2,1,1)*dky - t4ik(3,3,1,1)*dkz
c      elehessfieldi(2,1,1) = elehessfieldi(2,1,1) - t4ik(2,1,1,1)*dkx -
c     &     t4ik(2,2,1,1)*dky - t4ik(3,2,1,1)*dkz
      elehessfieldi(2,2,1) = elehessfieldi(2,2,1) - t4ik(2,2,1,1)*dkx -
     &     t4ik(2,2,2,1)*dky - t4ik(3,2,2,1)*dkz
      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t4ik(3,2,1,1)*dkx -
     &     t4ik(3,2,2,1)*dky - t4ik(3,3,2,1)*dkz
c      elehessfieldi(3,1,1) = elehessfieldi(3,1,1) - t4ik(3,1,1,1)*dkx -
c     &     t4ik(3,2,1,1)*dky - t4ik(3,3,1,1)*dkz
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t4ik(3,2,1,1)*dkx -
c     &     t4ik(3,2,2,1)*dky - t4ik(3,3,2,1)*dkz
      elehessfieldi(3,3,1) = elehessfieldi(3,3,1) - t4ik(3,3,1,1)*dkx -
     &     t4ik(3,3,2,1)*dky - t4ik(3,3,3,1)*dkz
c      elehessfieldi(2,1,1) = elehessfieldi(2,1,1) - t4ik(2,1,1,1)*dkx -
c     &     t4ik(2,2,1,1)*dky - t4ik(3,2,1,1)*dkz
c      elehessfieldi(2,2,1) = elehessfieldi(2,2,1) - t4ik(2,2,1,1)*dkx -
c     &     t4ik(2,2,2,1)*dky - t4ik(3,2,2,1)*dkz
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t4ik(3,2,1,1)*dkx -
c     &     t4ik(3,2,2,1)*dky - t4ik(3,3,2,1)*dkz
c      elehessfieldi(2,2,1) = elehessfieldi(2,2,1) - t4ik(2,2,1,1)*dkx -
c     &     t4ik(2,2,2,1)*dky - t4ik(3,2,2,1)*dkz
      elehessfieldi(2,2,2) = elehessfieldi(2,2,2) - t4ik(2,2,2,1)*dkx -
     &     t4ik(2,2,2,2)*dky - t4ik(3,2,2,2)*dkz
      elehessfieldi(3,2,2) = elehessfieldi(3,2,2) - t4ik(3,2,2,1)*dkx -
     &     t4ik(3,2,2,2)*dky - t4ik(3,3,2,2)*dkz
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t4ik(3,2,1,1)*dkx -
c     &     t4ik(3,2,2,1)*dky - t4ik(3,3,2,1)*dkz
c      elehessfieldi(3,2,2) = elehessfieldi(3,2,2) - t4ik(3,2,2,1)*dkx -
c     &     t4ik(3,2,2,2)*dky - t4ik(3,3,2,2)*dkz
      elehessfieldi(3,3,2) = elehessfieldi(3,3,2) - t4ik(3,3,2,1)*dkx -
     &     t4ik(3,3,2,2)*dky - t4ik(3,3,3,2)*dkz
c      elehessfieldi(3,1,1) = elehessfieldi(3,1,1) - t4ik(3,1,1,1)*dkx -
c     &     t4ik(3,2,1,1)*dky - t4ik(3,3,1,1)*dkz
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t4ik(3,2,1,1)*dkx -
c     &     t4ik(3,2,2,1)*dky - t4ik(3,3,2,1)*dkz
c      elehessfieldi(3,3,1) = elehessfieldi(3,3,1) - t4ik(3,3,1,1)*dkx -
c     &     t4ik(3,3,2,1)*dky - t4ik(3,3,3,1)*dkz
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) - t4ik(3,2,1,1)*dkx -
c     &     t4ik(3,2,2,1)*dky - t4ik(3,3,2,1)*dkz
c      elehessfieldi(3,2,2) = elehessfieldi(3,2,2) - t4ik(3,2,2,1)*dkx -
c     &     t4ik(3,2,2,2)*dky - t4ik(3,3,2,2)*dkz
c      elehessfieldi(3,3,2) = elehessfieldi(3,3,2) - t4ik(3,3,2,1)*dkx -
c     &     t4ik(3,3,2,2)*dky - t4ik(3,3,3,2)*dkz
c      elehessfieldi(3,3,1) = elehessfieldi(3,3,1) - t4ik(3,3,1,1)*dkx -
c     &     t4ik(3,3,2,1)*dky - t4ik(3,3,3,1)*dkz
c      elehessfieldi(3,3,2) = elehessfieldi(3,3,2) - t4ik(3,3,2,1)*dkx -
c     &     t4ik(3,3,2,2)*dky - t4ik(3,3,3,2)*dkz
      elehessfieldi(3,3,3) = elehessfieldi(3,3,3) - t4ik(3,3,3,1)*dkx -
     &     t4ik(3,3,3,2)*dky - t4ik(3,3,3,3)*dkz
c
      elehessfieldk(1,1,1) = elehessfieldk(1,1,1) - t4ik(1,1,1,1)*dix - 
     &     t4ik(2,1,1,1)*diy - t4ik(3,1,1,1)*diz
      elehessfieldk(2,1,1) = elehessfieldk(2,1,1) - t4ik(2,1,1,1)*dix -
     &     t4ik(2,2,1,1)*diy - t4ik(3,2,1,1)*diz
      elehessfieldk(3,1,1) = elehessfieldk(3,1,1) - t4ik(3,1,1,1)*dix -
     &     t4ik(3,2,1,1)*diy - t4ik(3,3,1,1)*diz
c      elehessfieldk(2,1,1) = elehessfieldk(2,1,1) - t4ik(2,1,1,1)*dix -
c     &     t4ik(2,2,1,1)*diy - t4ik(3,2,1,1)*diz
      elehessfieldk(2,2,1) = elehessfieldk(2,2,1) - t4ik(2,2,1,1)*dix -
     &     t4ik(2,2,2,1)*diy - t4ik(3,2,2,1)*diz
      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) - t4ik(3,2,1,1)*dix -
     &     t4ik(3,2,2,1)*diy - t4ik(3,3,2,1)*diz
c      elehessfieldk(3,1,1) = elehessfieldk(3,1,1) - t4ik(3,1,1,1)*dix -
c     &     t4ik(3,2,1,1)*diy - t4ik(3,3,1,1)*diz
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) - t4ik(3,2,1,1)*dix -
c     &     t4ik(3,2,2,1)*diy - t4ik(3,3,2,1)*diz
      elehessfieldk(3,3,1) = elehessfieldk(3,3,1) - t4ik(3,3,1,1)*dix -
     &     t4ik(3,3,2,1)*diy - t4ik(3,3,3,1)*diz
c      elehessfieldk(2,1,1) = elehessfieldk(2,1,1) - t4ik(2,1,1,1)*dix -
c     &     t4ik(2,2,1,1)*diy - t4ik(3,2,1,1)*diz
c      elehessfieldk(2,2,1) = elehessfieldk(2,2,1) - t4ik(2,2,1,1)*dix -
c     &     t4ik(2,2,2,1)*diy - t4ik(3,2,2,1)*diz
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) - t4ik(3,2,1,1)*dix -
c     &     t4ik(3,2,2,1)*diy - t4ik(3,3,2,1)*diz
c      elehessfieldk(2,2,1) = elehessfieldk(2,2,1) - t4ik(2,2,1,1)*dix -
c     &     t4ik(2,2,2,1)*diy - t4ik(3,2,2,1)*diz
      elehessfieldk(2,2,2) = elehessfieldk(2,2,2) - t4ik(2,2,2,1)*dix -
     &     t4ik(2,2,2,2)*diy - t4ik(3,2,2,2)*diz
      elehessfieldk(3,2,2) = elehessfieldk(3,2,2) - t4ik(3,2,2,1)*dix -
     &     t4ik(3,2,2,2)*diy - t4ik(3,3,2,2)*diz
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) - t4ik(3,2,1,1)*dix -
c     &     t4ik(3,2,2,1)*diy - t4ik(3,3,2,1)*diz
c      elehessfieldk(3,2,2) = elehessfieldk(3,2,2) - t4ik(3,2,2,1)*dix -
c     &     t4ik(3,2,2,2)*diy - t4ik(3,3,2,2)*diz
      elehessfieldk(3,3,2) = elehessfieldk(3,3,2) - t4ik(3,3,2,1)*dix -
     &     t4ik(3,3,2,2)*diy - t4ik(3,3,3,2)*diz
c      elehessfieldk(3,1,1) = elehessfieldk(3,1,1) - t4ik(3,1,1,1)*dix -
c     &     t4ik(3,2,1,1)*diy - t4ik(3,3,1,1)*diz
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) - t4ik(3,2,1,1)*dix -
c     &     t4ik(3,2,2,1)*diy - t4ik(3,3,2,1)*diz
c      elehessfieldk(3,3,1) = elehessfieldk(3,3,1) - t4ik(3,3,1,1)*dix -
c     &     t4ik(3,3,2,1)*diy - t4ik(3,3,3,1)*diz
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) - t4ik(3,2,1,1)*dix -
c     &     t4ik(3,2,2,1)*diy - t4ik(3,3,2,1)*diz
c      elehessfieldk(3,2,2) = elehessfieldk(3,2,2) - t4ik(3,2,2,1)*dix -
c     &     t4ik(3,2,2,2)*diy - t4ik(3,3,2,2)*diz
c      elehessfieldk(3,3,2) = elehessfieldk(3,3,2) - t4ik(3,3,2,1)*dix -
c     &     t4ik(3,3,2,2)*diy - t4ik(3,3,3,2)*diz
c      elehessfieldk(3,3,1) = elehessfieldk(3,3,1) - t4ik(3,3,1,1)*dix -
c     &     t4ik(3,3,2,1)*diy - t4ik(3,3,3,1)*diz
c      elehessfieldk(3,3,2) = elehessfieldk(3,3,2) - t4ik(3,3,2,1)*dix -
c     &     t4ik(3,3,2,2)*diy - t4ik(3,3,3,2)*diz
      elehessfieldk(3,3,3) = elehessfieldk(3,3,3) - t4ik(3,3,3,1)*dix -
     &     t4ik(3,3,3,2)*diy - t4ik(3,3,3,3)*diz
c
c     quadrupoles
c
      elehessfieldi(1,1,1) = elehessfieldi(1,1,1) -t5ik(1,1,1,1,1)*qkxx- 
     &     t5ik(2,1,1,1,1)*qkxy - t5ik(3,1,1,1,1)*qkxz - 
     &     t5ik(2,2,1,1,1)*qkyy-
     &     t5ik(3,2,1,1,1)*qkyz - t5ik(3,3,1,1,1)*qkzz
      elehessfieldi(2,1,1) = elehessfieldi(2,1,1) -t5ik(2,1,1,1,1)*qkxx-
     &     t5ik(2,2,1,1,1)*qkxy - t5ik(3,2,1,1,1)*qkxz - 
     &     t5ik(2,2,2,1,1)*qkyy-
     &     t5ik(3,2,2,1,1)*qkyz - t5ik(3,3,2,1,1)*qkzz
      elehessfieldi(3,1,1) = elehessfieldi(3,1,1) -t5ik(3,1,1,1,1)*qkxx-
     &     t5ik(3,2,1,1,1)*qkxy - t5ik(3,3,1,1,1)*qkxz - 
     &     t5ik(3,2,2,1,1)*qkyy-
     &     t5ik(3,3,2,1,1)*qkyz - t5ik(3,3,3,1,1)*qkzz
c      elehessfieldi(2,1,1) = elehessfieldi(2,1,1) -t5ik(2,1,1,1,1)*qkxx-
c     &     t5ik(2,2,1,1,1)*qkxy - t5ik(3,2,1,1,1)*qkxz - 
c     &     t5ik(2,2,2,1,1)*qkyy-
c     &     t5ik(3,2,2,1,1)*qkyz - t5ik(3,3,2,1,1)*qkzz
      elehessfieldi(2,2,1) = elehessfieldi(2,2,1) -t5ik(2,2,1,1,1)*qkxx-
     &     t5ik(2,2,2,1,1)*qkxy - t5ik(3,2,2,1,1)*qkxz - 
     &     t5ik(2,2,2,2,1)*qkyy-
     &     t5ik(3,2,2,2,1)*qkyz - t5ik(3,3,2,2,1)*qkzz
      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) -t5ik(3,2,1,1,1)*qkxx-
     &     t5ik(3,2,2,1,1)*qkxy - t5ik(3,3,2,1,1)*qkxz - 
     &     t5ik(3,2,2,2,1)*qkyy-
     &     t5ik(3,3,2,2,1)*qkyz - t5ik(3,3,3,2,1)*qkzz
c      elehessfieldi(3,1,1) = elehessfieldi(3,1,1) -t5ik(3,1,1,1,1)*qkxx-
c     &     t5ik(3,2,1,1,1)*qkxy - t5ik(3,3,1,1,1)*qkxz - 
c     &     t5ik(3,2,2,1,1)*qkyy-
c     &     t5ik(3,3,2,1,1)*qkyz - t5ik(3,3,3,1,1)*qkzz
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) -t5ik(3,2,1,1,1)*qkxx-
c     &     t5ik(3,2,2,1,1)*qkxy - t5ik(3,3,2,1,1)*qkxz - 
c     &     t5ik(3,2,2,2,1)*qkyy-
c     &     t5ik(3,3,2,2,1)*qkyz - t5ik(3,3,3,2,1)*qkzz
      elehessfieldi(3,3,1) = elehessfieldi(3,3,1) -t5ik(3,3,1,1,1)*qkxx-
     &     t5ik(3,3,2,1,1)*qkxy - t5ik(3,3,3,1,1)*qkxz - 
     &     t5ik(3,3,2,2,1)*qkyy-
     &     t5ik(3,3,3,2,1)*qkyz - t5ik(3,3,3,3,1)*qkzz
c      elehessfieldi(2,1,1) = elehessfieldi(2,1,1) -t5ik(2,1,1,1,1)*qkxx-
c     &     t5ik(2,2,1,1,1)*qkxy - t5ik(3,2,1,1,1)*qkxz - 
c     &     t5ik(2,2,2,1,1)*qkyy-
c     &     t5ik(3,2,2,1,1)*qkyz - t5ik(3,3,2,1,1)*qkzz
c      elehessfieldi(2,2,1) = elehessfieldi(2,2,1) -t5ik(2,2,1,1,1)*qkxx-
c     &     t5ik(2,2,2,1,1)*qkxy - t5ik(3,2,2,1,1)*qkxz - 
c     &     t5ik(2,2,2,2,1)*qkyy-
c     &     t5ik(3,2,2,2,1)*qkyz - t5ik(3,3,2,2,1)*qkzz
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) -t5ik(3,2,1,1,1)*qkxx-
c     &     t5ik(3,2,2,1,1)*qkxy - t5ik(3,3,2,1,1)*qkxz - 
c     &     t5ik(3,2,2,2,1)*qkyy-
c     &     t5ik(3,3,2,2,1)*qkyz - t5ik(3,3,3,2,1)*qkzz
c      elehessfieldi(2,2,1) = elehessfieldi(2,2,1) -t5ik(2,2,1,1,1)*qkxx-
c     &     t5ik(2,2,2,1,1)*qkxy - t5ik(3,2,2,1,1)*qkxz - 
c     &     t5ik(2,2,2,2,1)*qkyy-
c     &     t5ik(3,2,2,2,1)*qkyz - t5ik(3,3,2,2,1)*qkzz
      elehessfieldi(2,2,2) = elehessfieldi(2,2,2) -t5ik(2,2,2,1,1)*qkxx-
     &     t5ik(2,2,2,2,1)*qkxy - t5ik(3,2,2,2,1)*qkxz - 
     &     t5ik(2,2,2,2,2)*qkyy-
     &     t5ik(3,2,2,2,2)*qkyz - t5ik(3,3,2,2,2)*qkzz
      elehessfieldi(3,2,2) = elehessfieldi(3,2,2) -t5ik(3,2,2,1,1)*qkxx-
     &     t5ik(3,2,2,2,1)*qkxy - t5ik(3,3,2,2,1)*qkxz - 
     &     t5ik(3,2,2,2,2)*qkyy-
     &     t5ik(3,3,2,2,2)*qkyz - t5ik(3,3,3,2,2)*qkzz
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) -t5ik(3,2,1,1,1)*qkxx-
c     &     t5ik(3,2,2,1,1)*qkxy - t5ik(3,3,2,1,1)*qkxz - 
c     7     t5ik(3,2,2,2,1)*qkyy-
c     &     t5ik(3,3,2,2,1)*qkyz - t5ik(3,3,3,2,1)*qkzz
c      elehessfieldi(3,2,2) = elehessfieldi(3,2,2) -t5ik(3,2,2,1,1)*qkxx-
c     &     t5ik(3,2,2,2,1)*qkxy - t5ik(3,3,2,2,1)*qkxz - 
c     &     t5ik(3,2,2,2,2)*qkyy-
c     &     t5ik(3,3,2,2,2)*qkyz - t5ik(3,3,3,2,2)*qkzz
      elehessfieldi(3,3,2) = elehessfieldi(3,3,2) -t5ik(3,3,2,1,1)*qkxx-
     &     t5ik(3,3,2,2,1)*qkxy - t5ik(3,3,3,2,1)*qkxz - 
     &     t5ik(3,3,2,2,2)*qkyy-
     &     t5ik(3,3,3,2,2)*qkyz - t5ik(3,3,3,3,2)*qkzz
c      elehessfieldi(3,1,1) = elehessfieldi(3,1,1) -t5ik(3,1,1,1,1)*qkxx-
c     &     t5ik(3,2,1,1,1)*qkxy - t5ik(3,3,1,1,1)*qkxz - 
c     &     t5ik(3,2,2,1,1)*qkyy-
c     &     t5ik(3,3,2,1,1)*qkyz - t5ik(3,3,3,1,1)*qkzz
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) -t5ik(3,2,1,1,1)*qkxx-
c     &     t5ik(3,2,2,1,1)*qkxy - t5ik(3,3,2,1,1)*qkxz - 
c     &     t5ik(3,2,2,2,1)*qkyy-
c     &     t5ik(3,3,2,2,1)*qkyz - t5ik(3,3,3,2,1)*qkzz
c      elehessfieldi(3,3,1) = elehessfieldi(3,3,1) -t5ik(3,3,1,1,1)*qkxx-
c     &     t5ik(3,3,2,1,1)*qkxy - t5ik(3,3,3,1,1)*qkxz - 
c     &     t5ik(3,3,2,2,1)*qkyy-
c     &     t5ik(3,3,3,2,1)*qkyz - t5ik(3,3,3,3,1)*qkzz
c      elehessfieldi(3,2,1) = elehessfieldi(3,2,1) -t5ik(3,2,1,1,1)*qkxx-
c     &     t5ik(3,2,2,1,1)*qkxy - t5ik(3,3,2,1,1)*qkxz - 
c     &     t5ik(3,2,2,2,1)*qkyy-
c     &     t5ik(3,3,2,2,1)*qkyz - t5ik(3,3,3,2,1)*qkzz
c      elehessfieldi(3,2,2) = elehessfieldi(3,2,2) -t5ik(3,2,2,1,1)*qkxx-
c     &     t5ik(3,2,2,2,1)*qkxy - t5ik(3,3,2,2,1)*qkxz - 
c     &     t5ik(3,2,2,2,2)*qkyy-
c     &     t5ik(3,3,2,2,2)*qkyz - t5ik(3,3,3,2,2)*qkzz
c      elehessfieldi(3,3,2) = elehessfieldi(3,3,2) -t5ik(3,3,2,1,1)*qkxx-
c     &     t5ik(3,3,2,2,1)*qkxy - t5ik(3,3,3,2,1)*qkxz - 
c     &     t5ik(3,3,2,2,2)*qkyy-
c     &     t5ik(3,3,3,2,2)*qkyz - t5ik(3,3,3,3,2)*qkzz
c      elehessfieldi(3,3,1) = elehessfieldi(3,3,1) -t5ik(3,3,1,1,1)*qkxx-
c     &     t5ik(3,3,2,1,1)*qkxy - t5ik(3,3,3,1,1)*qkxz - 
c     &     t5ik(3,3,2,2,1)*qkyy-
c     &     t5ik(3,3,3,2,1)*qkyz - t5ik(3,3,3,3,1)*qkzz
c      elehessfieldi(3,3,2) = elehessfieldi(3,3,2) -t5ik(3,3,2,1,1)*qkxx-
c     &     t5ik(3,3,2,2,1)*qkxy - t5ik(3,3,3,2,1)*qkxz - 
c     &     t5ik(3,3,2,2,2)*qkyy-
c     &     t5ik(3,3,3,2,2)*qkyz - t5ik(3,3,3,3,2)*qkzz
      elehessfieldi(3,3,3) = elehessfieldi(3,3,3) -t5ik(3,3,3,1,1)*qkxx-
     &     t5ik(3,3,3,2,1)*qkxy - t5ik(3,3,3,3,1)*qkxz - 
     &     t5ik(3,3,3,2,2)*qkyy-
     &     t5ik(3,3,3,3,2)*qkyz - t5ik(3,3,3,3,3)*qkzz
c
      elehessfieldk(1,1,1) = elehessfieldk(1,1,1) +t5ik(1,1,1,1,1)*qixx+ 
     &     t5ik(2,1,1,1,1)*qixy + t5ik(3,1,1,1,1)*qixz + 
     &     t5ik(2,2,1,1,1)*qiyy+
     &     t5ik(3,2,1,1,1)*qiyz + t5ik(3,3,1,1,1)*qizz
      elehessfieldk(2,1,1) = elehessfieldk(2,1,1) +t5ik(2,1,1,1,1)*qixx+
     &     t5ik(2,2,1,1,1)*qixy + t5ik(3,2,1,1,1)*qixz + 
     &     t5ik(2,2,2,1,1)*qiyy+
     &     t5ik(3,2,2,1,1)*qiyz + t5ik(3,3,2,1,1)*qizz
      elehessfieldk(3,1,1) = elehessfieldk(3,1,1) +t5ik(3,1,1,1,1)*qixx+
     &     t5ik(3,2,1,1,1)*qixy + t5ik(3,3,1,1,1)*qixz + 
     &     t5ik(3,2,2,1,1)*qiyy+
     &     t5ik(3,3,2,1,1)*qiyz + t5ik(3,3,3,1,1)*qizz
c      elehessfieldk(2,1,1) = elehessfieldk(2,1,1) +t5ik(2,1,1,1,1)*qixx+
c     &     t5ik(2,2,1,1,1)*qixy + t5ik(3,2,1,1,1)*qixz + 
c     &     t5ik(2,2,2,1,1)*qiyy+
c     &     t5ik(3,2,2,1,1)*qiyz + t5ik(3,3,2,1,1)*qizz
      elehessfieldk(2,2,1) = elehessfieldk(2,2,1) +t5ik(2,2,1,1,1)*qixx+
     &     t5ik(2,2,2,1,1)*qixy + t5ik(3,2,2,1,1)*qixz + 
     &     t5ik(2,2,2,2,1)*qiyy+
     &     t5ik(3,2,2,2,1)*qiyz + t5ik(3,3,2,2,1)*qizz
      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) +t5ik(3,2,1,1,1)*qixx+
     &     t5ik(3,2,2,1,1)*qixy + t5ik(3,3,2,1,1)*qixz + 
     &     t5ik(3,2,2,2,1)*qiyy+
     &     t5ik(3,3,2,2,1)*qiyz + t5ik(3,3,3,2,1)*qizz
c      elehessfieldk(3,1,1) = elehessfieldk(3,1,1) +t5ik(3,1,1,1,1)*qixx+
c     &     t5ik(3,2,1,1,1)*qixy + t5ik(3,3,1,1,1)*qixz + 
c     &     t5ik(3,2,2,1,1)*qiyy+
c     &     t5ik(3,3,2,1,1)*qiyz + t5ik(3,3,3,1,1)*qizz
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) +t5ik(3,2,1,1,1)*qixx+
c     &     t5ik(3,2,2,1,1)*qixy + t5ik(3,3,2,1,1)*qixz + 
c     &     t5ik(3,2,2,2,1)*qiyy+
c     &     t5ik(3,3,2,2,1)*qiyz + t5ik(3,3,3,2,1)*qizz
      elehessfieldk(3,3,1) = elehessfieldk(3,3,1) +t5ik(3,3,1,1,1)*qixx+
     &     t5ik(3,3,2,1,1)*qixy + t5ik(3,3,3,1,1)*qixz + 
     &     t5ik(3,3,2,2,1)*qiyy+
     &     t5ik(3,3,3,2,1)*qiyz + t5ik(3,3,3,3,1)*qizz
c      elehessfieldk(2,1,1) = elehessfieldk(2,1,1) +t5ik(2,1,1,1,1)*qixx+
c     &     t5ik(2,2,1,1,1)*qixy + t5ik(3,2,1,1,1)*qixz + 
c     &     t5ik(2,2,2,1,1)*qiyy+
c     &     t5ik(3,2,2,1,1)*qiyz + t5ik(3,3,2,1,1)*qizz
c      elehessfieldk(2,2,1) = elehessfieldk(2,2,1) +t5ik(2,2,1,1,1)*qixx+
c     &     t5ik(2,2,2,1,1)*qixy + t5ik(3,2,2,1,1)*qixz + 
c     &     t5ik(2,2,2,2,1)*qiyy+
c     &     t5ik(3,2,2,2,1)*qiyz + t5ik(3,3,2,2,1)*qizz
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) +t5ik(3,2,1,1,1)*qixx+
c     &     t5ik(3,2,2,1,1)*qixy + t5ik(3,3,2,1,1)*qixz + 
c     &     t5ik(3,2,2,2,1)*qiyy+
c     &     t5ik(3,3,2,2,1)*qiyz + t5ik(3,3,3,2,1)*qizz
c      elehessfieldk(2,2,1) = elehessfieldk(2,2,1) +t5ik(2,2,1,1,1)*qixx+
c     &     t5ik(2,2,2,1,1)*qixy + t5ik(3,2,2,1,1)*qixz + 
c     &     t5ik(2,2,2,2,1)*qiyy+
c     &     t5ik(3,2,2,2,1)*qiyz + t5ik(3,3,2,2,1)*qizz
      elehessfieldk(2,2,2) = elehessfieldk(2,2,2) +t5ik(2,2,2,1,1)*qixx+
     &     t5ik(2,2,2,2,1)*qixy + t5ik(3,2,2,2,1)*qixz + 
     &     t5ik(2,2,2,2,2)*qiyy+
     &     t5ik(3,2,2,2,2)*qiyz + t5ik(3,3,2,2,2)*qizz
      elehessfieldk(3,2,2) = elehessfieldk(3,2,2) +t5ik(3,2,2,1,1)*qixx+
     &     t5ik(3,2,2,2,1)*qixy + t5ik(3,3,2,2,1)*qixz + 
     &     t5ik(3,2,2,2,2)*qiyy+
     &     t5ik(3,3,2,2,2)*qiyz + t5ik(3,3,3,2,2)*qizz
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) +t5ik(3,2,1,1,1)*qixx+
c     &     t5ik(3,2,2,1,1)*qixy + t5ik(3,3,2,1,1)*qixz + 
c     &     t5ik(3,2,2,2,1)*qiyy+
c     &     t5ik(3,3,2,2,1)*qiyz + t5ik(3,3,3,2,1)*qizz
c      elehessfieldk(3,2,2) = elehessfieldk(3,2,2) +t5ik(3,2,2,1,1)*qixx+
c     &     t5ik(3,2,2,2,1)*qixy + t5ik(3,3,2,2,1)*qixz + 
c     &     t5ik(3,2,2,2,2)*qiyy+
c     &     t5ik(3,3,2,2,2)*qiyz + t5ik(3,3,3,2,2)*qizz
      elehessfieldk(3,3,2) = elehessfieldk(3,3,2) +t5ik(3,3,2,1,1)*qixx+
     &     t5ik(3,3,2,2,1)*qixy + t5ik(3,3,3,2,1)*qixz + 
     &     t5ik(3,3,2,2,2)*qiyy+
     &     t5ik(3,3,3,2,2)*qiyz + t5ik(3,3,3,3,2)*qizz
c      elehessfieldk(3,1,1) = elehessfieldk(3,1,1) +t5ik(3,1,1,1,1)*qixx+
c     &     t5ik(3,2,1,1,1)*qixy + t5ik(3,3,1,1,1)*qixz + 
c     &     t5ik(3,2,2,1,1)*qiyy+
c     &     t5ik(3,3,2,1,1)*qiyz + t5ik(3,3,3,1,1)*qizz
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) +t5ik(3,2,1,1,1)*qixx+
c     &     t5ik(3,2,2,1,1)*qixy + t5ik(3,3,2,1,1)*qixz + 
c     &     t5ik(3,2,2,2,1)*qiyy+
c     &     t5ik(3,3,2,2,1)*qiyz + t5ik(3,3,3,2,1)*qizz
c      elehessfieldk(3,3,1) = elehessfieldk(3,3,1) +t5ik(3,3,1,1,1)*qixx+
c     &     t5ik(3,3,2,1,1)*qixy + t5ik(3,3,3,1,1)*qixz + 
c     &     t5ik(3,3,2,2,1)*qiyy+
c     &     t5ik(3,3,3,2,1)*qiyz + t5ik(3,3,3,3,1)*qizz
c      elehessfieldk(3,2,1) = elehessfieldk(3,2,1) +t5ik(3,2,1,1,1)*qixx+
c     &     t5ik(3,2,2,1,1)*qixy + t5ik(3,3,2,1,1)*qixz + 
c     &     t5ik(3,2,2,2,1)*qiyy+
c     &     t5ik(3,3,2,2,1)*qiyz + t5ik(3,3,3,2,1)*qizz
c      elehessfieldk(3,2,2) = elehessfieldk(3,2,2) +t5ik(3,2,2,1,1)*qixx+
c     &     t5ik(3,2,2,2,1)*qixy + t5ik(3,3,2,2,1)*qixz + 
c     &     t5ik(3,2,2,2,2)*qiyy+
c     &     t5ik(3,3,2,2,2)*qiyz + t5ik(3,3,3,2,2)*qizz
c      elehessfieldk(3,3,2) = elehessfieldk(3,3,2) +t5ik(3,3,2,1,1)*qixx+
c     &     t5ik(3,3,2,2,1)*qixy + t5ik(3,3,3,2,1)*qixz + 
c     &     t5ik(3,3,2,2,2)*qiyy+
c     &     t5ik(3,3,3,2,2)*qiyz + t5ik(3,3,3,3,2)*qizz
c      elehessfieldk(3,3,1) = elehessfieldk(3,3,1) +t5ik(3,3,1,1,1)*qixx+
c     &     t5ik(3,3,2,1,1)*qixy + t5ik(3,3,3,1,1)*qixz + 
c     &     t5ik(3,3,2,2,1)*qiyy+
c     &     t5ik(3,3,3,2,1)*qiyz + t5ik(3,3,3,3,1)*qizz
c      elehessfieldk(3,3,2) = elehessfieldk(3,3,2) +t5ik(3,3,2,1,1)*qixx+
c     &     t5ik(3,3,2,2,1)*qixy + t5ik(3,3,3,2,1)*qixz + 
c     &     t5ik(3,3,2,2,2)*qiyy+
c     &     t5ik(3,3,3,2,2)*qiyz + t5ik(3,3,3,3,2)*qizz
      elehessfieldk(3,3,3) = elehessfieldk(3,3,3) +t5ik(3,3,3,1,1)*qixx+
     &     t5ik(3,3,3,2,1)*qixy + t5ik(3,3,3,3,1)*qixz + 
     &     t5ik(3,3,3,2,2)*qiyy+
     &     t5ik(3,3,3,3,2)*qiyz + t5ik(3,3,3,3,3)*qizz
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine uhessfieldik  --  electric field hessian  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "uhessfieldik" computes the electric field hessian given two multipoles
c     and damped or undamped tmatrix
c
      subroutine uhessfieldik(i,k,t4,hessfieldid,hessfieldkd,
     &     hessfieldip,hessfieldkp)
      use polar
      implicit none
      integer i,k
      real*8 duix,duiy,duiz
      real*8 dukx,duky,dukz
      real*8 puix,puiy,puiz
      real*8 pukx,puky,pukz
      real*8 t4(3,3,3,3)
      real*8 hessfieldid(3,3,3),hessfieldkd(3,3,3)
      real*8 hessfieldip(3,3,3),hessfieldkp(3,3,3)
c
c     read in induced dipole values
c
      duix = uind(1,i)
      duiy = uind(2,i)
      duiz = uind(3,i)
      dukx = uind(1,k)
      duky = uind(2,k)
      dukz = uind(3,k)
      puix = uinp(1,i)
      puiy = uinp(2,i)
      puiz = uinp(3,i)
      pukx = uinp(1,k)
      puky = uinp(2,k)
      pukz = uinp(3,k)
c
c     calculate electric field hessian
c
      hessfieldid(1,1,1) = -t4(1,1,1,1)*dukx - 
     &     t4(2,1,1,1)*duky - t4(3,1,1,1)*dukz
      hessfieldid(2,1,1) = -t4(2,1,1,1)*dukx -
     &     t4(2,2,1,1)*duky - t4(3,2,1,1)*dukz
      hessfieldid(3,1,1) = -t4(3,1,1,1)*dukx -
     &     t4(3,2,1,1)*duky - t4(3,3,1,1)*dukz
c      hessfieldid(2,1,1) = -t4(2,1,1,1)*dukx -
c     &     t4(2,2,1,1)*duky - t4(3,2,1,1)*dukz
      hessfieldid(2,2,1) = -t4(2,2,1,1)*dukx -
     &     t4(2,2,2,1)*duky - t4(3,2,2,1)*dukz
      hessfieldid(3,2,1) = -t4(3,2,1,1)*dukx -
     &     t4(3,2,2,1)*duky - t4(3,3,2,1)*dukz
c      hessfieldid(3,1,1) = -t4(3,1,1,1)*dukx -
c     &     t4(3,2,1,1)*duky - t4(3,3,1,1)*dukz
c      hessfieldid(3,2,1) = -t4(3,2,1,1)*dukx -
c     &     t4(3,2,2,1)*duky - t4(3,3,2,1)*dukz
      hessfieldid(3,3,1) = -t4(3,3,1,1)*dukx -
     &     t4(3,3,2,1)*duky - t4(3,3,3,1)*dukz
c      hessfieldid(2,1,1) = -t4(2,1,1,1)*dukx -
c     &     t4(2,2,1,1)*duky - t4(3,2,1,1)*dukz
c      hessfieldid(2,2,1) = -t4(2,2,1,1)*dukx -
c     &     t4(2,2,2,1)*duky - t4(3,2,2,1)*dukz
c      hessfieldid(3,2,1) = -t4(3,2,1,1)*dukx -
c     &     t4(3,2,2,1)*duky - t4(3,3,2,1)*dukz
c      hessfieldid(2,2,1) = -t4(2,2,1,1)*dukx -
c     &     t4(2,2,2,1)*duky - t4(3,2,2,1)*dukz
      hessfieldid(2,2,2) = -t4(2,2,2,1)*dukx -
     &     t4(2,2,2,2)*duky - t4(3,2,2,2)*dukz
      hessfieldid(3,2,2) = -t4(3,2,2,1)*dukx -
     &     t4(3,2,2,2)*duky - t4(3,3,2,2)*dukz
c      hessfieldid(3,2,1) = -t4(3,2,1,1)*dukx -
c     &     t4(3,2,2,1)*duky - t4(3,3,2,1)*dukz
c      hessfieldid(3,2,2) = -t4(3,2,2,1)*dukx -
c     &     t4(3,2,2,2)*duky - t4(3,3,2,2)*dukz
      hessfieldid(3,3,2) = -t4(3,3,2,1)*dukx -
     &     t4(3,3,2,2)*duky - t4(3,3,3,2)*dukz
c      hessfieldid(3,1,1) = -t4(3,1,1,1)*dukx -
c     &     t4(3,2,1,1)*duky - t4(3,3,1,1)*dukz
c      hessfieldid(3,2,1) = -t4(3,2,1,1)*dukx -
c     &     t4(3,2,2,1)*duky - t4(3,3,2,1)*dukz
c      hessfieldid(3,3,1) = -t4(3,3,1,1)*dukx -
c     &     t4(3,3,2,1)*duky - t4(3,3,3,1)*dukz
c      hessfieldid(3,2,1) = -t4(3,2,1,1)*dukx -
c     &     t4(3,2,2,1)*duky - t4(3,3,2,1)*dukz
c      hessfieldid(3,2,2) = -t4(3,2,2,1)*dukx -
c     &     t4(3,2,2,2)*duky - t4(3,3,2,2)*dukz
c      hessfieldid(3,3,2) = -t4(3,3,2,1)*dukx -
c     &     t4(3,3,2,2)*duky - t4(3,3,3,2)*dukz
c      hessfieldid(3,3,1) = -t4(3,3,1,1)*dukx -
c     &     t4(3,3,2,1)*duky - t4(3,3,3,1)*dukz
c      hessfieldid(3,3,2) = -t4(3,3,2,1)*dukx -
c     &     t4(3,3,2,2)*duky - t4(3,3,3,2)*dukz
      hessfieldid(3,3,3) = -t4(3,3,3,1)*dukx -
     &     t4(3,3,3,2)*duky - t4(3,3,3,3)*dukz
c
      hessfieldkd(1,1,1) = -t4(1,1,1,1)*duix - 
     &     t4(2,1,1,1)*duiy - t4(3,1,1,1)*duiz
      hessfieldkd(2,1,1) = -t4(2,1,1,1)*duix -
     &     t4(2,2,1,1)*duiy - t4(3,2,1,1)*duiz
      hessfieldkd(3,1,1) = -t4(3,1,1,1)*duix -
     &     t4(3,2,1,1)*duiy - t4(3,3,1,1)*duiz
c      hessfieldkd(2,1,1) = -t4(2,1,1,1)*duix -
c     &     t4(2,2,1,1)*duiy - t4(3,2,1,1)*duiz
      hessfieldkd(2,2,1) = -t4(2,2,1,1)*duix -
     &     t4(2,2,2,1)*duiy - t4(3,2,2,1)*duiz
      hessfieldkd(3,2,1) = -t4(3,2,1,1)*duix -
     &     t4(3,2,2,1)*duiy - t4(3,3,2,1)*duiz
c      hessfieldkd(3,1,1) = -t4(3,1,1,1)*duix -
c     &     t4(3,2,1,1)*duiy - t4(3,3,1,1)*duiz
c      hessfieldkd(3,2,1) = -t4(3,2,1,1)*duix -
c     &     t4(3,2,2,1)*duiy - t4(3,3,2,1)*duiz
      hessfieldkd(3,3,1) = -t4(3,3,1,1)*duix -
     &     t4(3,3,2,1)*duiy - t4(3,3,3,1)*duiz
c      hessfieldkd(2,1,1) = -t4(2,1,1,1)*duix -
c     &     t4(2,2,1,1)*duiy - t4(3,2,1,1)*duiz
c      hessfieldkd(2,2,1) = -t4(2,2,1,1)*duix -
c     &     t4(2,2,2,1)*duiy - t4(3,2,2,1)*duiz
c      hessfieldkd(3,2,1) = -t4(3,2,1,1)*duix -
c     &     t4(3,2,2,1)*duiy - t4(3,3,2,1)*duiz
c      hessfieldkd(2,2,1) = -t4(2,2,1,1)*duix -
c     &     t4(2,2,2,1)*duiy - t4(3,2,2,1)*duiz
      hessfieldkd(2,2,2) = -t4(2,2,2,1)*duix -
     &     t4(2,2,2,2)*duiy - t4(3,2,2,2)*duiz
      hessfieldkd(3,2,2) = -t4(3,2,2,1)*duix -
     &     t4(3,2,2,2)*duiy - t4(3,3,2,2)*duiz
c      hessfieldkd(3,2,1) = -t4(3,2,1,1)*duix -
c     &     t4(3,2,2,1)*duiy - t4(3,3,2,1)*duiz
c      hessfieldkd(3,2,2) = -t4(3,2,2,1)*duix -
c     &     t4(3,2,2,2)*duiy - t4(3,3,2,2)*duiz
      hessfieldkd(3,3,2) = -t4(3,3,2,1)*duix -
     &     t4(3,3,2,2)*duiy - t4(3,3,3,2)*duiz
c      hessfieldkd(3,1,1) = -t4(3,1,1,1)*duix -
c     &     t4(3,2,1,1)*duiy - t4(3,3,1,1)*duiz
c      hessfieldkd(3,2,1) = -t4(3,2,1,1)*duix -
c     &     t4(3,2,2,1)*duiy - t4(3,3,2,1)*duiz
c      hessfieldkd(3,3,1) = -t4(3,3,1,1)*duix -
c     &     t4(3,3,2,1)*duiy - t4(3,3,3,1)*duiz
c      hessfieldkd(3,2,1) = -t4(3,2,1,1)*duix -
c     &     t4(3,2,2,1)*duiy - t4(3,3,2,1)*duiz
c      hessfieldkd(3,2,2) = -t4(3,2,2,1)*duix -
c     &     t4(3,2,2,2)*duiy - t4(3,3,2,2)*duiz
c      hessfieldkd(3,3,2) = -t4(3,3,2,1)*duix -
c     &     t4(3,3,2,2)*duiy - t4(3,3,3,2)*duiz
c      hessfieldkd(3,3,1) = -t4(3,3,1,1)*duix -
c     &     t4(3,3,2,1)*duiy - t4(3,3,3,1)*duiz
c      hessfieldkd(3,3,2) = -t4(3,3,2,1)*duix -
c     &     t4(3,3,2,2)*duiy - t4(3,3,3,2)*duiz
      hessfieldkd(3,3,3) = -t4(3,3,3,1)*duix -
     &     t4(3,3,3,2)*duiy - t4(3,3,3,3)*duiz
c
c     p induced dipoles
c
      hessfieldip(1,1,1) = -t4(1,1,1,1)*pukx - 
     &     t4(2,1,1,1)*puky - t4(3,1,1,1)*pukz
      hessfieldip(2,1,1) = -t4(2,1,1,1)*pukx -
     &     t4(2,2,1,1)*puky - t4(3,2,1,1)*pukz
      hessfieldip(3,1,1) = -t4(3,1,1,1)*pukx -
     &     t4(3,2,1,1)*puky - t4(3,3,1,1)*pukz
c      hessfieldip(2,1,1) = -t4(2,1,1,1)*pukx -
c     &     t4(2,2,1,1)*puky - t4(3,2,1,1)*pukz
      hessfieldip(2,2,1) = -t4(2,2,1,1)*pukx -
     &     t4(2,2,2,1)*puky - t4(3,2,2,1)*pukz
      hessfieldip(3,2,1) = -t4(3,2,1,1)*pukx -
     &     t4(3,2,2,1)*puky - t4(3,3,2,1)*pukz
c      hessfieldip(3,1,1) = -t4(3,1,1,1)*pukx -
c     &     t4(3,2,1,1)*puky - t4(3,3,1,1)*pukz
c      hessfieldip(3,2,1) = -t4(3,2,1,1)*pukx -
c     &     t4(3,2,2,1)*puky - t4(3,3,2,1)*pukz
      hessfieldip(3,3,1) = -t4(3,3,1,1)*pukx -
     &     t4(3,3,2,1)*puky - t4(3,3,3,1)*pukz
c      hessfieldip(2,1,1) = -t4(2,1,1,1)*pukx -
c     &     t4(2,2,1,1)*puky - t4(3,2,1,1)*pukz
c      hessfieldip(2,2,1) = -t4(2,2,1,1)*pukx -
c     &     t4(2,2,2,1)*puky - t4(3,2,2,1)*pukz
c      hessfieldip(3,2,1) = -t4(3,2,1,1)*pukx -
c     &     t4(3,2,2,1)*puky - t4(3,3,2,1)*pukz
c      hessfieldip(2,2,1) = -t4(2,2,1,1)*pukx -
c     &     t4(2,2,2,1)*puky - t4(3,2,2,1)*pukz
      hessfieldip(2,2,2) = -t4(2,2,2,1)*pukx -
     &     t4(2,2,2,2)*puky - t4(3,2,2,2)*pukz
      hessfieldip(3,2,2) = -t4(3,2,2,1)*pukx -
     &     t4(3,2,2,2)*puky - t4(3,3,2,2)*pukz
c      hessfieldip(3,2,1) = -t4(3,2,1,1)*pukx -
c     &     t4(3,2,2,1)*puky - t4(3,3,2,1)*pukz
c      hessfieldip(3,2,2) = -t4(3,2,2,1)*pukx -
c     &     t4(3,2,2,2)*puky - t4(3,3,2,2)*pukz
      hessfieldip(3,3,2) = -t4(3,3,2,1)*pukx -
     &     t4(3,3,2,2)*puky - t4(3,3,3,2)*pukz
c      hessfieldip(3,1,1) = -t4(3,1,1,1)*pukx -
c     &     t4(3,2,1,1)*puky - t4(3,3,1,1)*pukz
c      hessfieldip(3,2,1) = -t4(3,2,1,1)*pukx -
c     &     t4(3,2,2,1)*puky - t4(3,3,2,1)*pukz
c      hessfieldip(3,3,1) = -t4(3,3,1,1)*pukx -
c     &     t4(3,3,2,1)*puky - t4(3,3,3,1)*pukz
c      hessfieldip(3,2,1) = -t4(3,2,1,1)*pukx -
c     &     t4(3,2,2,1)*puky - t4(3,3,2,1)*pukz
c      hessfieldip(3,2,2) = -t4(3,2,2,1)*pukx -
c     &     t4(3,2,2,2)*puky - t4(3,3,2,2)*pukz
c      hessfieldip(3,3,2) = -t4(3,3,2,1)*pukx -
c     &     t4(3,3,2,2)*puky - t4(3,3,3,2)*pukz
c      hessfieldip(3,3,1) = -t4(3,3,1,1)*pukx -
c     &     t4(3,3,2,1)*puky - t4(3,3,3,1)*pukz
c      hessfieldip(3,3,2) = -t4(3,3,2,1)*pukx -
c     &     t4(3,3,2,2)*puky - t4(3,3,3,2)*pukz
      hessfieldip(3,3,3) = -t4(3,3,3,1)*pukx -
     &     t4(3,3,3,2)*puky - t4(3,3,3,3)*pukz
c
      hessfieldkp(1,1,1) = -t4(1,1,1,1)*puix - 
     &     t4(2,1,1,1)*puiy - t4(3,1,1,1)*puiz
      hessfieldkp(2,1,1) = -t4(2,1,1,1)*puix -
     &     t4(2,2,1,1)*puiy - t4(3,2,1,1)*puiz
      hessfieldkp(3,1,1) = -t4(3,1,1,1)*puix -
     &     t4(3,2,1,1)*puiy - t4(3,3,1,1)*puiz
c      hessfieldkp(2,1,1) = -t4(2,1,1,1)*puix -
c     &     t4(2,2,1,1)*puiy - t4(3,2,1,1)*puiz
      hessfieldkp(2,2,1) = -t4(2,2,1,1)*puix -
     &     t4(2,2,2,1)*puiy - t4(3,2,2,1)*puiz
      hessfieldkp(3,2,1) = -t4(3,2,1,1)*puix -
     &     t4(3,2,2,1)*puiy - t4(3,3,2,1)*puiz
c      hessfieldkp(3,1,1) = -t4(3,1,1,1)*puix -
c     &     t4(3,2,1,1)*puiy - t4(3,3,1,1)*puiz
c      hessfieldkp(3,2,1) = -t4(3,2,1,1)*puix -
c     &     t4(3,2,2,1)*puiy - t4(3,3,2,1)*puiz
      hessfieldkp(3,3,1) = -t4(3,3,1,1)*puix -
     &     t4(3,3,2,1)*puiy - t4(3,3,3,1)*puiz
c      hessfieldkp(2,1,1) = -t4(2,1,1,1)*puix -
c     &     t4(2,2,1,1)*puiy - t4(3,2,1,1)*puiz
c      hessfieldkp(2,2,1) = -t4(2,2,1,1)*puix -
c     &     t4(2,2,2,1)*puiy - t4(3,2,2,1)*puiz
c      hessfieldkp(3,2,1) = -t4(3,2,1,1)*puix -
c     &     t4(3,2,2,1)*puiy - t4(3,3,2,1)*puiz
c      hessfieldkp(2,2,1) = -t4(2,2,1,1)*puix -
c     &     t4(2,2,2,1)*puiy - t4(3,2,2,1)*puiz
      hessfieldkp(2,2,2) = -t4(2,2,2,1)*puix -
     &     t4(2,2,2,2)*puiy - t4(3,2,2,2)*puiz
      hessfieldkp(3,2,2) = -t4(3,2,2,1)*puix -
     &     t4(3,2,2,2)*puiy - t4(3,3,2,2)*puiz
c      hessfieldkp(3,2,1) = -t4(3,2,1,1)*puix -
c     &     t4(3,2,2,1)*puiy - t4(3,3,2,1)*puiz
c      hessfieldkp(3,2,2) = -t4(3,2,2,1)*puix -
c     &     t4(3,2,2,2)*puiy - t4(3,3,2,2)*puiz
      hessfieldkp(3,3,2) = -t4(3,3,2,1)*puix -
     &     t4(3,3,2,2)*puiy - t4(3,3,3,2)*puiz
c      hessfieldkp(3,1,1) = -t4(3,1,1,1)*puix -
c     &     t4(3,2,1,1)*puiy - t4(3,3,1,1)*puiz
c      hessfieldkp(3,2,1) = -t4(3,2,1,1)*puix -
c     &     t4(3,2,2,1)*puiy - t4(3,3,2,1)*puiz
c      hessfieldkp(3,3,1) = -t4(3,3,1,1)*puix -
c     &     t4(3,3,2,1)*puiy - t4(3,3,3,1)*puiz
c      hessfieldkp(3,2,1) = -t4(3,2,1,1)*puix -
c     &     t4(3,2,2,1)*puiy - t4(3,3,2,1)*puiz
c      hessfieldkp(3,2,2) = -t4(3,2,2,1)*puix -
c     &     t4(3,2,2,2)*puiy - t4(3,3,2,2)*puiz
c      hessfieldkp(3,3,2) = -t4(3,3,2,1)*puix -
c     &     t4(3,3,2,2)*puiy - t4(3,3,3,2)*puiz
c      hessfieldkp(3,3,1) = -t4(3,3,1,1)*puix -
c     &     t4(3,3,2,1)*puiy - t4(3,3,3,1)*puiz
c      hessfieldkp(3,3,2) = -t4(3,3,2,1)*puix -
c     &     t4(3,3,2,2)*puiy - t4(3,3,3,2)*puiz
      hessfieldkp(3,3,3) = -t4(3,3,3,1)*puix -
     &     t4(3,3,3,2)*puiy - t4(3,3,3,3)*puiz
      return
      end
