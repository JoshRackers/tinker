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
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
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
c     #################################################
c     ##                                             ##
c     ##  subroutine ufieldik  --  electric field     ##
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
c     WHERE SHOULD THE NEGATIVES BE???
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
cccccccccccccc THESE ARE RIGHT. THIS MEANS THAT GRADFIELD SHOULD BE THE
cccccccccccccc NEGATIVE DERIVATIVE OF FIELD.  (GRADFIELD = --dV/dxdx)
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
c     
c     everything here is negative (hessfield = ---dV/drdrdr)
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
