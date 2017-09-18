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
c      subroutine extra3
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
      real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9
      real*8 t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
      real*8 t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
      real*8 t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
      real*8 t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
      real*8 t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
      real*8 t60,t61,t62,t63,t64,t65,t66,t67,t68,t69
      real*8 t70,t71,t72,t73,t74,t75,t76,t77,t78,t79
      real*8 t80,t81,t82,t83,t84,t85,t86,t87,t88,t89
      real*8 t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
      real*8 t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
      real*8 t110,t111,t112,t113,t114,t115,t116,t117,t118,t119
      real*8 t120,t121,t122,t123,t124,t125,t126,t127,t128,t129
      real*8 t130,t131,t132,t133,t134,t135,t136,t137,t138,t139
      real*8 t140,t141,t142,t143,t144,t145,t146,t147,t148,t149
      real*8 t150,t151,t152,t153,t154,t155,t156,t157,t158,t159
      real*8 t160,t161,t162,t163,t164,t165,t166,t167,t168,t169
      real*8 t170,t171,t172,t173,t174,t175,t176,t177,t178,t179
      real*8 t180,t181,t182,t183,t184,t185,t186,t187,t188,t189
      real*8 t190,t191,t192,t193,t194,t195,t196,t197,t198,t199
      real*8 t200,t201,t202,t203,t204,t205,t206,t207,t208,t209
      real*8 t210,t211,t212,t213,t214,t215,t216,t217,t218,t219
      real*8 t220,t221,t222,t223,t224,t225,t226,t227,t228,t229
      real*8 t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
      real*8 t240,t241,t242,t243,t244,t245,t246,t247,t248,t249
      real*8 t250,t251,t252,t253,t254,t255,t256,t257,t258,t259
      real*8 t260,t261,t262,t263,t264,t265,t266,t267,t268,t269
      real*8 t270,t271,t272,t273,t274,t275,t276,t277,t278,t279
      real*8 t280,t281,t282,t283,t284,t285,t286,t287,t288,t289
      real*8 t290,t291,t292,t293,t294,t295,t296,t297,t298,t299
      real*8 t300,t301,t302,t303,t304,t305,t306,t307,t308,t309
      real*8 t310,t311,t312,t313,t314,t315,t316,t317,t318,t319
      real*8 t320,t321,t322,t323,t324,t325,t326,t327,t328,t329
      real*8 t330,t331,t332,t333,t334,t335,t336,t337,t338,t339
      real*8 t340,t341,t342,t343,t344,t345,t346,t347,t348,t349
      real*8 t350,t351,t352,t353,t354,t355,t356,t357,t358,t359
      real*8 t360,t361,t362,t363,t364,t365,t366,t367,t368,t369
      real*8 t370,t371,t372,t373,t374,t375,t376,t377,t378,t379
      real*8 t380,t381,t382,t383,t384,t385,t386,t387,t388,t389
      real*8 t390,t391,t392,t393,t394,t395,t396,t397,t398,t399
c
      real*8 ex_rpole(13)
      real*8 sc(10)
      real*8 gl(0:4),gln(0:4)
      real*8 iv,rdn
      real*8 diff
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      character*20 list(6)
ccccccccccccc
c      real*8, allocatable :: vscale(:)
c
c
c     zero out energy and partitioning due to extra potential terms
c
c      nex = 0
      ex = 0.0d0
c      do i = 1, n
c         aex(i) = 0.0d0
c      end do
      do i = 1, 13
         ex_rpole(i) = 0.0d0
      end do
      f = electric / dielec
      if (.not. use_mpole) then
         call kmpole
         call chkpole
         call rotpole
      end if
      if (.not. use_vdw) then 
         call kvdw
      end if
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     perform dynamic allocation of some local arrays
c
c      allocate (vscale(n))
c      do i = 1, n
c         vscale(i) = 1.0d0
c      end do
c
c     add any user-defined extra potentials and partitioning
c
      do i = 1, nvdw-1
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
         if (num_ele .eq. "VALENCE") then
            if (atomic(i) .gt. 2)  zi = zi - 2.0d0
            if (atomic(i) .gt. 10)  zi = zi - 8.0d0
            if (atomic(i) .gt. 18)  zi = zi - 8.0d0
            if (atomic(i) .gt. 20)  zi = zi - 10.0d0
         else if (num_ele .eq. "TWO") then
            if (atomic(i) .ne. 1) zi = ci + 2.0d0
         else if (num_ele .eq. "TWOPLUS") then
            if (atomic(i) .ne. 1) zi = 2.0d0
         else if (num_ele .eq. "TWO-VARIABLE") then
            if (atomic(i) .ne. 1) zi = ci + 2.0d0
         else if (num_ele .eq. "TWO-VARIABLE-H") then
            if (atomic(i) .ne. 1) zi = ci + 2.0d0
         else if (num_ele .eq. "ZERO") then
            zi = ci
         else if (num_ele .eq. "VARIABLE") then
            if (use_vdwclass) then
               zi = ci + val_ele(vdwclass(i))
            else
               zi = ci + val_ele(cpclass(i))
            end if
         else if (num_ele .eq. "VARIABLE-H") then
            if (use_vdwclass) then
               if (atomic(i) .ne. 1) zi = ci + val_ele(vdwclass(i))
            else
               if (atomic(i) .ne. 1) zi = ci + val_ele(cpclass(i))
            end if
         end if
         qi = ci - zi
c
         if (use_vdwclass) then
            overlapi = overlap(vdwclass(ii))
            boverlapi = boverlap(vdwclass(ii))
            overlapri = overlapr(vdwclass(ii))
            coverlapi = coverlap(vdwclass(ii))
            alphai = alpha(vdwclass(ii))
         else
            overlapi = overlap(cpclass(ii))
            boverlapi = boverlap(cpclass(ii))
            overlapri = overlapr(cpclass(ii))
            coverlapi = coverlap(cpclass(ii))
            alphai = alpha(cpclass(ii))
         end if
         if ((cptype .eq. "EX-DAMPING").or.
     &        (cptype .eq. "EX-DAMPING1")) then
            if (boverlapi .eq. 0.0d0) boverlapi = overlapi
            if (coverlapi .eq. 0.0d0) coverlapi = overlapi
c            qi = sqrt(overlapi)*qi
c            dix = sqrt(boverlapi)*dix
c            diy = sqrt(boverlapi)*diy
c            diz = sqrt(boverlapi)*diz
c            qixx = sqrt(coverlapi)*qixx
c            qixy = sqrt(coverlapi)*qixy
c            qixz = sqrt(coverlapi)*qixz
c            qiyy = sqrt(coverlapi)*qiyy
c            qiyz = sqrt(coverlapi)*qiyz
c            qizz = sqrt(coverlapi)*qizz
            qi = overlapi*qi
            dix = boverlapi*dix
            diy = boverlapi*diy
            diz = boverlapi*diz
            qixx = coverlapi*qixx
            qixy = coverlapi*qixy
            qixz = coverlapi*qixz
            qiyy = coverlapi*qiyy
            qiyz = coverlapi*qiyz
            qizz = coverlapi*qizz
         else if (cptype .eq. "EX-ORBITAL1-SQRT") then
c     flip the signs to make electron density positive
c     take the sqrt of the magnitude, but keep direction
c     (you can't do that)
c     you can if you normalize, then multipy each component
c     by the sqrt(old magnitude)
            qi = sqrt(-qi)
c     get magnitudes of full tensors
            dip_mag = sqrt(rpole(2,i)**2 + rpole(3,i)**2 + 
     &           rpole(4,i)**2)
            quad_mag = sqrt(rpole(5,i)**2 + rpole(6,i)**2 + 
     &           rpole(7,i)**2 + rpole(8,i)**2 + rpole(9,i)**2 +
     &           rpole(10,i)**2 + rpole(11,i)**2 + rpole(12,i)**2 +
     &           rpole(13,i)**2)
c     normalize dipole vector and flip sign
            do m = 2, 4
               ex_rpole(m) = -rpole(m,i)/dip_mag
            end do
c     normalize quadrupole tensor and flip sign
            do m = 5, 13
               ex_rpole(m) = -rpole(m,i)/quad_mag
            end do
c     multiple by sqrt(mag) to get new tensors
            dix = ex_rpole(2)*sqrt(dip_mag)
            diy = ex_rpole(3)*sqrt(dip_mag)
            diz = ex_rpole(4)*sqrt(dip_mag)
            qixx = ex_rpole(5)*sqrt(quad_mag)
            qixy = ex_rpole(6)*sqrt(quad_mag)
            qixz = ex_rpole(7)*sqrt(quad_mag)
            qiyy = ex_rpole(9)*sqrt(quad_mag)
            qiyz = ex_rpole(10)*sqrt(quad_mag)
            qizz = ex_rpole(13)*sqrt(quad_mag)
         else if (cptype .eq. "EX-ORBITAL1-ABS") then
            qi = sqrt(abs(qi))
            dip_mag = sqrt(rpole(2,i)**2 + rpole(3,i)**2 +
     &           rpole(4,i)**2)
            quad_mag = sqrt(rpole(5,i)**2 + rpole(6,i)**2 +
     &           rpole(7,i)**2 + rpole(8,i)**2 + rpole(9,i)**2 +
     &           rpole(10,i)**2 + rpole(11,i)**2 + rpole(12,i)**2+
     &           rpole(13,i)**2)
            do m = 2, 4
               ex_rpole(m) = abs(rpole(m,i))/dip_mag
            end do
            do m = 5, 13
               ex_rpole(m) = abs(rpole(m,i))/quad_mag
            end do
            dix = ex_rpole(2)*sqrt(dip_mag)
            diy = ex_rpole(3)*sqrt(dip_mag)
            diz = ex_rpole(4)*sqrt(dip_mag)
            qixx = ex_rpole(5)*sqrt(quad_mag)
            qixy = ex_rpole(6)*sqrt(quad_mag)
            qixz = ex_rpole(7)*sqrt(quad_mag)
            qiyy = ex_rpole(9)*sqrt(quad_mag)
            qiyz = ex_rpole(10)*sqrt(quad_mag)
            qizz = ex_rpole(13)*sqrt(quad_mag)
         end if
c         do j = 1, n12(i)
c            vscale(i12(j,i)) = v2scale
c         end do
c         do j = 1, n13(i)
c            vscale(i13(j,i)) = v3scale
c         end do
c         do j = 1, n14(i)
c            vscale(i14(j,i)) = v4scale
c            iv14(i14(j,i)) = i
c         end do
c         do j = 1, n15(i)
c            vscale(i15(j,i)) = v5scale
c         end do
         do k = i+1, nvdw
            kk = ivdw(k)
            ck = rpole(1,k)
c            xr = x(i) - x(k)
c            yr = y(i) - y(k)
c            zr = z(i) - z(k)
c
c            xr = x(k) - x(i)
c            yr = y(k) - y(i)
c            zr = z(k) - z(i)
            xr = xred(k) - xred(i)
            yr = yred(k) - yred(i)
            zr = zred(k) - zred(i)
c
            r2 = xr*xr + yr*yr + zr*zr
            r = sqrt(r2)
c            overlapi = overlap(vdwclass(ii))
c            overlapk = overlap(vdwclass(kk))
c            boverlapi = boverlap(vdwclass(ii))
c            boverlapk = boverlap(vdwclass(kk))
            if (use_vdwclass) then
               overlapk = overlap(vdwclass(kk))
               boverlapk = boverlap(vdwclass(kk))
               overlaprk = overlapr(vdwclass(kk))
               coverlapk = coverlap(vdwclass(kk))
c     comment this out if you don't want to fit alpha's
               alphak = alpha(vdwclass(kk))
            else
               overlapk = overlap(cpclass(kk))
               boverlapk = boverlap(cpclass(kk))
               overlaprk = overlapr(cpclass(kk))
               coverlapk = coverlap(cpclass(kk))
               alphak = alpha(cpclass(kk))
            end if
            oik = overlapi * overlapk
c            boik = boverlapi * boverlapk
            boik = (boverlapi + boverlapk) / 2.0d0
            coik = (coverlapi + coverlapk) / 2.0d0
            roik = (overlapri + overlaprk) / 2.0d0
            sigma = (overlapri + overlaprk) / 2.0d0
            if (singoverlap) then
c               boik = soverlap
               sigma = soverlap
            end if
            if (exmix .eq. 'ARITHMETIC') then
               oik = 0.5d0 * (overlapi + overlapk)
            else if (exmix .eq. 'GEOMETRIC') then
c               oik = sqrt(overlapi) * sqrt(overlapk)
               boik = sqrt(boverlapi*boverlapk)
               oik = sqrt(overlapi*overlapk)
ccccccccccccccccccccccccccccccc HACK!!!!!!!!!
c               alphamix = 2.0d0 * alphai * alphak / (alphai + alphak)
c               termi = overlapi**(1.0d0/alphai)
c               termk = overlapk**(1.0d0/alphak)
c               oik = (overlapi*overlapk)**(alphamix/2.0d0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            else if (exmix .eq. 'HARMONIC') then
               oik = 2.0d0 * (overlapi*overlapk)/(overlapi + overlapk)
            else if (exmix .eq. 'HHG') then
               oik = 4.0d0 * (overlapi*overlapk) / 
     &              (sqrt(overlapi) + sqrt(overlapk))**2
            else if (exmix .eq. "W-H") then
               oik2 = overlapi * overlapk
               oik = sqrt(oik2)
            else if (exmix .eq."W-H1") then
               termi = 2.0d0*(alphai**3)*(alphak**3)
               termk = (alphai**6) + (alphak**6)
               oik = sqrt(overlapi*overlapk)*termi/termk
c               print *,"wh1",termi,termk,oik
            else if (exmix .eq."W-H2") then
               invalphai = 1.0d0/alphai
               invalphak = 1.0d0/alphak
               termi = 2.0d0*(invalphai**3)*(invalphak**3)
               termk = (invalphai**6) + (invalphak**6)
               oik = sqrt(overlapi*overlapk)*termi/termk
c               print *,"wh2",termi,termk,oik
c     bohm ahlrichs
            else if (exmix .eq. "BA") then
               alphamix = 2.0d0 * alphai * alphak / (alphai + alphak)
               termi = overlapi**(1.0d0/alphai)
               termk = overlapk**(1.0d0/alphak)
c               oik = (overlapi*overlapk)**(alphamix/2.0d0)
               oik = (termi*termk)**(alphamix/2.0d0)
c
c               alphamix = 2.0d0 * overlapri * overlaprk / 
c     &              (overlapri + overlaprk)
c               termi = overlapi**(1.0d0/overlapri)
c               termk = overlapk**(1.0d0/overlaprk)
c               boik = (boverlapi*boverlapk)**(alphamix/2.0d0)
c               boik = sqrt(boverlapi*boverlapk)
            end if
            zk = dble(atomic(k))
            if (num_ele .eq. "VALENCE") then
               if (atomic(k) .gt. 2)  zk = zk - 2.0d0
               if (atomic(k) .gt. 10)  zk = zk - 8.0d0
               if (atomic(k) .gt. 18)  zk = zk - 8.0d0
               if (atomic(k) .gt. 20)  zk = zk - 10.0d0
            else if (num_ele .eq. "TWO") then
               if (atomic(k) .ne. 1) zk = ck + 2.0d0
            else if (num_ele .eq. "TWOPLUS") then
               if (atomic(k) .ne. 1) zk = 2.0d0
            else if (num_ele .eq. "TWO-VARIABLE") then
               if (atomic(k) .ne. 1) zk = ck + 2.0d0
            else if (num_ele .eq. "TWO-VARIABLE-H") then
               if (atomic(k) .ne. 1) zk = ck + 2.0d0
            else if (num_ele .eq. "ZERO") then
               zk = ck
            else if (num_ele .eq. "VARIABLE") then
               if (use_vdwclass) then
                  zk = ck + val_ele(vdwclass(k))
               else
                  zk = ck + val_ele(cpclass(k))
               end if
            else if (num_ele .eq. "VARIABLE-H") then
               if (use_vdwclass) then
                  if (atomic(k) .ne. 1) zk = ck + val_ele(vdwclass(k))
               else
                  if (atomic(k) .ne. 1) zk = ck + val_ele(cpclass(k))
               end if
            end if
            qk = ck - zk
            dampi = alphai*r
            dampk = alphak*r
            expdampi = exp(-dampi)
            expdampk = exp(-dampk)
            pre = qi*qk*(alphai**2)*(alphak**2)
            diff = abs(alphai - alphak)
c            if (alphai.ne.alphak) then
            if (diff .gt. 0.0001d0) then
               term = alphai**2 - alphak**2
c               e = (pre * oik * 4.0d0 * pi / (r*term))*
c     &              (expdampk - expdampi)
c               e = (oik + boik*r)*(pre * 4.0d0 * pi / (r*term))*
c     &              (expdampk - expdampi)
c               e = oik*(pre * 4.0d0 * pi / (r*term))*
c     &                 (expdampk - expdampi)
               e = oik*(2.0*pre / (r*term))*
     &                 (expdampk - expdampi)
               if (exmodel .eq. "RNEW") then
                  e = (oik + boik*r)*(pre * 4.0d0 * pi / (r*term))*
     &                 (expdampk - expdampi)
               else if (exmodel .eq. "R2") then
c                  e = oik*(1.0d0+boik/r)*(pre * 4.0d0 * pi / (r*term))*
c     &                 (expdampk - expdampi)
                  e = (oik + boik/r + coik/r2)*
     &                 (pre * 4.0d0 * pi / (r*term))*
     &                 (expdampk - expdampi)
               else if (exmodel .eq. "EXP") then
                  e = oik*(1.0d0 + boik*exp(-sigma*r))*
     &                 (pre * 4.0d0 * pi / (r*term))*
     &                 (expdampk - expdampi)
               else if (exmodel .eq. "EXP2") then
                  e = oik*(pre * 4.0d0 * pi / (r*term))*
     &                 (expdampk - expdampi) + 
     &                 boik*pre*pi*exp(-sigma*r)
               else if (exmodel .eq. "GORDON1") then
c                  e = oik*pre*(8.0d0/(r*term**3))*
c     &                 (alphai*(r*term - 4.0d0*alphak)*expdampk +
c     &                 alphak*(r*term + 4.0d0*alphai)*expdampi)
                  e = oik*pre*(8.0d0/(r*term**2))*
     &                 (alphai*(r - 4.0d0*alphak/term)*expdampk +
     &                 alphak*(r + 4.0d0*alphai/term)*expdampi)
               end if
c               e = (oik * 8.0d0 / (r * term**3))*
c     &              ((alphai*(r*term - 4.0d0*alphak))*expdampk +
c     &              (alphak*(r*term + 4.0d0*alphai))*expdampi)
            else
c               e = (oik / alphai**3)*(1.0d0 + dampi +
c     &              (1.0d0/3.0d0)*dampi**2)*expdampi
c               e = (pre * oik * 2.0d0 * pi / alphai)*expdampi
c               e = (oik + boik*r)*(pre * 2.0d0 * pi / alphai)*expdampi
c               e = oik*(pre * 2.0d0 * pi / alphai)
c     &                 *expdampi
               e = oik*(pre / alphai)
     &                 *expdampi
               if (exmodel .eq. "RNEW") then
                  e = (oik + boik*r)*(pre * 2.0d0 * pi / alphai)
     &                 *expdampi
               else if (exmodel .eq. "R2") then
c                  e = oik*(1.0d0 + boik/r)*(pre * 2.0d0 * pi / alphai)
c     &                 *expdampi
                  e = (oik + boik/r + coik/r2)*
     &                 (pre * 2.0d0 * pi / alphai)
     &                 *expdampi
               else if (exmodel .eq. "EXP") then
                  e = oik*(1.0d0 + boik*exp(-sigma*r))*
     &                 (pre * 2.0d0 * pi / alphai)*expdampi
               else if (exmodel .eq. "EXP2") then
                  e = oik*(pre * 2.0d0 * pi / alphai)
     &                 *expdampi +
     &                 boik*pre*pi*exp(-sigma*r)
               else if (exmodel .eq. "GORDON1") then
                  e = oik*pre*(1.0d0/(alphai**3))*(1.0d0 + dampi + 
     &                 (1.0d0/3.0d0)*(dampi)**2)*expdampi
               end if
            end if
            if (exmodel .eq. "SLATER") then
               roik = sqrt(overlapri*overlaprk)
               e = oik*(1.0d0 + roik*r + (1.0d0/3.0d0)*(roik*r)**2)*
     &              exp(-roik*r)
            end if
            if (exmodel .eq. "BORN-MAYER") then
               roik = (overlapri*overlaprk*(overlapri + overlaprk))/
     &              (overlapri**2 + overlaprk**2)
               e = oik*exp(-roik*r)
            end if
            if (exmodel .eq. "R7") then
               e = 12345.6789d0
            end if
c           e = e*vscale(kk)
            list = ( /"CP","COMBO","COMBO2","COMBO3","COMBO4",
     &           "COMBO5" /)
            if (any(list .eq. exmodel)) then
c            if (((exmodel .eq. "CP").or.(exmodel .eq. "COMBO")).or.
c     &           ((exmodel.eq."COMBO2").or.(exmodel.eq."COMBO3"))) then
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               if ((cptype .eq. "EX-DAMPING").or.
     &              (cptype .eq. "EX-DAMPING1")) then
                  if (boverlapk .eq. 0.0d0) boverlapk = overlapk
                  if (coverlapk .eq. 0.0d0) coverlapk = overlapk
c                  qk = sqrt(overlapk)*qk
c                  dkx = sqrt(boverlapk)*dkx
c                  dky = sqrt(boverlapk)*dky
c                  dkz = sqrt(boverlapk)*dkz
c                  qkxx = sqrt(coverlapk)*qkxx
c                  qkxy = sqrt(coverlapk)*qkxy
c                  qkxz = sqrt(coverlapk)*qkxz
c                  qkyy = sqrt(coverlapk)*qkyy
c                  qkyz = sqrt(coverlapk)*qkyz
c                  qkzz = sqrt(coverlapk)*qkzz
                  qk = overlapk*qk
                  dkx = boverlapk*dkx
                  dky = boverlapk*dky
                  dkz = boverlapk*dkz
                  qkxx = coverlapk*qkxx
                  qkxy = coverlapk*qkxy
                  qkxz = coverlapk*qkxz
                  qkyy = coverlapk*qkyy
                  qkyz = coverlapk*qkyz
                  qkzz = coverlapk*qkzz
               else if (cptype .eq. "EX-ORBITAL1-SQRT") then
                  qk = sqrt(-qk)
c     get magnitudes of full tensors
                  dip_mag = sqrt(rpole(2,k)**2 + rpole(3,k)**2 + 
     &                 rpole(4,k)**2)
                  quad_mag = sqrt(rpole(5,k)**2 + rpole(6,k)**2 + 
     &                 rpole(7,k)**2 + rpole(8,k)**2 + rpole(9,k)**2 +
     &                 rpole(10,k)**2 + rpole(11,k)**2 + rpole(12,k)**2+
     &                 rpole(13,k)**2)
c     normalize dipole vector and flip sign
                  do m = 2, 4
                     ex_rpole(m) = -rpole(m,k)/dip_mag
                  end do
c     normalize quadrupole tensor and flip sign
                  do m = 5, 13
                     ex_rpole(m) = -rpole(m,k)/quad_mag
                  end do
c     multiple by sqrt(mag) to get new tensors
                  dkx = ex_rpole(2)*sqrt(dip_mag)
                  dky = ex_rpole(3)*sqrt(dip_mag)
                  dkz = ex_rpole(4)*sqrt(dip_mag)
                  qkxx = ex_rpole(5)*sqrt(quad_mag)
                  qkxy = ex_rpole(6)*sqrt(quad_mag)
                  qkxz = ex_rpole(7)*sqrt(quad_mag)
                  qkyy = ex_rpole(9)*sqrt(quad_mag)
                  qkyz = ex_rpole(10)*sqrt(quad_mag)
                  qkzz = ex_rpole(13)*sqrt(quad_mag)
c
c                  print *,"qk",ck,qk
c                  print *,"k",ex_rpole
c                  print *,"pole",pole(2,k),pole(3,k),pole(4,k)
c                  print *,"real dkpole mag k",sqrt(rpole(2,k)**2 + 
c     &                 rpole(3,k)**2 + rpole(4,k)**2)
c                  print *,"dipole magnitude k",sqrt(dkx**2 + dky**2 + 
c     &                 dkz**2)
               else if (cptype .eq. "EX-ORBITAL1-ABS") then
                  qk = sqrt(abs(qk))
                  dip_mag = sqrt(rpole(2,k)**2 + rpole(3,k)**2 +
     &                 rpole(4,k)**2)
                  quad_mag = sqrt(rpole(5,k)**2 + rpole(6,k)**2 +
     &                 rpole(7,k)**2 + rpole(8,k)**2 + rpole(9,k)**2 +
     &                 rpole(10,k)**2 + rpole(11,k)**2 + rpole(12,k)**2+
     &                 rpole(13,k)**2)
                  do m = 2, 4
                     ex_rpole(m) = abs(rpole(m,k))/dip_mag
                  end do
                  do m = 5, 13
                     ex_rpole(m) = abs(rpole(m,k))/quad_mag
                  end do
                  dkx = ex_rpole(2)*sqrt(dip_mag)
                  dky = ex_rpole(3)*sqrt(dip_mag)
                  dkz = ex_rpole(4)*sqrt(dip_mag)
                  qkxx = ex_rpole(5)*sqrt(quad_mag)
                  qkxy = ex_rpole(6)*sqrt(quad_mag)
                  qkxz = ex_rpole(7)*sqrt(quad_mag)
                  qkyy = ex_rpole(9)*sqrt(quad_mag)
                  qkyz = ex_rpole(10)*sqrt(quad_mag)
                  qkzz = ex_rpole(13)*sqrt(quad_mag)
               end if
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
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
               if (gtype .eq. "NEW") then
                  scale1i = 1.0d0 - (1.0d0 + dampi/2.0d0)*expdampi
                  scale1k = 1.0d0 - (1.0d0 + dampk/2.0d0)*expdampk
                  scale3i = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2)
     &                 *expdampi
                  scale3k = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk**2)
     &                 *expdampk
                  scale5i = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &                 (1.0d0/6.0d0)*dampi**3)*expdampi
                  scale5k = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk**2 + 
     &                 (1.0d0/6.0d0)*dampk**3)*expdampk
               else
                  scale1i = 1.0d0 - expdampi
                  scale3i = 1.0d0 - (1.0d0 + dampi)*expdampi
                  scale5i = 1.0d0 - (1.0d0 + dampi + 
     &                 (1.0d0/3.0d0)*dampi**2)*expdampi
                  scale1k = 1.0d0 - expdampk
                  scale3k = 1.0d0 - (1.0d0 +dampk)*expdampk
                  scale5k = 1.0d0 - (1.0d0 +dampk +
     &                 (1.0d0/3.0d0)*dampk**2)*expdampk
               end if
               if (penetration .eq. "GORDON") then
                  if (alphai .ne. alphak) then
                     termi = alphak**2/(alphak**2 - alphai**2)
                     termk = alphai**2/(alphai**2 - alphak**2)
                     if (gtype .eq. "NEW") then
                        scale1ik = 1.0d0 - (termi**2)*
     &                       (1.0d0+2.0d0*termk+0.5d0*dampi)*expdampi -
     &                       (termk**2)*(1.0d0 + 2.0d0*termi + 
     &                       0.5d0*dampk)*expdampk
                        scale3ik = 1.0d0 - (termi**2)*(1.0d0 + dampi + 
     &                       0.5d0*dampi**2)*expdampi - 
     &                       (termk**2)*(1.0d0 + dampk + 0.5d0*dampk**2)
     &                       *expdampk
     &                       - 2.0d0*(termi**2)*termk*(1.0d0 + dampi)
     &                       *expdampi
     &                       - 2.0d0*(termk**2)*termi*(1.0d0 + dampk)
     &                       *expdampk
                        scale5ik = 1.0d0 - (termi**2)*
     &                       (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &                       (1.0d0/6.0d0)*dampi**3)*expdampi - 
     &                       (termk**2)*
     &                       (1.0d0 + dampk + 0.5d0*dampk**2 + 
     &                       (1.0d0/6.0d0)*dampk**3)*expdampk - 
     &                       2.0d0*(termi**2)*termk*(1.0 + dampi + 
     &                       (1.0d0/3.0d0)*dampi**2)*expdampi - 
     &                       2.0d0*(termk**2)*termi*(1.0 + dampk +
     &                       (1.0d0/3.0d0)*dampk**2)*expdampk
                        scale7ik = 1.0d0 - (termi**2)*
     &                       (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &                       (1.0d0/6.0d0)*dampi**3 + 
     &                       (1.0d0/30.0d0)*dampi**4)*
     &                       expdampi - 
     &                       (termk**2)*(1.0d0 + dampk + 0.5d0*dampk**2+ 
     &                       (1.0d0/6.0d0)*dampk**3 + 
     &                       (1.0d0/30.0d0)*dampk**4)*
     &                       expdampk - 
     &                       2.0d0*(termi**2)*termk*(1.0d0 + dampi +
     &                       (2.0d0/5.0d0)*dampi**2 + 
     &                       (1.0d0/15.0d0)*dampi**3)*
     &                       expdampi - 
     &                       2.0d0*(termk**2)*termi*(1.0d0 + dampk + 
     &                       (2.0d0/5.0d0)*dampk**2 + 
     &                       (1.0d0/15.0d0)*dampk**3)*
     &                       expdampk
                        scale9ik = 1.0d0 - (termi**2)*
     &                       (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &                       (1.0d0/6.0d0)*dampi**3 + 
     &                       (4.0d0/105.0d0)*dampi**4 +
     &                       (1.0d0/210.0d0)*dampi**5)*expdampi - 
     &                       (termk**2)*
     &                       (1.0d0 + dampk + 0.5d0*dampk**2 +
     &                       (1.0d0/6.0d0)*dampk**3 + 
     &                       (4.0d0/105.0d0)*dampk**4 +
     &                       (1.0d0/210.0d0)*dampk**5)*expdampk -
     &                       2.0d0*(termi**2)*termk*
     &                       (1.0d0 + dampi + (3.0d0/7.0d0)*dampi**2 + 
     &                       (2.0d0/21.0d0)*dampi**3 + 
     &                       (1.0d0/105.0d0)*dampi**4)*
     &                       expdampi - 
     &                       2.0d0*(termk**2)*termi*
     &                       (1.0d0 + dampk + (3.0d0/7.0d0)*dampk**2 +
     &                       (2.0d0/21.0d0)*dampk**3 + 
     &                       (1.0d0/105.0d0)*dampk**4)*
     &                       expdampk
c                     else if (cptype .eq. "EX-DAMPING") then
c                        term = 2.0d0*(alphai**2)*(alphak**2)
c     &                       /(alphai**2 - alphak**2)
c                        scale1ik = term*(expdampk - expdampi)
c                        scale3ik = term*(1.0d0 + dampk)*expdampk - 
c     &                       term*(1.0d0 + dampi)*expdampi
c                        scale5ik = term*(1.0d0 + dampk + 
c     &                       (1.0d0/3.0d0)*dampk**2)*expdampk -
c     &                       term*(1.0d0 + dampi + 
c     &                       (1.0d0/3.0d0)*dampi**2)*expdampi
c                        scale7ik = term*(1.0d0 + dampk +
c     &                       (2.0d0/5.0d0)*dampk**2 +
c     &                       (1.0d0/15.0d0)*dampk**3)*expdampk -
c     &                       term*(1.0d0 + dampi +
c     &                       (2.0d0/5.0d0)*dampi**2 +
c     &                       (1.0d0/15.0d0)*dampi**3)*expdampi
c                        scale9ik = term*(1.0d0 + dampk +
c     &                       (3.0d0/7.0d0)*dampk**2 +
c     &                       (2.0d0/21.0d0)*dampk**3 + 
c     &                       (1.0d0/105.0d0)*dampk**4)*expdampk -
c     &                       term*(1.0d0 + dampi +
c     &                       (3.0d0/7.0d0)*dampi**2 +
c     &                       (2.0d0/21.0d0)*dampi**3 +
c     &                       (1.0d0/105.0d0)*dampi**4)*expdampi
                     else
                        scale1ik = 1.0d0 -termi*expdampi -termk*expdampk
                        scale3ik = 1.0d0 - termi*(1.0d0 +dampi)*expdampi 
     &                       - termk*(1.0d0 + dampk)*expdampk
                        scale5ik = 1.0d0 - termi*(1.0d0 + dampi + 
     &                       (1.0d0/3.0d0)*dampi**2)*expdampi - 
     &                       termk*(1.0d0 + dampk + 
     &                       (1.0d0/3.0d0)*dampk**2)*expdampk
                        scale7ik = 1.0d0 - termi*(1.0d0 + dampi +
     &                       0.4d0*dampi**2 + (1.0d0/15.0d0)*dampi**3)*
     &                       expdampi -
     &                       termk*(1.0d0 + dampk +
     &                       0.4d0*dampk**2 + (1.0d0/15.0d0)*dampk**3)*
     &                       expdampk
                        scale9ik = 1.0d0 - termi*(1.0d0 + dampi + 
     &                       (3.0d0/7.0d0)*dampi**2 + 
     &                       (2.0d0/21.0d0)*dampi**3 + 
     &                       (1.0d0/105.0d0)*dampi**4)*expdampi -
     &                       termk*(1.0d0 + dampk +
     &                       (3.0d0/7.0d0)*dampk**2 +
     &                       (2.0d0/21.0d0)*dampk**3 +
     &                       (1.0d0/105.0d0)*dampk**4)*expdampk
                     end if
                  else
                     if (gtype .eq. "NEW") then
                        scale1ik = 1.0d0 - (1.0d0 + 
     &                       (11.0d0/16.0d0)*dampi + 
     &                       (3.0d0/16.0d0)*dampi**2 + 
     &                       (1.0d0/48.0d0)*dampi**3)
     &                       *expdampi
                        scale3ik = 1.0d0 - (1.0d0 + dampi + 
     &                       0.5d0*dampi**2 + 
     &                       (7.0d0/48.0d0)*dampi**3 + 
     &                       (1.0d0/48.0d0)*dampi**4)
     &                       *expdampi
                        scale5ik = 1.0d0 - (1.0d0 + dampi + 
     &                       0.5d0*dampi**2 + 
     &                       (1.0d0/6.0d0)*dampi**3 + 
     &                       (1.0d0/24.0d0)*dampi**4 +
     &                       (1.0d0/144.0d0)*dampi**5)*expdampi
                        scale7ik = 1.0d0 - (1.0d0 + dampi + 
     &                       0.5d0*dampi**2 + 
     &                       (1.0d0/6.0d0)*dampi**3 + 
     &                       (1.0d0/24.0d0)*dampi**4 + 
     &                       (1.0d0/120.0d0)*dampi**5 + 
     &                       (1.0d0/720.0d0)*dampi**6)
     &                       *expdampi
                        scale9ik = 1.0d0 - (1.0d0 + dampi + 
     &                       0.5d0*dampi**2 + 
     &                       (1.0d0/6.0d0)*dampi**3 + 
     &                       (1.0d0/24.0d0)*dampi**4 + 
     &                       (1.0d0/120.0d0)*dampi**5 + 
     &                       (1.0d0/720.0d0)*dampi**6 + 
     &                       (1.0d0/5040.0d0)*dampi**7)*expdampi
c                     else if (cptype .eq. "EX-DAMPING") then
c     insert a=b terms for ex-damping
c                        scale1ik = (alphai**3)*r*expdampi
c                        scale3ik = (alphai**4)*r2*expdampi
c                        scale5ik = (alphai**4)*((1.0d0/3.0d0)*r2 + 
c     &                       (1.0d0/3.0d0)*alphai*r**3)*expdampi
c                        scale7ik = (alphai**4)*((1.0d0/5.0d0)*r2 + 
c     &                       (1.0d0/5.0d0)*alphai*r**3 +
c     &                       (1.0d0/15.0d0)*(alphai**2)*r**4)*expdampi
c                        scale9ik = (alphai**4)*((1.0d0/7.0d0)*r2 +
c     &                       (1.0d0/7.0d0)*alphai*r**3 +
c     &                       (2.0d0/35.0d0)*(alphai**2)*r**4 +
c     &                       (1.0d0/105.0d0)*(alphai**3)*r**5)*expdampi
                     else
                        scale1ik = 1.0d0 - (1.0d0+0.5d0*dampi)*expdampi
                        scale3ik = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2)
     &                       *expdampi
                        scale5ik = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                       + (1.0d0/6.0d0)*dampi**3)*expdampi
                        scale7ik = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                       + (1.0d0/6.0d0)*dampi**3  
     &                       + (1.0d0/30.0d0)*dampi**4)*expdampi
                        scale9ik = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                       + (1.0d0/6.0d0)*dampi**3 
     &                       + (4.0d0/105.0d0)*dampi**4 
     &                       + (1.0d0/210.0d0)*dampi**5)*expdampi
                     end if
                  end if
               end if
c
c     exchange overlap damping
c
               if (cptype .eq. "EX-DAMPING") then
                  if (alphai .ne. alphak) then
                     term = 2.0d0*(alphai**2)*(alphak**2)
     &                    /(alphai**2 - alphak**2)
                     oscale1ik = term*(expdampk - expdampi)
                     oscale3ik = term*(1.0d0 + dampk)*expdampk -
     &                    term*(1.0d0 + dampi)*expdampi
                     oscale5ik = term*(1.0d0 + dampk +
     &                    (1.0d0/3.0d0)*dampk**2)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (1.0d0/3.0d0)*dampi**2)*expdampi
                     oscale7ik = term*(1.0d0 + dampk +
     &                    (2.0d0/5.0d0)*dampk**2 +
     &                    (1.0d0/15.0d0)*dampk**3)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (2.0d0/5.0d0)*dampi**2 +
     &                    (1.0d0/15.0d0)*dampi**3)*expdampi
                     oscale9ik = term*(1.0d0 + dampk +
     &                    (3.0d0/7.0d0)*dampk**2 +
     &                    (2.0d0/21.0d0)*dampk**3 +
     &                    (1.0d0/105.0d0)*dampk**4)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (3.0d0/7.0d0)*dampi**2 +
     &                    (2.0d0/21.0d0)*dampi**3 +
     &                    (1.0d0/105.0d0)*dampi**4)*expdampi
                  else
                     oscale1ik = (alphai**3)*r*expdampi
                     oscale3ik = (alphai**4)*r2*expdampi
                     oscale5ik = (alphai**4)*((1.0d0/3.0d0)*r2 +
     &                    (1.0d0/3.0d0)*alphai*r**3)*expdampi
                     oscale7ik = (alphai**4)*((1.0d0/5.0d0)*r2 +
     &                    (1.0d0/5.0d0)*alphai*r**3 +
     &                    (1.0d0/15.0d0)*(alphai**2)*r**4)*expdampi
                     oscale9ik = (alphai**4)*((1.0d0/7.0d0)*r2 +
     &                    (1.0d0/7.0d0)*alphai*r**3 +
     &                    (2.0d0/35.0d0)*(alphai**2)*r**4 +
     &                    (1.0d0/105.0d0)*(alphai**3)*r**5)*expdampi
                  end if
               else if (cptype .eq. "EX-DAMPING-EXP") then
c
c     this implements the exp augmented overlap damping scheme
c     of Par Soderhjelm, E = A*overlap*(1 + exp(-kr))
c
                  if (alphai .ne. alphak) then
                     term = 2.0d0*(alphai**2)*(alphak**2)
     &                    /(alphai**2 - alphak**2)
                     oscale1ik = term*(expdampk - expdampi)
                     oscale3ik = term*(1.0d0 + dampk)*expdampk -
     &                    term*(1.0d0 + dampi)*expdampi
                     oscale5ik = term*(1.0d0 + dampk +
     &                    (1.0d0/3.0d0)*dampk**2)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (1.0d0/3.0d0)*dampi**2)*expdampi
                     oscale7ik = term*(1.0d0 + dampk +
     &                    (2.0d0/5.0d0)*dampk**2 +
     &                    (1.0d0/15.0d0)*dampk**3)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (2.0d0/5.0d0)*dampi**2 +
     &                    (1.0d0/15.0d0)*dampi**3)*expdampi
                     oscale9ik = term*(1.0d0 + dampk +
     &                    (3.0d0/7.0d0)*dampk**2 +
     &                    (2.0d0/21.0d0)*dampk**3 +
     &                    (1.0d0/105.0d0)*dampk**4)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (3.0d0/7.0d0)*dampi**2 +
     &                    (2.0d0/21.0d0)*dampi**3 +
     &                    (1.0d0/105.0d0)*dampi**4)*expdampi
                  else
c     save for later
                     oldalpha3 = alphai**3
                     oldalpha4 = alphai**4
c
                     oscale1ik = (alphai**3)*r*expdampi
                     oscale3ik = (alphai**4)*r2*expdampi
                     oscale5ik = (alphai**4)*((1.0d0/3.0d0)*r2 +
     &                    (1.0d0/3.0d0)*alphai*r**3)*expdampi
                     oscale7ik = (alphai**4)*((1.0d0/5.0d0)*r2 +
     &                    (1.0d0/5.0d0)*alphai*r**3 +
     &                    (1.0d0/15.0d0)*(alphai**2)*r**4)*expdampi
                     oscale9ik = (alphai**4)*((1.0d0/7.0d0)*r2 +
     &                    (1.0d0/7.0d0)*alphai*r**3 +
     &                    (2.0d0/35.0d0)*(alphai**2)*r**4 +
     &                    (1.0d0/105.0d0)*(alphai**3)*r**5)*expdampi
                  end if
c
c     now compute exp augmented terms
c
                  alphai = alphai + overlapri
                  alphak = alphak + overlaprk
                  dampi = alphai*r
                  dampk = alphak*r
                  expdampi = exp(-dampi)
                  expdampk = exp(-dampk)
                  if (alphai .ne. alphak) then
c     this used "old" term with orginal alphas
                     expscale1ik = term*(expdampk - expdampi)
                     expscale3ik = term*(1.0d0 + dampk)*expdampk -
     &                    term*(1.0d0 + dampi)*expdampi
                     expscale5ik = term*(1.0d0 + dampk +
     &                    (1.0d0/3.0d0)*dampk**2)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (1.0d0/3.0d0)*dampi**2)*expdampi
                     expscale7ik = term*(1.0d0 + dampk +
     &                    (2.0d0/5.0d0)*dampk**2 +
     &                    (1.0d0/15.0d0)*dampk**3)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (2.0d0/5.0d0)*dampi**2 +
     &                    (1.0d0/15.0d0)*dampi**3)*expdampi
                     expscale9ik = term*(1.0d0 + dampk +
     &                    (3.0d0/7.0d0)*dampk**2 +
     &                    (2.0d0/21.0d0)*dampk**3 +
     &                    (1.0d0/105.0d0)*dampk**4)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (3.0d0/7.0d0)*dampi**2 +
     &                    (2.0d0/21.0d0)*dampi**3 +
     &                    (1.0d0/105.0d0)*dampi**4)*expdampi
                  else
c     use old alphas for prefactor
                     expscale1ik = oldalpha3*r*expdampi
                     expscale3ik = oldalpha3*alphai*r2*expdampi
                     expscale5ik = oldalpha3*alphai*((1.0d0/3.0d0)*r2 +
     &                    (1.0d0/3.0d0)*alphai*r**3)*expdampi
                     expscale7ik = oldalpha3*alphai*((1.0d0/5.0d0)*r2 +
     &                    (1.0d0/5.0d0)*alphai*r**3 +
     &                    (1.0d0/15.0d0)*(alphai**2)*r**4)*expdampi
                     expscale9ik = oldalpha3*alphai*((1.0d0/7.0d0)*r2 +
     &                    (1.0d0/7.0d0)*alphai*r**3 +
     &                    (2.0d0/35.0d0)*(alphai**2)*r**4 +
     &                    (1.0d0/105.0d0)*(alphai**3)*r**5)*expdampi
                  end if
               else if (cptype .eq. "EX-DAMPING1") then
                  if (alphai .ne. alphak) then
                     term = alphai**2 - alphak**2
c                     pre = 8.0d0*(alphai**3)*(alphak**3)/(r*(term**2))
                     pre = 8.0d0*(alphai**3)*(alphak**3)/(term**2) 
                     oscale1ik = pre*(alphai*(r - 4.0d0*alphak/term)
     &                    *expdampk + alphak*(r + 4.0d0*alphai/term)
     &                    *expdampi)
                     oscale3ik = pre*alphai*alphak*((r2 - 4.0d0/term -
     &                    4.0d0*dampk/term)*expdampk + (r2 + 
     &                    4.0d0/term + 4.0d0*dampi/term)*expdampi)
                     oscale5ik = pre*alphai*alphak*(((1.0d0/3.0d0)*r2 + 
     &                    (1.0d0/3.0d0)*dampk*r2 - 4.0d0/term -
     &                    4.0d0*dampk/term - 
     &                    4.0d0*(dampk**2)/(3.0d0*term))*expdampk + 
     &                    ((1.0d0/3.0d0)*r2 + 
     &                    (1.0d0/3.0d0)*dampi*r2 + 4.0d0/term + 
     &                    4.0d0*dampi/term + 
     &                    4.0d0*(dampi**2)/(3.0d0*term))*expdampi)
                     oscale7ik = pre*alphai*alphak*(((1.0d0/5.0d0)*r2 + 
     &                    (1.0d0/5.0d0)*dampk*r2 + 
     &                    (1.0d0/15.0d0)*(dampk**2)*r2 - 4.0d0/term
     &                    - 4.0d0*dampk/term  
     &                    - (8.0d0/5.0d0)*(dampk**2)/term
     &                    - (4.0d0/15.0d0)*(dampk**3)/term)
     &                    *expdampk
     &                    + ((1.0d0/5.0d0)*r2 + (1.0d0/5.0d0)*dampi*r2 +
     &                    (1.0d0/15.0d0)*(dampi**2)*r2 + 4.0d0/term +
     &                    4.0d0*dampi/term + 
     &                    (8.0d0/5.0d0)*(dampi**2)/term +
     &                    (4.0d0/15.0d0)*(dampi**3)/term)*expdampi)
                     oscale9ik = pre*alphai*alphak*(((1.0d0/7.0d0)*r2 + 
     &                    (1.0d0/7.0d0)*dampk*r2 + 
     &                    (2.0d0/35.0d0)*(dampk**2)*r2 + 
     &                    (1.0d0/105.0d0)*(dampk**3)*r2 - 
     &                    4.0d0/term - 4.0d0*dampk/term - 
     &                    (12.0d0/7.0d0)*(dampk**2)/term -
     &                    (8.0d0/21.0d0)*(dampk**3)/term - 
     &                    (4.0d0/105.0d0)*(dampk**4)/term)*expdampk + 
     &                    ((1.0d0/7.0d0)*r2 + (1.0d0/7.0d0)*dampi*r2 + 
     &                    (2.0d0/35.0d0)*(dampi**2)*r2 + 
     &                    (1.0d0/105.0d0)*(dampi**3)*r2 +
     &                    4.0d0/term + 4.0d0*dampi/term +
     &                    (12.0d0/7.0d0)*(dampi**2)/term + 
     &                    (8.0d0/21.0d0)*(dampi**3)/term +
     &                    (4.0d0/105.0d0)*(dampi**4)/term)*expdampi)
                  else
                     pre = alphai**3
                     oscale1ik = pre*r*(1.0d0 + dampi +
     &                    (1.0d0/3.0d0)*(dampi**2))*expdampi
c                     oscale1ik = pre*(1.0d0/r)*(1.0d0 + dampi + 
c     &                    (1.0d0/3.0d0)*dampi**2)*expdampi
                     oscale3ik = pre*(1.0d0/3.0d0)*((alphai**2)*r**3 + 
     &                    (alphai**3)*r**4)*expdampi
                     oscale5ik = pre*(1.0d0/9.0d0)*(alphai**4)*(r**5)*
     &                    expdampi
                     oscale7ik = pre*(1.0d0/45.0d0)*(alphai**5)*(r**6)*
     &                    expdampi
                     oscale9ik = pre*(1.0d0/315.0d0)*((alphai**5)*r**6 + 
     &                    (alphai**6)*r**7)*expdampi
                  end if
               else if (cptype .eq. "EX-DAMPING-R2") then
                  if (alphai .ne. alphak) then
                     term = 2.0d0*(alphai**2)*(alphak**2)
     &                    /(alphai**2 - alphak**2)
                     oscale1ik = term*(expdampk - expdampi)
                     oscale3ik = term*(1.0d0 + dampk)*expdampk -
     &                    term*(1.0d0 + dampi)*expdampi
                     oscale5ik = term*(1.0d0 + dampk +
     &                    (1.0d0/3.0d0)*dampk**2)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (1.0d0/3.0d0)*dampi**2)*expdampi
                     oscale7ik = term*(1.0d0 + dampk +
     &                    (2.0d0/5.0d0)*dampk**2 +
     &                    (1.0d0/15.0d0)*dampk**3)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (2.0d0/5.0d0)*dampi**2 +
     &                    (1.0d0/15.0d0)*dampi**3)*expdampi
                     oscale9ik = term*(1.0d0 + dampk +
     &                    (3.0d0/7.0d0)*dampk**2 +
     &                    (2.0d0/21.0d0)*dampk**3 +
     &                    (1.0d0/105.0d0)*dampk**4)*expdampk -
     &                    term*(1.0d0 + dampi +
     &                    (3.0d0/7.0d0)*dampi**2 +
     &                    (2.0d0/21.0d0)*dampi**3 +
     &                    (1.0d0/105.0d0)*dampi**4)*expdampi
                  else
                     oscale1ik = (alphai**3)*r*expdampi
                     oscale3ik = (alphai**4)*r2*expdampi
                     oscale5ik = (alphai**4)*((1.0d0/3.0d0)*r2 +
     &                    (1.0d0/3.0d0)*alphai*r**3)*expdampi
                     oscale7ik = (alphai**4)*((1.0d0/5.0d0)*r2 +
     &                    (1.0d0/5.0d0)*alphai*r**3 +
     &                    (1.0d0/15.0d0)*(alphai**2)*r**4)*expdampi
                     oscale9ik = (alphai**4)*((1.0d0/7.0d0)*r2 +
     &                    (1.0d0/7.0d0)*alphai*r**3 +
     &                    (2.0d0/35.0d0)*(alphai**2)*r**4 +
     &                    (1.0d0/105.0d0)*(alphai**3)*r**5)*expdampi
                  end if
               else if ((cptype .eq. "EX-ORBITAL1").or.(cptype .eq. 
     &                 "EX-ORBITAL1-ABS")) then
c     need to add in if ex-orbital1-sqrt
                  if (use_vdwclass) then
                     alphai2 = 0.5d0 * alpha(vdwclass(ii))
                     alphak2 = 0.5d0 * alpha(vdwclass(kk))
                  else
                     alphai2 = 0.5d0 * alpha(cpclass(ii))
                     alphak2 = 0.5d0 * alpha(cpclass(kk))
                  end if
                  dampi = alphai2*r
                  dampk = alphak2*r
                  expdampi = exp(-dampi)
                  expdampk = exp(-dampk)
                  if (alphai .ne. alphak) then
                     term = alphai2**2 - alphak2**2
                     pre = 8.0d0*((sqrt(alphai))**3)*((sqrt(alphak))**3)
     &                    /(term**2)
                     oscale1ik = pre*(alphai2*(r - 4.0d0*alphak2/term)
     &                    *expdampk + alphak2*(r + 4.0d0*alphai2/term)
     &                    *expdampi)
                     oscale3ik = pre*alphai2*alphak2*((r2 - 4.0d0/term -
     &                    4.0d0*dampk/term)*expdampk + (r2 +
     &                    4.0d0/term + 4.0d0*dampi/term)*expdampi)
                     oscale5ik = pre*alphai2*alphak2*(((1.0d0/3.0d0)*r2+
     &                    (1.0d0/3.0d0)*dampk*r2 - 4.0d0/term -
     &                    4.0d0*dampk/term -
     &                    4.0d0*(dampk**2)/(3.0d0*term))*expdampk +
     &                    ((1.0d0/3.0d0)*r2 +
     &                    (1.0d0/3.0d0)*dampi*r2 + 4.0d0/term +
     &                    4.0d0*dampi/term +
     &                    4.0d0*(dampi**2)/(3.0d0*term))*expdampi)
                     oscale7ik = pre*alphai2*alphak2*(((1.0d0/5.0d0)*r2+
     &                    (1.0d0/5.0d0)*dampk*r2 +
     &                    (1.0d0/15.0d0)*(dampk**2)*r2 - 4.0d0/term
     &                    - 4.0d0*dampk/term
     &                    - (8.0d0/5.0d0)*(dampk**2)/term
     &                    - (4.0d0/15.0d0)*(dampk**3)/term)
     &                    *expdampk
     &                    + ((1.0d0/5.0d0)*r2 + (1.0d0/5.0d0)*dampi*r2 +
     &                    (1.0d0/15.0d0)*(dampi**2)*r2 + 4.0d0/term +
     &                    4.0d0*dampi/term +
     &                    (8.0d0/5.0d0)*(dampi**2)/term +
     &                    (4.0d0/15.0d0)*(dampi**3)/term)*expdampi)
                     oscale9ik = pre*alphai2*alphak2*(((1.0d0/7.0d0)*r2+
     &                    (1.0d0/7.0d0)*dampk*r2 +
     &                    (2.0d0/35.0d0)*(dampk**2)*r2 +
     &                    (1.0d0/105.0d0)*(dampk**3)*r2 -
     &                    4.0d0/term - 4.0d0*dampk/term -
     &                    (12.0d0/7.0d0)*(dampk**2)/term -
     &                    (8.0d0/21.0d0)*(dampk**3)/term -
     &                    (4.0d0/105.0d0)*(dampk**4)/term)*expdampk +
     &                    ((1.0d0/7.0d0)*r2 + (1.0d0/7.0d0)*dampi*r2 +
     &                    (2.0d0/35.0d0)*(dampi**2)*r2 +
     &                    (1.0d0/105.0d0)*(dampi**3)*r2 +
     &                    4.0d0/term + 4.0d0*dampi/term +
     &                    (12.0d0/7.0d0)*(dampi**2)/term +
     &                    (8.0d0/21.0d0)*(dampi**3)/term +
     &                    (4.0d0/105.0d0)*(dampi**4)/term)*expdampi)
c
                  else
                     pre = (alphai**3)/(alphai2**3)
                     oscale1ik = pre*r*(1.0d0 + dampi +
     &                    (1.0d0/3.0d0)*(dampi**2))*expdampi
                     oscale3ik = pre*(1.0d0/3.0d0)*((alphai2**2)*r**3 +
     &                    (alphai2**3)*r**4)*expdampi
                     oscale5ik = pre*(1.0d0/9.0d0)*(alphai2**4)*(r**5)*
     &                    expdampi
                     oscale7ik = pre*(1.0d0/45.0d0)*(alphai2**5)*(r**6)*
     &                    expdampi
                     oscale9ik = pre*(1.0d0/315.0d0)*((alphai2**5)*r**6+ 
     &                    (alphai2**6)*r**7)*expdampi
                  end if
c
c     we need ((1/r)*fdamp)**2
c
c                  oscale1ik = (oscale1ik**2)
c                  oscale3ik = (oscale3ik**2)
c                  oscale5ik = (oscale5ik**2)
c                  oscale7ik = (oscale7ik**2)
c                  oscale9ik = (oscale9ik**2)
c
               else if (cptype .eq. "SIMPLE") then
                  damp = roik*r
                  expdamp = exp(-damp)
c
                  oscale1ik = r*expdamp
                  oscale3ik = r*damp*expdamp
                  oscale5ik = r*((1.0d0/3.0d0) + (1.0d0/3.0d0)*damp)
     &                 *damp*expdamp
                  oscale7ik = r*((1.0d0/5.0d0) + (1.0d0/5.0d0)*damp + 
     &                 (1.0d0/15.0d0)*damp**2)*damp*expdamp
                  oscale9ik = r*((1.0d0/7.0d0) + (1.0d0/7.0d0)*damp +
     &                 (2.0d0/35.0d0)*damp**2 + (1.0d0/105.0d0)*damp**3)
     &                 *damp*expdamp
c     
               else if (cptype .eq. "ORBITAL2R") then
c
c     this implements S^2/R for exchange repulsion
c     
                  if (use_vdwclass) then
                     alphai2 = alpha(vdwclass(ii))
                     alphak2 = alpha(vdwclass(kk))
                     alphai = 0.5d0 * alpha(vdwclass(ii))
                     alphak = 0.5d0 * alpha(vdwclass(kk))
                  else
                     alphai2 = alpha(cpclass(ii))
                     alphak2 = alpha(cpclass(kk))
                     alphai = 0.5d0 * alpha(cpclass(ii))
                     alphak = 0.5d0 * alpha(cpclass(kk))
                  end if
                  dampi = alphai*r
                  dampk = alphak*r
                  expdampi = exp(-dampi)
                  expdampk = exp(-dampk)
                  diff = alphai2 - alphak2
c                  if (alphai .ne. alphak) then
                  if (abs(diff) .ge. 0.01d0) then
                     term = alphai**2 - alphak**2
                     pre = 64.0d0*(alphai2**3)*(alphak2**3)/(term**4)
c
c     Maple generated code
c
      oscale1ik = r2 ** (-0.3D1 / 0.2D1) * (alphai * (r - 0.40D1 * alpha
     &k / term) * expdampk + alphak * (r + 0.40D1 * alphai / term) * exp
     &dampi) ** 2
c
      oscale3ik = (-0.240D2 * r2 ** (-0.5D1 / 0.2D1) * alphai ** 2 * r *
     & alphak / term + 0.80D1 / r2 ** 2 * alphai ** 2 * alphak / term - 
     &0.2D1 / r2 ** 2 * alphai ** 2 * r + 0.3200D2 / r2 ** 2 * alphai **
     & 2 * alphak ** 3 / term ** 2 + 0.3D1 * r2 ** (-0.5D1 / 0.2D1) * al
     &phai ** 2 * r ** 2 - 0.160D2 / r2 ** 2 * alphak ** 2 * r * alphai 
     &** 2 / term + 0.4800D2 * r2 ** (-0.5D1 / 0.2D1) * alphak ** 2 * al
     &phai ** 2 / term ** 2 + 0.2D1 / r2 ** 2 * alphai ** 2 * r ** 2 * a
     &lphak) * expdampk ** 2 + (0.2D1 / r2 ** 2 * alphai ** 2 * r ** 2 *
     & alphak + 0.2D1 / r2 ** 2 * alphak ** 2 * r ** 2 * alphai + 0.80D1
     & / r2 ** 2 * alphai ** 3 * r * alphak / term - 0.3200D2 / r2 ** 2 
     &* alphak ** 2 * alphai ** 3 / term ** 2 - 0.240D2 * r2 ** (-0.5D1 
     &/ 0.2D1) * alphak ** 2 * r * alphai / term - 0.80D1 / r2 ** 2 * al
     &phai ** 2 * alphak / term - 0.4D1 / r2 ** 2 * alphai * r * alphak 
     &- 0.9600D2 * r2 ** (-0.5D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / 
     &term ** 2 + 0.6D1 * r2 ** (-0.5D1 / 0.2D1) * alphai * r ** 2 * alp
     &hak - 0.3200D2 / r2 ** 2 * alphai ** 2 * alphak ** 3 / term ** 2 -
     & 0.80D1 / r2 ** 2 * alphak ** 3 * r * alphai / term + 0.240D2 * r2
     & ** (-0.5D1 / 0.2D1) * alphai ** 2 * r * alphak / term + 0.80D1 / 
     &r2 ** 2 * alphak ** 2 * alphai / term) * expdampk * expdampi + (0.
     &240D2 * r2 ** (-0.5D1 / 0.2D1) * alphak ** 2 * r * alphai / term +
     & 0.160D2 / r2 ** 2 * alphak ** 2 * r * alphai ** 2 / term + 0.4800
     &D2 * r2 ** (-0.5D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term ** 
     &2 + 0.2D1 / r2 ** 2 * alphak ** 2 * r ** 2 * alphai + 0.3200D2 / r
     &2 ** 2 * alphak ** 2 * alphai ** 3 / term ** 2 + 0.3D1 * r2 ** (-0
     &.5D1 / 0.2D1) * alphak ** 2 * r ** 2 - 0.2D1 / r2 ** 2 * alphak **
     & 2 * r - 0.80D1 / r2 ** 2 * alphak ** 2 * alphai / term) * expdamp
     &i ** 2
c
      oscale5ik = (0.6400D2 * r2 ** (-0.5D1 / 0.2D1) * alphai ** 2 * alp
     &hak ** 4 / term ** 2 + 0.2D1 * r2 ** (-0.5D1 / 0.2D1) * alphai ** 
     &2 + 0.560D2 / r2 ** 3 * alphai ** 2 * alphak / term + 0.22400D3 / 
     &r2 ** 3 * alphai ** 2 * alphak ** 3 / term ** 2 + 0.15D2 * r2 ** (
     &-0.7D1 / 0.2D1) * alphai ** 2 * r ** 2 + 0.4D1 * r2 ** (-0.5D1 / 0
     &.2D1) * alphak ** 2 * alphai ** 2 * r ** 2 - 0.320D2 * r2 ** (-0.5
     &D1 / 0.2D1) * alphai ** 2 * alphak ** 3 * r / term + 0.24000D3 * r
     &2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term ** 2 - 0.
     &1200D3 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * r * alphak / term 
     &+ 0.14D2 / r2 ** 3 * alphai ** 2 * r ** 2 * alphak - 0.1120D3 / r2
     & ** 3 * alphak ** 2 * r * alphai ** 2 / term - 0.14D2 / r2 ** 3 * 
     &alphai ** 2 * r - 0.8D1 * r2 ** (-0.5D1 / 0.2D1) * alphai ** 2 * a
     &lphak * r + 0.320D2 * r2 ** (-0.5D1 / 0.2D1) * alphak ** 2 * alpha
     &i ** 2 / term) * expdampk ** 2 + (0.4D1 * r2 ** (-0.5D1 / 0.2D1) *
     & alphak ** 2 * alphai ** 2 * r ** 2 + 0.14D2 / r2 ** 3 * alphai **
     & 2 * r ** 2 * alphak - 0.8D1 * r2 ** (-0.5D1 / 0.2D1) * alphai ** 
     &2 * alphak * r + 0.14D2 / r2 ** 3 * alphak ** 2 * r ** 2 * alphai 
     &- 0.8D1 * r2 ** (-0.5D1 / 0.2D1) * alphak ** 2 * alphai * r - 0.28
     &D2 / r2 ** 3 * alphai * r * alphak + 0.2D1 * r2 ** (-0.5D1 / 0.2D1
     &) * alphai ** 3 * r ** 2 * alphak + 0.2D1 * r2 ** (-0.5D1 / 0.2D1)
     & * alphak ** 3 * r ** 2 * alphai + 0.30D2 * r2 ** (-0.7D1 / 0.2D1)
     & * alphai * r ** 2 * alphak - 0.22400D3 / r2 ** 3 * alphai ** 2 * 
     &alphak ** 3 / term ** 2 - 0.48000D3 * r2 ** (-0.7D1 / 0.2D1) * alp
     &hak ** 2 * alphai ** 2 / term ** 2 - 0.3200D2 * r2 ** (-0.5D1 / 0.
     &2D1) * alphai ** 2 * alphak ** 4 / term ** 2 - 0.560D2 / r2 ** 3 *
     & alphai ** 2 * alphak / term - 0.6400D2 * r2 ** (-0.5D1 / 0.2D1) *
     & alphai ** 3 * alphak ** 3 / term ** 2 + 0.560D2 / r2 ** 3 * alpha
     &k ** 2 * alphai / term - 0.160D2 * r2 ** (-0.5D1 / 0.2D1) * alphai
     & ** 3 * alphak / term + 0.160D2 * r2 ** (-0.5D1 / 0.2D1) * alphak 
     &** 3 * alphai / term - 0.3200D2 * r2 ** (-0.5D1 / 0.2D1) * alphak 
     &** 2 * alphai ** 4 / term ** 2 - 0.22400D3 / r2 ** 3 * alphak ** 2
     & * alphai ** 3 / term ** 2 - 0.80D1 * r2 ** (-0.5D1 / 0.2D1) * alp
     &hai ** 2 * alphak ** 3 * r / term + 0.1200D3 * r2 ** (-0.7D1 / 0.2
     &D1) * alphai ** 2 * r * alphak / term - 0.80D1 * r2 ** (-0.5D1 / 0
     &.2D1) * alphak ** 4 * r * alphai / term + 0.560D2 / r2 ** 3 * alph
     &ai ** 3 * r * alphak / term - 0.560D2 / r2 ** 3 * alphak ** 3 * r 
     &* alphai / term + 0.80D1 * r2 ** (-0.5D1 / 0.2D1) * alphai ** 4 * 
     &r * alphak / term + 0.80D1 * r2 ** (-0.5D1 / 0.2D1) * alphak ** 2 
     &* alphai ** 3 * r / term - 0.1200D3 * r2 ** (-0.7D1 / 0.2D1) * alp
     &hak ** 2 * r * alphai / term + 0.4D1 * r2 ** (-0.5D1 / 0.2D1) * al
     &phai * alphak) * expdampk * expdampi + (0.22400D3 / r2 ** 3 * alph
     &ak ** 2 * alphai ** 3 / term ** 2 - 0.14D2 / r2 ** 3 * alphak ** 2
     & * r + 0.320D2 * r2 ** (-0.5D1 / 0.2D1) * alphak ** 2 * alphai ** 
     &3 * r / term + 0.15D2 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * r *
     &* 2 + 0.4D1 * r2 ** (-0.5D1 / 0.2D1) * alphak ** 2 * alphai ** 2 *
     & r ** 2 + 0.2D1 * r2 ** (-0.5D1 / 0.2D1) * alphak ** 2 - 0.8D1 * r
     &2 ** (-0.5D1 / 0.2D1) * alphak ** 2 * alphai * r - 0.320D2 * r2 **
     & (-0.5D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term - 0.560D2 / r
     &2 ** 3 * alphak ** 2 * alphai / term + 0.6400D2 * r2 ** (-0.5D1 / 
     &0.2D1) * alphak ** 2 * alphai ** 4 / term ** 2 + 0.1120D3 / r2 ** 
     &3 * alphak ** 2 * r * alphai ** 2 / term + 0.24000D3 * r2 ** (-0.7
     &D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term ** 2 + 0.14D2 / r2 
     &** 3 * alphak ** 2 * r ** 2 * alphai + 0.1200D3 * r2 ** (-0.7D1 / 
     &0.2D1) * alphak ** 2 * r * alphai / term) * expdampi ** 2
c
      oscale7ik = (0.168000D4 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * a
     &lphai ** 2 / term ** 2 - 0.114D3 / r2 ** 4 * alphai ** 2 * r + 0.1
     &82400D4 / r2 ** 4 * alphai ** 2 * alphak ** 3 / term ** 2 + 0.105D
     &3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * r ** 2 + 0.48D2 * r2 **
     & (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 * r ** 2 + 0.12D2 / 
     &r2 ** 3 * alphai ** 2 * alphak + 0.8D1 / r2 ** 3 * alphai ** 2 * a
     &lphak ** 3 * r ** 2 + 0.24D2 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 
     &2 - 0.3840D3 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alphak ** 3 
     &* r / term + 0.114D3 / r2 ** 4 * alphai ** 2 * r ** 2 * alphak + 0
     &.4560D3 / r2 ** 4 * alphai ** 2 * alphak / term - 0.96D2 * r2 ** (
     &-0.7D1 / 0.2D1) * alphai ** 2 * alphak * r + 0.3840D3 * r2 ** (-0.
     &7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term - 0.8400D3 * r2 **
     & (-0.9D1 / 0.2D1) * alphai ** 2 * r * alphak / term - 0.640D2 / r2
     & ** 3 * alphai ** 2 * alphak ** 4 * r / term + 0.76800D3 * r2 ** (
     &-0.7D1 / 0.2D1) * alphai ** 2 * alphak ** 4 / term ** 2 + 0.12800D
     &3 / r2 ** 3 * alphai ** 2 * alphak ** 5 / term ** 2 - 0.9120D3 / r
     &2 ** 4 * alphak ** 2 * r * alphai ** 2 / term - 0.24D2 / r2 ** 3 *
     & alphak ** 2 * alphai ** 2 * r + 0.960D2 / r2 ** 3 * alphai ** 2 *
     & alphak ** 3 / term) * expdampk ** 2 + (-0.182400D4 / r2 ** 4 * al
     &phai ** 2 * alphak ** 3 / term ** 2 + 0.48D2 * r2 ** (-0.7D1 / 0.2
     &D1) * alphak ** 2 * alphai ** 2 * r ** 2 - 0.336000D4 * r2 ** (-0.
     &9D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term ** 2 + 0.12D2 / r2
     & ** 3 * alphai ** 2 * alphak + 0.12D2 / r2 ** 3 * alphak ** 2 * al
     &phai + 0.48D2 * r2 ** (-0.7D1 / 0.2D1) * alphai * alphak - 0.960D2
     & * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alphak ** 3 * r / term +
     & 0.8400D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * r * alphak / te
     &rm - 0.160D2 / r2 ** 3 * alphai ** 2 * alphak ** 4 * r / term + 0.
     &160D2 / r2 ** 3 * alphak ** 2 * alphai ** 4 * r / term + 0.960D2 *
     & r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 3 * r / term - 0
     &.8400D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * r * alphai / term
     & - 0.960D2 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 4 * r * alphai / t
     &erm + 0.4560D3 / r2 ** 4 * alphai ** 3 * r * alphak / term - 0.456
     &0D3 / r2 ** 4 * alphak ** 3 * r * alphai / term - 0.80D1 / r2 ** 3
     & * alphak ** 5 * r * alphai / term + 0.80D1 / r2 ** 3 * alphai ** 
     &5 * r * alphak / term + 0.960D2 * r2 ** (-0.7D1 / 0.2D1) * alphai 
     &** 4 * r * alphak / term + 0.6D1 / r2 ** 3 * alphai ** 2 * alphak 
     &** 3 * r ** 2 + 0.114D3 / r2 ** 4 * alphai ** 2 * r ** 2 * alphak 
     &- 0.4560D3 / r2 ** 4 * alphai ** 2 * alphak / term - 0.96D2 * r2 *
     &* (-0.7D1 / 0.2D1) * alphai ** 2 * alphak * r - 0.38400D3 * r2 ** 
     &(-0.7D1 / 0.2D1) * alphai ** 2 * alphak ** 4 / term ** 2 - 0.3200D
     &2 / r2 ** 3 * alphai ** 2 * alphak ** 5 / term ** 2 - 0.24D2 / r2 
     &** 3 * alphak ** 2 * alphai ** 2 * r + 0.240D2 / r2 ** 3 * alphai 
     &** 2 * alphak ** 3 / term + 0.2D1 / r2 ** 3 * alphai ** 4 * r ** 2
     & * alphak + 0.2D1 / r2 ** 3 * alphak ** 4 * r ** 2 * alphai - 0.12
     &D2 / r2 ** 3 * alphai * alphak ** 3 * r + 0.240D2 / r2 ** 3 * alph
     &ai * alphak ** 4 / term - 0.9600D2 / r2 ** 3 * alphai ** 3 * alpha
     &k ** 4 / term ** 2 - 0.12D2 / r2 ** 3 * alphak * alphai ** 3 * r -
     & 0.240D2 / r2 ** 3 * alphak * alphai ** 4 / term - 0.9600D2 / r2 *
     &* 3 * alphak ** 3 * alphai ** 4 / term ** 2 - 0.228D3 / r2 ** 4 * 
     &alphai * r * alphak + 0.24D2 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 
     &3 * r ** 2 * alphak - 0.1920D3 * r2 ** (-0.7D1 / 0.2D1) * alphai *
     &* 3 * alphak / term + 0.1920D3 * r2 ** (-0.7D1 / 0.2D1) * alphai *
     & alphak ** 3 / term - 0.76800D3 * r2 ** (-0.7D1 / 0.2D1) * alphai 
     &** 3 * alphak ** 3 / term ** 2 + 0.24D2 * r2 ** (-0.7D1 / 0.2D1) *
     & alphak ** 3 * r ** 2 * alphai + 0.210D3 * r2 ** (-0.9D1 / 0.2D1) 
     &* alphai * r ** 2 * alphak - 0.240D2 / r2 ** 3 * alphak ** 2 * alp
     &hai ** 3 / term + 0.6D1 / r2 ** 3 * alphak ** 2 * alphai ** 3 * r 
     &** 2 - 0.3200D2 / r2 ** 3 * alphak ** 2 * alphai ** 5 / term ** 2 
     &+ 0.114D3 / r2 ** 4 * alphak ** 2 * r ** 2 * alphai + 0.4560D3 / r
     &2 ** 4 * alphak ** 2 * alphai / term - 0.38400D3 * r2 ** (-0.7D1 /
     & 0.2D1) * alphak ** 2 * alphai ** 4 / term ** 2 - 0.96D2 * r2 ** (
     &-0.7D1 / 0.2D1) * alphak ** 2 * alphai * r - 0.182400D4 / r2 ** 4 
     &* alphak ** 2 * alphai ** 3 / term ** 2) * expdampk * expdampi + (
     &0.3840D3 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 3 * r 
     &/ term + 0.640D2 / r2 ** 3 * alphak ** 2 * alphai ** 4 * r / term 
     &+ 0.114D3 / r2 ** 4 * alphak ** 2 * r ** 2 * alphai + 0.9120D3 / r
     &2 ** 4 * alphak ** 2 * r * alphai ** 2 / term + 0.182400D4 / r2 **
     & 4 * alphak ** 2 * alphai ** 3 / term ** 2 + 0.8D1 / r2 ** 3 * alp
     &hak ** 2 * alphai ** 3 * r ** 2 + 0.12800D3 / r2 ** 3 * alphak ** 
     &2 * alphai ** 5 / term ** 2 + 0.105D3 * r2 ** (-0.9D1 / 0.2D1) * a
     &lphak ** 2 * r ** 2 - 0.24D2 / r2 ** 3 * alphak ** 2 * alphai ** 2
     & * r + 0.12D2 / r2 ** 3 * alphak ** 2 * alphai - 0.3840D3 * r2 ** 
     &(-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term + 0.168000D4 *
     & r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term ** 2 + 
     &0.76800D3 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 4 / t
     &erm ** 2 - 0.96D2 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai 
     &* r - 0.114D3 / r2 ** 4 * alphak ** 2 * r + 0.8400D3 * r2 ** (-0.9
     &D1 / 0.2D1) * alphak ** 2 * r * alphai / term + 0.24D2 * r2 ** (-0
     &.7D1 / 0.2D1) * alphak ** 2 - 0.4560D3 / r2 ** 4 * alphak ** 2 * a
     &lphai / term - 0.960D2 / r2 ** 3 * alphak ** 2 * alphai ** 3 / ter
     &m + 0.48D2 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 * 
     &r ** 2) * expdampi ** 2
c
      oscale9ik = (0.945D3 * r2 ** (-0.11D2 / 0.2D1) * alphai ** 2 * r *
     &* 2 + 0.48D2 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 
     &- 0.1122D4 / r2 ** 5 * alphai ** 2 * r + 0.216D3 / r2 ** 4 * alpha
     &i ** 2 * alphak + 0.45120D4 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2
     & * alphai ** 2 / term + 0.144D3 / r2 ** 4 * alphai ** 2 * alphak *
     &* 3 * r ** 2 + 0.230400D4 / r2 ** 4 * alphai ** 2 * alphak ** 5 / 
     &term ** 2 + 0.1122D4 / r2 ** 5 * alphai ** 2 * r ** 2 * alphak + 0
     &.44880D4 / r2 ** 5 * alphai ** 2 * alphak / term + 0.1795200D5 / r
     &2 ** 5 * alphai ** 2 * alphak ** 3 / term ** 2 + 0.17280D4 / r2 **
     & 4 * alphai ** 2 * alphak ** 3 / term - 0.64D2 * r2 ** (-0.7D1 / 0
     &.2D1) * alphai ** 2 * alphak ** 3 * r + 0.2560D3 * r2 ** (-0.7D1 /
     & 0.2D1) * alphai ** 2 * alphak ** 4 / term + 0.16D2 * r2 ** (-0.7D
     &1 / 0.2D1) * alphai ** 2 * alphak ** 4 * r ** 2 + 0.25600D3 * r2 *
     &* (-0.7D1 / 0.2D1) * alphai ** 2 * alphak ** 6 / term ** 2 - 0.112
     &8D4 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak * r + 0.902400
     &D4 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak ** 4 / term ** 
     &2 - 0.89760D4 / r2 ** 5 * alphak ** 2 * r * alphai ** 2 / term - 0
     &.1280D3 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alphak ** 5 * r /
     & term - 0.75600D4 * r2 ** (-0.11D2 / 0.2D1) * alphai ** 2 * r * al
     &phak / term - 0.45120D4 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * a
     &lphak ** 3 * r / term - 0.11520D4 / r2 ** 4 * alphai ** 2 * alphak
     & ** 4 * r / term - 0.432D3 / r2 ** 4 * alphak ** 2 * alphai ** 2 *
     & r + 0.564D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 2 
     &* r ** 2 + 0.1512000D5 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 2 * a
     &lphai ** 2 / term ** 2 + 0.282D3 * r2 ** (-0.9D1 / 0.2D1) * alphai
     & ** 2) * expdampk ** 2 + (0.48D2 * r2 ** (-0.7D1 / 0.2D1) * alphak
     & ** 2 * alphai ** 2 + 0.216D3 / r2 ** 4 * alphai ** 2 * alphak + 0
     &.24D2 * r2 ** (-0.7D1 / 0.2D1) * alphai * alphak ** 3 + 0.24D2 * r
     &2 ** (-0.7D1 / 0.2D1) * alphak * alphai ** 3 + 0.564D3 * r2 ** (-0
     &.9D1 / 0.2D1) * alphai * alphak + 0.216D3 / r2 ** 4 * alphak ** 2 
     &* alphai + 0.44880D4 / r2 ** 5 * alphak ** 2 * alphai / term - 0.1
     &795200D5 / r2 ** 5 * alphak ** 2 * alphai ** 3 / term ** 2 + 0.112
     &2D4 / r2 ** 5 * alphak ** 2 * r ** 2 * alphai - 0.48D2 * r2 ** (-0
     &.7D1 / 0.2D1) * alphak ** 2 * alphai ** 3 * r - 0.640D2 * r2 ** (-
     &0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 4 / term + 0.8D1 * r2 ** 
     &(-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 4 * r ** 2 - 0.3200D2 *
     & r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 6 / term ** 2 - 
     &0.451200D4 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 4 / 
     &term ** 2 - 0.1128D4 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alph
     &ai * r + 0.36D2 / r2 ** 4 * alphai ** 4 * r ** 2 * alphak - 0.216D
     &3 / r2 ** 4 * alphak * alphai ** 3 * r - 0.4320D3 / r2 ** 4 * alph
     &ak * alphai ** 4 / term - 0.172800D4 / r2 ** 4 * alphak ** 3 * alp
     &hai ** 4 / term ** 2 - 0.2244D4 / r2 ** 5 * alphai * r * alphak - 
     &0.216D3 / r2 ** 4 * alphai * alphak ** 3 * r + 0.4320D3 / r2 ** 4 
     &* alphai * alphak ** 4 / term - 0.172800D4 / r2 ** 4 * alphai ** 3
     & * alphak ** 4 / term ** 2 + 0.36D2 / r2 ** 4 * alphak ** 4 * r **
     & 2 * alphai + 0.22560D4 * r2 ** (-0.9D1 / 0.2D1) * alphai * alphak
     & ** 3 / term - 0.902400D4 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 3 *
     & alphak ** 3 / term ** 2 - 0.19200D3 * r2 ** (-0.7D1 / 0.2D1) * al
     &phai ** 4 * alphak ** 4 / term ** 2 + 0.2D1 * r2 ** (-0.7D1 / 0.2D
     &1) * alphai ** 5 * r ** 2 * alphak - 0.16D2 * r2 ** (-0.7D1 / 0.2D
     &1) * alphak * alphai ** 4 * r - 0.320D2 * r2 ** (-0.7D1 / 0.2D1) *
     & alphak * alphai ** 5 / term - 0.12800D3 * r2 ** (-0.7D1 / 0.2D1) 
     &* alphak ** 3 * alphai ** 5 / term ** 2 + 0.282D3 * r2 ** (-0.9D1 
     &/ 0.2D1) * alphai ** 3 * r ** 2 * alphak + 0.282D3 * r2 ** (-0.9D1
     & / 0.2D1) * alphak ** 3 * r ** 2 * alphai - 0.16D2 * r2 ** (-0.7D1
     & / 0.2D1) * alphai * alphak ** 4 * r + 0.320D2 * r2 ** (-0.7D1 / 0
     &.2D1) * alphai * alphak ** 5 / term - 0.12800D3 * r2 ** (-0.7D1 / 
     &0.2D1) * alphai ** 3 * alphak ** 5 / term ** 2 + 0.1890D4 * r2 ** 
     &(-0.11D2 / 0.2D1) * alphai * r ** 2 * alphak - 0.22560D4 * r2 ** (
     &-0.9D1 / 0.2D1) * alphai ** 3 * alphak / term + 0.2D1 * r2 ** (-0.
     &7D1 / 0.2D1) * alphak ** 5 * r ** 2 * alphai + 0.12D2 * r2 ** (-0.
     &7D1 / 0.2D1) * alphai ** 3 * alphak ** 3 * r ** 2 + 0.108D3 / r2 *
     &* 4 * alphai ** 2 * alphak ** 3 * r ** 2 - 0.57600D3 / r2 ** 4 * a
     &lphai ** 2 * alphak ** 5 / term ** 2 + 0.1122D4 / r2 ** 5 * alphai
     & ** 2 * r ** 2 * alphak - 0.44880D4 / r2 ** 5 * alphai ** 2 * alph
     &ak / term - 0.1795200D5 / r2 ** 5 * alphai ** 2 * alphak ** 3 / te
     &rm ** 2 + 0.4320D3 / r2 ** 4 * alphai ** 2 * alphak ** 3 / term - 
     &0.48D2 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alphak ** 3 * r + 
     &0.640D2 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alphak ** 4 / ter
     &m + 0.8D1 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alphak ** 4 * r
     & ** 2 - 0.3200D2 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alphak *
     &* 6 / term ** 2 - 0.1128D4 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 
     &* alphak * r - 0.451200D4 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 *
     & alphak ** 4 / term ** 2 - 0.4320D3 / r2 ** 4 * alphak ** 2 * alph
     &ai ** 3 / term + 0.108D3 / r2 ** 4 * alphak ** 2 * alphai ** 3 * r
     & ** 2 - 0.57600D3 / r2 ** 4 * alphak ** 2 * alphai ** 5 / term ** 
     &2 - 0.160D2 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 3 * alphak ** 4 *
     & r / term + 0.160D2 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 3 * alpha
     &i ** 4 * r / term + 0.80D1 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 6 
     &* r * alphak / term + 0.44880D4 / r2 ** 5 * alphai ** 3 * r * alph
     &ak / term - 0.44880D4 / r2 ** 5 * alphak ** 3 * r * alphai / term 
     &- 0.80D1 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 6 * r * alphai / ter
     &m + 0.11280D4 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 4 * r * alphak 
     &/ term - 0.11280D4 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 4 * r * al
     &phai / term - 0.240D2 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alp
     &hak ** 5 * r / term + 0.75600D4 * r2 ** (-0.11D2 / 0.2D1) * alphai
     & ** 2 * r * alphak / term - 0.11280D4 * r2 ** (-0.9D1 / 0.2D1) * a
     &lphai ** 2 * alphak ** 3 * r / term - 0.2880D3 / r2 ** 4 * alphai 
     &** 2 * alphak ** 4 * r / term + 0.11280D4 * r2 ** (-0.9D1 / 0.2D1)
     & * alphak ** 2 * alphai ** 3 * r / term + 0.240D2 * r2 ** (-0.7D1 
     &/ 0.2D1) * alphak ** 2 * alphai ** 5 * r / term - 0.75600D4 * r2 *
     &* (-0.11D2 / 0.2D1) * alphak ** 2 * r * alphai / term + 0.2880D3 /
     & r2 ** 4 * alphak ** 2 * alphai ** 4 * r / term + 0.1440D3 / r2 **
     & 4 * alphai ** 5 * r * alphak / term - 0.1440D3 / r2 ** 4 * alphak
     & ** 5 * r * alphai / term - 0.432D3 / r2 ** 4 * alphak ** 2 * alph
     &ai ** 2 * r + 0.564D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alp
     &hai ** 2 * r ** 2 - 0.3024000D5 * r2 ** (-0.11D2 / 0.2D1) * alphak
     & ** 2 * alphai ** 2 / term ** 2) * expdampk * expdampi + (0.48D2 *
     & r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 + 0.216D3 / r2
     & ** 4 * alphak ** 2 * alphai - 0.1122D4 / r2 ** 5 * alphak ** 2 * 
     &r + 0.945D3 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 2 * r ** 2 - 0.4
     &4880D4 / r2 ** 5 * alphak ** 2 * alphai / term + 0.1795200D5 / r2 
     &** 5 * alphak ** 2 * alphai ** 3 / term ** 2 + 0.1122D4 / r2 ** 5 
     &* alphak ** 2 * r ** 2 * alphai - 0.64D2 * r2 ** (-0.7D1 / 0.2D1) 
     &* alphak ** 2 * alphai ** 3 * r - 0.2560D3 * r2 ** (-0.7D1 / 0.2D1
     &) * alphak ** 2 * alphai ** 4 / term + 0.16D2 * r2 ** (-0.7D1 / 0.
     &2D1) * alphak ** 2 * alphai ** 4 * r ** 2 + 0.25600D3 * r2 ** (-0.
     &7D1 / 0.2D1) * alphak ** 2 * alphai ** 6 / term ** 2 + 0.902400D4 
     &* r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 4 / term ** 2 -
     & 0.1128D4 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai * r - 0.
     &45120D4 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / ter
     &m - 0.17280D4 / r2 ** 4 * alphak ** 2 * alphai ** 3 / term + 0.144
     &D3 / r2 ** 4 * alphak ** 2 * alphai ** 3 * r ** 2 + 0.230400D4 / r
     &2 ** 4 * alphak ** 2 * alphai ** 5 / term ** 2 + 0.89760D4 / r2 **
     & 5 * alphak ** 2 * r * alphai ** 2 / term + 0.45120D4 * r2 ** (-0.
     &9D1 / 0.2D1) * alphak ** 2 * alphai ** 3 * r / term + 0.1280D3 * r
     &2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 5 * r / term + 0.7
     &5600D4 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 2 * r * alphai / term
     & + 0.11520D4 / r2 ** 4 * alphak ** 2 * alphai ** 4 * r / term - 0.
     &432D3 / r2 ** 4 * alphak ** 2 * alphai ** 2 * r + 0.564D3 * r2 ** 
     &(-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 2 * r ** 2 + 0.1512000D
     &5 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 2 * alphai ** 2 / term ** 
     &2 + 0.282D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2) * expdampi ** 
     &2
c
c     end Maple generated code
c
      oscale1ik = pre*oscale1ik
      oscale3ik = pre*oscale3ik
      oscale5ik = pre*oscale5ik
      oscale7ik = pre*oscale7ik
      oscale9ik = pre*oscale9ik
c
                  else
                     pre = (alphai2**6)/(alphai**6)
c     insert alphai .ne. alphak terms here
c
c     start Maple generated code
c
      oscale1ik = (0.1D1 + alphai * r + alphai ** 2 * r2 / 0.3D1) ** 2 *
     & expdampi ** 2 * r2 ** (-0.1D1 / 0.2D1)
c
      oscale3ik = (-0.2D1 / 0.3D1 * r2 ** (-0.1D1 / 0.2D1) * alphai ** 2
     & + 0.2D1 / r2 * alphai ** 2 * r - 0.2D1 / 0.3D1 * r2 ** (-0.1D1 / 
     &0.2D1) * alphai ** 3 * r + 0.2D1 / 0.3D1 * alphai ** 3 - sqrt(r2) 
     &* alphai ** 4 / 0.3D1 + 0.2D1 / r2 * alphai ** 3 * r ** 2 + 0.4D1 
     &/ 0.3D1 * alphai ** 4 * r + 0.2D1 / 0.9D1 * r2 * alphai ** 5 + r2 
     &** (-0.3D1 / 0.2D1) + 0.2D1 * r2 ** (-0.3D1 / 0.2D1) * alphai * r 
     &+ r2 ** (-0.3D1 / 0.2D1) * alphai ** 2 * r ** 2) * expdampi ** 2
c
      oscale5ik = (dble(4 * r2 ** (-0.3D1 / 0.2D1) * alphai ** 4 * r ** 
     &2) + 0.8D1 / 0.3D1 * dble(r2 ** (-0.1D1 / 0.2D1)) * dble(alphai **
     & 5) * dble(r) + dble(6 / r2 ** 2 * alphai ** 3 * r ** 2) + dble(6 
     &* r2 ** (-0.5D1 / 0.2D1) * alphai * r) + dble(3 * r2 ** (-0.5D1 / 
     &0.2D1) * alphai ** 2 * r ** 2) - 0.2D1 / 0.3D1 * dble(r2 ** (-0.3D
     &1 / 0.2D1)) * dble(alphai ** 3) * dble(r) - 0.4D1 / 0.3D1 / dble(r
     &2) * dble(alphai ** 4) * dble(r) + dble(6 / r2 ** 2 * alphai ** 2 
     &* r) + 0.4D1 / 0.9D1 * sqrt(dble(r2)) * dble(alphai ** 6) + dble(3
     & * r2 ** (-0.5D1 / 0.2D1)) - 0.8D1 / 0.3D1 * dble(r2 ** (-0.3D1 / 
     &0.2D1)) * dble(alphai ** 2) - 0.2D1 / 0.3D1 / dble(r2) * dble(alph
     &ai ** 3) + dble(r2 ** (-0.1D1 / 0.2D1) * alphai ** 4) / 0.3D1 - 0.
     &10D2 / 0.9D1 * dble(alphai ** 5)) * expdampi ** 2
c
      oscale7ik = (dble(30 / r2 ** 3 * alphai ** 2 * r) + dble(30 / r2 *
     &* 3 * alphai ** 3 * r ** 2) - dble(12 / r2 ** 2 * alphai ** 4 * r)
     & - dble(2 * r2 ** (-0.5D1 / 0.2D1) * alphai ** 3 * r) + dble(15 * 
     &r2 ** (-0.7D1 / 0.2D1)) + 0.8D1 / 0.9D1 * dble(alphai ** 7) - dble
     &(6 / r2 ** 2 * alphai ** 3) - dble(2 / r2 * alphai ** 5) + dble(r2
     & ** (-0.3D1 / 0.2D1) * alphai ** 4) / 0.3D1 - 0.8D1 / 0.3D1 * dble
     &(r2 ** (-0.1D1 / 0.2D1)) * dble(alphai ** 6) - dble(14 * r2 ** (-0
     &.5D1 / 0.2D1) * alphai ** 2) + dble(8 / r2 ** 2 * alphai ** 5 * r 
     &** 2) + 0.16D2 / 0.3D1 / dble(r2) * dble(alphai ** 6) * dble(r) + 
     &dble(24 * r2 ** (-0.5D1 / 0.2D1) * alphai ** 4 * r ** 2) + dble(30
     & * r2 ** (-0.7D1 / 0.2D1) * alphai * r) + dble(15 * r2 ** (-0.7D1 
     &/ 0.2D1) * alphai ** 2 * r ** 2)) * expdampi ** 2
c
      oscale9ik = (0.16D2 / 0.9D1 * r2 ** (-0.1D1 / 0.2D1) * alphai ** 8
     & + 0.105D3 * r2 ** (-0.9D1 / 0.2D1) - 0.100D3 * r2 ** (-0.7D1 / 0.
     &2D1) * alphai ** 2 - 0.50D2 / r2 ** 3 * alphai ** 3 + r2 ** (-0.5D
     &1 / 0.2D1) * alphai ** 4 - 0.10D2 / 0.3D1 / r2 ** 2 * alphai ** 5 
     &- 0.12D2 * r2 ** (-0.3D1 / 0.2D1) * alphai ** 6 - 0.16D2 / 0.3D1 /
     & r2 * alphai ** 7 + 0.210D3 / r2 ** 4 * alphai ** 3 * r ** 2 + 0.1
     &6D2 * r2 ** (-0.5D1 / 0.2D1) * alphai ** 6 * r ** 2 + 0.32D2 / 0.3
     &D1 * r2 ** (-0.3D1 / 0.2D1) * alphai ** 7 * r + 0.80D2 / r2 ** 3 *
     & alphai ** 5 * r ** 2 - 0.100D3 / r2 ** 3 * alphai ** 4 * r - 0.40
     &D2 * r2 ** (-0.5D1 / 0.2D1) * alphai ** 5 * r + 0.32D2 / 0.3D1 / r
     &2 ** 2 * alphai ** 6 * r - 0.10D2 * r2 ** (-0.7D1 / 0.2D1) * alpha
     &i ** 3 * r + 0.210D3 / r2 ** 4 * alphai ** 2 * r + 0.210D3 * r2 **
     & (-0.9D1 / 0.2D1) * alphai * r + 0.105D3 * r2 ** (-0.9D1 / 0.2D1) 
     &* alphai ** 2 * r ** 2 + 0.180D3 * r2 ** (-0.7D1 / 0.2D1) * alphai
     & ** 4 * r ** 2) * expdampi ** 2
c
c     end Maple generated code
c
c
      oscale1ik = pre*oscale1ik
      oscale3ik = pre*oscale3ik
      oscale5ik = pre*oscale5ik
      oscale7ik = pre*oscale7ik
      oscale9ik = pre*oscale9ik
c
                  end if
c
                  print *,"i,k",i,k,molcule(i),molcule(k)
                  print *,"alphai,alphak",alphai,alphak
                  print *,"oscale1ik",oscale1ik,gl(0)*oscale1ik
                  print *,"oscale3ik",oscale3ik,gl(1)*oscale3ik
                  print *,"oscale5ik",oscale5ik,gl(2)*oscale5ik
                  print *,"oscale7ik",oscale7ik,gl(3)*oscale7ik
                  print *,"oscale9ik",oscale9ik,gl(4)*oscale9ik
c                  print *,"gl",gl
c     
                  es2r = gl(0)*oscale1ik
     &                 + gl(1)*oscale3ik
     &                 + gl(2)*oscale5ik
     &                 + gl(3)*oscale7ik
     &                 + gl(4)*oscale9ik
c
c                  print *,"es2r",es2r
c
               else if (cptype .eq. "ORBITAL2ONLYR") then
c
c     this implements S^2 for exchange repulsion
c     
                  if (use_vdwclass) then
                     alphai2 = alpha(vdwclass(ii))
                     alphak2 = alpha(vdwclass(kk))
                     alphai = 0.5d0 * alpha(vdwclass(ii))
                     alphak = 0.5d0 * alpha(vdwclass(kk))
                  else
                     alphai2 = alpha(cpclass(ii))
                     alphak2 = alpha(cpclass(kk))
                     alphai = 0.5d0 * alpha(cpclass(ii))
                     alphak = 0.5d0 * alpha(cpclass(kk))
                  end if
                  dampi = alphai*r
                  dampk = alphak*r
                  expdampi = exp(-dampi)
                  expdampk = exp(-dampk)
                  diff = alphai2 - alphak2
                  term = alphai**2 - alphak**2
c                  if (alphai .ne. alphak) then
                  if (abs(diff) .ge. 0.001d0) then
c                  if ((term .ne. 0.0d0).and.(term .le. something) then
                     term = alphai**2 - alphak**2
                     pre = 64.0d0*(alphai2**3)*(alphak2**3)/(term**4)
c
c     Maple generated code
c
c
c     chain rule
c
c
      s = (alphai * r - 0.4D1 / term * alphak * alphai) / exp(al
     #phak * r) + (alphak * r + 0.4D1 / term * alphak * alphai) / exp(al
     #phai * r)
      ds = (alphai * alphak * r ** 2 - 0.4D1 * alphai * alphak **
     # 2 / term * r - 0.4D1 / term * alphak * alphai) / exp(alphak * r) 
     #+ (alphai * alphak * r ** 2 + 0.4D1 / term * alphak * alphai + 0.4
     #D1 * alphai ** 2 * alphak / term * r) / exp(alphai * r)
      dds = (alphai * alphak * r ** 2 / 0.3D1 - 0.4D1 * alphai * a
     #lphak ** 2 / term * r + alphai * alphak ** 2 * r ** 3 / 0.3D1 - 0.
     #4D1 / 0.3D1 * alphai * alphak ** 3 / term * r ** 2 - 0.4D1 / term 
     #* alphak * alphai) / exp(alphak * r) + (alphai * alphak * r ** 2 /
     # 0.3D1 + 0.4D1 / term * alphak * alphai + alphai ** 2 * alphak * r
     # ** 3 / 0.3D1 + 0.4D1 * alphai ** 2 * alphak / term * r + 0.4D1 / 
     #0.3D1 * alphai ** 3 * alphak / term * r ** 2) / exp(alphai * r)
      ddds = (-0.8D1 / 0.5D1 * alphai * alphak ** 3 / term * r ** 2
     # + alphai * alphak * r ** 2 / 0.5D1 + alphai * alphak ** 2 * r ** 
     #3 / 0.5D1 + alphai * alphak ** 3 * r ** 4 / 0.15D2 - 0.4D1 / 0.15D
     #2 * alphai * alphak ** 4 / term * r ** 3 - 0.4D1 / term * alphak *
     # alphai - 0.4D1 * alphai * alphak ** 2 / term * r) / exp(alphak * 
     #r) + (0.4D1 / term * alphak * alphai + 0.4D1 * alphai ** 2 * alpha
     #k / term * r + 0.8D1 / 0.5D1 * alphai ** 3 * alphak / term * r ** 
     #2 + alphai * alphak * r ** 2 / 0.5D1 + alphai ** 2 * alphak * r **
     # 3 / 0.5D1 + alphai ** 3 * alphak * r ** 4 / 0.15D2 + 0.4D1 / 0.15
     #D2 * alphai ** 4 * alphak / term * r ** 3) / exp(alphai * r)
      dddds = (-0.12D2 / 0.7D1 * alphai * alphak ** 3 / term * r **
     #2 - 0.4D1 / term * alphak * alphai + alphai * alphak * r ** 2 / 0.
     #7D1 + alphai * alphak ** 2 * r ** 3 / 0.7D1 - 0.8D1 / 0.21D2 * alp
     #hai * alphak ** 4 / term * r ** 3 + 0.2D1 / 0.35D2 * alphai * alph
     #ak ** 3 * r ** 4 - 0.4D1 / 0.105D3 * alphai * alphak ** 5 / term *
     # r ** 4 - 0.4D1 * alphai * alphak ** 2 / term * r + alphai * alpha
     #k ** 4 * r ** 5 / 0.105D3) / exp(alphak * r) + (0.8D1 / 0.21D2 * a
     #lphai ** 4 * alphak / term * r ** 3 + 0.4D1 / term * alphak * alph
     #ai + 0.4D1 / 0.105D3 * alphai ** 5 * alphak / term * r ** 4 + 0.12
     #D2 / 0.7D1 * alphai ** 3 * alphak / term * r ** 2 + alphai ** 4 * 
     #alphak * r ** 5 / 0.105D3 + alphai ** 2 * alphak * r ** 3 / 0.7D1 
     #+ alphai * alphak * r ** 2 / 0.7D1 + 0.4D1 * alphai ** 2 * alphak 
     #/ term * r + 0.2D1 / 0.35D2 * alphai ** 3 * alphak * r ** 4) / exp
     #(alphai * r)
c
      s = s*rr1
      ds = ds*rr3
      dds = dds*rr5
      ddds = ddds*rr7
      dddds = dddds*rr9
c
      oscale1ik = s**2
      oscale3ik = 2.0d0*s*ds
      oscale5ik = 2.0d0*(ds*ds + s*dds)
      oscale7ik = 2.0d0*(dds*ds + ds*dds + ds*dds + s*ddds)
      oscale9ik = 2.0d0*(ddds*ds + dds*dds + dds*dds + ds*ddds +
     &     dds*dds + ds*ddds + ds*ddds + s*dddds)
c
c     end Maple generated code
c
      oscale1ik = pre*oscale1ik
      oscale3ik = pre*oscale3ik
      oscale5ik = pre*oscale5ik
      oscale7ik = pre*oscale7ik
      oscale9ik = pre*oscale9ik
c
                  else
                     pre = (alphai2**6)/(alphai**6)
c     insert alphai .ne. alphak terms here
c
c     start Maple generated code
c
      dampi = alphai*r
      expdampi2 = exp(-2.0d0*dampi)
      dampi2 = dampi**2
      dampi3 = dampi**3
      dampi4 = dampi**4
      oscale1ik = (0.9D1 + 0.18D2 * dampi + 0.15D2 * dampi2 + 0.6D1 * da
     &mpi3 + dampi4) * expdampi2 / 0.9D1
c
      oscale3ik = 0.2D1 / 0.9D1 * alphai ** 2 * (0.3D1 + 0.6D1 * dampi +
     & 0.4D1 * dampi2 + dampi3) * expdampi2
c
      oscale5ik = 0.2D1 / 0.9D1 * alphai ** 4 * (0.4D1 + 0.5D1 * dampi +
     & 0.2D1 * dampi2) * expdampi2
c
      oscale7ik = 0.2D1 / 0.9D1 * alphai ** 5 * (0.4D1 * dampi2 + 0.3D1 
     &+ 0.6D1 * dampi) / r * expdampi2
c
      oscale9ik = 0.2D1 / 0.9D1 * alphai ** 5 * (0.6D1 * dampi + 0.3D1 +
     & 0.8D1 * dampi3 + 0.8D1 * dampi2) / r ** 3 * expdampi2
c
c     end Maple generated code
c
c
      oscale1ik = pre*oscale1ik
      oscale3ik = pre*oscale3ik
      oscale5ik = pre*oscale5ik
      oscale7ik = pre*oscale7ik
      oscale9ik = pre*oscale9ik
c
                  end if
c
c                  print *,"orbital^2onlyr"
c                  print *,"i,k",i,k,molcule(i),molcule(k)
c                  print *,"alphai,alphak",alphai,alphak
c                  print *,"oscale1ik",oscale1ik,gl(0)*oscale1ik
c                  print *,"oscale1ik/r",oscale1ik/r,gl(0)*oscale1ik/r
c                  print *,"oscale3ik",oscale3ik,gl(1)*oscale3ik
c                  print *,"oscale5ik",oscale5ik,gl(2)*oscale5ik
c                  print *,"oscale7ik",oscale7ik,gl(3)*oscale7ik
c                  print *,"oscale9ik",oscale9ik,gl(4)*oscale9ik
c                  print *,"i,k",i,k,molcule(i),molcule(k)
c                  print *,"alphai,alphak",alphai,alphak
c                  print *,"oscale1ik",oscale1ik
c                  print *,"oscale3ik",oscale3ik
c                  print *,"oscale5ik",oscale5ik
c                  print *,"oscale7ik",oscale7ik
c                  print *,"oscale9ik",oscale9ik
c                  print *,"gl",gl
c     
                  es2r = gl(0)*oscale1ik
     &                 + gl(1)*oscale3ik
     &                 + gl(2)*oscale5ik
     &                 + gl(3)*oscale7ik
     &                 + gl(4)*oscale9ik
c
c                  print *,"es2r",es2r
c
               else if (cptype .eq. "ORBITAL2THENR") then
c
c     this implements (del(S))^2/R
c     this is problematic because then we end up with rx^2 instead of rx?
c     to solve this, we must factor out the multipole moments and the 
c     geometry (rx,ry,rz,etc) before squaring.
c     
                  if (use_vdwclass) then
                     alphai2 = alpha(vdwclass(ii))
                     alphak2 = alpha(vdwclass(kk))
                     alphai = 0.5d0 * alpha(vdwclass(ii))
                     alphak = 0.5d0 * alpha(vdwclass(kk))
                  else
                     alphai2 = alpha(cpclass(ii))
                     alphak2 = alpha(cpclass(kk))
                     alphai = 0.5d0 * alpha(cpclass(ii))
                     alphak = 0.5d0 * alpha(cpclass(kk))
                  end if
                  dampi = alphai*r
                  dampk = alphak*r
                  expdampi = exp(-dampi)
                  expdampk = exp(-dampk)
                  if (alphai .ne. alphak) then
                     term = alphai**2 - alphak**2
                     pre = 64.0d0*(alphai2**3)*(alphak2**3)/(term**4)
c
c     Maple generated code
c
      oscale1ik = r2 ** (-0.1D1 / 0.2D1) * (alphai * (r - 0.40D1 * alpha
     &k / term) * expdampk + alphak * (r + 0.40D1 * alphai / term) * exp
     &dampi)
c
      oscale3ik = (r2 ** (-0.3D1 / 0.2D1) * alphai * r - 0.40D1 * r2 ** 
     &(-0.3D1 / 0.2D1) * alphak * alphai / term - 0.40D1 / r2 * alphai *
     & alphak ** 2 / term - 0.1D1 / r2 * alphai + 0.1D1 / r2 * alphak * 
     &alphai * r) * expdampk + (0.1D1 / r2 * alphak * alphai * r + 0.40D
     &1 / r2 * alphak * alphai ** 2 / term + r2 ** (-0.3D1 / 0.2D1) * al
     &phak * r + 0.40D1 * r2 ** (-0.3D1 / 0.2D1) * alphak * alphai / ter
     &m - 0.1D1 / r2 * alphak) * expdampi
c
      oscale5ik = (r2 ** (-0.3D1 / 0.2D1) * alphai * alphak ** 2 * r - 0
     &.120D2 * r2 ** (-0.5D1 / 0.2D1) * alphak * alphai / term - 0.40D1 
     &* r2 ** (-0.3D1 / 0.2D1) * alphai * alphak ** 3 / term - 0.2D1 * a
     &lphak * alphai * r2 ** (-0.3D1 / 0.2D1) + 0.3D1 * r2 ** (-0.5D1 / 
     &0.2D1) * alphai * r - 0.120D2 / r2 ** 2 * alphai * alphak ** 2 / t
     &erm - 0.3D1 / r2 ** 2 * alphai + 0.3D1 / r2 ** 2 * alphak * alphai
     & * r) * expdampk + (0.3D1 * r2 ** (-0.5D1 / 0.2D1) * alphak * r + 
     &0.3D1 / r2 ** 2 * alphak * alphai * r + 0.120D2 / r2 ** 2 * alphak
     & * alphai ** 2 / term + r2 ** (-0.3D1 / 0.2D1) * alphak * alphai *
     &* 2 * r - 0.2D1 * alphak * alphai * r2 ** (-0.3D1 / 0.2D1) + 0.120
     &D2 * r2 ** (-0.5D1 / 0.2D1) * alphak * alphai / term - 0.3D1 / r2 
     &** 2 * alphak + 0.40D1 * r2 ** (-0.3D1 / 0.2D1) * alphak * alphai 
     &** 3 / term) * expdampi
c
      oscale7ik = (0.15D2 * r2 ** (-0.7D1 / 0.2D1) * alphai * r - 0.15D2
     & / r2 ** 3 * alphai - 0.40D1 / r2 ** 2 * alphai * alphak ** 4 / te
     &rm + 0.15D2 / r2 ** 3 * alphak * alphai * r - 0.12D2 * r2 ** (-0.5
     &D1 / 0.2D1) * alphak * alphai + 0.6D1 * r2 ** (-0.5D1 / 0.2D1) * a
     &lphai * alphak ** 2 * r - 0.600D2 / r2 ** 3 * alphai * alphak ** 2
     & / term - 0.3D1 * alphai * alphak ** 2 / r2 ** 2 + 0.1D1 / r2 ** 2
     & * alphai * alphak ** 3 * r - 0.600D2 * r2 ** (-0.7D1 / 0.2D1) * a
     &lphak * alphai / term - 0.240D2 * r2 ** (-0.5D1 / 0.2D1) * alphai 
     &* alphak ** 3 / term) * expdampk + (0.1D1 / r2 ** 2 * alphak * alp
     &hai ** 3 * r + 0.15D2 * r2 ** (-0.7D1 / 0.2D1) * alphak * r + 0.60
     &0D2 / r2 ** 3 * alphak * alphai ** 2 / term - 0.12D2 * r2 ** (-0.5
     &D1 / 0.2D1) * alphak * alphai - 0.15D2 / r2 ** 3 * alphak + 0.240D
     &2 * r2 ** (-0.5D1 / 0.2D1) * alphak * alphai ** 3 / term - 0.3D1 *
     & alphak * alphai ** 2 / r2 ** 2 + 0.40D1 / r2 ** 2 * alphak * alph
     &ai ** 4 / term + 0.15D2 / r2 ** 3 * alphak * alphai * r + 0.600D2 
     &* r2 ** (-0.7D1 / 0.2D1) * alphak * alphai / term + 0.6D1 * r2 ** 
     &(-0.5D1 / 0.2D1) * alphak * alphai ** 2 * r) * expdampi
c
      oscale9ik = (r2 ** (-0.5D1 / 0.2D1) * alphai * alphak ** 4 * r - 0
     &.400D2 / r2 ** 3 * alphai * alphak ** 4 / term + 0.10D2 / r2 ** 3 
     &* alphai * alphak ** 3 * r + 0.105D3 / r2 ** 4 * alphak * alphai *
     & r - 0.4200D3 / r2 ** 4 * alphai * alphak ** 2 / term - 0.40D1 * r
     &2 ** (-0.5D1 / 0.2D1) * alphai * alphak ** 5 / term + 0.45D2 * r2 
     &** (-0.7D1 / 0.2D1) * alphai * alphak ** 2 * r - 0.30D2 / r2 ** 3 
     &* alphai * alphak ** 2 + 0.105D3 * r2 ** (-0.9D1 / 0.2D1) * alphai
     & * r - 0.1800D3 * r2 ** (-0.7D1 / 0.2D1) * alphai * alphak ** 3 / 
     &term - 0.90D2 * r2 ** (-0.7D1 / 0.2D1) * alphak * alphai - 0.105D3
     & / r2 ** 4 * alphai - 0.4200D3 * r2 ** (-0.9D1 / 0.2D1) * alphak *
     & alphai / term - 0.4D1 * alphai * alphak ** 3 * r2 ** (-0.5D1 / 0.
     &2D1)) * expdampk + (-0.4D1 * alphak * alphai ** 3 * r2 ** (-0.5D1 
     &/ 0.2D1) - 0.105D3 / r2 ** 4 * alphak + r2 ** (-0.5D1 / 0.2D1) * a
     &lphak * alphai ** 4 * r + 0.45D2 * r2 ** (-0.7D1 / 0.2D1) * alphak
     & * alphai ** 2 * r + 0.40D1 * r2 ** (-0.5D1 / 0.2D1) * alphak * al
     &phai ** 5 / term + 0.4200D3 / r2 ** 4 * alphak * alphai ** 2 / ter
     &m + 0.105D3 * r2 ** (-0.9D1 / 0.2D1) * alphak * r + 0.10D2 / r2 **
     & 3 * alphak * alphai ** 3 * r + 0.400D2 / r2 ** 3 * alphak * alpha
     &i ** 4 / term + 0.4200D3 * r2 ** (-0.9D1 / 0.2D1) * alphak * alpha
     &i / term + 0.1800D3 * r2 ** (-0.7D1 / 0.2D1) * alphak * alphai ** 
     &3 / term + 0.105D3 / r2 ** 4 * alphak * alphai * r - 0.30D2 / r2 *
     &* 3 * alphak * alphai ** 2 - 0.90D2 * r2 ** (-0.7D1 / 0.2D1) * alp
     &hak * alphai) * expdampi
c
c     end Maple generated code
c
      oscale1ik = pre*(oscale1ik**2)/r
      oscale3ik = pre*(oscale3ik**2)/r
      oscale5ik = pre*(oscale5ik**2)/r
      oscale7ik = pre*(oscale7ik**2)/r
      oscale9ik = pre*(oscale9ik**2)/r
c
                  else
                     pre = (alphai2**6)/(alphai**6)
c     insert alphai .ne. alphak terms here
c
c     start Maple generated code
c
      oscale1ik = (0.1D1 + alphai * r + alphai ** 2 * r2 / 0.3D1) * expd
     &ampi
c
      oscale3ik = (-0.2D1 / 0.3D1 * alphai ** 2 + alphai ** 2 * r2 ** (-
     &0.1D1 / 0.2D1) * r + alphai ** 3 * sqrt(r2) / 0.3D1) * expdampi
c
      oscale5ik = (-alphai ** 2 / r2 - alphai ** 3 * r2 ** (-0.1D1 / 0.2
     &D1) + alphai ** 2 * r2 ** (-0.3D1 / 0.2D1) * r + alphai ** 3 / r2 
     &* r + alphai ** 4 / 0.3D1) * expdampi
c
      oscale7ik = (-dble(3 * alphai ** 2 / r2 ** 2) - dble(3 * alphai **
     & 3 * r2 ** (-0.3D1 / 0.2D1)) - dble(alphai ** 4 / r2) + dble(3 * a
     &lphai ** 2 * r2 ** (-0.5D1 / 0.2D1) * r) + dble(3 * alphai ** 3 / 
     &r2 ** 2 * r) + dble(alphai ** 4 * r2 ** (-0.3D1 / 0.2D1) * r) + db
     &le(alphai ** 5 * r2 ** (-0.1D1 / 0.2D1)) / 0.3D1) * expdampi
c
      oscale9ik = (-dble(6 * alphai ** 4 / r2 ** 2) + dble(alphai ** 5 /
     & r2 ** 2 * r) + dble(15 * alphai ** 2 * r2 ** (-0.7D1 / 0.2D1) * r
     &) - dble(15 * alphai ** 2 / r2 ** 3) + dble(15 * alphai ** 3 / r2 
     &** 3 * r) + dble(6 * alphai ** 4 * r2 ** (-0.5D1 / 0.2D1) * r) + d
     &ble(alphai ** 6 / r2) / 0.3D1 - 0.2D1 / 0.3D1 * dble(alphai ** 5) 
     &* dble(r2 ** (-0.3D1 / 0.2D1)) - dble(15 * alphai ** 3 * r2 ** (-0
     &.5D1 / 0.2D1))) * expdampi
c
c     end Maple generated code
c
c
      oscale1ik = pre*(oscale1ik**2)/r
      oscale3ik = pre*(oscale3ik**2)/r
      oscale5ik = pre*(oscale5ik**2)/r
      oscale7ik = pre*(oscale7ik**2)/r
      oscale9ik = pre*(oscale9ik**2)/r
c
                  end if
c
                  print *,"orbital^2thenR"
                  print *,"i,k",i,k,molcule(i),molcule(k)
                  print *,"alphai,alphak",alphai,alphak
                  print *,"oscale1ik",oscale1ik
                  print *,"oscale3ik",oscale3ik
                  print *,"oscale5ik",oscale5ik
                  print *,"oscale7ik",oscale7ik
                  print *,"oscale9ik",oscale9ik
                  print *,"gl",gl
c     
                  es2r = gl(0)*oscale1ik
     &                 + gl(1)*oscale3ik
     &                 + gl(2)*oscale5ik
     &                 + gl(3)*oscale7ik
     &                 + gl(4)*oscale9ik
c
                  print *,"es2r",es2r
c
               else if (cptype .eq. "ORBITALR") then
c     this implements S/R
c     this method has a problem because the units of S are Coulombs
c     not Coulomb^2
                  if (use_vdwclass) then
                     alphai2 = alpha(vdwclass(ii))
                     alphak2 = alpha(vdwclass(kk))
                     alphai = 0.5d0 * alpha(vdwclass(ii))
                     alphak = 0.5d0 * alpha(vdwclass(kk))
                  else
                     alphai2 = alpha(cpclass(ii))
                     alphak2 = alpha(cpclass(kk))
                     alphai = 0.5d0 * alpha(cpclass(ii))
                     alphak = 0.5d0 * alpha(cpclass(kk))
                  end if
                  dampi = alphai*r
                  dampk = alphak*r
                  expdampi = exp(-dampi)
                  expdampk = exp(-dampk)
                  if (alphai .ne. alphak) then
                     term = alphai**2 - alphak**2
                     pre = 8.0d0*sqrt(alphai2**3)*sqrt(alphak2**3)
     &                    /(term**2)
c
      oscale1ik = 0.1D1 / r2 * (alphai * (r - 0.40D1 * alphak / term) * 
     &expdampk + alphak * (r + 0.40D1 * alphai / term) * expdampi)
c
      oscale3ik = (0.2D1 / r2 ** 2 * alphai * r - 0.80D1 / r2 ** 2 * alp
     &hak * alphai / term - 0.40D1 * r2 ** (-0.3D1 / 0.2D1) * alphai * a
     &lphak ** 2 / term - 0.1D1 * r2 ** (-0.3D1 / 0.2D1) * alphai + r2 *
     &* (-0.3D1 / 0.2D1) * alphak * alphai * r) * expdampk + (r2 ** (-0.
     &3D1 / 0.2D1) * alphak * alphai * r + 0.40D1 * r2 ** (-0.3D1 / 0.2D
     &1) * alphak * alphai ** 2 / term + 0.2D1 / r2 ** 2 * alphak * r + 
     &0.80D1 / r2 ** 2 * alphak * alphai / term - 0.1D1 * r2 ** (-0.3D1 
     &/ 0.2D1) * alphak) * expdampi
c
      oscale5ik = (0.5D1 * r2 ** (-0.5D1 / 0.2D1) * alphak * alphai * r 
     &- 0.320D2 / r2 ** 3 * alphak * alphai / term - 0.40D1 / r2 ** 2 * 
     &alphai * alphak ** 3 / term + 0.8D1 / r2 ** 3 * alphai * r - 0.2D1
     & / r2 ** 2 * alphak * alphai + 0.1D1 / r2 ** 2 * alphai * alphak *
     &* 2 * r - 0.5D1 * r2 ** (-0.5D1 / 0.2D1) * alphai - 0.200D2 * r2 *
     &* (-0.5D1 / 0.2D1) * alphai * alphak ** 2 / term) * expdampk + (-0
     &.2D1 / r2 ** 2 * alphak * alphai + 0.320D2 / r2 ** 3 * alphak * al
     &phai / term + 0.1D1 / r2 ** 2 * alphak * alphai ** 2 * r + 0.40D1 
     &/ r2 ** 2 * alphak * alphai ** 3 / term + 0.5D1 * r2 ** (-0.5D1 / 
     &0.2D1) * alphak * alphai * r + 0.200D2 * r2 ** (-0.5D1 / 0.2D1) * 
     &alphak * alphai ** 2 / term + 0.8D1 / r2 ** 3 * alphak * r - 0.5D1
     & * r2 ** (-0.5D1 / 0.2D1) * alphak) * expdampi
c
      oscale7ik = (-0.40D1 * r2 ** (-0.5D1 / 0.2D1) * alphai * alphak **
     & 4 / term - 0.360D2 / r2 ** 3 * alphai * alphak ** 3 / term - 0.18
     &D2 / r2 ** 3 * alphak * alphai + 0.9D1 / r2 ** 3 * alphai * alphak
     & ** 2 * r + 0.48D2 / r2 ** 4 * alphai * r + 0.33D2 * r2 ** (-0.7D1
     & / 0.2D1) * alphak * alphai * r - 0.3D1 * r2 ** (-0.5D1 / 0.2D1) *
     & alphai * alphak ** 2 - 0.1920D3 / r2 ** 4 * alphak * alphai / ter
     &m + r2 ** (-0.5D1 / 0.2D1) * alphai * alphak ** 3 * r - 0.33D2 * r
     &2 ** (-0.7D1 / 0.2D1) * alphai - 0.1320D3 * r2 ** (-0.7D1 / 0.2D1)
     & * alphai * alphak ** 2 / term) * expdampk + (-0.3D1 * r2 ** (-0.5
     &D1 / 0.2D1) * alphak * alphai ** 2 - 0.18D2 / r2 ** 3 * alphak * a
     &lphai + 0.1920D3 / r2 ** 4 * alphak * alphai / term + 0.48D2 / r2 
     &** 4 * alphak * r + r2 ** (-0.5D1 / 0.2D1) * alphak * alphai ** 3 
     &* r + 0.40D1 * r2 ** (-0.5D1 / 0.2D1) * alphak * alphai ** 4 / ter
     &m + 0.1320D3 * r2 ** (-0.7D1 / 0.2D1) * alphak * alphai ** 2 / ter
     &m + 0.9D1 / r2 ** 3 * alphak * alphai ** 2 * r + 0.360D2 / r2 ** 3
     & * alphak * alphai ** 3 / term + 0.33D2 * r2 ** (-0.7D1 / 0.2D1) *
     & alphak * alphai * r - 0.33D2 * r2 ** (-0.7D1 / 0.2D1) * alphak) *
     & expdampi
c
      oscale9ik = (0.384D3 / r2 ** 5 * alphai * r - 0.42D2 * r2 ** (-0.7
     &D1 / 0.2D1) * alphai * alphak ** 2 + 0.1D1 / r2 ** 3 * alphai * al
     &phak ** 4 * r - 0.560D2 * r2 ** (-0.7D1 / 0.2D1) * alphai * alphak
     & ** 4 / term + 0.14D2 * r2 ** (-0.7D1 / 0.2D1) * alphai * alphak *
     &* 3 * r + 0.279D3 * r2 ** (-0.9D1 / 0.2D1) * alphak * alphai * r -
     & 0.11160D4 * r2 ** (-0.9D1 / 0.2D1) * alphai * alphak ** 2 / term 
     &- 0.279D3 * r2 ** (-0.9D1 / 0.2D1) * alphai - 0.4D1 / r2 ** 3 * al
     &phai * alphak ** 3 + 0.87D2 / r2 ** 4 * alphai * alphak ** 2 * r -
     & 0.3480D3 / r2 ** 4 * alphai * alphak ** 3 / term - 0.174D3 / r2 *
     &* 4 * alphak * alphai - 0.40D1 / r2 ** 3 * alphai * alphak ** 5 / 
     &term - 0.15360D4 / r2 ** 5 * alphak * alphai / term) * expdampk + 
     &(0.3480D3 / r2 ** 4 * alphak * alphai ** 3 / term - 0.279D3 * r2 *
     &* (-0.9D1 / 0.2D1) * alphak + 0.11160D4 * r2 ** (-0.9D1 / 0.2D1) *
     & alphak * alphai ** 2 / term + 0.384D3 / r2 ** 5 * alphak * r - 0.
     &4D1 / r2 ** 3 * alphak * alphai ** 3 + 0.15360D4 / r2 ** 5 * alpha
     &k * alphai / term - 0.42D2 * r2 ** (-0.7D1 / 0.2D1) * alphak * alp
     &hai ** 2 + 0.279D3 * r2 ** (-0.9D1 / 0.2D1) * alphak * alphai * r 
     &- 0.174D3 / r2 ** 4 * alphak * alphai + 0.1D1 / r2 ** 3 * alphak *
     & alphai ** 4 * r + 0.40D1 / r2 ** 3 * alphak * alphai ** 5 / term 
     &+ 0.560D2 * r2 ** (-0.7D1 / 0.2D1) * alphak * alphai ** 4 / term +
     & 0.14D2 * r2 ** (-0.7D1 / 0.2D1) * alphak * alphai ** 3 * r + 0.87
     &D2 / r2 ** 4 * alphak * alphai ** 2 * r) * expdampi
c
      oscale1ik = pre*oscale1ik
      oscale3ik = pre*oscale3ik
      oscale5ik = pre*oscale5ik
      oscale7ik = pre*oscale7ik
      oscale9ik = pre*oscale9ik
c
                  else
                     pre = (alphai2**3)/(alphai**3)
c
      oscale1ik = (0.1D1 + alphai * r + alphai ** 2 * r2 / 0.3D1) * expd
     &ampi * r2 ** (-0.1D1 / 0.2D1)
c
      oscale3ik = (-r2 ** (-0.1D1 / 0.2D1) * alphai ** 2 / 0.3D1 + 0.1D1
     & / r2 * alphai ** 2 * r + alphai ** 3 / 0.3D1 + r2 ** (-0.3D1 / 0.
     &2D1) + r2 ** (-0.3D1 / 0.2D1) * alphai * r) * expdampi
c
      oscale5ik = (-0.4D1 / 0.3D1 * alphai ** 2 * r2 ** (-0.3D1 / 0.2D1)
     & - alphai ** 3 / r2 / 0.3D1 + 0.3D1 * alphai ** 2 / r2 ** 2 * r + 
     &alphai ** 3 * r2 ** (-0.3D1 / 0.2D1) * r + alphai ** 4 * r2 ** (-0
     &.1D1 / 0.2D1) / 0.3D1 + 0.3D1 * r2 ** (-0.5D1 / 0.2D1) + 0.3D1 * r
     &2 ** (-0.5D1 / 0.2D1) * alphai * r) * expdampi
c
      oscale7ik = (dble(15 * alphai ** 2 / r2 ** 3 * r) + dble(6 * alpha
     &i ** 3 * r2 ** (-0.5D1 / 0.2D1) * r) + dble(alphai ** 4 / r2 ** 2 
     &* r) + dble(15 * r2 ** (-0.7D1 / 0.2D1) * alphai * r) - dble(7 * a
     &lphai ** 2 * r2 ** (-0.5D1 / 0.2D1)) - dble(3 * alphai ** 3 / r2 *
     &* 2) + dble(alphai ** 5 / r2) / 0.3D1 + dble(15 * r2 ** (-0.7D1 / 
     &0.2D1))) * expdampi
c
      oscale9ik = (dble(105 * r2 ** (-0.9D1 / 0.2D1) * alphai * r) - dbl
     &e(4 * alphai ** 4 * r2 ** (-0.5D1 / 0.2D1)) + dble(alphai ** 5 * r
     &2 ** (-0.5D1 / 0.2D1) * r) + dble(105 * alphai ** 2 / r2 ** 4 * r)
     & + dble(45 * alphai ** 3 * r2 ** (-0.7D1 / 0.2D1) * r) + dble(10 *
     & alphai ** 4 / r2 ** 3 * r) + dble(alphai ** 6 * r2 ** (-0.3D1 / 0
     &.2D1)) / 0.3D1 - dble(50 * alphai ** 2 * r2 ** (-0.7D1 / 0.2D1)) -
     & dble(25 * alphai ** 3 / r2 ** 3) + 0.2D1 / 0.3D1 * dble(alphai **
     & 5) / dble(r2 ** 2) + dble(105 * r2 ** (-0.9D1 / 0.2D1))) * expdam
     &pi
c
      oscale1ik = pre*oscale1ik
      oscale3ik = pre*oscale3ik
      oscale5ik = pre*oscale5ik
      oscale7ik = pre*oscale7ik
      oscale9ik = pre*oscale9ik
c
                  end if
c
                  print *,"orbital_over_R"
                  print *,"i,k",i,k,molcule(i),molcule(k)
                  print *,"alphai,alphak",alphai,alphak
                  print *,"oscale1ik",oscale1ik
                  print *,"oscale3ik",oscale3ik
                  print *,"oscale5ik",oscale5ik
                  print *,"oscale7ik",oscale7ik
                  print *,"oscale9ik",oscale9ik
                  print *,"gl",gl
c
                  es2r = gl(0)*oscale1ik
     &                 + gl(1)*oscale3ik
     &                 + gl(2)*oscale5ik
     &                 + gl(3)*oscale7ik
     &                 + gl(4)*oscale9ik
c
                  print *,"es2r",es2r
c
c     
               else if (cptype .eq. "ORBITAL2R3") then
c
c     this implements S^2/R for exchange repulsion
c     
                  if (use_vdwclass) then
                     alphai2 = alpha(vdwclass(ii))
                     alphak2 = alpha(vdwclass(kk))
                     alphai = 0.5d0 * alpha(vdwclass(ii))
                     alphak = 0.5d0 * alpha(vdwclass(kk))
                  else
                     alphai2 = alpha(cpclass(ii))
                     alphak2 = alpha(cpclass(kk))
                     alphai = 0.5d0 * alpha(cpclass(ii))
                     alphak = 0.5d0 * alpha(cpclass(kk))
                  end if
                  dampi = alphai*r
                  dampk = alphak*r
                  expdampi = exp(-dampi)
                  expdampk = exp(-dampk)
                  if (alphai .ne. alphak) then
                     term = alphai**2 - alphak**2
                     pre = 64.0d0*(alphai2**3)*(alphak2**3)/(term**4)
c
c     Maple generated code
c
      oscale1ik = r2 ** (-0.5D1 / 0.2D1) * (alphai * (r - 0.40D1 * alpha
     &k / term) * expdampk + alphak * (r + 0.40D1 * alphai / term) * exp
     &dampi) ** 2
c
      oscale3ik = (-0.400D2 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * r *
     & alphak / term + 0.2D1 / r2 ** 3 * alphai ** 2 * r ** 2 * alphak -
     & 0.2D1 / r2 ** 3 * alphai ** 2 * r + 0.3200D2 / r2 ** 3 * alphai *
     &* 2 * alphak ** 3 / term ** 2 - 0.160D2 / r2 ** 3 * alphak ** 2 * 
     &r * alphai ** 2 / term + 0.5D1 * r2 ** (-0.7D1 / 0.2D1) * alphai *
     &* 2 * r ** 2 + 0.80D1 / r2 ** 3 * alphai ** 2 * alphak / term + 0.
     &8000D2 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term
     & ** 2) * expdampk ** 2 + (-0.400D2 * r2 ** (-0.7D1 / 0.2D1) * alph
     &ak ** 2 * r * alphai / term - 0.80D1 / r2 ** 3 * alphai ** 2 * alp
     &hak / term - 0.4D1 / r2 ** 3 * alphai * r * alphak + 0.80D1 / r2 *
     &* 3 * alphai ** 3 * r * alphak / term + 0.400D2 * r2 ** (-0.7D1 / 
     &0.2D1) * alphai ** 2 * r * alphak / term + 0.2D1 / r2 ** 3 * alpha
     &k ** 2 * r ** 2 * alphai + 0.2D1 / r2 ** 3 * alphai ** 2 * r ** 2 
     &* alphak + 0.80D1 / r2 ** 3 * alphak ** 2 * alphai / term - 0.3200
     &D2 / r2 ** 3 * alphai ** 2 * alphak ** 3 / term ** 2 - 0.3200D2 / 
     &r2 ** 3 * alphak ** 2 * alphai ** 3 / term ** 2 + 0.10D2 * r2 ** (
     &-0.7D1 / 0.2D1) * alphai * r ** 2 * alphak - 0.16000D3 * r2 ** (-0
     &.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term ** 2 - 0.80D1 / r
     &2 ** 3 * alphak ** 3 * r * alphai / term) * expdampk * expdampi + 
     &(-0.80D1 / r2 ** 3 * alphak ** 2 * alphai / term + 0.5D1 * r2 ** (
     &-0.7D1 / 0.2D1) * alphak ** 2 * r ** 2 + 0.160D2 / r2 ** 3 * alpha
     &k ** 2 * r * alphai ** 2 / term - 0.2D1 / r2 ** 3 * alphak ** 2 * 
     &r + 0.8000D2 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 
     &/ term ** 2 + 0.400D2 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * r *
     & alphai / term + 0.2D1 / r2 ** 3 * alphak ** 2 * r ** 2 * alphai +
     & 0.3200D2 / r2 ** 3 * alphak ** 2 * alphai ** 3 / term ** 2) * exp
     &dampi ** 2
c
      oscale5ik = (0.6400D2 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alp
     &hak ** 4 / term ** 2 + 0.2D1 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 
     &2 + 0.35D2 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * r ** 2 + 0.560
     &00D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term *
     &* 2 - 0.2800D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * r * alphak
     & / term + 0.22D2 / r2 ** 4 * alphai ** 2 * r ** 2 * alphak + 0.880
     &D2 / r2 ** 4 * alphai ** 2 * alphak / term + 0.35200D3 / r2 ** 4 *
     & alphai ** 2 * alphak ** 3 / term ** 2 - 0.22D2 / r2 ** 4 * alphai
     & ** 2 * r - 0.1760D3 / r2 ** 4 * alphak ** 2 * r * alphai ** 2 / t
     &erm + 0.4D1 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 *
     & r ** 2 - 0.8D1 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alphak * 
     &r + 0.320D2 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 /
     & term - 0.320D2 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alphak **
     & 3 * r / term) * expdampk ** 2 + (-0.112000D4 * r2 ** (-0.9D1 / 0.
     &2D1) * alphak ** 2 * alphai ** 2 / term ** 2 + 0.22D2 / r2 ** 4 * 
     &alphai ** 2 * r ** 2 * alphak - 0.880D2 / r2 ** 4 * alphai ** 2 * 
     &alphak / term - 0.35200D3 / r2 ** 4 * alphai ** 2 * alphak ** 3 / 
     &term ** 2 + 0.4D1 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai 
     &** 2 * r ** 2 - 0.8D1 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * alp
     &hak * r - 0.44D2 / r2 ** 4 * alphai * r * alphak + 0.2D1 * r2 ** (
     &-0.7D1 / 0.2D1) * alphai ** 3 * r ** 2 * alphak + 0.2D1 * r2 ** (-
     &0.7D1 / 0.2D1) * alphak ** 3 * r ** 2 * alphai - 0.160D2 * r2 ** (
     &-0.7D1 / 0.2D1) * alphai ** 3 * alphak / term + 0.160D2 * r2 ** (-
     &0.7D1 / 0.2D1) * alphai * alphak ** 3 / term - 0.6400D2 * r2 ** (-
     &0.7D1 / 0.2D1) * alphai ** 3 * alphak ** 3 / term ** 2 + 0.70D2 * 
     &r2 ** (-0.9D1 / 0.2D1) * alphai * r ** 2 * alphak + 0.22D2 / r2 **
     & 4 * alphak ** 2 * r ** 2 * alphai + 0.880D2 / r2 ** 4 * alphak **
     & 2 * alphai / term - 0.35200D3 / r2 ** 4 * alphak ** 2 * alphai **
     & 3 / term ** 2 - 0.3200D2 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 *
     & alphai ** 4 / term ** 2 - 0.8D1 * r2 ** (-0.7D1 / 0.2D1) * alphak
     & ** 2 * alphai * r - 0.3200D2 * r2 ** (-0.7D1 / 0.2D1) * alphai **
     & 2 * alphak ** 4 / term ** 2 + 0.4D1 * r2 ** (-0.7D1 / 0.2D1) * al
     &phai * alphak - 0.80D1 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 2 * al
     &phak ** 3 * r / term + 0.880D2 / r2 ** 4 * alphai ** 3 * r * alpha
     &k / term - 0.880D2 / r2 ** 4 * alphak ** 3 * r * alphai / term + 0
     &.80D1 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 4 * r * alphak / term -
     & 0.80D1 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 4 * r * alphai / term
     & + 0.80D1 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 3 * r
     & / term - 0.2800D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * r * al
     &phai / term + 0.2800D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * r 
     &* alphak / term) * expdampk * expdampi + (0.6400D2 * r2 ** (-0.7D1
     & / 0.2D1) * alphak ** 2 * alphai ** 4 / term ** 2 - 0.8D1 * r2 ** 
     &(-0.7D1 / 0.2D1) * alphak ** 2 * alphai * r - 0.320D2 * r2 ** (-0.
     &7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 / term + 0.4D1 * r2 ** (-
     &0.7D1 / 0.2D1) * alphak ** 2 * alphai ** 2 * r ** 2 - 0.22D2 / r2 
     &** 4 * alphak ** 2 * r + 0.56000D3 * r2 ** (-0.9D1 / 0.2D1) * alph
     &ak ** 2 * alphai ** 2 / term ** 2 + 0.320D2 * r2 ** (-0.7D1 / 0.2D
     &1) * alphak ** 2 * alphai ** 3 * r / term + 0.35200D3 / r2 ** 4 * 
     &alphak ** 2 * alphai ** 3 / term ** 2 + 0.22D2 / r2 ** 4 * alphak 
     &** 2 * r ** 2 * alphai - 0.880D2 / r2 ** 4 * alphak ** 2 * alphai 
     &/ term + 0.2800D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * r * alp
     &hai / term + 0.1760D3 / r2 ** 4 * alphak ** 2 * r * alphai ** 2 / 
     &term + 0.35D2 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * r ** 2 + 0.
     &2D1 * r2 ** (-0.7D1 / 0.2D1) * alphak ** 2) * expdampi ** 2
c
      oscale7ik = (-0.5760D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * al
     &phak ** 3 * r / term + 0.393600D4 / r2 ** 5 * alphai ** 2 * alphak
     & ** 3 / term ** 2 + 0.36D2 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 
     &+ 0.315D3 * r2 ** (-0.11D2 / 0.2D1) * alphai ** 2 * r ** 2 + 0.12D
     &2 / r2 ** 4 * alphai ** 2 * alphak - 0.246D3 / r2 ** 5 * alphai **
     & 2 * r - 0.24D2 / r2 ** 4 * alphak ** 2 * alphai ** 2 * r + 0.960D
     &2 / r2 ** 4 * alphai ** 2 * alphak ** 3 / term + 0.9840D3 / r2 ** 
     &5 * alphai ** 2 * alphak / term + 0.8D1 / r2 ** 4 * alphai ** 2 * 
     &alphak ** 3 * r ** 2 + 0.12800D3 / r2 ** 4 * alphai ** 2 * alphak 
     &** 5 / term ** 2 - 0.144D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 
     &* alphak * r + 0.5760D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * a
     &lphai ** 2 / term - 0.25200D4 * r2 ** (-0.11D2 / 0.2D1) * alphai *
     &* 2 * r * alphak / term + 0.72D2 * r2 ** (-0.9D1 / 0.2D1) * alphak
     & ** 2 * alphai ** 2 * r ** 2 + 0.115200D4 * r2 ** (-0.9D1 / 0.2D1)
     & * alphai ** 2 * alphak ** 4 / term ** 2 + 0.504000D4 * r2 ** (-0.
     &11D2 / 0.2D1) * alphak ** 2 * alphai ** 2 / term ** 2 - 0.19680D4 
     &/ r2 ** 5 * alphak ** 2 * r * alphai ** 2 / term - 0.640D2 / r2 **
     & 4 * alphai ** 2 * alphak ** 4 * r / term + 0.246D3 / r2 ** 5 * al
     &phai ** 2 * r ** 2 * alphak) * expdampk ** 2 + (-0.393600D4 / r2 *
     &* 5 * alphai ** 2 * alphak ** 3 / term ** 2 + 0.2D1 / r2 ** 4 * al
     &phak ** 4 * r ** 2 * alphai - 0.240D2 / r2 ** 4 * alphak * alphai 
     &** 4 / term - 0.9600D2 / r2 ** 4 * alphak ** 3 * alphai ** 4 / ter
     &m ** 2 - 0.12D2 / r2 ** 4 * alphai * alphak ** 3 * r + 0.240D2 / r
     &2 ** 4 * alphai * alphak ** 4 / term - 0.9600D2 / r2 ** 4 * alphai
     & ** 3 * alphak ** 4 / term ** 2 - 0.12D2 / r2 ** 4 * alphak * alph
     &ai ** 3 * r - 0.492D3 / r2 ** 5 * alphai * r * alphak - 0.24D2 / r
     &2 ** 4 * alphak ** 2 * alphai ** 2 * r + 0.240D2 / r2 ** 4 * alpha
     &i ** 2 * alphak ** 3 / term - 0.9840D3 / r2 ** 5 * alphai ** 2 * a
     &lphak / term + 0.6D1 / r2 ** 4 * alphai ** 2 * alphak ** 3 * r ** 
     &2 - 0.3200D2 / r2 ** 4 * alphai ** 2 * alphak ** 5 / term ** 2 - 0
     &.144D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak * r + 0.72D
     &2 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 2 * r ** 2 - 
     &0.57600D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak ** 4 / t
     &erm ** 2 - 0.1008000D5 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 2 * a
     &lphai ** 2 / term ** 2 + 0.246D3 / r2 ** 5 * alphai ** 2 * r ** 2 
     &* alphak + 0.246D3 / r2 ** 5 * alphak ** 2 * r ** 2 * alphai + 0.9
     &840D3 / r2 ** 5 * alphak ** 2 * alphai / term - 0.393600D4 / r2 **
     & 5 * alphak ** 2 * alphai ** 3 / term ** 2 - 0.240D2 / r2 ** 4 * a
     &lphak ** 2 * alphai ** 3 / term + 0.6D1 / r2 ** 4 * alphak ** 2 * 
     &alphai ** 3 * r ** 2 - 0.3200D2 / r2 ** 4 * alphak ** 2 * alphai *
     &* 5 / term ** 2 - 0.57600D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2
     & * alphai ** 4 / term ** 2 - 0.144D3 * r2 ** (-0.9D1 / 0.2D1) * al
     &phak ** 2 * alphai * r + 0.36D2 * r2 ** (-0.9D1 / 0.2D1) * alphai 
     &** 3 * r ** 2 * alphak + 0.36D2 * r2 ** (-0.9D1 / 0.2D1) * alphak 
     &** 3 * r ** 2 * alphai - 0.2880D3 * r2 ** (-0.9D1 / 0.2D1) * alpha
     &i ** 3 * alphak / term + 0.2880D3 * r2 ** (-0.9D1 / 0.2D1) * alpha
     &i * alphak ** 3 / term - 0.115200D4 * r2 ** (-0.9D1 / 0.2D1) * alp
     &hai ** 3 * alphak ** 3 / term ** 2 + 0.630D3 * r2 ** (-0.11D2 / 0.
     &2D1) * alphai * r ** 2 * alphak + 0.2D1 / r2 ** 4 * alphai ** 4 * 
     &r ** 2 * alphak - 0.1440D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 
     &* alphak ** 3 * r / term + 0.25200D4 * r2 ** (-0.11D2 / 0.2D1) * a
     &lphai ** 2 * r * alphak / term - 0.160D2 / r2 ** 4 * alphai ** 2 *
     & alphak ** 4 * r / term - 0.80D1 / r2 ** 4 * alphak ** 5 * r * alp
     &hai / term + 0.80D1 / r2 ** 4 * alphai ** 5 * r * alphak / term + 
     &0.9840D3 / r2 ** 5 * alphai ** 3 * r * alphak / term - 0.9840D3 / 
     &r2 ** 5 * alphak ** 3 * r * alphai / term + 0.1440D3 * r2 ** (-0.9
     &D1 / 0.2D1) * alphai ** 4 * r * alphak / term - 0.1440D3 * r2 ** (
     &-0.9D1 / 0.2D1) * alphak ** 4 * r * alphai / term - 0.25200D4 * r2
     & ** (-0.11D2 / 0.2D1) * alphak ** 2 * r * alphai / term + 0.1440D3
     & * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 3 * r / term +
     & 0.160D2 / r2 ** 4 * alphak ** 2 * alphai ** 4 * r / term + 0.12D2
     & / r2 ** 4 * alphai ** 2 * alphak + 0.12D2 / r2 ** 4 * alphak ** 2
     & * alphai + 0.72D2 * r2 ** (-0.9D1 / 0.2D1) * alphai * alphak) * e
     &xpdampk * expdampi + (0.115200D4 * r2 ** (-0.9D1 / 0.2D1) * alphak
     & ** 2 * alphai ** 4 / term ** 2 + 0.12D2 / r2 ** 4 * alphak ** 2 *
     & alphai - 0.9840D3 / r2 ** 5 * alphak ** 2 * alphai / term + 0.315
     &D3 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 2 * r ** 2 + 0.640D2 / r2
     & ** 4 * alphak ** 2 * alphai ** 4 * r / term - 0.246D3 / r2 ** 5 *
     & alphak ** 2 * r - 0.5760D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2
     & * alphai ** 2 / term + 0.504000D4 * r2 ** (-0.11D2 / 0.2D1) * alp
     &hak ** 2 * alphai ** 2 / term ** 2 + 0.8D1 / r2 ** 4 * alphak ** 2
     & * alphai ** 3 * r ** 2 + 0.12800D3 / r2 ** 4 * alphak ** 2 * alph
     &ai ** 5 / term ** 2 + 0.19680D4 / r2 ** 5 * alphak ** 2 * r * alph
     &ai ** 2 / term + 0.36D2 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 + 0
     &.25200D4 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 2 * r * alphai / te
     &rm + 0.393600D4 / r2 ** 5 * alphak ** 2 * alphai ** 3 / term ** 2 
     &+ 0.246D3 / r2 ** 5 * alphak ** 2 * r ** 2 * alphai + 0.5760D3 * r
     &2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 3 * r / term - 0.1
     &44D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai * r - 0.24D2 
     &/ r2 ** 4 * alphak ** 2 * alphai ** 2 * r - 0.960D2 / r2 ** 4 * al
     &phak ** 2 * alphai ** 3 / term + 0.72D2 * r2 ** (-0.9D1 / 0.2D1) *
     & alphak ** 2 * alphai ** 2 * r ** 2) * expdampi ** 2
c
      oscale9ik = (0.48D2 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alpha
     &i ** 2 + 0.312D3 / r2 ** 5 * alphai ** 2 * alphak - 0.3090D4 / r2 
     &** 6 * alphai ** 2 * r + 0.3465D4 * r2 ** (-0.13D2 / 0.2D1) * alph
     &ai ** 2 * r ** 2 + 0.91200D4 * r2 ** (-0.11D2 / 0.2D1) * alphak **
     & 2 * alphai ** 2 / term + 0.1140D4 * r2 ** (-0.11D2 / 0.2D1) * alp
     &hak ** 2 * alphai ** 2 * r ** 2 + 0.24960D4 / r2 ** 5 * alphai ** 
     &2 * alphak ** 3 / term + 0.208D3 / r2 ** 5 * alphai ** 2 * alphak 
     &** 3 * r ** 2 + 0.4944000D5 / r2 ** 6 * alphai ** 2 * alphak ** 3 
     &/ term ** 2 + 0.3090D4 / r2 ** 6 * alphai ** 2 * r ** 2 * alphak +
     & 0.123600D5 / r2 ** 6 * alphai ** 2 * alphak / term + 0.332800D4 /
     & r2 ** 5 * alphai ** 2 * alphak ** 5 / term ** 2 + 0.1824000D5 * r
     &2 ** (-0.11D2 / 0.2D1) * alphai ** 2 * alphak ** 4 / term ** 2 - 0
     &.64D2 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak ** 3 * r + 0
     &.2560D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak ** 4 / ter
     &m + 0.16D2 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak ** 4 * 
     &r ** 2 + 0.25600D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak
     & ** 6 / term ** 2 - 0.2280D4 * r2 ** (-0.11D2 / 0.2D1) * alphai **
     & 2 * alphak * r - 0.624D3 / r2 ** 5 * alphak ** 2 * alphai ** 2 * 
     &r + 0.570D3 * r2 ** (-0.11D2 / 0.2D1) * alphai ** 2 + 0.5544000D5 
     &* r2 ** (-0.13D2 / 0.2D1) * alphak ** 2 * alphai ** 2 / term ** 2 
     &- 0.16640D4 / r2 ** 5 * alphai ** 2 * alphak ** 4 * r / term - 0.1
     &280D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak ** 5 * r / t
     &erm - 0.91200D4 * r2 ** (-0.11D2 / 0.2D1) * alphai ** 2 * alphak *
     &* 3 * r / term - 0.277200D5 * r2 ** (-0.13D2 / 0.2D1) * alphai ** 
     &2 * r * alphak / term - 0.247200D5 / r2 ** 6 * alphak ** 2 * r * a
     &lphai ** 2 / term) * expdampk ** 2 + (0.48D2 * r2 ** (-0.9D1 / 0.2
     &D1) * alphak ** 2 * alphai ** 2 + 0.312D3 / r2 ** 5 * alphai ** 2 
     &* alphak + 0.312D3 / r2 ** 5 * alphak ** 2 * alphai + 0.24D2 * r2 
     &** (-0.9D1 / 0.2D1) * alphak ** 3 * alphai + 0.24D2 * r2 ** (-0.9D
     &1 / 0.2D1) * alphai ** 3 * alphak + 0.1140D4 * r2 ** (-0.11D2 / 0.
     &2D1) * alphai * alphak - 0.312D3 / r2 ** 5 * alphai * alphak ** 3 
     &* r - 0.249600D4 / r2 ** 5 * alphak ** 3 * alphai ** 4 / term ** 2
     & - 0.6180D4 / r2 ** 6 * alphai * r * alphak + 0.52D2 / r2 ** 5 * a
     &lphai ** 4 * r ** 2 * alphak + 0.52D2 / r2 ** 5 * alphak ** 4 * r 
     &** 2 * alphai + 0.6240D3 / r2 ** 5 * alphai * alphak ** 4 / term -
     & 0.249600D4 / r2 ** 5 * alphai ** 3 * alphak ** 4 / term ** 2 - 0.
     &312D3 / r2 ** 5 * alphak * alphai ** 3 * r - 0.12800D3 * r2 ** (-0
     &.9D1 / 0.2D1) * alphai ** 5 * alphak ** 3 / term ** 2 - 0.16D2 * r
     &2 ** (-0.9D1 / 0.2D1) * alphak ** 4 * alphai * r + 0.320D2 * r2 **
     & (-0.9D1 / 0.2D1) * alphak ** 5 * alphai / term + 0.45600D4 * r2 *
     &* (-0.11D2 / 0.2D1) * alphai * alphak ** 3 / term - 0.1824000D5 * 
     &r2 ** (-0.11D2 / 0.2D1) * alphai ** 3 * alphak ** 3 / term ** 2 - 
     &0.3200D2 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 6 / te
     &rm ** 2 - 0.48D2 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai *
     &* 3 * r - 0.640D2 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alphai 
     &** 4 / term + 0.8D1 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * alpha
     &i ** 4 * r ** 2 - 0.2280D4 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 2
     & * alphai * r - 0.912000D4 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 2
     & * alphai ** 4 / term ** 2 - 0.6240D3 / r2 ** 5 * alphak * alphai 
     &** 4 / term + 0.3090D4 / r2 ** 6 * alphak ** 2 * r ** 2 * alphai +
     & 0.123600D5 / r2 ** 6 * alphak ** 2 * alphai / term - 0.4944000D5 
     &/ r2 ** 6 * alphak ** 2 * alphai ** 3 / term ** 2 - 0.6240D3 / r2 
     &** 5 * alphak ** 2 * alphai ** 3 / term + 0.156D3 / r2 ** 5 * alph
     &ak ** 2 * alphai ** 3 * r ** 2 - 0.83200D3 / r2 ** 5 * alphak ** 2
     & * alphai ** 5 / term ** 2 + 0.1140D4 * r2 ** (-0.11D2 / 0.2D1) * 
     &alphak ** 2 * alphai ** 2 * r ** 2 + 0.6240D3 / r2 ** 5 * alphai *
     &* 2 * alphak ** 3 / term + 0.156D3 / r2 ** 5 * alphai ** 2 * alpha
     &k ** 3 * r ** 2 - 0.4944000D5 / r2 ** 6 * alphai ** 2 * alphak ** 
     &3 / term ** 2 + 0.3090D4 / r2 ** 6 * alphai ** 2 * r ** 2 * alphak
     & - 0.123600D5 / r2 ** 6 * alphai ** 2 * alphak / term - 0.83200D3 
     &/ r2 ** 5 * alphai ** 2 * alphak ** 5 / term ** 2 - 0.912000D4 * r
     &2 ** (-0.11D2 / 0.2D1) * alphai ** 2 * alphak ** 4 / term ** 2 - 0
     &.48D2 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak ** 3 * r + 0
     &.640D2 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak ** 4 / term
     & + 0.8D1 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak ** 4 * r 
     &** 2 - 0.3200D2 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * alphak **
     & 6 / term ** 2 - 0.2280D4 * r2 ** (-0.11D2 / 0.2D1) * alphai ** 2 
     &* alphak * r - 0.624D3 / r2 ** 5 * alphak ** 2 * alphai ** 2 * r +
     & 0.6930D4 * r2 ** (-0.13D2 / 0.2D1) * alphai * r ** 2 * alphak - 0
     &.45600D4 * r2 ** (-0.11D2 / 0.2D1) * alphai ** 3 * alphak / term +
     & 0.570D3 * r2 ** (-0.11D2 / 0.2D1) * alphai ** 3 * r ** 2 * alphak
     & + 0.2D1 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 5 * r ** 2 * alphai 
     &+ 0.2D1 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 5 * r ** 2 * alphak -
     & 0.19200D3 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 4 * alphak ** 4 / 
     &term ** 2 + 0.12D2 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 3 * alphak
     & ** 3 * r ** 2 - 0.16D2 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 4 * a
     &lphak * r - 0.320D2 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 5 * alpha
     &k / term + 0.570D3 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 3 * r ** 
     &2 * alphai - 0.12800D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 5 * al
     &phai ** 3 / term ** 2 - 0.11088000D6 * r2 ** (-0.13D2 / 0.2D1) * a
     &lphak ** 2 * alphai ** 2 / term ** 2 - 0.4160D3 / r2 ** 5 * alphai
     & ** 2 * alphak ** 4 * r / term - 0.240D2 * r2 ** (-0.9D1 / 0.2D1) 
     &* alphai ** 2 * alphak ** 5 * r / term - 0.22800D4 * r2 ** (-0.11D
     &2 / 0.2D1) * alphai ** 2 * alphak ** 3 * r / term + 0.277200D5 * r
     &2 ** (-0.13D2 / 0.2D1) * alphai ** 2 * r * alphak / term - 0.2080D
     &3 / r2 ** 5 * alphak ** 5 * r * alphai / term - 0.22800D4 * r2 ** 
     &(-0.11D2 / 0.2D1) * alphak ** 4 * r * alphai / term - 0.80D1 * r2 
     &** (-0.9D1 / 0.2D1) * alphak ** 6 * r * alphai / term + 0.80D1 * r
     &2 ** (-0.9D1 / 0.2D1) * alphai ** 6 * r * alphak / term + 0.160D2 
     &* r2 ** (-0.9D1 / 0.2D1) * alphai ** 4 * alphak ** 3 / term * r - 
     &0.160D2 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 4 * alphai ** 3 / ter
     &m * r + 0.22800D4 * r2 ** (-0.11D2 / 0.2D1) * alphai ** 4 * r * al
     &phak / term + 0.123600D5 / r2 ** 6 * alphai ** 3 * r * alphak / te
     &rm + 0.2080D3 / r2 ** 5 * alphai ** 5 * r * alphak / term - 0.1236
     &00D5 / r2 ** 6 * alphak ** 3 * r * alphai / term + 0.240D2 * r2 **
     & (-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 5 * r / term - 0.27720
     &0D5 * r2 ** (-0.13D2 / 0.2D1) * alphak ** 2 * r * alphai / term + 
     &0.22800D4 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 2 * alphai ** 3 * 
     &r / term + 0.4160D3 / r2 ** 5 * alphak ** 2 * alphai ** 4 * r / te
     &rm) * expdampk * expdampi + (0.48D2 * r2 ** (-0.9D1 / 0.2D1) * alp
     &hak ** 2 * alphai ** 2 + 0.312D3 / r2 ** 5 * alphak ** 2 * alphai 
     &- 0.3090D4 / r2 ** 6 * alphak ** 2 * r + 0.3465D4 * r2 ** (-0.13D2
     & / 0.2D1) * alphak ** 2 * r ** 2 + 0.25600D3 * r2 ** (-0.9D1 / 0.2
     &D1) * alphak ** 2 * alphai ** 6 / term ** 2 - 0.64D2 * r2 ** (-0.9
     &D1 / 0.2D1) * alphak ** 2 * alphai ** 3 * r - 0.2560D3 * r2 ** (-0
     &.9D1 / 0.2D1) * alphak ** 2 * alphai ** 4 / term + 0.16D2 * r2 ** 
     &(-0.9D1 / 0.2D1) * alphak ** 2 * alphai ** 4 * r ** 2 - 0.2280D4 *
     & r2 ** (-0.11D2 / 0.2D1) * alphak ** 2 * alphai * r + 0.1824000D5 
     &* r2 ** (-0.11D2 / 0.2D1) * alphak ** 2 * alphai ** 4 / term ** 2 
     &+ 0.3090D4 / r2 ** 6 * alphak ** 2 * r ** 2 * alphai - 0.123600D5 
     &/ r2 ** 6 * alphak ** 2 * alphai / term + 0.4944000D5 / r2 ** 6 * 
     &alphak ** 2 * alphai ** 3 / term ** 2 - 0.24960D4 / r2 ** 5 * alph
     &ak ** 2 * alphai ** 3 / term + 0.208D3 / r2 ** 5 * alphak ** 2 * a
     &lphai ** 3 * r ** 2 + 0.332800D4 / r2 ** 5 * alphak ** 2 * alphai 
     &** 5 / term ** 2 - 0.91200D4 * r2 ** (-0.11D2 / 0.2D1) * alphak **
     & 2 * alphai ** 2 / term + 0.1140D4 * r2 ** (-0.11D2 / 0.2D1) * alp
     &hak ** 2 * alphai ** 2 * r ** 2 - 0.624D3 / r2 ** 5 * alphak ** 2 
     &* alphai ** 2 * r + 0.570D3 * r2 ** (-0.11D2 / 0.2D1) * alphak ** 
     &2 + 0.5544000D5 * r2 ** (-0.13D2 / 0.2D1) * alphak ** 2 * alphai *
     &* 2 / term ** 2 + 0.247200D5 / r2 ** 6 * alphak ** 2 * r * alphai 
     &** 2 / term + 0.1280D3 * r2 ** (-0.9D1 / 0.2D1) * alphak ** 2 * al
     &phai ** 5 * r / term + 0.277200D5 * r2 ** (-0.13D2 / 0.2D1) * alph
     &ak ** 2 * r * alphai / term + 0.91200D4 * r2 ** (-0.11D2 / 0.2D1) 
     &* alphak ** 2 * alphai ** 3 * r / term + 0.16640D4 / r2 ** 5 * alp
     &hak ** 2 * alphai ** 4 * r / term) * expdampi ** 2
c
c     end Maple generated code
c
      oscale1ik = pre*oscale1ik
      oscale3ik = pre*oscale3ik
      oscale5ik = pre*oscale5ik
      oscale7ik = pre*oscale7ik
      oscale9ik = pre*oscale9ik
c
                  else
                     pre = (alphai2**6)/(alphai**6)
c     insert alphai .ne. alphak terms here
c
c     start Maple generated code
c
      oscale1ik = (0.1D1 + alphai * r + alphai ** 2 * r2 / 0.3D1) ** 2 *
     & expdampi ** 2 * r2 ** (-0.3D1 / 0.2D1)
c
      oscale3ik = (0.2D1 / 0.3D1 * r2 ** (-0.3D1 / 0.2D1) * alphai ** 2 
     &+ 0.2D1 / r2 ** 2 * alphai ** 2 * r + 0.2D1 / 0.3D1 * r2 ** (-0.3D
     &1 / 0.2D1) * alphai ** 3 * r + 0.2D1 / 0.3D1 / r2 * alphai ** 3 - 
     &r2 ** (-0.1D1 / 0.2D1) * alphai ** 4 / 0.9D1 + 0.2D1 / r2 ** 2 * a
     &lphai ** 3 * r ** 2 + 0.4D1 / 0.3D1 / r2 * alphai ** 4 * r + 0.2D1
     & / 0.9D1 * alphai ** 5 + 0.3D1 * r2 ** (-0.5D1 / 0.2D1) + 0.6D1 * 
     &r2 ** (-0.5D1 / 0.2D1) * alphai * r + 0.3D1 * r2 ** (-0.5D1 / 0.2D
     &1) * alphai ** 2 * r ** 2) * expdampi ** 2
c
      oscale5ik = (dble(14 / r2 ** 3 * alphai ** 3 * r ** 2) + dble(30 *
     & r2 ** (-0.7D1 / 0.2D1) * alphai * r) + dble(15 * r2 ** (-0.7D1 / 
     &0.2D1) * alphai ** 2 * r ** 2) + dble(2 * r2 ** (-0.5D1 / 0.2D1) *
     & alphai ** 3 * r) + dble(4 / r2 ** 2 * alphai ** 4 * r) + dble(14 
     &/ r2 ** 3 * alphai ** 2 * r) + dble(4 * r2 ** (-0.5D1 / 0.2D1) * a
     &lphai ** 4 * r ** 2) + 0.8D1 / 0.3D1 * dble(r2 ** (-0.3D1 / 0.2D1)
     &) * dble(alphai ** 5) * dble(r) + dble(2 / r2 ** 2 * alphai ** 3) 
     &- dble(r2 ** (-0.3D1 / 0.2D1) * alphai ** 4) / 0.9D1 - 0.2D1 / 0.9
     &D1 / dble(r2) * dble(alphai ** 5) + 0.4D1 / 0.9D1 * dble(r2 ** (-0
     &.1D1 / 0.2D1)) * dble(alphai ** 6) + dble(15 * r2 ** (-0.7D1 / 0.2
     &D1))) * expdampi ** 2
c
      oscale7ik = (dble(114 / r2 ** 4 * alphai ** 2 * r) + dble(114 / r2
     & ** 4 * alphai ** 3 * r ** 2) + dble(210 * r2 ** (-0.9D1 / 0.2D1) 
     &* alphai * r) + dble(105 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 2 * 
     &r ** 2) + dble(6 / r2 ** 3 * alphai ** 3) + dble(12 / r2 ** 3 * al
     &phai ** 4 * r) + dble(8 / r2 ** 3 * alphai ** 5 * r ** 2) + 0.16D2
     & / 0.3D1 / dble(r2 ** 2) * dble(alphai ** 6) * dble(r) - dble(14 *
     & r2 ** (-0.7D1 / 0.2D1) * alphai ** 2) + dble(10 * r2 ** (-0.7D1 /
     & 0.2D1) * alphai ** 3 * r) + dble(48 * r2 ** (-0.7D1 / 0.2D1) * al
     &phai ** 4 * r ** 2) + dble(16 * r2 ** (-0.5D1 / 0.2D1) * alphai **
     & 5 * r) - 0.10D2 / 0.3D1 / dble(r2 ** 2) * dble(alphai ** 5) + 0.8
     &D1 / 0.9D1 / dble(r2) * dble(alphai ** 7) - dble(r2 ** (-0.5D1 / 0
     &.2D1) * alphai ** 4) / 0.3D1 + dble(105 * r2 ** (-0.9D1 / 0.2D1)))
     & * expdampi ** 2
c
      oscale9ik = (dble(945 * r2 ** (-0.11D2 / 0.2D1)) + 0.16D2 / 0.9D1 
     &* dble(r2 ** (-0.3D1 / 0.2D1)) * dble(alphai ** 8) + 0.16D2 / 0.9D
     &1 / dble(r2 ** 2) * dble(alphai ** 7) - dble(212 * r2 ** (-0.9D1 /
     & 0.2D1) * alphai ** 2) - dble(2 / r2 ** 4 * alphai ** 3) - 0.5D1 /
     & 0.3D1 * dble(r2 ** (-0.7D1 / 0.2D1)) * dble(alphai ** 4) - dble(3
     &0 / r2 ** 3 * alphai ** 5) - dble(12 * r2 ** (-0.5D1 / 0.2D1) * al
     &phai ** 6) + dble(1890 * r2 ** (-0.11D2 / 0.2D1) * alphai * r) + d
     &ble(945 * r2 ** (-0.11D2 / 0.2D1) * alphai ** 2 * r ** 2) + dble(1
     &122 / r2 ** 5 * alphai ** 2 * r) + dble(1122 / r2 ** 5 * alphai **
     & 3 * r ** 2) + dble(16 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 6 * r 
     &** 2) + 0.32D2 / 0.3D1 * dble(r2 ** (-0.5D1 / 0.2D1)) * dble(alpha
     &i ** 7) * dble(r) + dble(144 / r2 ** 4 * alphai ** 5 * r ** 2) + 0
     &.160D3 / 0.3D1 / dble(r2 ** 3) * dble(alphai ** 6) * dble(r) + dbl
     &e(70 * r2 ** (-0.9D1 / 0.2D1) * alphai ** 3 * r) - dble(4 / r2 ** 
     &4 * alphai ** 4 * r) + dble(564 * r2 ** (-0.9D1 / 0.2D1) * alphai 
     &** 4 * r ** 2) + dble(88 * r2 ** (-0.7D1 / 0.2D1) * alphai ** 5 * 
     &r)) * expdampi ** 2
c
c     end Maple generated code
c
c
      oscale1ik = pre*oscale1ik
      oscale3ik = pre*oscale3ik
      oscale5ik = pre*oscale5ik
      oscale7ik = pre*oscale7ik
      oscale9ik = pre*oscale9ik
c
                  end if
c
c                  print *,"i,k",i,k,molcule(i),molcule(k)
c                  print *,"alphai,alphak",alphai,alphak
c                  print *,"oscale1ik",oscale1ik
c                  print *,"oscale3ik",oscale3ik
c                  print *,"oscale5ik",oscale5ik
c                  print *,"oscale7ik",oscale7ik
c                  print *,"oscale9ik",oscale9ik
c                  print *,"gl",gl
c     
                  es2r = gl(0)*oscale1ik
     &                 + gl(1)*oscale3ik
     &                 + gl(2)*oscale5ik
     &                 + gl(3)*oscale7ik
     &                 + gl(4)*oscale9ik
c
c                  print *,"es2r",es2r
c
               else if (cptype .eq. "DENSITY") then
c     this branch is dead
                  if (alphai .ne. alphak) then
c                     term = alphai**2 - alphak**2
c                     pre = 2.0d0*qi*qk*(alphai**3)*(alphak**3)/
c     &                    (r*term**3)
c                     chgchg = pre*(alphai*(r*term - 
c     &                    4.0d0*alphak)*expdampk + alphak*(r*term +
c     &                    4.0d0*alphai)*expdampi)
c     check these for the direction of r vector
c                     unormi = sqrt(uix**2 + uiy**2 + uiz**2)
c                     unormk = sqrt(ukx**2 + uky**2 + ukz**2)
c                     pre = -2.0d0*unormi*qk*(alphai**4)*(alphak**3)*
c     &                    alphak/(r2*(term**3))
c                     dipchg = pre*((-4.0d0 - 4.0d0*alphak*r + 
c     &                    term*r2)*expdampk + (4.0d0 + 4.0d0*alphai*r
c     &                    + term*r2)*expdampi)
c                     pre = -2.0d0*unormk*qi*(alphak**4)*(alphai**3)*
c     &                    alphai/(r2*(term**3))
c                     chgdip = pre*((-4.0d0 - 4.0d0*alphai*r +
c     &                    term*r2)*expdampi + (4.0d0 + 4.0d0*alphak*r
c     &                    + term*r2)*expdampk)
c     check check check
c     
c     check these two - get all signs of dot products
c                     pre = 2.0d0*unormi*unormk*
c     &                    (alphai**4)*(alphak**4)/((r**3)*(term**3))
c                     dipdip = pre*((8.0d0*(1.0d0 + dampi) + 
c     &                    4.0d0*dampi**2 + alphai*term*r**3)*expdampi
c     &                    - (8.0d0*(1.0d0 + dampk) + 4.0d0*dampk**2 -
c     &                    alphak*term*r**3)*expdampk)
                  else
c                     
                  end if
c
c     multipole overlap exchange energy
c
               end if
               eover = gl(0)*oscale1ik*rr1
     &              + gl(1)*oscale3ik*rr3
     &              + gl(2)*oscale5ik*rr5
     &              + gl(3)*oscale7ik*rr7
     &              + gl(4)*oscale9ik*rr9
c      
c     undamped energy
c
               ee = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &              + gl(3)*rr7 + gl(4)*rr9
               en = gln(0)*rr1 + gln(1)*rr3 + gln(2)*rr5
               eundamp = ee + en
               eundamp_ele = ee
c     
c     damped energy
c
               ee = gl(0)*scale1ik*rr1 + gl(1)*scale3ik*rr3 
     &              + gl(2)*scale5ik*rr5
     &              + gl(3)*scale7ik*rr7  
     &              + gl(4)*scale9ik*rr9
c     
c     recompute nuclear - damped electron interactions
c
               gln(0) = zi*zk + zi*scale1k*qk + zk*scale1i*qi
               gln(1) = zk*scale3i*sc(3) - zi*scale3k*sc(4)
               gln(2) = zi*scale5k*sc(6) + zk*scale5i*sc(5)
c     
c     total damped electrostatic energy
c     
               en = gln(0)*rr1 + gln(1)*rr3 + gln(2)*rr5
               edamp = ee + en
               edamp_ele = ee
c     
c     total charge penetration energy
c
               if (cptype .eq. "ELE") then
                  e = edamp_ele - eundamp_ele
                  if (exmodel .eq. "COMBO5") then
                     e = (-oik*f*e) + (-boik*r*f*e)
                  else
                     e = -oik*f*e
                  end if
               else if ((cptype .eq. "EX-DAMPING").or.
     &                 (cptype .eq. "EX-DAMPING1")) then
c                  print *,qi,qk
c                  print *,gl(0)*oscale1ik*rr1
c                  print *,gl(1)*oscale3ik*rr3
c                  print *,gl(2)*oscale5ik*rr5
c                  print *,gl(3)*oscale7ik*rr7
c                  print *,gl(4)*oscale9ik*rr9
                  print *,"1",oscale1ik
                  print *,"3",oscale3ik
                  print *,"5",oscale5ik
                  print *,"7",oscale7ik
                  print *,"9",oscale9ik
                  e = eover
                  if (exmodel .eq. "COMBO4") then
                     e2 = edamp_ele - eundamp_ele
                     e = -boik*f*e2 + oik*f*eover
                  end if
               else if (cptype .eq. "EX-DAMPING-EXP") then
                  eover = gl(0)*oscale1ik*rr1 
     &              + gl(1)*oscale3ik*rr3
     &              + gl(2)*oscale5ik*rr5
     &              + gl(3)*oscale7ik*rr7
     &              + gl(4)*oscale9ik*rr9
                  eoverexp = gl(0)*expscale1ik*rr1
     &              + gl(1)*expscale3ik*rr3
     &              + gl(2)*expscale5ik*rr5
     &              + gl(3)*expscale7ik*rr7
     &              + gl(4)*expscale9ik*rr9
                  e = oik*(eover - boik*eoverexp)
               else if (cptype .eq. "EX-ORBITAL1") then
c                  print *,"i,k,r",i,k,r
c                  print *,"1",oscale1ik
c                  print *,"3",oscale3ik
c                  print *,"5",oscale5ik
c                  print *,"7",oscale7ik
c                  print *,"9",oscale9ik
c     this is (sort of) S^2/R
c                  eover = gl(0)*(oscale1ik*rr1)**2
c     &                 + gl(1)*(oscale3ik*rr3)**2
c     &                 + gl(2)*(oscale5ik*rr5)**2
c     &                 + gl(3)*(oscale7ik*rr7)**2
c     &                 + gl(4)*(oscale9ik*rr9)**2
c
c                  e = oik*(eover**2)/r
c     this is S then /R
                  e = oik*eover/r
               else if (cptype .eq. "EX-ORBITAL1-ABS") then
c     what units is eover?
                  print *,"eover",eover
                  e = oik*(eover**2)/r
c     need to add in ex-orbital1-sqrt
               else if (cptype .eq. "EX-DAMPING-R2") then
                  e = oik*eover*r2
               else if (cptype .eq. "SIMPLE") then
                  e = oik*f*eover
               else if (cptype .eq. "ORBITAL2R") then
c                  e = oik*es2r
                  e = oik*zi*zk*es2r
                  if (molcule(i) .ne. molcule(k)) then 
                     print *,"S^2/R",i,k,es2r
                  end if
               else if (cptype .eq. "ORBITAL2R3") then
                  e = oik*es2r
c                  print *,"e",e
               else if (cptype .eq. "ORBITAL2ONLYR") then
c                  print *,"dividing by r"
c                  e = oik*es2r
                  e = oik*zi*zk*es2r/r
c                  e = oik*atomic(ii)*atomic(kk)*es2r/r
c     zi and zk are the "core" charge, not atomic number
c                  if (molcule(i) .ne. molcule(k)) then
c                     print *,"(S^2_tot)/R",i,k,e
c                  end if
               else if (cptype .eq. "ORBITAL2THENR") then
                  e = oik*es2r
                  print *,"e",e
               else if (cptype .eq. "ORBITALR") then
                  e = oik*es2r
                  print *,"e",e
               else
                  e = edamp - eundamp
                  e = -oik * f * e
               end if
ccccccccccccccccc TOTAL ELECTROSTATIC ENERGY WITH DAMPING cccccccc
c                  e = edamp
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     apply the energy adjustments for scaled interactions
c     
               if (exmodel .eq. "COMBO") then
                  if (alphai.ne.alphak) then
                     term = alphai**2 - alphak**2
                     if (gtype .eq. "NEW") then
                        e2 = pre*(8.0d0/(r*term**2))*
     &                       (alphai*(r - 4.0d0*alphak/term)*expdampk +
     &                       alphak*(r + 4.0d0*alphai/term)*expdampi)
                     else
                        e2 = (pre * 4.0d0 * pi / (r*term))*
     &                       (expdampk - expdampi)
                     end if
                  else
                     if (gtype .eq. "NEW") then
                        e2 = pre*(1.0d0/(alphai**3))*(1.0d0 + dampi +
     &                       (1.0d0/3.0d0)*(dampi)**2)*expdampi
                     else
                        e2 = (pre * 2.0d0 * pi / alphai)
     &                       *expdampi
                     end if
                  end if
                  e = e + boik*e2
               else if (exmodel .eq. "COMBO2") then
                  roik = (overlapri*overlaprk*(overlapri + overlaprk))/
     &                 (overlapri**2 + overlaprk**2)
c                  roik = (overlapri + overlaprk) / 2.0d0
                  e2 = -boik*exp(-roik*r)
                  e = e + e2
               else if (exmodel .eq. "COMBO3") then
                  e2 = r*e/oik
                  e = e + boik*e2
               end if               
            end if
c
            if (molcule(i).ne.molcule(k)) then
c               print *,"energy",i,k,r,e
               ex = ex + e
c               nex = nex + 1
c               aex(i) = aex(i) + e
c               aex(k) = aex(k) + e
c               einter = einter + e
               if (debug) then
                  print *,atomic(ii),atomic(kk),r,e
                  if (e .lt. 0.0d0) then
                     print *,"Exchange less than zero!!!",e
                  end if
               end if
            end if
         end do
      end do
c      deallocate (vscale)
      return
      end
