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
      use math
      use molcul
      use mpole
      use couple
      use chgpen
      use vdw
      use potent
c      use vdwpot
      implicit none
      integer i,k,j
      integer ii,kk
      integer rorder
      real*8 r,r2
      real*8 xr,yr,zr
      real*8 ci,qi,zi
      real*8 ck,qk,zk
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 overlapi,overlapk
      real*8 oik,oik2,term,pre
      real*8 e
      real*8, allocatable :: vscale(:)
c
c
c     zero out energy and partitioning due to extra potential terms
c
      nex = 0
      ex = 0.0d0
      do i = 1, n
         aex(i) = 0.0d0
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
      if (.not.use_mpole) then
         call kmpole
         call rotpole
      end if
      do i = 1, nvdw-1
         ii = ivdw(i)
         ci = rpole(1,i)
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
            xr = x(i) - x(k)
            yr = y(i) - y(k)
            zr = z(i) - z(k)
            r2 = xr*xr + yr*yr + zr*zr
            r = sqrt(r2)
            overlapi = overlap(vdwclass(ii))
            overlapk = overlap(vdwclass(kk))
            if (singoverlap) then
               overlapi = soverlap
               overlapk = soverlap
            end if
            oik = overlapi * overlapk
            if (exmix .eq. "W-H") then
               oik2 = overlapi * overlapk
               oik = sqrt(oik2)
            end if
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
            alphai = alpha(cpclass(ii))
            alphak = alpha(cpclass(kk))
            dampi = alphai*r
            dampk = alphak*r
            expdampi = exp(-dampi)
            expdampk = exp(-dampk)
            pre = qi*qk*(alphai**2)*(alphak**2)
            if (alphai.ne.alphak) then
               term = alphai**2 - alphak**2
               e = (pre * oik * 4.0d0 * pi / (r*term))*
     &              (expdampk - expdampi)
c               e = (oik * 8.0d0 / (r * term**3))*
c     &              ((alphai*(r*term - 4.0d0*alphak))*expdampk +
c     &              (alphak*(r*term + 4.0d0*alphai))*expdampi)
            else
c               e = (oik / alphai**3)*(1.0d0 + dampi +
c     &              (1.0d0/3.0d0)*dampi**2)*expdampi
               e = (pre * oik * 2.0d0 * pi / alphai)*expdampi
            end if
c            e = e*vscale(kk)
            if (molcule(i).ne.molcule(k)) then
               ex = ex + e
               nex = nex + 1
               aex(i) = aex(i) + e
               aex(k) = aex(k) + e
            end if
         end do
      end do
c      deallocate (vscale)
      return
      end
