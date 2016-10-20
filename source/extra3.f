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
      use couple
      use chgpen
      use vdw
      use vdwpot
      implicit none
      integer i,k,j
      integer ii,kk
      integer rorder
      real*8 r,r2
      real*8 xr,yr,zr
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 overlapi,overlapk
      real*8 oik,term
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
      allocate (vscale(n))
      do i = 1, n
         vscale(i) = 1.0d0
      end do
c
c     add any user-defined extra potentials and partitioning
c
      do i = 1, nvdw-1
         ii = ivdw(i)
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
c            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
         do k = i+1, nvdw
            kk = ivdw(k)
            xr = x(i) - x(k)
            yr = y(i) - y(k)
            zr = z(i) - z(k)
            r2 = xr*xr + yr*yr + zr*zr
            r = sqrt(r2)
            overlapi = overlap(cpclass(ii))
            overlapk = overlap(cpclass(kk))
            oik = overlapi * overlapk
c     put other mixing rules here
            alphai = alpha(cpclass(ii))
            alphak = alpha(cpclass(kk))
            dampi = alphai*r
            dampk = alphak*r
            expdampi = exp(-dampi)
            expdampk = exp(-dampk)
            if (alphai.ne.alphak) then
               term = alphai**2 - alphak**2
               e = (oik * 8.0d0 / (r * term**3))*
     &              ((alphai*(r*term - 4.0d0*alphak))*expdampk +
     &              (alphak*(r*term + 4.0d0*alphai))*expdampi)
            else
               e = (oik / alphai**3)*(1.0d0 + dampi +
     &              (1.0d0/3.0d0)*dampi**2)*expdampi
            end if
            e = e*vscale(kk)
            ex = ex + e
            nex = nex + 1
         end do
      end do
      return
      end
