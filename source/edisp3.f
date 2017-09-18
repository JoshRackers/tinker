c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine edisp3  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "edisp3" calculates any additional user defined potential
c     contribution and also partitions the energy among the atoms
c
c
      subroutine edisp3
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
      use inter
c      use vdwpot
      implicit none
      integer i,k,j
      integer ii,kk
      integer rorder
      real*8 r,r2,r6
      real*8 xr,yr,zr
      real*8 ci,qi,zi
      real*8 ck,qk,zk
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 damping
      real*8 csixi,csixk
      real*8 cik
      real*8 term
      real*8 termi,termk
      real*8 damp3,damp5
      real*8 e
      real*8 iv,rdn
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
c      real*8, allocatable :: vscale(:)
c
c
c     zero out energy and partitioning due to extra potential terms
c
c      nedis = 0
      edis = 0.0d0
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
c      if (.not.use_mpole) then
c         call kmpole
c         call rotpole
c      end if
cccccccccccccccccccccc
      if (.not.allocated(xred)) allocate (xred(n))
      if (.not.allocated(yred)) allocate (yred(n))
      if (.not.allocated(zred)) allocate (zred(n))
      do k = 1, npole
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
ccccccccccccccccccccccc
      do i = 1, npole-1
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
         do k = i+1, npole
            kk = ivdw(k)
            ck = rpole(1,k)
c            xr = x(i) - x(k)
c            yr = y(i) - y(k)
c            zr = z(i) - z(k)
c
            xr = xred(k) - xred(i)
            yr = yred(k) - yred(i)
            zr = zred(k) - zred(i)
c
            r2 = xr*xr + yr*yr + zr*zr
            r = sqrt(r2)
            r6 = r2**3
            if (use_vdwclass) then
               csixi = csix(vdwclass(ii))
               csixk = csix(vdwclass(kk))
               alphai = alpha(vdwclass(ii))
               alphak = alpha(vdwclass(kk))
            else
               csixi = csix(cpclass(ii))
               csixk = csix(cpclass(kk))
               alphai = alpha(cpclass(ii))
               alphak = alpha(cpclass(kk))
            end if
            cik = csixi*csixk
c
            dampi = alphai*r
            dampk = alphak*r
            expdampi = exp(-dampi)
            expdampk = exp(-dampk)            
            if (dispersiondamping .eq. "YES") then
               if (alphai.ne.alphak) then
                  termi = alphak**2/(alphak**2 - alphai**2)
                  termk = alphai**2/(alphai**2 - alphak**2)
c     this damping term can't be right...
c     what SHOULD it be?
c     this should be the rr3 overlap damping term (squared)
                  damping = 1.0d0 - termi*(1.0d0 +dampi)*expdampi
     &                 - termk*(1.0d0 + dampk)*expdampk
c                  print *,i,k,damping
c                  damping = 1.0d0 - (4.0d0 * pi / (r*term))*
c       &                 (expdampk - expdampi)
                  e = -(cik * damping**2) / r6
               else
                  damping = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2)
     &                 *expdampi
c                  print *,i,k,damping
c                  damping = 1.0d0 - (2.0d0 * pi / alphai)
c     &                 *expdampi
                  e = -(cik * damping**2) / r6
               end if
            else if (dispersiondamping .eq. "FANCY") then
               if (alphai.ne.alphak) then
                  termi = alphak**2/(alphak**2 - alphai**2)
                  termk = alphai**2/(alphai**2 - alphak**2)
                  damp3 = 1.0d0 - termi*(1.0d0 +dampi)*expdampi
     &                 - termk*(1.0d0 + dampk)*expdampk
                  damp5 = 1.0d0 - termi*(1.0d0 + dampi +
     &                 (1.0d0/3.0d0)*dampi**2)*expdampi -
     &                 termk*(1.0d0 + dampk +
     &                 (1.0d0/3.0d0)*dampk**2)*expdampk
                  damping = (3.0d0*damp5 - damp3)/2.0d0
c                  print *,"damping",i,k,r,damp3,damp5,damping
                  e = -(cik * damping**2) / r6
               else
                  damp3 = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2)
     &                 *expdampi
                  damp5 = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3)*expdampi
                  damping = (3.0d0*damp5 - damp3)/2.0d0
c                  print*,"damping",i,k,r,damp3,damp5,damping
                  e = -(cik * damping**2) / r6
               end if
            else if (dispersiondamping .eq. "FANCY1") then
               if (alphai .ne. alphak) then
                  termi = alphak**2/(alphak**2 - alphai**2)
                  termk = alphai**2/(alphai**2 - alphak**2)
                  damp3 = 1.0d0 - (termi**2)*(1.0d0 + dampi + 
     &                 0.5d0*dampi**2)*expdampi - 
     &             (termk**2)*(1.0d0 + dampk + 0.5d0*dampk**2)*expdampk
     &                 - 2.0d0*(termi**2)*termk*(1.0d0 + dampi)*expdampi
     &                 - 2.0d0*(termk**2)*termi*(1.0d0 + dampk)*expdampk
                  damp5 = 1.0d0 - (termi**2)*
     &                 (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &                 (1.0d0/6.0d0)*dampi**3)*expdampi - (termk**2)*
     &                 (1.0d0 + dampk + 0.5d0*dampk**2 + 
     &                 (1.0d0/6.0d0)*dampk**3)*expdampk - 
     &                 2.0d0*(termi**2)*termk*(1.0 + dampi + 
     &                 (1.0d0/3.0d0)*dampi**2)*expdampi - 
     &                 2.0d0*(termk**2)*termi*(1.0 + dampk +
     &                 (1.0d0/3.0d0)*dampk**2)*expdampk
                  damping = (3.0d0*damp5 - damp3)/2.0d0
                  e = -(cik * damping**2) / r6
               else
                  damp3 = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &                (7.0d0/48.0d0)*dampi**3 + (1.0d0/48.0d0)*dampi**4)
     &                 *expdampi
                  damp5 = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &                (1.0d0/6.0d0)*dampi**3 + (1.0d0/24.0d0)*dampi**4 +
     &                 (1.0d0/144.0d0)*dampi**5)*expdampi
                  damping = (3.0d0*damp5 - damp3)/2.0d0
                  e = -(cik * damping**2) / r6
               end if
            else
               e = -cik / r6
            end if
            if (molcule(i).ne.molcule(k)) then
               edis = edis + e
c               nedis = nedis + 1
               einter = einter + e
            end if
         end do
      end do
c      deallocate (vscale)
      return
      end
