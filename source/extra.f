c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine extra  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "extra" calculates any additional user defined potential
c     contribution and also partitions the energy among the atoms
c
c
      subroutine extra
      use limits
      implicit none
      integer i
c
c     choose method for summing over charge transfer interactions
c
      if (use_mlist) then
         call extra0b
      else
         call extra0a
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine extra0a  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "extra0a" calculates any additional user defined potential
c     contribution and also partitions the energy among the atoms
c     double loop
c
c
      subroutine extra0a
      use sizes
      use atoms
      use atoms
      use bound
      use cell
      use couple
      use energi
      use group
      use shunt
      use usage
      use xtrapot
      use mpole
      use mplpot
      implicit none
      integer i,j,k
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,fgrp,r,r2
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 alphai,alphak
      real*8 transferi,transferk
      real*8, allocatable :: mscale(:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     zero out energy and partitioning due to extra potential terms
c
c      nex = 0
      ex = 0.0d0
c      do i = 1, n
c         aex(i) = 0.0d0
c      end do
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
c     set cutoff and switching coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     calculate the charge transfer energy term
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         transferi = chgct(i)
         alphai = alphact(i)
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
                  transferk = chgct(k)
                  alphak = alphact(k)
c
c     compute charge transfer energy for the pair
c
                  e = -transferi * exp(-alphak*r) - 
     &                 transferk * exp(-alphai*r) 
                  ex = ex + e*mscale(kk)
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
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine extra0b  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "extra0b" calculates any additional user defined potential
c     contribution and also partitions the energy among the atoms
c     neighbor list
c
c
      subroutine extra0b
      use sizes
      use atoms
      use atoms
      use bound
      use cell
      use couple
      use energi
      use group
      use neigh
      use shunt
      use usage
      use xtrapot
      use mpole
      use mplpot
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,fgrp,r,r2
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 alphai,alphak
      real*8 transferi,transferk
      real*8, allocatable :: mscale(:)
      logical proceed,usei,usek
      character*6 mode
c
c     zero out energy
c
      ex = 0.0d0
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
c     set cutoff and switching coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,chgct,alphact,n12,i12,n13,i13,
!$OMP& n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,
!$OMP& nelst,elst,use_bounds,off2)
!$OMP& firstprivate(mscale) shared (ex)
!$OMP DO reduction(+:ex) schedule(guided)
c
c     compute the charge transfer energy
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         transferi = chgct(i)
         alphai = alphact(i)
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
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               transferk = chgct(k)
               alphak = alphact(k)
c     
c     compute charge transfer energy for the pair
c
               e = -transferi * exp(-alphak*r) - 
     &              transferk * exp(-alphai*r) 
               ex = ex + e*mscale(kk)
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
