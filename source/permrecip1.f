c
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine permrecip1  --  Ewald recip permanent induced field  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "permrecip1" computes the reciprocal space contribution of the
c     permanent atomic multipole moments to the field
c
c
      subroutine permrecip1
      use sizes
      use bound
      use boxes
      use ewald
      use math
      use mpole
      use pme
      use potderivs
      implicit none
      integer i,j,k,l,ntot
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff,nf1,nf2,nf3
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: fphi(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(pot_recip)) allocate (pot_recip(npole))
      if (.not. allocated(field_recip)) allocate (field_recip(3,npole))
      if (.not. allocated(gradfield_recip)) 
     &     allocate (gradfield_recip(3,3,npole))
      if (.not. allocated(hessfield_recip)) 
     &     allocate (hessfield_recip(3,3,3,npole))
c
c     zero out the value of the potential, field, etc. at each site
c
      do i = 1, npole
         pot_recip(i) = 0.0d0
         do j = 1, 3
            field_recip(j,i) = 0.0d0
            do k = 1, 3
               gradfield_recip(k,j,i) = 0.0d0
               do l = 1, 3
                  hessfield_recip(l,k,j,i) = 0.0d0
               end do
            end do
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
      allocate (cphi(20,npole))
      allocate (fphi(20,npole))
c
c     copy multipole moments and coordinates to local storage
c
      do i = 1, npole
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
         cmp(5,i) = rpole(5,i)
         cmp(6,i) = rpole(9,i)
         cmp(7,i) = rpole(13,i)
         cmp(8,i) = 2.0d0 * rpole(6,i)
         cmp(9,i) = 2.0d0 * rpole(7,i)
         cmp(10,i) = 2.0d0 * rpole(10,i)
      end do
c
c     compute B-spline coefficients and spatial decomposition
c
      call bspline_fill
      call table_fill
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp (cmp,fmp)
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_mpole (fmp)
      call fftfront
c
c     make the scalar summation over reciprocal lattice
c
      qfac(1,1,1) = 0.0d0
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      ntot = nfft1 * nfft2 * nfft3
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the PME grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get field
c
      call fftback
      call fphi_mpole (fphi)
c
c     convert the field from fractional to Cartesian
c
      call fphi_to_cphi (fphi,cphi)
c
c     increment the field at each multipole site
c
      do i = 1, npole
         pot_recip(i) = pot_recip(i) + cphi(1,i)
c
         field_recip(1,i) = field_recip(1,i) + cphi(2,i)
         field_recip(2,i) = field_recip(2,i) + cphi(3,i)
         field_recip(3,i) = field_recip(3,i) + cphi(4,i)
c
         gradfield_recip(1,1,i) = gradfield_recip(1,1,i) + cphi(5,i)
         gradfield_recip(2,2,i) = gradfield_recip(2,2,i) + cphi(6,i)
         gradfield_recip(3,3,i) = gradfield_recip(3,3,i) + cphi(7,i)
         gradfield_recip(1,2,i) = gradfield_recip(1,2,i) + cphi(8,i)
         gradfield_recip(1,3,i) = gradfield_recip(1,3,i) + cphi(9,i)
         gradfield_recip(2,3,i) = gradfield_recip(2,3,i) + cphi(10,i)
c
         hessfield_recip(1,1,1,i) = hessfield_recip(1,1,1,i) +cphi(11,i)
         hessfield_recip(2,2,2,i) = hessfield_recip(2,2,2,i) +cphi(12,i)
         hessfield_recip(3,3,3,i) = hessfield_recip(3,3,3,i) +cphi(13,i)
         hessfield_recip(1,1,2,i) = hessfield_recip(1,1,2,i) +cphi(14,i)
         hessfield_recip(1,1,3,i) = hessfield_recip(1,1,3,i) +cphi(15,i)
         hessfield_recip(1,2,2,i) = hessfield_recip(1,2,2,i) +cphi(16,i)
         hessfield_recip(2,2,3,i) = hessfield_recip(2,2,3,i) +cphi(17,i)
         hessfield_recip(1,3,3,i) = hessfield_recip(1,3,3,i) +cphi(18,i)
         hessfield_recip(2,3,3,i) = hessfield_recip(2,3,3,i) +cphi(19,i)
         hessfield_recip(1,2,3,i) = hessfield_recip(1,2,3,i) +cphi(20,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cmp)
      deallocate (fmp)
      deallocate (cphi)
      deallocate (fphi)
      return
      end
