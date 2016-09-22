c
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine mutualrecip1  --  ewald recip mutual induced field  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "mutualrecip1" computes the reciprocal space contribution of the
c     induced atomic dipole moments to the field
c
c
      subroutine mutualrecip1
      use sizes
      use boxes
      use ewald
      use math
      use mpole
      use pme
      use polar
      use potderivs
      implicit none
      integer i,j,k,l
      real*8 term
      real*8 a(3,3)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: fdip_phi1(:,:)
      real*8, allocatable :: fdip_phi2(:,:)
      real*8, allocatable :: fdip_sum_phi(:,:)
      real*8, allocatable :: dipfield1(:,:)
      real*8, allocatable :: dipfield2(:,:)
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some global arrays
c
      if (.not.allocated(udfield_recip)) 
     &     allocate (udfield_recip(3,npole))
      if (.not.allocated(upfield_recip)) 
     &     allocate (upfield_recip(3,npole))
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            udfield_recip(j,i) = 0.0d0
            upfield_recip(j,i) = 0.0d0
c            do k = 1, 3
c               udgradfield_recip(k,j,i) = 0.0d0
c               upgradfield_recip(k,j,i) = 0.0d0
c               do l = 1, 3
c                  udhessfield_recip(l,k,j,i) = 0.0d0
c                  uphessfield_recip(l,k,j,i) = 0.0d0
c               end do
c            end do
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,npole))
      allocate (fuinp(3,npole))
      allocate (fdip_phi1(20,npole))
      allocate (fdip_phi2(20,npole))
      allocate (fdip_sum_phi(20,npole))
      allocate (dipfield1(3,npole))
      allocate (dipfield2(3,npole))
c
c     convert Cartesian dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
      do i = 1, npole
         do k = 1, 3
            fuind(k,i) = a(k,1)*uind(1,i) + a(k,2)*uind(2,i)
     &                      + a(k,3)*uind(3,i)
            fuinp(k,i) = a(k,1)*uinp(1,i) + a(k,2)*uinp(2,i)
     &                      + a(k,3)*uinp(3,i)
         end do
      end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_uind (fuind,fuinp)
      call fftfront
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
      call fphi_uind (fdip_phi1,fdip_phi2,fdip_sum_phi)
c
c     convert the dipole fields from fractional to Cartesian
c
      do i = 1, 3
         a(i,1) = dble(nfft1) * recip(i,1)
         a(i,2) = dble(nfft2) * recip(i,2)
         a(i,3) = dble(nfft3) * recip(i,3)
      end do
      do i = 1, npole
         do k = 1, 3
            dipfield1(k,i) = a(k,1)*fdip_phi1(2,i)
     &                          + a(k,2)*fdip_phi1(3,i)
     &                          + a(k,3)*fdip_phi1(4,i)
            dipfield2(k,i) = a(k,1)*fdip_phi2(2,i)
     &                          + a(k,2)*fdip_phi2(3,i)
     &                          + a(k,3)*fdip_phi2(4,i)
         end do
      end do
c
c     increment the field at each multipole site
c
      do i = 1, npole
         do k = 1, 3
            udfield_recip(k,i) = udfield_recip(k,i) + dipfield1(k,i)
            upfield_recip(k,i) = upfield_recip(k,i) + dipfield2(k,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (fdip_phi1)
      deallocate (fdip_phi2)
      deallocate (fdip_sum_phi)
      deallocate (dipfield1)
      deallocate (dipfield2)
      return
      end
