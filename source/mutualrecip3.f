c
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine mutualrecip3  --  ewald recip mutual induced field  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "mutualrecip3" computes the reciprocal space contribution of the
c     induced atomic dipole moments to the field
c
c
      subroutine mutualrecip3
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
      real*8, allocatable :: fphid(:,:)
      real*8, allocatable :: fphip(:,:)
      real*8, allocatable :: fphidp(:,:)
      real*8, allocatable :: cphid(:,:)
      real*8, allocatable :: cphip(:,:)
      real*8, allocatable :: cphidp(:,:)
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
      if (.not.allocated(udgradfield_recip))
     &     allocate (udgradfield_recip(3,3,npole))
      if (.not.allocated(upgradfield_recip))
     &     allocate (upgradfield_recip(3,3,npole))
      if (.not.allocated(udhessfield_recip))
     &     allocate (udhessfield_recip(3,3,3,npole))
      if (.not.allocated(uphessfield_recip))
     &     allocate (uphessfield_recip(3,3,3,npole))
      if (.not.allocated(udphessfield_recip))
     &     allocate (udphessfield_recip(3,3,3,npole))
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            udfield_recip(j,i) = 0.0d0
            upfield_recip(j,i) = 0.0d0
            do k = 1, 3
               udgradfield_recip(k,j,i) = 0.0d0
               upgradfield_recip(k,j,i) = 0.0d0
               do l = 1, 3
                  udhessfield_recip(l,k,j,i) = 0.0d0
                  uphessfield_recip(l,k,j,i) = 0.0d0
                  udphessfield_recip(l,k,j,i) = 0.0d0
               end do
            end do
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,npole))
      allocate (fuinp(3,npole))
c     d dipoles
      allocate (fphid(20,npole))
c     p dipoles
      allocate (fphip(20,npole))
c     p+d dipoles
      allocate (fphidp(20,npole))
      allocate (cphid(20,npole))
      allocate (cphip(20,npole))
      allocate (cphidp(20,npole))
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
      call fphi_uind (fphid,fphip,fphidp)
c
c     convert the dipole fields from fractional to Cartesian
c
      call fphi_to_cphi (fphid,cphid)
      call fphi_to_cphi (fphip,cphip)
      call fphi_to_cphi (fphidp,cphidp)
c
c     accumulate mutual cartesian reciprocal fields
c
      do i = 1, npole
         udfield_recip(1,i) = udfield_recip(1,i) + cphid(2,i)
         udfield_recip(2,i) = udfield_recip(2,i) + cphid(3,i)
         udfield_recip(3,i) = udfield_recip(3,i) + cphid(4,i)
c
         upfield_recip(1,i) = upfield_recip(1,i) + cphip(2,i)
         upfield_recip(2,i) = upfield_recip(2,i) + cphip(3,i)
         upfield_recip(3,i) = upfield_recip(3,i) + cphip(4,i)
c     FIELD GRADIENT
         udgradfield_recip(1,1,i) = udgradfield_recip(1,1,i) + 
     &        cphid(5,i)
         udgradfield_recip(2,2,i) = udgradfield_recip(2,2,i) + 
     &        cphid(6,i)
         udgradfield_recip(3,3,i) = udgradfield_recip(3,3,i) + 
     &        cphid(7,i)
         udgradfield_recip(2,1,i) = udgradfield_recip(2,1,i) + 
     &        cphid(8,i)
         udgradfield_recip(3,1,i) = udgradfield_recip(3,1,i) + 
     &        cphid(9,i)
         udgradfield_recip(3,2,i) = udgradfield_recip(3,2,i) + 
     &        cphid(10,i)
c
         upgradfield_recip(1,1,i) = upgradfield_recip(1,1,i) +
     &        cphip(5,i)
         upgradfield_recip(2,2,i) = upgradfield_recip(2,2,i) +
     &        cphip(6,i)
         upgradfield_recip(3,3,i) = upgradfield_recip(3,3,i) +
     &        cphip(7,i)
         upgradfield_recip(2,1,i) = upgradfield_recip(2,1,i) +
     &        cphip(8,i)
         upgradfield_recip(3,1,i) = upgradfield_recip(3,1,i) +
     &        cphip(9,i)
         upgradfield_recip(3,2,i) = upgradfield_recip(3,2,i) +
     &        cphip(10,i)
c
c         upgradfield_recip(1,1,i) = upgradfield_recip(1,1,i) +
c     &        0.5d0*cphidp(5,i)
c         upgradfield_recip(2,2,i) = upgradfield_recip(2,2,i) +
c     &        0.5d0*cphidp(6,i)
c         upgradfield_recip(3,3,i) = upgradfield_recip(3,3,i) +
c     &        0.5d0*cphidp(7,i)
c         upgradfield_recip(2,1,i) = upgradfield_recip(2,1,i) +
c     &        0.5d0*cphidp(8,i)
c         upgradfield_recip(3,1,i) = upgradfield_recip(3,1,i) +
c     &        0.5d0*cphidp(9,i)
c         upgradfield_recip(3,2,i) = upgradfield_recip(3,2,i) +
c     &        0.5d0*cphidp(10,i)
c
c         udhessfield_recip(1,1,1,i) = udhessfield_recip(1,1,1,i) +
c     &        cphid(11,i)
c         udhessfield_recip(2,2,2,i) = udhessfield_recip(2,2,2,i) +
c     &        cphid(12,i)
c         udhessfield_recip(3,3,3,i) = udhessfield_recip(3,3,3,i) +
c     &        cphid(13,i)
c         udhessfield_recip(2,1,1,i) = udhessfield_recip(2,1,1,i) +
c     &        cphid(14,i)
c         udhessfield_recip(3,1,1,i) = udhessfield_recip(3,1,1,i) +
c     &        cphid(15,i)
c         udhessfield_recip(2,2,1,i) = udhessfield_recip(2,2,1,i) +
c     &        cphid(16,i)
c         udhessfield_recip(3,2,2,i) = udhessfield_recip(3,2,2,i) +
c     &        cphid(17,i)
c         udhessfield_recip(3,3,1,i) = udhessfield_recip(3,3,1,i) +
c     &        cphid(18,i)
c         udhessfield_recip(3,3,2,i) = udhessfield_recip(3,3,2,i) +
c     &        cphid(19,i)
c         udhessfield_recip(3,2,1,i) = udhessfield_recip(3,2,1,i) +
c     &        cphid(20,i)
c
c         uphessfield_recip(1,1,1,i) = uphessfield_recip(1,1,1,i) +
c     &        cphip(11,i)
c         uphessfield_recip(2,2,2,i) = uphessfield_recip(2,2,2,i) +
c     &        cphip(12,i)
c         uphessfield_recip(3,3,3,i) = uphessfield_recip(3,3,3,i) +
c     &        cphip(13,i)
c         uphessfield_recip(2,1,1,i) = uphessfield_recip(2,1,1,i) +
c     &        cphip(14,i)
c         uphessfield_recip(3,1,1,i) = uphessfield_recip(3,1,1,i) +
c     &        cphip(15,i)
c         uphessfield_recip(2,2,1,i) = uphessfield_recip(2,2,1,i) +
c     &        cphip(16,i)
c         uphessfield_recip(3,2,2,i) = uphessfield_recip(3,2,2,i) +
c     &        cphip(17,i)
c         uphessfield_recip(3,3,1,i) = uphessfield_recip(3,3,1,i) +
c     &        cphip(18,i)
c         uphessfield_recip(3,3,2,i) = uphessfield_recip(3,3,2,i) +
c     &        cphip(19,i)
c         uphessfield_recip(3,2,1,i) = uphessfield_recip(3,2,1,i) +
c     &        cphip(20,i)
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c     this is a hack. must fix fphi_und to give up to 20 for d & p
         udphessfield_recip(1,1,1,i) = udphessfield_recip(1,1,1,i) +
     &        cphidp(11,i)
         udphessfield_recip(2,2,2,i) = udphessfield_recip(2,2,2,i) +
     &        cphidp(12,i)
         udphessfield_recip(3,3,3,i) = udphessfield_recip(3,3,3,i) +
     &        cphidp(13,i)
         udphessfield_recip(2,1,1,i) = udphessfield_recip(2,1,1,i) +
     &        cphidp(14,i)
         udphessfield_recip(3,1,1,i) = udphessfield_recip(3,1,1,i) +
     &        cphidp(15,i)
         udphessfield_recip(2,2,1,i) = udphessfield_recip(2,2,1,i) +
     &        cphidp(16,i)
         udphessfield_recip(3,2,2,i) = udphessfield_recip(3,2,2,i) +
     &        cphidp(17,i)
         udphessfield_recip(3,3,1,i) = udphessfield_recip(3,3,1,i) +
     &        cphidp(18,i)
         udphessfield_recip(3,3,2,i) = udphessfield_recip(3,3,2,i) +
     &        cphidp(19,i)
         udphessfield_recip(3,2,1,i) = udphessfield_recip(3,2,1,i) +
     &        cphidp(20,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (fphid)
      deallocate (fphip)
      deallocate (fphidp)
      deallocate (cphid)
      deallocate (cphip)
      deallocate (cphidp)
      return
      end
