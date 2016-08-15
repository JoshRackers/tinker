c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine mutualself  --  ewald induced dipole self field  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "mutualself" computes the self mutual electrostatic field due to
c     induced dipole moments and the cell dipole boundary correction if needed
c
c
      subroutine mutualself
      use sizes
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
      use polar
      use potderivs
      implicit none
      integer i,j
      real*8 term
      real*8 ucell(3)
      real*8 ucellp(3)
c
c     perform dynamic allocation of some global arrays
c
      if (.not.allocated(udfield_self)) allocate (udfield_self(3,npole))
      if (.not.allocated(upfield_self)) allocate (upfield_self(3,npole))
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            udfield_self(j,i) = 0.0d0
            upfield_self(j,i) = 0.0d0
         end do
      end do
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
c     changed + to -
            udfield_self(j,i) = udfield_self(j,i) - term*uind(j,i)
            upfield_self(j,i) = upfield_self(j,i) - term*uinp(j,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to the field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
            ucellp(i) = 0.0d0
         end do
         do i = 1, npole
            do j = 1, 3
               ucell(j) = ucell(j) + uind(j,i)
               ucellp(j) = ucellp(j) + uinp(j,i)
            end do
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
c     changed - to +
               udfield_self(j,i) = udfield_self(j,i) + term*ucell(j)
               upfield_self(j,i) = upfield_self(j,i) + term*ucellp(j)
            end do
         end do
      end if
      return
      end
