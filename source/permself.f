c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine permself  --  ewald permanent self fields     ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "permself" computes the electrostatic self field
c     due to permanent monents and computes the 
c     cell dipole boundary correction if needed
c
c
      subroutine permself
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
      integer i,j,ii
      real*8 term
      real*8 ucell(3)
c
c
c     perform dynamic allocation of global arrays
c
      if (.not.allocated(field_self)) allocate (field_self(3,npole))
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field_self(j,i) = 0.0d0
         end do
      end do
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
            field_self(j,i) = field_self(j,i) - term*rpole(j+1,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do i = 1, npole
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field_self(j,i) = field_self(j,i) + term*ucell(j)
            end do
         end do
      end if
      return
      end
