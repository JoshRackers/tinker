c
c
c     #############################################################
c     ## COPYRIGHT (C) 2016 by Josh Rackers & Jay William Ponder ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mutualfield3  --  electric field gradient     ##
c     ##                       and hessian from induce dipoles     ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mutualfield3" computes the electric field gradient and hessian
c     from the induced dipoles
c
c     (the 3 stands for field hessian)
c
c     assumes multipole components have already been rotated into
c     the global coordinate frame
c
c
      subroutine mutualfield3
      use sizes
      use atoms
      use limits
      use mpole
      use mplpot
      use polar
      use potent
      use potderivs
      implicit none
c
c     perform dynamic allocation of global arrays
c
      if (.not.allocated(udfield)) allocate (udfield(3,npole))
      if (.not.allocated(upfield)) allocate (upfield(3,npole))
      if (.not.allocated(udgradfield)) allocate (udgradfield(3,3,npole))
      if (.not.allocated(upgradfield)) allocate (upgradfield(3,3,npole))
      if (.not.allocated(udhessfield)) 
     &     allocate (udhessfield(3,3,3,npole))
      if (.not.allocated(uphessfield)) 
     &     allocate (uphessfield(3,3,3,npole))
c
      if (.not.allocated(udfield_ewald))
     &     allocate (udfield_ewald(3,npole))
      if (.not.allocated(upfield_ewald))
     &     allocate (upfield_ewald(3,npole))
      if (.not.allocated(udgradfield_ewald)) 
     &     allocate (udgradfield_ewald(3,3,npole))
      if (.not.allocated(upgradfield_ewald)) 
     &     allocate (upgradfield_ewald(3,3,npole))
      if (.not.allocated(udhessfield_ewald))
     &     allocate (udhessfield_ewald(3,3,3,npole))
      if (.not.allocated(uphessfield_ewald))
     &     allocate (uphessfield_ewald(3,3,3,npole))
c
      if (.not.allocated(udfield_thole))
     &     allocate (udfield_thole(3,npole))
      if (.not.allocated(upfield_thole))
     &     allocate (upfield_thole(3,npole))
      if (.not.allocated(udgradfield_thole)) 
     &     allocate (udgradfield_thole(3,3,npole))
      if (.not.allocated(upgradfield_thole)) 
     &     allocate (upgradfield_thole(3,3,npole))
      if (.not.allocated(udhessfield_thole))
     &     allocate (udhessfield_thole(3,3,3,npole))
      if (.not.allocated(uphessfield_thole))
     &     allocate (uphessfield_thole(3,3,3,npole))
c
      if (.not.allocated(udfieldd_thole))
     &     allocate (udfieldd_thole(3,npole))
      if (.not.allocated(udgradfieldd_thole))
     &     allocate (udgradfieldd_thole(3,3,npole))
      if (.not.allocated(udhessfieldd_thole))
     &     allocate (udhessfieldd_thole(3,3,3,npole))
c
      if (.not.allocated(upfieldp_thole))
     &     allocate (upfieldp_thole(3,npole))
      if (.not.allocated(upgradfieldp_thole))
     &     allocate (upgradfieldp_thole(3,3,npole))
      if (.not.allocated(uphessfieldp_thole))
     &     allocate (uphessfieldp_thole(3,3,3,npole))


c
c     get the electrostatic potential, field and field gradient
c     due to the induced dipoles
c
      if (use_mlist) then
         call mutualfield3b
      else
         call mutualfield3a
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mutualfield3a  --  mutual field                ##
c     ##                                gradient via double loop    ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mutualfield3a" computes the mutual 
c     field gradient due to induced dipoles via a double loop
c
c
      subroutine mutualfield3a 
      use sizes
      use atoms
      use bound
      use cell
      use couple
      use group
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use potderivs
      use shunt
      use usage
      implicit none
      integer i,j,k,l,m,h
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      integer order,rorder
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr1,rr3,rr5,rr7,rr9
      real*8 t2(3,3),t3(3,3,3),t4(3,3,3,3)
      real*8 t2rr3(3,3),t2rr5(3,3)
      real*8 t3rr5(3,3,3),t3rr7(3,3,3)
      real*8 t4rr5(3,3,3,3),t4rr7(3,3,3,3),t4rr9(3,3,3,3)
      real*8 fieldid(3),fieldkd(3)
      real*8 fieldip(3),fieldkp(3)
      real*8 gradfieldid(3,3),gradfieldkd(3,3)
      real*8 gradfieldip(3,3),gradfieldkp(3,3)
      real*8 hessfieldid(3,3,3),hessfieldkd(3,3,3)
      real*8 hessfieldip(3,3,3),hessfieldkp(3,3,3)
      real*8, allocatable :: scale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      logical proceed
      logical usei,usek
      character*6 mode
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            udfield(j,i) = 0.0d0
            upfield(j,i) = 0.0d0
            udfield_ewald(j,i) = 0.0d0
            upfield_ewald(j,i) = 0.0d0
            udfield_thole(j,i) = 0.0d0
            upfield_thole(j,i) = 0.0d0
            udfieldd_thole(j,i) = 0.0d0
            upfieldp_thole(j,i) = 0.0d0
            do k = 1, 3
               udgradfield(k,j,i) = 0.0d0
               upgradfield(k,j,i) = 0.0d0
               udgradfield_ewald(k,j,i) = 0.0d0
               upgradfield_ewald(k,j,i) = 0.0d0
               udgradfield_thole(k,j,i) = 0.0d0
               upgradfield_thole(k,j,i) = 0.0d0
               udgradfieldd_thole(k,j,i) = 0.0d0
               upgradfieldp_thole(k,j,i) = 0.0d0
               do l = 1, 3
                  udhessfield(l,k,j,i) = 0.0d0
                  uphessfield(l,k,j,i) = 0.0d0
                  udhessfield_ewald(l,k,j,i) = 0.0d0
                  uphessfield_ewald(l,k,j,i) = 0.0d0
                  udhessfield_thole(l,k,j,i) = 0.0d0
                  uphessfield_thole(l,k,j,i) = 0.0d0
                  udhessfieldd_thole(l,k,j,i) = 0.0d0
                  uphessfieldp_thole(l,k,j,i) = 0.0d0
               end do
            end do
         end do
      end do
c
c     set the switching function coefficients
c
      if (damp_ewald) then
         mode = 'EWALD'
         call switch (mode)
      else
         mode = 'MPOLE'
         call switch (mode)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set highest order rr and damping terms needed
c     3 = up to field hessian (rr9 for induced dipoles)
c
      order = 3
      rorder = order*2 + 3
      allocate (scale(rorder))
      do i = 1,rorder
          scale(i) = 0.0d0
      end do
c
c     find the electrostatic potential, field and field gradient 
c     due to induced dipoles
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
c
c     set d and p exclusion rules
c
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  call t2matrixrr3(xr,yr,zr,rr3,t2rr3)
                  call t2matrixrr5(xr,yr,zr,rr5,t2rr5)
                  call t3matrixrr5(xr,yr,zr,rr5,t3rr5)
                  call t3matrixrr7(xr,yr,zr,rr7,t3rr7)
                  call t4matrixrr5(xr,yr,zr,rr5,t4rr5)
                  call t4matrixrr7(xr,yr,zr,rr7,t4rr7)
                  call t4matrixrr9(xr,yr,zr,rr9,t4rr9)
c
c     call routines that produce potential, field, field gradient
c     for types of damping
c
                  if (damp_none) then
                     t2 = t2rr3 + t2rr5
                     t3 = t3rr5 + t3rr7
                     t4 = t4rr5 + t4rr7 + t4rr9
                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                    fieldkp)
                     call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                    gradfieldip,gradfieldkp)
                     call uhessfieldik(i,k,t4,hessfieldid,hessfieldkd,
     &                    hessfieldip,hessfieldkp)
                     do j = 1, 3
                        udfield(j,i) = udfield(j,i) + fieldid(j)
                        udfield(j,k) = udfield(j,k) + fieldkd(j)
                        upfield(j,i) = upfield(j,i) + fieldip(j)
                        upfield(j,k) = upfield(j,k) + fieldkp(j)
                        do l = 1, 3
                           udgradfield(l,j,i) = udgradfield(l,j,i) + 
     &                          gradfieldid(l,j)
                           udgradfield(l,j,k) = udgradfield(l,j,k) + 
     &                          gradfieldkd(l,j)
                           upgradfield(l,j,i) = upgradfield(l,j,i) + 
     &                          gradfieldip(l,j)
                           upgradfield(l,j,k) = upgradfield(l,j,k) + 
     &                          gradfieldkp(l,j)
                           do h = 1, 3
                              udhessfield(h,l,j,i) = 
     &                             udhessfield(h,l,j,i) +
     &                             hessfieldid(h,l,j)
                              udhessfield(h,l,j,k) = 
     &                             udhessfield(h,l,j,k) +
     &                             hessfieldkd(h,l,j)
                              uphessfield(h,l,j,i) = 
     &                             uphessfield(h,l,j,i) +
     &                             hessfieldip(h,l,j)
                              uphessfield(h,l,j,k) = 
     &                             uphessfield(h,l,j,k) +
     &                             hessfieldkp(h,l,j)
                           end do
                        end do
                     end do
                  end if
c
c     error function damping for ewald
c
                  if (damp_ewald) then
                     call dampewald(i,k,rorder,r,r2,scale)
c
c     the ewald damping factors already contain their powers of r (rrx)
c
                     t2 = t2rr3*scale(3)/rr3 + t2rr5*scale(5)/rr5
                     t3 = t3rr5*scale(5)/rr5 + t3rr7*scale(7)/rr7
                     t4 = t4rr5*scale(5)/rr5 + t4rr7*scale(7)/rr7 +
     &                    t4rr9*scale(9)/rr9
                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                    fieldkp)
                     call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                    gradfieldip,gradfieldkp)
                     call uhessfieldik(i,k,t4,hessfieldid,hessfieldkd,
     &                    hessfieldip,hessfieldkp)
                     do j = 1, 3
                        udfield_ewald(j,i) = udfield_ewald(j,i) +
     &                       fieldid(j)
                        udfield_ewald(j,k) = udfield_ewald(j,k) +
     &                       fieldkd(j)
                        upfield_ewald(j,i) = upfield_ewald(j,i) +
     &                       fieldip(j)
                        upfield_ewald(j,k) = upfield_ewald(j,k) +
     &                       fieldkp(j)
                        do l = 1, 3
                           udgradfield_ewald(l,j,i) = 
     &                          udgradfield_ewald(l,j,i) + 
     &                          gradfieldid(l,j)
                           udgradfield_ewald(l,j,k) = 
     &                          udgradfield_ewald(l,j,k) + 
     &                          gradfieldkd(l,j)
                           upgradfield_ewald(l,j,i) = 
     &                          upgradfield_ewald(l,j,i) +
     &                          gradfieldip(l,j)
                           upgradfield_ewald(l,j,k) = 
     &                          upgradfield_ewald(l,j,k) +
     &                          gradfieldkp(l,j)
                           do h = 1, 3
                              udhessfield_ewald(h,l,j,i) =
     &                             udhessfield_ewald(h,l,j,i) +
     &                             hessfieldid(h,l,j)
                              udhessfield_ewald(h,l,j,k) =
     &                             udhessfield_ewald(h,l,j,k) +
     &                             hessfieldkd(h,l,j)
                              uphessfield_ewald(h,l,j,i) =
     &                             uphessfield_ewald(h,l,j,i) +
     &                             hessfieldip(h,l,j)
                              uphessfield_ewald(h,l,j,k) =
     &                             uphessfield_ewald(h,l,j,k) +
     &                             hessfieldkp(h,l,j)
                           end do
                        end do
                     end do
                  end if
c
c     thole damping
c
                  if (damp_thole) then
                     call dampthole(i,k,rorder,r,scale)
                     t2 = t2rr3*scale(3) + t2rr5*scale(5)
                     t3 = t3rr5*scale(5) + t3rr7*scale(7)
                     t4 = t4rr5*scale(5) + t4rr7*scale(7) +
     &                    t4rr9*scale(9)
                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                    fieldkp)
                     call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                    gradfieldip,gradfieldkp)
                     call uhessfieldik(i,k,t4,hessfieldid,hessfieldkd,
     &                    hessfieldip,hessfieldkp)
                     do j = 1, 3
                        udfield_thole(j,i) = udfield_thole(j,i) +
     &                       fieldid(j)
                        udfield_thole(j,k) = udfield_thole(j,k) +
     &                       fieldkd(j)
                        upfield_thole(j,i) = upfield_thole(j,i) +
     &                       fieldip(j)
                        upfield_thole(j,k) = upfield_thole(j,k) +
     &                       fieldkp(j)
                        udfieldd_thole(j,i) = udfieldd_thole(j,i) +
     &                       fieldid(j)*dscale(kk)
                        udfieldd_thole(j,k) = udfieldd_thole(j,k) +
     &                       fieldkd(j)*dscale(kk)
                        upfieldp_thole(j,i) = upfieldp_thole(j,i) +
     &                       fieldip(j)*pscale(kk)
                        upfieldp_thole(j,k) = upfieldp_thole(j,k) +
     &                       fieldkp(j)*pscale(kk)
                        do l = 1, 3
                           udgradfield_thole(l,j,i) = 
     &                          udgradfield_thole(l,j,i) +
     &                          gradfieldid(l,j)
                           udgradfield_thole(l,j,k) = 
     &                          udgradfield_thole(l,j,k) +
     &                          gradfieldkd(l,j)
                           upgradfield_thole(l,j,i) = 
     &                          upgradfield_thole(l,j,i) +
     &                          gradfieldip(l,j)
                           upgradfield_thole(l,j,k) = 
     &                          upgradfield_thole(l,j,k) +
     &                          gradfieldkp(l,j)
                           udgradfieldd_thole(l,j,i) =
     &                          udgradfieldd_thole(l,j,i) +
     &                          gradfieldid(l,j)*dscale(kk)
                           udgradfieldd_thole(l,j,k) =
     &                          udgradfieldd_thole(l,j,k) +
     &                          gradfieldkd(l,j)*dscale(kk)
                           upgradfieldp_thole(l,j,i) =
     &                          upgradfieldp_thole(l,j,i) +
     &                          gradfieldip(l,j)*pscale(kk)
                           upgradfieldp_thole(l,j,k) =
     &                          upgradfieldp_thole(l,j,k) +
     &                          gradfieldkp(l,j)*pscale(kk)
                           do h = 1, 3
                              udhessfield_thole(h,l,j,i) =
     &                             udhessfield_thole(h,l,j,i) +
     &                             hessfieldid(h,l,j)
                              udhessfield_thole(h,l,j,k) =
     &                             udhessfield_thole(h,l,j,k) +
     &                             hessfieldkd(h,l,j)
                              uphessfield_thole(h,l,j,i) =
     &                             uphessfield_thole(h,l,j,i) +
     &                             hessfieldip(h,l,j)
                              uphessfield_thole(h,l,j,k) =
     &                             uphessfield_thole(h,l,j,k) +
     &                             hessfieldkp(h,l,j)
                              udhessfieldd_thole(h,l,j,i) =
     &                             udhessfieldd_thole(h,l,j,i) +
     &                             hessfieldid(h,l,j)*dscale(kk)
                              udhessfieldd_thole(h,l,j,k) =
     &                             udhessfieldd_thole(h,l,j,k) +
     &                             hessfieldkd(h,l,j)*dscale(kk)
                              uphessfieldp_thole(h,l,j,i) =
     &                             uphessfieldp_thole(h,l,j,i) +
     &                             hessfieldip(h,l,j)*pscale(kk)
                              uphessfieldp_thole(h,l,j,k) =
     &                             uphessfieldp_thole(h,l,j,k) +
     &                             hessfieldkp(h,l,j)*pscale(kk)
                           end do
                        end do
                     end do
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c     
c     calculate interaction with other unit cells
c
      do i = 1, npole
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         ii = ipole(i)
c
c     set d and p exclusion rules
c
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
c
c     decide whether to compute the current interaction
c
         do k = i, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
c
c     compute the contributions for this interaction
c
            if (proceed) then
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     rr1 = 1.0d0 / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
                     rr9 = 7.0d0 * rr7 / r2
                     call t2matrixrr3(xr,yr,zr,rr3,t2rr3)
                     call t2matrixrr5(xr,yr,zr,rr5,t2rr5)
                     call t3matrixrr5(xr,yr,zr,rr5,t3rr5)
                     call t3matrixrr7(xr,yr,zr,rr7,t3rr7)
                     call t4matrixrr5(xr,yr,zr,rr5,t4rr5)
                     call t4matrixrr7(xr,yr,zr,rr7,t4rr7)
                     call t4matrixrr9(xr,yr,zr,rr9,t4rr9)
c
c     call routines that produce potential, field, field gradient
c     for types of damping
c
c
c     no damping
c
                     if (damp_none) then
                        t2 = t2rr3 + t2rr5
                        t3 = t3rr5 + t3rr7
                        t4 = t4rr5 + t4rr7 + t4rr9
                        call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                       fieldkp)
                        call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd
     &                       ,gradfieldip,gradfieldkp)
                        call uhessfieldik(i,k,t4,hessfieldid,hessfieldkd
     &                       ,hessfieldip,hessfieldkp)
                        do j = 1, 3
                           udfield(j,i) = udfield(j,i) + fieldid(j)
                           udfield(j,k) = udfield(j,k) + fieldkd(j)
                           upfield(j,i) = upfield(j,i) + fieldip(j)
                           upfield(j,k) = upfield(j,k) + fieldkp(j)
                           do l = 1, 3
                              udgradfield(l,j,i) = udgradfield(l,j,i) + 
     &                             gradfieldid(l,j)
                              udgradfield(l,j,k) = udgradfield(l,j,k) + 
     &                             gradfieldkd(l,j)
                              upgradfield(l,j,i) = upgradfield(l,j,i) + 
     &                             gradfieldip(l,j)
                              upgradfield(l,j,k) = upgradfield(l,j,k) + 
     &                             gradfieldkp(l,j)
                              do h = 1, 3
                                 udhessfield(h,l,j,i) = 
     &                                udhessfield(h,l,j,i) +
     &                                hessfieldid(h,l,j)
                                 udhessfield(h,l,j,k) = 
     &                                udhessfield(h,l,j,k) +
     &                                hessfieldkd(h,l,j)
                                 uphessfield(h,l,j,i) = 
     &                                uphessfield(h,l,j,i) +
     &                                hessfieldip(h,l,j)
                                 uphessfield(h,l,j,k) = 
     &                                uphessfield(h,l,j,k) +
     &                                hessfieldkp(h,l,j)
                              end do
                           end do
                        end do
                     end if
c     
c     error function damping for ewald
c     
                     if (damp_ewald) then
                        call dampewald(i,k,rorder,r,r2,scale)
c     
c     the ewald damping factors already contain their powers of r (rrx)
c     
                        t2 = t2rr3*scale(3)/rr3 + t2rr5*scale(5)/rr5
                        t3 = t3rr5*scale(5)/rr5 + t3rr7*scale(7)/rr7
                        t4 = t4rr5*scale(5)/rr5 + t4rr7*scale(7)/rr7 +
     &                       t4rr9*scale(9)/rr9
                        call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                       fieldkp)
                        call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd
     &                       ,gradfieldip,gradfieldkp)
                        call uhessfieldik(i,k,t4,hessfieldid,hessfieldkd
     &                       ,hessfieldip,hessfieldkp)
                        do j = 1, 3
                           udfield_ewald(j,i) = udfield_ewald(j,i) +
     &                          fieldid(j)
                           udfield_ewald(j,k) = udfield_ewald(j,k) +
     &                          fieldkd(j)
                           upfield_ewald(j,i) = upfield_ewald(j,i) +
     &                          fieldip(j)
                           upfield_ewald(j,k) = upfield_ewald(j,k) +
     &                          fieldkp(j)
                           do l = 1, 3
                              udgradfield_ewald(l,j,i) = 
     &                             udgradfield_ewald(l,j,i) + 
     &                             gradfieldid(l,j)
                              udgradfield_ewald(l,j,k) = 
     &                             udgradfield_ewald(l,j,k) + 
     &                             gradfieldkd(l,j)
                              upgradfield_ewald(l,j,i) = 
     &                             upgradfield_ewald(l,j,i) +
     &                             gradfieldip(l,j)
                              upgradfield_ewald(l,j,k) = 
     &                             upgradfield_ewald(l,j,k) +
     &                             gradfieldkp(l,j)
                              do h = 1, 3
                                 udhessfield_ewald(h,l,j,i) =
     &                                udhessfield_ewald(h,l,j,i) +
     &                                hessfieldid(h,l,j)
                                 udhessfield_ewald(h,l,j,k) =
     &                                udhessfield_ewald(h,l,j,k) +
     &                                hessfieldkd(h,l,j)
                                 uphessfield_ewald(h,l,j,i) =
     &                                uphessfield_ewald(h,l,j,i) +
     &                                hessfieldip(h,l,j)
                                 uphessfield_ewald(h,l,j,k) =
     &                                uphessfield_ewald(h,l,j,k) +
     &                                hessfieldkp(h,l,j)
                              end do
                           end do
                        end do
                     end if
c     
c     thole damping
c     
                     if (damp_thole) then
                        call dampthole(i,k,rorder,r,scale)
                        t2 = t2rr3*scale(3) + t2rr5*scale(5)
                        t3 = t3rr5*scale(5) + t3rr7*scale(7)
                        t4 = t4rr5*scale(5) + t4rr7*scale(7) +
     &                       t4rr9*scale(9)
                        call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                       fieldkp)
                        call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd
     &                       ,gradfieldip,gradfieldkp)
                        call uhessfieldik(i,k,t4,hessfieldid,hessfieldkd
     &                       ,hessfieldip,hessfieldkp)
                        do j = 1, 3
                           udfield_thole(j,i) = udfield_thole(j,i) +
     &                          fieldid(j)
                           udfield_thole(j,k) = udfield_thole(j,k) +
     &                          fieldkd(j)
                           upfield_thole(j,i) = upfield_thole(j,i) +
     &                          fieldip(j)
                           upfield_thole(j,k) = upfield_thole(j,k) +
     &                          fieldkp(j)
                           udfieldd_thole(j,i) = udfieldd_thole(j,i) +
     &                          fieldid(j)*dscale(kk)
                           udfieldd_thole(j,k) = udfieldd_thole(j,k) +
     &                          fieldkd(j)*dscale(kk)
                           upfieldp_thole(j,i) = upfieldp_thole(j,i) +
     &                          fieldip(j)*pscale(kk)
                           upfieldp_thole(j,k) = upfieldp_thole(j,k) +
     &                          fieldkp(j)*pscale(kk)
                           do l = 1, 3
                              udgradfield_thole(l,j,i) = 
     &                             udgradfield_thole(l,j,i) +
     &                             gradfieldid(l,j)
                              udgradfield_thole(l,j,k) = 
     &                             udgradfield_thole(l,j,k) +
     &                             gradfieldkd(l,j)
                              upgradfield_thole(l,j,i) = 
     &                             upgradfield_thole(l,j,i) +
     &                             gradfieldip(l,j)
                              upgradfield_thole(l,j,k) = 
     &                             upgradfield_thole(l,j,k) +
     &                             gradfieldkp(l,j)
                              udgradfieldd_thole(l,j,i) =
     &                             udgradfieldd_thole(l,j,i) +
     &                             gradfieldid(l,j)*dscale(kk)
                              udgradfieldd_thole(l,j,k) =
     &                             udgradfieldd_thole(l,j,k) +
     &                             gradfieldkd(l,j)*dscale(kk)
                              upgradfieldp_thole(l,j,i) =
     &                             upgradfieldp_thole(l,j,i) +
     &                             gradfieldip(l,j)*pscale(kk)
                              upgradfieldp_thole(l,j,k) =
     &                             upgradfieldp_thole(l,j,k) +
     &                             gradfieldkp(l,j)*pscale(kk)
                              do h = 1, 3
                                 udhessfield_thole(h,l,j,i) =
     &                                udhessfield_thole(h,l,j,i) +
     &                                hessfieldid(h,l,j)
                                 udhessfield_thole(h,l,j,k) =
     &                                udhessfield_thole(h,l,j,k) +
     &                                hessfieldkd(h,l,j)
                                 uphessfield_thole(h,l,j,i) =
     &                                uphessfield_thole(h,l,j,i) +
     &                                hessfieldip(h,l,j)
                                 uphessfield_thole(h,l,j,k) =
     &                                uphessfield_thole(h,l,j,k) +
     &                                hessfieldkp(h,l,j)
                                 udhessfieldd_thole(h,l,j,i) =
     &                                udhessfieldd_thole(h,l,j,i) +
     &                                hessfieldid(h,l,j)*dscale(kk)
                                 udhessfieldd_thole(h,l,j,k) =
     &                                udhessfieldd_thole(h,l,j,k) +
     &                                hessfieldkd(h,l,j)*dscale(kk)
                                 uphessfieldp_thole(h,l,j,i) =
     &                                uphessfieldp_thole(h,l,j,i) +
     &                                hessfieldip(h,l,j)*pscale(kk)
                                 uphessfieldp_thole(h,l,j,k) =
     &                                uphessfieldp_thole(h,l,j,k) +
     &                                hessfieldkp(h,l,j)*pscale(kk)
                              end do
                           end do
                        end do
                     end if
                  end if
               end do
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (scale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mutualfield3b  --  potential, field and field  ##
c     ##                                gradient via neighbor list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mutualfield3b" computes the direct electrostatic potential, field and 
c     field gradient due to induced dipoles using a neighbor list
c
c
      subroutine mutualfield3b
      use sizes
      use atoms
      use bound
      use cell
      use couple
      use group
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potderivs
      use shunt
      use usage
      implicit none
      integer i,j,k,l,m,h
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      integer order,rorder
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr1,rr3,rr5,rr7,rr9
      real*8 t2(3,3),t3(3,3,3),t4(3,3,3,3)
      real*8 t2rr3(3,3),t2rr5(3,3)
      real*8 t3rr5(3,3,3),t3rr7(3,3,3)
      real*8 t4rr5(3,3,3,3),t4rr7(3,3,3,3),t4rr9(3,3,3,3)
      real*8 fieldid(3),fieldkd(3)
      real*8 fieldip(3),fieldkp(3)
      real*8 gradfieldid(3,3),gradfieldkd(3,3)
      real*8 gradfieldip(3,3),gradfieldkp(3,3)
      real*8 hessfieldid(3,3,3),hessfieldkd(3,3,3)
      real*8 hessfieldip(3,3,3),hessfieldkp(3,3,3)
      real*8, allocatable :: scale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: udfieldo(:,:)
      real*8, allocatable :: upfieldo(:,:)
      real*8, allocatable :: udgradfieldo(:,:,:)
      real*8, allocatable :: upgradfieldo(:,:,:)
      real*8, allocatable :: udhessfieldo(:,:,:,:)
      real*8, allocatable :: uphessfieldo(:,:,:,:)
      real*8, allocatable :: udfield_ewaldo(:,:)
      real*8, allocatable :: upfield_ewaldo(:,:)
      real*8, allocatable :: udgradfield_ewaldo(:,:,:)
      real*8, allocatable :: upgradfield_ewaldo(:,:,:)
      real*8, allocatable :: udhessfield_ewaldo(:,:,:,:)
      real*8, allocatable :: uphessfield_ewaldo(:,:,:,:)
      real*8, allocatable :: udfield_tholeo(:,:)
      real*8, allocatable :: upfield_tholeo(:,:)
      real*8, allocatable :: udgradfield_tholeo(:,:,:)
      real*8, allocatable :: upgradfield_tholeo(:,:,:)
      real*8, allocatable :: udhessfield_tholeo(:,:,:,:)
      real*8, allocatable :: uphessfield_tholeo(:,:,:,:)
      real*8, allocatable :: udfieldd_tholeo(:,:)
      real*8, allocatable :: udgradfieldd_tholeo(:,:,:)
      real*8, allocatable :: udhessfieldd_tholeo(:,:,:,:)
      real*8, allocatable :: upfieldp_tholeo(:,:)
      real*8, allocatable :: upgradfieldp_tholeo(:,:,:)
      real*8, allocatable :: uphessfieldp_tholeo(:,:,:,:)
      logical proceed
      logical usei,usek
      character*6 mode
c
c     perform dynamic allocation of some local arrays
c
      allocate (udfieldo(3,npole))
      allocate (upfieldo(3,npole))
      allocate (udgradfieldo(3,3,npole))
      allocate (upgradfieldo(3,3,npole))
      allocate (udhessfieldo(3,3,3,npole))
      allocate (uphessfieldo(3,3,3,npole))
c
      allocate (udfield_ewaldo(3,npole))
      allocate (upfield_ewaldo(3,npole))
      allocate (udgradfield_ewaldo(3,3,npole))
      allocate (upgradfield_ewaldo(3,3,npole))
      allocate (udhessfield_ewaldo(3,3,3,npole))
      allocate (uphessfield_ewaldo(3,3,3,npole))
c
      allocate (udfield_tholeo(3,npole))
      allocate (upfield_tholeo(3,npole))
      allocate (udgradfield_tholeo(3,3,npole))
      allocate (upgradfield_tholeo(3,3,npole))
      allocate (udhessfield_tholeo(3,3,3,npole))
      allocate (uphessfield_tholeo(3,3,3,npole))
c
      allocate (udfieldd_tholeo(3,npole))
      allocate (udgradfieldd_tholeo(3,3,npole))
      allocate (udhessfieldd_tholeo(3,3,3,npole))
c
      allocate (upfieldp_tholeo(3,npole))
      allocate (upgradfieldp_tholeo(3,3,npole))
      allocate (uphessfieldp_tholeo(3,3,3,npole))
c
c     zero out the value of the gradfield at each site
c
      do i = 1, npole
         do j = 1, 3
            udfield(j,i) = 0.0d0
            upfield(j,i) = 0.0d0
            udfield_ewald(j,i) = 0.0d0
            upfield_ewald(j,i) = 0.0d0
            udfield_thole(j,i) = 0.0d0
            upfield_thole(j,i) = 0.0d0
            udfieldd_thole(j,i) = 0.0d0
            upfieldp_thole(j,i) = 0.0d0
c
            udfieldo(j,i) = 0.0d0
            upfieldo(j,i) = 0.0d0
            udfield_ewaldo(j,i) = 0.0d0
            upfield_ewaldo(j,i) = 0.0d0
            udfield_tholeo(j,i) = 0.0d0
            upfield_tholeo(j,i) = 0.0d0
            udfieldd_tholeo(j,i) = 0.0d0
            upfieldp_tholeo(j,i) = 0.0d0
            do k = 1, 3
               udgradfield(k,j,i) = 0.0d0
               upgradfield(k,j,i) = 0.0d0
               udgradfield_ewald(k,j,i) = 0.0d0
               upgradfield_ewald(k,j,i) = 0.0d0
               udgradfield_thole(k,j,i) = 0.0d0
               upgradfield_thole(k,j,i) = 0.0d0
               udgradfieldd_thole(k,j,i) = 0.0d0
               upgradfieldp_thole(k,j,i) = 0.0d0
c
               udgradfieldo(k,j,i) = 0.0d0
               upgradfieldo(k,j,i) = 0.0d0
               udgradfield_ewaldo(k,j,i) = 0.0d0
               upgradfield_ewaldo(k,j,i) = 0.0d0
               udgradfield_tholeo(k,j,i) = 0.0d0
               upgradfield_tholeo(k,j,i) = 0.0d0
               udgradfieldd_tholeo(k,j,i) = 0.0d0
               upgradfieldp_tholeo(k,j,i) = 0.0d0
               do l = 1, 3
                  udhessfield(l,k,j,i) = 0.0d0
                  uphessfield(l,k,j,i) = 0.0d0
                  udhessfield_ewald(l,k,j,i) = 0.0d0
                  uphessfield_ewald(l,k,j,i) = 0.0d0
                  udhessfield_thole(l,k,j,i) = 0.0d0
                  uphessfield_thole(l,k,j,i) = 0.0d0
                  udhessfieldd_thole(l,k,j,i) = 0.0d0
                  uphessfieldp_thole(l,k,j,i) = 0.0d0
c
                  udhessfieldo(l,k,j,i) = 0.0d0
                  uphessfieldo(l,k,j,i) = 0.0d0
                  udhessfield_ewaldo(l,k,j,i) = 0.0d0
                  uphessfield_ewaldo(l,k,j,i) = 0.0d0
                  udhessfield_tholeo(l,k,j,i) = 0.0d0
                  uphessfield_tholeo(l,k,j,i) = 0.0d0
                  udhessfieldd_tholeo(l,k,j,i) = 0.0d0
                  uphessfieldp_tholeo(l,k,j,i) = 0.0d0
               end do
            end do
         end do
      end do
c
c     set the switching function coefficients
c
      if (damp_ewald) then
         mode = 'EWALD'
         call switch (mode)
      else
         mode = 'MPOLE'
         call switch (mode)
      end if
c
c     set highest order rr and damping terms needed
c     2 = up to field gradient (rr9)
c
      order = 3
      rorder = order*2 + 3
      allocate (scale(rorder))
      do i = 1,rorder
          scale(i) = 0.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared)
!$OMP& private(i,j,k,ii,ix,iy,iz,usei,kk,kx,ky,kz,usek,kkk,proceed,
!$OMP& xr,yr,zr,r,r2,rr1,rr3,rr5,rr7,rr9,fgrp,
!$OMP& fieldid,fieldkd,fieldip,fieldkp,
!$OMP& gradfieldid,gradfieldkd,gradfieldip,gradfieldkp,
!$OMP& hessfieldid,hessfieldkd,hessfieldip,hessfieldkp,
!$OMP& t2rr3,t2rr5,t3rr5,t3rr7,t4rr5,t4rr7,t4rr9,
!$OMP& t2,t3,t4,scale)
!$OMP& firstprivate(dscale,pscale)
!$OMP DO reduction(+:udfieldo,upfieldo,
!$OMP& udgradfieldo,upgradfieldo,
!$OMP& udhessfieldo,uphessfieldo,
!$OMP& udfield_ewaldo,upfield_ewaldo,
!$OMP& udgradfield_ewaldo,upgradfield_ewaldo,
!$OMP& udhessfield_ewaldo,uphessfield_ewaldo,
!$OMP& udfield_tholeo,upfield_tholeo,
!$OMP& udgradfield_tholeo,upgradfield_tholeo,
!$OMP& udhessfield_tholeo,uphessfield_tholeo,
!$OMP& udfieldd_tholeo,upfieldp_tholeo,
!$OMP& udgradfieldd_tholeo,upgradfieldp_tholeo,
!$OMP& udhessfieldd_tholeo,uphessfieldp_tholeo)
!$OMP& schedule(guided)
c
c     calculate the multipole interaction
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set d and p exclusion rules
c
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
c
c     decide whether to compute the current interaction
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               call t2matrixrr3(xr,yr,zr,rr3,t2rr3)
               call t2matrixrr5(xr,yr,zr,rr5,t2rr5)
               call t3matrixrr5(xr,yr,zr,rr5,t3rr5)
               call t3matrixrr7(xr,yr,zr,rr7,t3rr7)
               call t4matrixrr5(xr,yr,zr,rr5,t4rr5)
               call t4matrixrr7(xr,yr,zr,rr7,t4rr7)
               call t4matrixrr9(xr,yr,zr,rr9,t4rr9)
c     
c     call routines that produce potential, field, field gradient
c     for types of damping
c     
               if (damp_none) then
                  t2 = t2rr3 + t2rr5
                  t3 = t3rr5 + t3rr7
                  t4 = t4rr5 + t4rr7 + t4rr9
                  call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                 fieldkp)
                  call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                 gradfieldip,gradfieldkp)
                  call uhessfieldik(i,k,t4,hessfieldid,hessfieldkd,
     &                 hessfieldip,hessfieldkp)
                  do j = 1, 3
                     udfieldo(j,i) = udfieldo(j,i) + fieldid(j)
                     udfieldo(j,k) = udfieldo(j,k) + fieldkd(j)
                     upfieldo(j,i) = upfieldo(j,i) + fieldip(j)
                     upfieldo(j,k) = upfieldo(j,k) + fieldkp(j)
                     do l = 1, 3
                        udgradfieldo(l,j,i) = udgradfieldo(l,j,i) + 
     &                       gradfieldid(l,j)
                        udgradfieldo(l,j,k) = udgradfieldo(l,j,k) + 
     &                       gradfieldkd(l,j)
                        upgradfieldo(l,j,i) = upgradfieldo(l,j,i) + 
     &                       gradfieldip(l,j)
                        upgradfieldo(l,j,k) = upgradfieldo(l,j,k) + 
     &                       gradfieldkp(l,j)
                        do h = 1, 3
                           udhessfieldo(h,l,j,i) = 
     &                          udhessfieldo(h,l,j,i) +
     &                          hessfieldid(h,l,j)
                           udhessfieldo(h,l,j,k) = 
     &                          udhessfieldo(h,l,j,k) +
     &                          hessfieldkd(h,l,j)
                           uphessfieldo(h,l,j,i) = 
     &                          uphessfieldo(h,l,j,i) +
     &                          hessfieldip(h,l,j)
                           uphessfieldo(h,l,j,k) = 
     &                          uphessfieldo(h,l,j,k) +
     &                          hessfieldkp(h,l,j)
                        end do
                     end do
                  end do
               end if
c     
c     error function damping for ewald
c     
               if (damp_ewald) then
                  call dampewald(i,k,rorder,r,r2,scale)
c     
c     the ewald damping factors already contain their powers of r (rrx)
c     
                  t2 = t2rr3*scale(3)/rr3 + t2rr5*scale(5)/rr5
                  t3 = t3rr5*scale(5)/rr5 + t3rr7*scale(7)/rr7
                  t4 = t4rr5*scale(5)/rr5 + t4rr7*scale(7)/rr7 +
     &                 t4rr9*scale(9)/rr9
                  call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                 fieldkp)
                  call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                 gradfieldip,gradfieldkp)
                  call uhessfieldik(i,k,t4,hessfieldid,hessfieldkd,
     &                 hessfieldip,hessfieldkp)
                  do j = 1, 3
                     udfield_ewaldo(j,i) = udfield_ewaldo(j,i) +
     &                    fieldid(j)
                     udfield_ewaldo(j,k) = udfield_ewaldo(j,k) +
     &                    fieldkd(j)
                     upfield_ewaldo(j,i) = upfield_ewaldo(j,i) +
     &                    fieldip(j)
                     upfield_ewaldo(j,k) = upfield_ewaldo(j,k) +
     &                    fieldkp(j)
                     do l = 1, 3
                        udgradfield_ewaldo(l,j,i) = 
     &                       udgradfield_ewaldo(l,j,i) + 
     &                       gradfieldid(l,j)
                        udgradfield_ewaldo(l,j,k) = 
     &                       udgradfield_ewaldo(l,j,k) + 
     &                       gradfieldkd(l,j)
                        upgradfield_ewaldo(l,j,i) = 
     &                       upgradfield_ewaldo(l,j,i) +
     &                       gradfieldip(l,j)
                        upgradfield_ewaldo(l,j,k) = 
     &                       upgradfield_ewaldo(l,j,k) +
     &                       gradfieldkp(l,j)
                        do h = 1, 3
                           udhessfield_ewaldo(h,l,j,i) =
     &                          udhessfield_ewaldo(h,l,j,i) +
     &                          hessfieldid(h,l,j)
                           udhessfield_ewaldo(h,l,j,k) =
     &                          udhessfield_ewaldo(h,l,j,k) +
     &                          hessfieldkd(h,l,j)
                           uphessfield_ewaldo(h,l,j,i) =
     &                          uphessfield_ewaldo(h,l,j,i) +
     &                          hessfieldip(h,l,j)
                           uphessfield_ewaldo(h,l,j,k) =
     &                          uphessfield_ewaldo(h,l,j,k) +
     &                          hessfieldkp(h,l,j)
                        end do
                     end do
                  end do
               end if
c
c     thole damping
c
               if (damp_thole) then
                  call dampthole(i,k,rorder,r,scale)
                  t2 = t2rr3*scale(3) + t2rr5*scale(5)
                  t3 = t3rr5*scale(5) + t3rr7*scale(7)
                  t4 = t4rr5*scale(5) + t4rr7*scale(7) +
     &                 t4rr9*scale(9)
                  call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                 fieldkp)
                  call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                 gradfieldip,gradfieldkp)
                  call uhessfieldik(i,k,t4,hessfieldid,hessfieldkd,
     &                 hessfieldip,hessfieldkp)
                  do j = 1, 3
                     udfield_tholeo(j,i) = udfield_tholeo(j,i) +
     &                    fieldid(j)
                     udfield_tholeo(j,k) = udfield_tholeo(j,k) +
     &                    fieldkd(j)
                     upfield_tholeo(j,i) = upfield_tholeo(j,i) +
     &                    fieldip(j)
                     upfield_tholeo(j,k) = upfield_tholeo(j,k) +
     &                    fieldkp(j)
                     udfieldd_tholeo(j,i) = udfieldd_tholeo(j,i) +
     &                    fieldid(j)*dscale(kk)
                     udfieldd_tholeo(j,k) = udfieldd_tholeo(j,k) +
     &                    fieldkd(j)*dscale(kk)
                     upfieldp_tholeo(j,i) = upfieldp_tholeo(j,i) +
     &                    fieldip(j)*pscale(kk)
                     upfieldp_tholeo(j,k) = upfieldp_tholeo(j,k) +
     &                    fieldkp(j)*pscale(kk)
                     do l = 1, 3
                        udgradfield_tholeo(l,j,i) = 
     &                       udgradfield_tholeo(l,j,i) +
     &                       gradfieldid(l,j)
                        udgradfield_tholeo(l,j,k) = 
     &                       udgradfield_tholeo(l,j,k) +
     &                       gradfieldkd(l,j)
                        upgradfield_tholeo(l,j,i) = 
     &                       upgradfield_tholeo(l,j,i) +
     &                       gradfieldip(l,j)
                        upgradfield_tholeo(l,j,k) = 
     &                       upgradfield_tholeo(l,j,k) +
     &                       gradfieldkp(l,j)
                        udgradfieldd_tholeo(l,j,i) =
     &                       udgradfieldd_tholeo(l,j,i) +
     &                       gradfieldid(l,j)*dscale(kk)
                        udgradfieldd_tholeo(l,j,k) =
     &                       udgradfieldd_tholeo(l,j,k) +
     &                       gradfieldkd(l,j)*dscale(kk)
                        upgradfieldp_tholeo(l,j,i) =
     &                       upgradfieldp_tholeo(l,j,i) +
     &                       gradfieldip(l,j)*pscale(kk)
                        upgradfieldp_tholeo(l,j,k) =
     &                       upgradfieldp_tholeo(l,j,k) +
     &                       gradfieldkp(l,j)*pscale(kk)
                        do h = 1, 3
                           udhessfield_tholeo(h,l,j,i) =
     &                          udhessfield_tholeo(h,l,j,i) +
     &                          hessfieldid(h,l,j)
                           udhessfield_tholeo(h,l,j,k) =
     &                          udhessfield_tholeo(h,l,j,k) +
     &                          hessfieldkd(h,l,j)
                           uphessfield_tholeo(h,l,j,i) =
     &                          uphessfield_tholeo(h,l,j,i) +
     &                          hessfieldip(h,l,j)
                           uphessfield_tholeo(h,l,j,k) =
     &                          uphessfield_tholeo(h,l,j,k) +
     &                          hessfieldkp(h,l,j)
                           udhessfieldd_tholeo(h,l,j,i) =
     &                          udhessfieldd_tholeo(h,l,j,i) +
     &                          hessfieldid(h,l,j)*dscale(kk)
                           udhessfieldd_tholeo(h,l,j,k) =
     &                          udhessfieldd_tholeo(h,l,j,k) +
     &                          hessfieldkd(h,l,j)*dscale(kk)
                           uphessfieldp_tholeo(h,l,j,i) =
     &                          uphessfieldp_tholeo(h,l,j,i) +
     &                          hessfieldip(h,l,j)*pscale(kk)
                           uphessfieldp_tholeo(h,l,j,k) =
     &                          uphessfieldp_tholeo(h,l,j,k) +
     &                          hessfieldkp(h,l,j)*pscale(kk)
                        end do
                     end do
                  end do
               end if
            end if
         end do
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c
      udfield = udfieldo
      upfield = upfieldo
      udgradfield = udgradfieldo
      upgradfield = upgradfieldo
      udhessfield = udhessfieldo
      uphessfield = uphessfieldo
c
      udfield_ewald = udfield_ewaldo
      upfield_ewald = upfield_ewaldo
      udgradfield_ewald = udgradfield_ewaldo
      upgradfield_ewald = upgradfield_ewaldo
      udhessfield_ewald = udhessfield_ewaldo
      uphessfield_ewald = uphessfield_ewaldo
c
      udfield_thole = udfield_tholeo
      upfield_thole = upfield_tholeo
      udgradfield_thole = udgradfield_tholeo
      upgradfield_thole = upgradfield_tholeo
      udhessfield_thole = udhessfield_tholeo
      uphessfield_thole = uphessfield_tholeo
c
      udfieldd_thole = udfieldd_tholeo
      upfieldp_thole = upfieldp_tholeo
      udgradfieldd_thole = udgradfieldd_tholeo
      upgradfieldp_thole = upgradfieldp_tholeo
      udhessfieldd_thole = udhessfieldd_tholeo
      uphessfieldp_thole = uphessfieldp_tholeo
c
c     perform deallocation of some local arrays
c
      deallocate (udfieldo)
      deallocate (upfieldo)
      deallocate (udgradfieldo)
      deallocate (upgradfieldo)
      deallocate (udhessfieldo)
      deallocate (uphessfieldo)
c
      deallocate (udfield_ewaldo)
      deallocate (upfield_ewaldo)
      deallocate (udgradfield_ewaldo)
      deallocate (upgradfield_ewaldo)
      deallocate (udhessfield_ewaldo)
      deallocate (uphessfield_ewaldo)
c
      deallocate (udfield_tholeo)
      deallocate (upfield_tholeo)
      deallocate (udgradfield_tholeo)
      deallocate (upgradfield_tholeo)
      deallocate (udhessfield_tholeo)
      deallocate (uphessfield_tholeo)
c
      deallocate (udfieldd_tholeo)
      deallocate (udgradfieldd_tholeo)
      deallocate (udhessfieldd_tholeo)
c
      deallocate (upfieldp_tholeo)
      deallocate (upgradfieldp_tholeo)
      deallocate (uphessfieldp_tholeo)
      return
      end
