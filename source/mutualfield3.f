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
      use chgpen
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
      if ((directdamp.eq."THOLE").or.(mutualdamp.eq."THOLE")) then
         if (.not.allocated(udfield_thole))
     &        allocate (udfield_thole(3,npole))
         if (.not.allocated(upfield_thole))
     &        allocate (upfield_thole(3,npole))
         if (.not.allocated(udgradfield_thole)) 
     &        allocate (udgradfield_thole(3,3,npole))
         if (.not.allocated(upgradfield_thole)) 
     &        allocate (upgradfield_thole(3,3,npole))
         if (.not.allocated(udhessfield_thole))
     &        allocate (udhessfield_thole(3,3,3,npole))
         if (.not.allocated(uphessfield_thole))
     &        allocate (uphessfield_thole(3,3,3,npole))
c
         if (.not.allocated(upfieldd_thole))
     &        allocate (upfieldd_thole(3,npole))
         if (.not.allocated(upgradfieldd_thole))
     &        allocate (upgradfieldd_thole(3,3,npole))
         if (.not.allocated(uphessfieldd_thole))
     &        allocate (uphessfieldd_thole(3,3,3,npole))
c
         if (.not.allocated(udfieldp_thole))
     &        allocate (udfieldp_thole(3,npole))
         if (.not.allocated(udgradfieldp_thole))
     &        allocate (udgradfieldp_thole(3,3,npole))
         if (.not.allocated(udhessfieldp_thole))
     &        allocate (udhessfieldp_thole(3,3,3,npole))
      end if
      if ((directdamp.eq."GORDON").or.(mutualdamp.eq."GORDON"))then
         if (.not.allocated(udfield_gordon))
     &        allocate (udfield_gordon(3,npole))
         if (.not.allocated(upfield_gordon))
     &        allocate (upfield_gordon(3,npole))
         if (.not.allocated(udgradfield_gordon))
     &        allocate (udgradfield_gordon(3,3,npole))
         if (.not.allocated(upgradfield_gordon))
     &        allocate (upgradfield_gordon(3,3,npole))
         if (.not.allocated(udhessfield_gordon))
     &        allocate (udhessfield_gordon(3,3,3,npole))
         if (.not.allocated(uphessfield_gordon))
     &        allocate (uphessfield_gordon(3,3,3,npole))
c
         if (.not.allocated(upnucfieldd_gordon))
     &        allocate (upnucfieldd_gordon(3,npole))
         if (.not.allocated(upfieldd_gordon))
     &        allocate (upfieldd_gordon(3,npole))
         if (.not.allocated(upgradfieldd_gordon))
     &        allocate (upgradfieldd_gordon(3,3,npole))
         if (.not.allocated(uphessfieldd_gordon))
     &        allocate (uphessfieldd_gordon(3,3,3,npole))
c
         if (.not.allocated(udnucfieldp_gordon))
     &        allocate (udnucfieldp_gordon(3,npole))
         if (.not.allocated(udfieldp_gordon))
     &        allocate (udfieldp_gordon(3,npole))
         if (.not.allocated(udgradfieldp_gordon))
     &        allocate (udgradfieldp_gordon(3,3,npole))
         if (.not.allocated(udhessfieldp_gordon))
     &        allocate (udhessfieldp_gordon(3,3,3,npole))
      end if
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
      use chgpen
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
      real*8 t2i(3,3),t3i(3,3,3),t4i(3,3,3,3)
      real*8 t2k(3,3),t3k(3,3,3),t4k(3,3,3,3)
      real*8 t2ik(3,3),t3ik(3,3,3),t4ik(3,3,3,3)
      real*8 t2rr3(3,3),t2rr5(3,3)
      real*8 t3rr5(3,3,3),t3rr7(3,3,3)
      real*8 t4rr5(3,3,3,3),t4rr7(3,3,3,3),t4rr9(3,3,3,3)
      real*8 fieldid(3),fieldkd(3)
      real*8 fieldip(3),fieldkp(3)
      real*8 nucfieldid(3),nucfieldkd(3)
      real*8 nucfieldip(3),nucfieldkp(3)
      real*8 elefieldid(3),elefieldkd(3)
      real*8 elefieldip(3),elefieldkp(3)
      real*8 gradfieldid(3,3),gradfieldkd(3,3)
      real*8 gradfieldip(3,3),gradfieldkp(3,3)
      real*8 hessfieldid(3,3,3),hessfieldkd(3,3,3)
      real*8 hessfieldip(3,3,3),hessfieldkp(3,3,3)
      real*8, allocatable :: scale(:)
      real*8, allocatable :: scalei(:)
      real*8, allocatable :: scalek(:)
      real*8, allocatable :: scaleik(:)
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
            if ((directdamp.eq."THOLE").or.(mutualdamp.eq."THOLE")) then
               udfield_thole(j,i) = 0.0d0
               upfield_thole(j,i) = 0.0d0
               upfieldd_thole(j,i) = 0.0d0
               udfieldp_thole(j,i) = 0.0d0
            end if
            if ((directdamp .eq. "GORDON").or.
     &              (mutualdamp .eq. "GORDON")) then
               udfield_gordon(j,i) = 0.0d0
               upfield_gordon(j,i) = 0.0d0
               upnucfieldd_gordon(j,i) = 0.0d0
               udnucfieldp_gordon(j,i) = 0.0d0
               upfieldd_gordon(j,i) = 0.0d0
               udfieldp_gordon(j,i) = 0.0d0
            end if
            do k = 1, 3
               udgradfield(k,j,i) = 0.0d0
               upgradfield(k,j,i) = 0.0d0
               udgradfield_ewald(k,j,i) = 0.0d0
               upgradfield_ewald(k,j,i) = 0.0d0
               if ((directdamp .eq. "THOLE").or.
     &              (mutualdamp .eq. "THOLE")) then
                  udgradfield_thole(k,j,i) = 0.0d0
                  upgradfield_thole(k,j,i) = 0.0d0
                  upgradfieldd_thole(k,j,i) = 0.0d0
                  udgradfieldp_thole(k,j,i) = 0.0d0
               end if
               if ((directdamp .eq. "GORDON").or.
     &                 (mutualdamp .eq. "GORDON")) then
                  udgradfield_gordon(k,j,i) = 0.0d0
                  upgradfield_gordon(k,j,i) = 0.0d0
                  upgradfieldd_gordon(k,j,i) = 0.0d0
                  udgradfieldp_gordon(k,j,i) = 0.0d0
               end if
               do l = 1, 3
                  udhessfield(l,k,j,i) = 0.0d0
                  uphessfield(l,k,j,i) = 0.0d0
                  udhessfield_ewald(l,k,j,i) = 0.0d0
                  uphessfield_ewald(l,k,j,i) = 0.0d0
                  if ((directdamp .eq. "THOLE").or.
     &                 (mutualdamp .eq. "THOLE")) then
                     udhessfield_thole(l,k,j,i) = 0.0d0
                     uphessfield_thole(l,k,j,i) = 0.0d0
                     uphessfieldd_thole(l,k,j,i) = 0.0d0
                     udhessfieldp_thole(l,k,j,i) = 0.0d0
                  end if
                  if ((directdamp.eq."GORDON").or.
     &                    (mutualdamp.eq."GORDON")) then
                     udhessfield_gordon(l,k,j,i) = 0.0d0
                     uphessfield_gordon(l,k,j,i) = 0.0d0
                     uphessfieldd_gordon(l,k,j,i) = 0.0d0
                     udhessfieldp_gordon(l,k,j,i) = 0.0d0
                  end if
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
      allocate (scalei(rorder))
      allocate (scalek(rorder))
      allocate (scaleik(rorder))
      do i = 1,rorder
          scale(i) = 0.0d0
          scalei(i) = 0.0d0
          scalek(i) = 0.0d0
          scaleik(i) = 0.0d0
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
                        udfieldp_thole(j,i) = udfieldp_thole(j,i) +
     &                       fieldid(j)*pscale(kk)
                        udfieldp_thole(j,k) = udfieldp_thole(j,k) +
     &                       fieldkd(j)*pscale(kk)
                        upfieldd_thole(j,i) = upfieldd_thole(j,i) +
     &                       fieldip(j)*dscale(kk)
                        upfieldd_thole(j,k) = upfieldd_thole(j,k) +
     &                       fieldkp(j)*dscale(kk)
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
                           udgradfieldp_thole(l,j,i) =
     &                          udgradfieldp_thole(l,j,i) +
     &                          gradfieldid(l,j)*pscale(kk)
                           udgradfieldp_thole(l,j,k) =
     &                          udgradfieldp_thole(l,j,k) +
     &                          gradfieldkd(l,j)*pscale(kk)
                           upgradfieldd_thole(l,j,i) =
     &                          upgradfieldd_thole(l,j,i) +
     &                          gradfieldip(l,j)*dscale(kk)
                           upgradfieldd_thole(l,j,k) =
     &                          upgradfieldd_thole(l,j,k) +
     &                          gradfieldkp(l,j)*dscale(kk)
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
                              udhessfieldp_thole(h,l,j,i) =
     &                             udhessfieldp_thole(h,l,j,i) +
     &                             hessfieldid(h,l,j)*pscale(kk)
                              udhessfieldp_thole(h,l,j,k) =
     &                             udhessfieldp_thole(h,l,j,k) +
     &                             hessfieldkd(h,l,j)*pscale(kk)
                              uphessfieldd_thole(h,l,j,i) =
     &                             uphessfieldd_thole(h,l,j,i) +
     &                             hessfieldip(h,l,j)*dscale(kk)
                              uphessfieldd_thole(h,l,j,k) =
     &                             uphessfieldd_thole(h,l,j,k) +
     &                             hessfieldkp(h,l,j)*dscale(kk)
                           end do
                        end do
                     end do
                  end if
c
c     gordon damping
c
                  if (damp_gordon) then
                     call dampgordon(i,k,rorder,r,scalei,scalek,scaleik)
                     t2i = t2rr3*scalei(3) + t2rr5*scalei(5)
                     t2k = t2rr3*scalek(3) + t2rr5*scalek(5)
                     t2ik = t2rr3*scaleik(3) + t2rr5*scaleik(5)
                     t3ik = t3rr5*scaleik(5) + t3rr7*scaleik(7)
                     t4ik = t4rr5*scaleik(5) + t4rr7*scaleik(7) +
     &                      t4rr9*scaleik(9)
c     need this for field at nuclei from induced dipoles
                     call cp_ufieldik(i,k,t2i,t2k,t2ik,
     &                    nucfieldid,nucfieldkd,elefieldid,elefieldkd,
     &                    nucfieldip,nucfieldkp,elefieldip,elefieldkp)
c                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
c     &                    fieldkp)
                     call ugradfieldik(i,k,t3ik,gradfieldid,gradfieldkd
     &                    ,gradfieldip,gradfieldkp)
                     call uhessfieldik(i,k,t4ik,hessfieldid,hessfieldkd
     &                    ,hessfieldip,hessfieldkp)
                     do j = 1, 3
                        udfield_gordon(j,i) = udfield_gordon(j,i) +
     &                       elefieldid(j)
                        udfield_gordon(j,k) = udfield_gordon(j,k) +
     &                       elefieldkd(j)
                        upfield_gordon(j,i) = upfield_gordon(j,i) +
     &                       elefieldip(j)
                        upfield_gordon(j,k) = upfield_gordon(j,k) +
     &                       elefieldkp(j)
                        udfieldp_gordon(j,i) = udfieldp_gordon(j,i) +
     &                       elefieldid(j)*pscale(kk)
                        udfieldp_gordon(j,k) = udfieldp_gordon(j,k) +
     &                       elefieldkd(j)*pscale(kk)
                        upfieldd_gordon(j,i) = upfieldd_gordon(j,i) +
     &                       elefieldip(j)*dscale(kk)
                        upfieldd_gordon(j,k) = upfieldd_gordon(j,k) +
     &                       elefieldkp(j)*dscale(kk)
c
                        udnucfieldp_gordon(j,i)=udnucfieldp_gordon(j,i)+
     &                       nucfieldid(j)*pscale(kk)
                        udnucfieldp_gordon(j,k)=udnucfieldp_gordon(j,k)+
     &                       nucfieldkd(j)*pscale(kk)
                        if (.not.damp_gordonreg) then
                           upnucfieldd_gordon(j,i) = 
     &                          upnucfieldd_gordon(j,i)+
     &                          nucfieldip(j)*dscale(kk)
                           upnucfieldd_gordon(j,k) = 
     &                          upnucfieldd_gordon(j,k)+
     &                          nucfieldkp(j)*dscale(kk)
                        end if
                        do l = 1, 3
                           udgradfield_gordon(l,j,i) = 
     &                          udgradfield_gordon(l,j,i) +
     &                          gradfieldid(l,j)
                           udgradfield_gordon(l,j,k) = 
     &                          udgradfield_gordon(l,j,k) +
     &                          gradfieldkd(l,j)
                           upgradfield_gordon(l,j,i) = 
     &                          upgradfield_gordon(l,j,i) +
     &                          gradfieldip(l,j)
                           upgradfield_gordon(l,j,k) = 
     &                          upgradfield_gordon(l,j,k) +
     &                          gradfieldkp(l,j)
                           udgradfieldp_gordon(l,j,i) =
     &                          udgradfieldp_gordon(l,j,i) +
     &                          gradfieldid(l,j)*pscale(kk)
                           udgradfieldp_gordon(l,j,k) =
     &                          udgradfieldp_gordon(l,j,k) +
     &                          gradfieldkd(l,j)*pscale(kk)
                           upgradfieldd_gordon(l,j,i) =
     &                          upgradfieldd_gordon(l,j,i) +
     &                          gradfieldip(l,j)*dscale(kk)
                           upgradfieldd_gordon(l,j,k) =
     &                          upgradfieldd_gordon(l,j,k) +
     &                          gradfieldkp(l,j)*dscale(kk)
                           do h = 1, 3
                              udhessfield_gordon(h,l,j,i) =
     &                             udhessfield_gordon(h,l,j,i) +
     &                             hessfieldid(h,l,j)
                              udhessfield_gordon(h,l,j,k) =
     &                             udhessfield_gordon(h,l,j,k) +
     &                             hessfieldkd(h,l,j)
                              uphessfield_gordon(h,l,j,i) =
     &                             uphessfield_gordon(h,l,j,i) +
     &                             hessfieldip(h,l,j)
                              uphessfield_gordon(h,l,j,k) =
     &                             uphessfield_gordon(h,l,j,k) +
     &                             hessfieldkp(h,l,j)
                              udhessfieldp_gordon(h,l,j,i) =
     &                             udhessfieldp_gordon(h,l,j,i) +
     &                             hessfieldid(h,l,j)*pscale(kk)
                              udhessfieldp_gordon(h,l,j,k) =
     &                             udhessfieldp_gordon(h,l,j,k) +
     &                             hessfieldkd(h,l,j)*pscale(kk)
                              uphessfieldd_gordon(h,l,j,i) =
     &                             uphessfieldd_gordon(h,l,j,i) +
     &                             hessfieldip(h,l,j)*dscale(kk)
                              uphessfieldd_gordon(h,l,j,k) =
     &                             uphessfieldd_gordon(h,l,j,k) +
     &                             hessfieldkp(h,l,j)*dscale(kk)
                           end do
                        end do
                     end do
                     if (damp_gordonreg) then
                        call dampgordonreg(i,k,rorder,r,scalei,scalek)
                        t2i = t2rr3*scalei(3) + t2rr5*scalei(5)
                        t2k = t2rr3*scalek(3) + t2rr5*scalek(5)
                        t2ik = t2rr3*scaleik(3) + t2rr5*scaleik(5)
                        call cp_ufieldik(i,k,t2i,t2k,t2ik,
     &                      nucfieldid,nucfieldkd,elefieldid,elefieldkd,
     &                      nucfieldip,nucfieldkp,elefieldip,elefieldkp)
                        do j = 1, 3
                           upnucfieldd_gordon(j,i) =
     &                          upnucfieldd_gordon(j,i)+
     &                          nucfieldip(j)*dscale(kk)
                           upnucfieldd_gordon(j,k) =
     &                          upnucfieldd_gordon(j,k)+
     &                          nucfieldkp(j)*dscale(kk)
                        end do
                     end if
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
      print *,"USING REPLICA 3"
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
                        udfieldp_thole(j,i) = udfieldp_thole(j,i) +
     &                       fieldid(j)*pscale(kk)
                        udfieldp_thole(j,k) = udfieldp_thole(j,k) +
     &                       fieldkd(j)*pscale(kk)
                        upfieldd_thole(j,i) = upfieldd_thole(j,i) +
     &                       fieldip(j)*dscale(kk)
                        upfieldd_thole(j,k) = upfieldd_thole(j,k) +
     &                       fieldkp(j)*dscale(kk)
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
                           udgradfieldp_thole(l,j,i) =
     &                          udgradfieldp_thole(l,j,i) +
     &                          gradfieldid(l,j)*pscale(kk)
                           udgradfieldp_thole(l,j,k) =
     &                          udgradfieldp_thole(l,j,k) +
     &                          gradfieldkd(l,j)*pscale(kk)
                           upgradfieldd_thole(l,j,i) =
     &                          upgradfieldd_thole(l,j,i) +
     &                          gradfieldip(l,j)*dscale(kk)
                           upgradfieldd_thole(l,j,k) =
     &                          upgradfieldd_thole(l,j,k) +
     &                          gradfieldkp(l,j)*dscale(kk)
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
                              udhessfieldp_thole(h,l,j,i) =
     &                             udhessfieldp_thole(h,l,j,i) +
     &                             hessfieldid(h,l,j)*pscale(kk)
                              udhessfieldp_thole(h,l,j,k) =
     &                             udhessfieldp_thole(h,l,j,k) +
     &                             hessfieldkd(h,l,j)*pscale(kk)
                              uphessfieldd_thole(h,l,j,i) =
     &                             uphessfieldd_thole(h,l,j,i) +
     &                             hessfieldip(h,l,j)*dscale(kk)
                              uphessfieldd_thole(h,l,j,k) =
     &                             uphessfieldd_thole(h,l,j,k) +
     &                             hessfieldkp(h,l,j)*dscale(kk)
                           end do
                        end do
                     end do
                  end if
c
c     gordon damping
c
                  if (damp_gordon) then
                     call dampgordon(i,k,rorder,r,scalei,scalek,scaleik)
                     t2i = t2rr3*scalei(3) + t2rr5*scalei(5)
                     t2k = t2rr3*scalek(3) + t2rr5*scalek(5)
                     t2ik = t2rr3*scaleik(3) + t2rr5*scaleik(5)
                     t3ik = t3rr5*scaleik(5) + t3rr7*scaleik(7)
                     t4ik = t4rr5*scaleik(5) + t4rr7*scaleik(7) +
     &                      t4rr9*scaleik(9)
c     need this for field at nuclei from induced dipoles
                     call cp_ufieldik(i,k,t2i,t2k,t2ik,
     &                    nucfieldid,nucfieldkd,elefieldid,elefieldkd,
     &                    nucfieldip,nucfieldkp,elefieldip,elefieldkp)
c                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
c     &                    fieldkp)
                     call ugradfieldik(i,k,t3ik,gradfieldid,gradfieldkd
     &                    ,gradfieldip,gradfieldkp)
                     call uhessfieldik(i,k,t4ik,hessfieldid,hessfieldkd
     &                    ,hessfieldip,hessfieldkp)
                     do j = 1, 3
                        udfield_gordon(j,i) = udfield_gordon(j,i) +
     &                       elefieldid(j)
                        udfield_gordon(j,k) = udfield_gordon(j,k) +
     &                       elefieldkd(j)
                        upfield_gordon(j,i) = upfield_gordon(j,i) +
     &                       elefieldip(j)
                        upfield_gordon(j,k) = upfield_gordon(j,k) +
     &                       elefieldkp(j)
                        udfieldp_gordon(j,i) = udfieldp_gordon(j,i) +
     &                       elefieldid(j)*pscale(kk)
                        udfieldp_gordon(j,k) = udfieldp_gordon(j,k) +
     &                       elefieldkd(j)*pscale(kk)
                        upfieldd_gordon(j,i) = upfieldd_gordon(j,i) +
     &                       elefieldip(j)*dscale(kk)
                        upfieldd_gordon(j,k) = upfieldd_gordon(j,k) +
     &                       elefieldkp(j)*dscale(kk)
c
                        udnucfieldp_gordon(j,i)=udnucfieldp_gordon(j,i)+
     &                       nucfieldid(j)*pscale(kk)
                        udnucfieldp_gordon(j,k)=udnucfieldp_gordon(j,k)+
     &                       nucfieldkd(j)*pscale(kk)
                        if (.not.damp_gordonreg) then
                           upnucfieldd_gordon(j,i) = 
     &                          upnucfieldd_gordon(j,i)+
     &                          nucfieldip(j)*dscale(kk)
                           upnucfieldd_gordon(j,k) = 
     &                          upnucfieldd_gordon(j,k)+
     &                          nucfieldkp(j)*dscale(kk)
                        end if
                        do l = 1, 3
                           udgradfield_gordon(l,j,i) = 
     &                          udgradfield_gordon(l,j,i) +
     &                          gradfieldid(l,j)
                           udgradfield_gordon(l,j,k) = 
     &                          udgradfield_gordon(l,j,k) +
     &                          gradfieldkd(l,j)
                           upgradfield_gordon(l,j,i) = 
     &                          upgradfield_gordon(l,j,i) +
     &                          gradfieldip(l,j)
                           upgradfield_gordon(l,j,k) = 
     &                          upgradfield_gordon(l,j,k) +
     &                          gradfieldkp(l,j)
                           udgradfieldp_gordon(l,j,i) =
     &                          udgradfieldp_gordon(l,j,i) +
     &                          gradfieldid(l,j)*pscale(kk)
                           udgradfieldp_gordon(l,j,k) =
     &                          udgradfieldp_gordon(l,j,k) +
     &                          gradfieldkd(l,j)*pscale(kk)
                           upgradfieldd_gordon(l,j,i) =
     &                          upgradfieldd_gordon(l,j,i) +
     &                          gradfieldip(l,j)*dscale(kk)
                           upgradfieldd_gordon(l,j,k) =
     &                          upgradfieldd_gordon(l,j,k) +
     &                          gradfieldkp(l,j)*dscale(kk)
                           do h = 1, 3
                              udhessfield_gordon(h,l,j,i) =
     &                             udhessfield_gordon(h,l,j,i) +
     &                             hessfieldid(h,l,j)
                              udhessfield_gordon(h,l,j,k) =
     &                             udhessfield_gordon(h,l,j,k) +
     &                             hessfieldkd(h,l,j)
                              uphessfield_gordon(h,l,j,i) =
     &                             uphessfield_gordon(h,l,j,i) +
     &                             hessfieldip(h,l,j)
                              uphessfield_gordon(h,l,j,k) =
     &                             uphessfield_gordon(h,l,j,k) +
     &                             hessfieldkp(h,l,j)
                              udhessfieldp_gordon(h,l,j,i) =
     &                             udhessfieldp_gordon(h,l,j,i) +
     &                             hessfieldid(h,l,j)*pscale(kk)
                              udhessfieldp_gordon(h,l,j,k) =
     &                             udhessfieldp_gordon(h,l,j,k) +
     &                             hessfieldkd(h,l,j)*pscale(kk)
                              uphessfieldd_gordon(h,l,j,i) =
     &                             uphessfieldd_gordon(h,l,j,i) +
     &                             hessfieldip(h,l,j)*dscale(kk)
                              uphessfieldd_gordon(h,l,j,k) =
     &                             uphessfieldd_gordon(h,l,j,k) +
     &                             hessfieldkp(h,l,j)*dscale(kk)
                           end do
                        end do
                     end do
                     if (damp_gordonreg) then
                        call dampgordonreg(i,k,rorder,r,scalei,scalek)
                        t2i = t2rr3*scalei(3) + t2rr5*scalei(5)
                        t2k = t2rr3*scalek(3) + t2rr5*scalek(5)
                        t2ik = t2rr3*scaleik(3) + t2rr5*scaleik(5)
                        call cp_ufieldik(i,k,t2i,t2k,t2ik,
     &                      nucfieldid,nucfieldkd,elefieldid,elefieldkd,
     &                      nucfieldip,nucfieldkp,elefieldip,elefieldkp)
                        do j = 1, 3
                           upnucfieldd_gordon(j,i) =
     &                          upnucfieldd_gordon(j,i)+
     &                          nucfieldip(j)*dscale(kk)
                           upnucfieldd_gordon(j,k) =
     &                          upnucfieldd_gordon(j,k)+
     &                          nucfieldkp(j)*dscale(kk)
                        end do
                     end if
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
      use chgpen
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
      real*8 t2i(3,3),t3i(3,3,3),t4i(3,3,3,3)
      real*8 t2k(3,3),t3k(3,3,3),t4k(3,3,3,3)
      real*8 t2ik(3,3),t3ik(3,3,3),t4ik(3,3,3,3)
      real*8 t2rr3(3,3),t2rr5(3,3)
      real*8 t3rr5(3,3,3),t3rr7(3,3,3)
      real*8 t4rr5(3,3,3,3),t4rr7(3,3,3,3),t4rr9(3,3,3,3)
      real*8 fieldid(3),fieldkd(3)
      real*8 fieldip(3),fieldkp(3)
      real*8 gradfieldid(3,3),gradfieldkd(3,3)
      real*8 gradfieldip(3,3),gradfieldkp(3,3)
      real*8 hessfieldid(3,3,3),hessfieldkd(3,3,3)
      real*8 hessfieldip(3,3,3),hessfieldkp(3,3,3)
      real*8 nucfieldid(3),nucfieldkd(3)
      real*8 nucfieldip(3),nucfieldkp(3)
      real*8 elefieldid(3),elefieldkd(3)
      real*8 elefieldip(3),elefieldkp(3)
      real*8, allocatable :: scale(:)
      real*8, allocatable :: scalei(:)
      real*8, allocatable :: scalek(:)
      real*8, allocatable :: scaleik(:)
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
      real*8, allocatable :: upfieldd_tholeo(:,:)
      real*8, allocatable :: upgradfieldd_tholeo(:,:,:)
      real*8, allocatable :: uphessfieldd_tholeo(:,:,:,:)
      real*8, allocatable :: udfieldp_tholeo(:,:)
      real*8, allocatable :: udgradfieldp_tholeo(:,:,:)
      real*8, allocatable :: udhessfieldp_tholeo(:,:,:,:)
      real*8, allocatable :: udfield_gordono(:,:)
      real*8, allocatable :: upfield_gordono(:,:)
      real*8, allocatable :: udgradfield_gordono(:,:,:)
      real*8, allocatable :: upgradfield_gordono(:,:,:)
      real*8, allocatable :: udhessfield_gordono(:,:,:,:)
      real*8, allocatable :: uphessfield_gordono(:,:,:,:)
      real*8, allocatable :: upnucfieldd_gordono(:,:)
      real*8, allocatable :: upfieldd_gordono(:,:)
      real*8, allocatable :: upgradfieldd_gordono(:,:,:)
      real*8, allocatable :: uphessfieldd_gordono(:,:,:,:)
      real*8, allocatable :: udnucfieldp_gordono(:,:)
      real*8, allocatable :: udfieldp_gordono(:,:)
      real*8, allocatable :: udgradfieldp_gordono(:,:,:)
      real*8, allocatable :: udhessfieldp_gordono(:,:,:,:)
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
      if ((directdamp.eq."THOLE").or.(mutualdamp.eq."THOLE")) then
         allocate (udfield_tholeo(3,npole))
         allocate (upfield_tholeo(3,npole))
         allocate (udgradfield_tholeo(3,3,npole))
         allocate (upgradfield_tholeo(3,3,npole))
         allocate (udhessfield_tholeo(3,3,3,npole))
         allocate (uphessfield_tholeo(3,3,3,npole))
c
         allocate (upfieldd_tholeo(3,npole))
         allocate (upgradfieldd_tholeo(3,3,npole))
         allocate (uphessfieldd_tholeo(3,3,3,npole))
c
         allocate (udfieldp_tholeo(3,npole))
         allocate (udgradfieldp_tholeo(3,3,npole))
         allocate (udhessfieldp_tholeo(3,3,3,npole))
      end if
      if ((directdamp.eq."GORDON").or.(mutualdamp.eq."GORDON"))then
         allocate (udfield_gordono(3,npole))
         allocate (upfield_gordono(3,npole))
         allocate (udgradfield_gordono(3,3,npole))
         allocate (upgradfield_gordono(3,3,npole))
         allocate (udhessfield_gordono(3,3,3,npole))
         allocate (uphessfield_gordono(3,3,3,npole))
c
         allocate (upnucfieldd_gordono(3,npole))
         allocate (upfieldd_gordono(3,npole))
         allocate (upgradfieldd_gordono(3,3,npole))
         allocate (uphessfieldd_gordono(3,3,3,npole))
c
         allocate (udnucfieldp_gordono(3,npole))
         allocate (udfieldp_gordono(3,npole))
         allocate (udgradfieldp_gordono(3,3,npole))
         allocate (udhessfieldp_gordono(3,3,3,npole))
      end if
c
c     zero out the value of the gradfield at each site
c
      do i = 1, npole
         do j = 1, 3
            udfield(j,i) = 0.0d0
            upfield(j,i) = 0.0d0
            udfield_ewald(j,i) = 0.0d0
            upfield_ewald(j,i) = 0.0d0
c
            udfieldo(j,i) = 0.0d0
            upfieldo(j,i) = 0.0d0
            udfield_ewaldo(j,i) = 0.0d0
            upfield_ewaldo(j,i) = 0.0d0
            if ((directdamp.eq."THOLE").or.(mutualdamp.eq."THOLE")) then
               udfield_thole(j,i) = 0.0d0
               upfield_thole(j,i) = 0.0d0
               upfieldd_thole(j,i) = 0.0d0
               udfieldp_thole(j,i) = 0.0d0
c
               udfield_tholeo(j,i) = 0.0d0
               upfield_tholeo(j,i) = 0.0d0
               upfieldd_tholeo(j,i) = 0.0d0
               udfieldp_tholeo(j,i) = 0.0d0
            end if
            if ((directdamp .eq. "GORDON").or.
     &              (mutualdamp .eq. "GORDON")) then
               udfield_gordon(j,i) = 0.0d0
               upfield_gordon(j,i) = 0.0d0
               upnucfieldd_gordon(j,i) = 0.0d0
               udnucfieldp_gordon(j,i) = 0.0d0
               upfieldd_gordon(j,i) = 0.0d0
               udfieldp_gordon(j,i) = 0.0d0
c
               udfield_gordono(j,i) = 0.0d0
               upfield_gordono(j,i) = 0.0d0
               upnucfieldd_gordono(j,i) = 0.0d0
               udnucfieldp_gordono(j,i) = 0.0d0
               upfieldd_gordono(j,i) = 0.0d0
               udfieldp_gordono(j,i) = 0.0d0
            end if
            do k = 1, 3
               udgradfield(k,j,i) = 0.0d0
               upgradfield(k,j,i) = 0.0d0
               udgradfield_ewald(k,j,i) = 0.0d0
               upgradfield_ewald(k,j,i) = 0.0d0
c
               udgradfieldo(k,j,i) = 0.0d0
               upgradfieldo(k,j,i) = 0.0d0
               udgradfield_ewaldo(k,j,i) = 0.0d0
               upgradfield_ewaldo(k,j,i) = 0.0d0
               if ((directdamp .eq. "THOLE").or.
     &              (mutualdamp .eq. "THOLE")) then
                  udgradfield_thole(k,j,i) = 0.0d0
                  upgradfield_thole(k,j,i) = 0.0d0
                  upgradfieldd_thole(k,j,i) = 0.0d0
                  udgradfieldp_thole(k,j,i) = 0.0d0
c
                  udgradfield_tholeo(k,j,i) = 0.0d0
                  upgradfield_tholeo(k,j,i) = 0.0d0
                  upgradfieldd_tholeo(k,j,i) = 0.0d0
                  udgradfieldp_tholeo(k,j,i) = 0.0d0
               end if
               if ((directdamp .eq. "GORDON").or.
     &                 (mutualdamp .eq. "GORDON")) then
                  udgradfield_gordon(k,j,i) = 0.0d0
                  upgradfield_gordon(k,j,i) = 0.0d0
                  upgradfieldd_gordon(k,j,i) = 0.0d0
                  udgradfieldp_gordon(k,j,i) = 0.0d0
c
                  udgradfield_gordono(k,j,i) = 0.0d0
                  upgradfield_gordono(k,j,i) = 0.0d0
                  upgradfieldd_gordono(k,j,i) = 0.0d0
                  udgradfieldp_gordono(k,j,i) = 0.0d0
               end if
               do l = 1, 3
                  udhessfield(l,k,j,i) = 0.0d0
                  uphessfield(l,k,j,i) = 0.0d0
                  udhessfield_ewald(l,k,j,i) = 0.0d0
                  uphessfield_ewald(l,k,j,i) = 0.0d0
c
                  udhessfieldo(l,k,j,i) = 0.0d0
                  uphessfieldo(l,k,j,i) = 0.0d0
                  udhessfield_ewaldo(l,k,j,i) = 0.0d0
                  uphessfield_ewaldo(l,k,j,i) = 0.0d0
                  if ((directdamp .eq. "THOLE").or.
     &                 (mutualdamp .eq. "THOLE")) then
                     udhessfield_thole(l,k,j,i) = 0.0d0
                     uphessfield_thole(l,k,j,i) = 0.0d0
                     uphessfieldd_thole(l,k,j,i) = 0.0d0
                     udhessfieldp_thole(l,k,j,i) = 0.0d0
c
                     udhessfield_tholeo(l,k,j,i) = 0.0d0
                     uphessfield_tholeo(l,k,j,i) = 0.0d0
                     uphessfieldd_tholeo(l,k,j,i) = 0.0d0
                     udhessfieldp_tholeo(l,k,j,i) = 0.0d0
                  end if
                  if ((directdamp.eq."GORDON").or.
     &                    (mutualdamp.eq."GORDON")) then
                     udhessfield_gordon(l,k,j,i) = 0.0d0
                     uphessfield_gordon(l,k,j,i) = 0.0d0
                     uphessfieldd_gordon(l,k,j,i) = 0.0d0
                     udhessfieldp_gordon(l,k,j,i) = 0.0d0
c
                     udhessfield_gordono(l,k,j,i) = 0.0d0
                     uphessfield_gordono(l,k,j,i) = 0.0d0
                     uphessfieldd_gordono(l,k,j,i) = 0.0d0
                     udhessfieldp_gordono(l,k,j,i) = 0.0d0
                  end if
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
      allocate (scalei(rorder))
      allocate (scalek(rorder))
      allocate (scaleik(rorder))
      do i = 1,rorder
          scale(i) = 0.0d0
          scalei(i) = 0.0d0
          scalek(i) = 0.0d0
          scaleik(i) = 0.0d0
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
!$OMP& nucfieldid,nucfieldkd,nucfieldip,nucfieldkp,
!$OMP& elefieldid,elefieldkd,elefieldip,elefieldkp, 
!$OMP& t2rr3,t2rr5,t3rr5,t3rr7,t4rr5,t4rr7,t4rr9,
!$OMP& t2,t3,t4,t2i,t3i,t4i,t2k,t3k,t4k,t2ik,t3ik,t4ik,
!$OMP& scale,scalei,scalek,scaleik)
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
!$OMP& udfieldp_tholeo,upfieldd_tholeo,
!$OMP& udgradfieldp_tholeo,upgradfieldd_tholeo,
!$OMP& udhessfieldp_tholeo,uphessfieldd_tholeo,
!$OMP& udfield_gordono,upfield_gordono, 
!$OMP& udgradfield_gordono,upgradfield_gordono,
!$OMP& udhessfield_gordono,uphessfield_gordono,
!$OMP& udnucfieldp_gordono,upnucfieldd_gordono,
!$OMP& udfieldp_gordono,upfieldd_gordono,
!$OMP& udgradfieldp_gordono,upgradfieldd_gordono,
!$OMP& udhessfieldp_gordono,uphessfieldd_gordono)
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
                     udfieldp_tholeo(j,i) = udfieldp_tholeo(j,i) +
     &                    fieldid(j)*pscale(kk)
                     udfieldp_tholeo(j,k) = udfieldp_tholeo(j,k) +
     &                    fieldkd(j)*pscale(kk)
                     upfieldd_tholeo(j,i) = upfieldd_tholeo(j,i) +
     &                    fieldip(j)*dscale(kk)
                     upfieldd_tholeo(j,k) = upfieldd_tholeo(j,k) +
     &                    fieldkp(j)*dscale(kk)
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
                        udgradfieldp_tholeo(l,j,i) =
     &                       udgradfieldp_tholeo(l,j,i) +
     &                       gradfieldid(l,j)*pscale(kk)
                        udgradfieldp_tholeo(l,j,k) =
     &                       udgradfieldp_tholeo(l,j,k) +
     &                       gradfieldkd(l,j)*pscale(kk)
                        upgradfieldd_tholeo(l,j,i) =
     &                       upgradfieldd_tholeo(l,j,i) +
     &                       gradfieldip(l,j)*dscale(kk)
                        upgradfieldd_tholeo(l,j,k) =
     &                       upgradfieldd_tholeo(l,j,k) +
     &                       gradfieldkp(l,j)*dscale(kk)
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
                           udhessfieldp_tholeo(h,l,j,i) =
     &                          udhessfieldp_tholeo(h,l,j,i) +
     &                          hessfieldid(h,l,j)*pscale(kk)
                           udhessfieldp_tholeo(h,l,j,k) =
     &                          udhessfieldp_tholeo(h,l,j,k) +
     &                          hessfieldkd(h,l,j)*pscale(kk)
                           uphessfieldd_tholeo(h,l,j,i) =
     &                          uphessfieldd_tholeo(h,l,j,i) +
     &                          hessfieldip(h,l,j)*dscale(kk)
                           uphessfieldd_tholeo(h,l,j,k) =
     &                          uphessfieldd_tholeo(h,l,j,k) +
     &                          hessfieldkp(h,l,j)*dscale(kk)
                        end do
                     end do
                  end do
               end if
c
c     gordon damping
c
               if (damp_gordon) then
                  call dampgordon(i,k,rorder,r,scalei,scalek,scaleik)
                  t2i = t2rr3*scalei(3) + t2rr5*scalei(5)
                  t2k = t2rr3*scalek(3) + t2rr5*scalek(5)
                  t2ik = t2rr3*scaleik(3) + t2rr5*scaleik(5)
                  t3ik = t3rr5*scaleik(5) + t3rr7*scaleik(7)
                  t4ik = t4rr5*scaleik(5) + t4rr7*scaleik(7) +
     &                 t4rr9*scaleik(9)
c     need this for field at nuclei from induced dipoles
                  call cp_ufieldik(i,k,t2i,t2k,t2ik,
     &                 nucfieldid,nucfieldkd,elefieldid,elefieldkd,
     &                 nucfieldip,nucfieldkp,elefieldip,elefieldkp)
c     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
c     &                    fieldkp)
                  call ugradfieldik(i,k,t3ik,gradfieldid,gradfieldkd
     &                 ,gradfieldip,gradfieldkp)
                  call uhessfieldik(i,k,t4ik,hessfieldid,hessfieldkd
     &                 ,hessfieldip,hessfieldkp)
                  do j = 1, 3
                     udfield_gordono(j,i) = udfield_gordono(j,i) +
     &                    elefieldid(j)
                     udfield_gordono(j,k) = udfield_gordono(j,k) +
     &                    elefieldkd(j)
                     upfield_gordono(j,i) = upfield_gordono(j,i) +
     &                    elefieldip(j)
                     upfield_gordono(j,k) = upfield_gordono(j,k) +
     &                    elefieldkp(j)
                     udfieldp_gordono(j,i) = udfieldp_gordono(j,i) +
     &                    elefieldid(j)*pscale(kk)
                     udfieldp_gordono(j,k) = udfieldp_gordono(j,k) +
     &                    elefieldkd(j)*pscale(kk)
                     upfieldd_gordono(j,i) = upfieldd_gordono(j,i) +
     &                    elefieldip(j)*dscale(kk)
                     upfieldd_gordono(j,k) = upfieldd_gordono(j,k) +
     &                    elefieldkp(j)*dscale(kk)
c     
                     udnucfieldp_gordono(j,i)=udnucfieldp_gordono(j,i)+
     &                    nucfieldid(j)*pscale(kk)
                     udnucfieldp_gordono(j,k)=udnucfieldp_gordono(j,k)+
     &                    nucfieldkd(j)*pscale(kk)
                     if (.not.damp_gordonreg) then
                        upnucfieldd_gordono(j,i) = 
     &                       upnucfieldd_gordono(j,i)+
     &                       nucfieldip(j)*dscale(kk)
                        upnucfieldd_gordono(j,k) = 
     &                       upnucfieldd_gordono(j,k)+
     &                       nucfieldkp(j)*dscale(kk)
                     end if
                     do l = 1, 3
                        udgradfield_gordono(l,j,i) = 
     &                       udgradfield_gordono(l,j,i) +
     &                       gradfieldid(l,j)
                        udgradfield_gordono(l,j,k) = 
     &                       udgradfield_gordono(l,j,k) +
     &                       gradfieldkd(l,j)
                        upgradfield_gordono(l,j,i) = 
     &                       upgradfield_gordono(l,j,i) +
     &                       gradfieldip(l,j)
                        upgradfield_gordono(l,j,k) = 
     &                       upgradfield_gordono(l,j,k) +
     &                       gradfieldkp(l,j)
                        udgradfieldp_gordono(l,j,i) =
     &                       udgradfieldp_gordono(l,j,i) +
     &                       gradfieldid(l,j)*pscale(kk)
                        udgradfieldp_gordono(l,j,k) =
     &                       udgradfieldp_gordono(l,j,k) +
     &                       gradfieldkd(l,j)*pscale(kk)
                        upgradfieldd_gordono(l,j,i) =
     &                       upgradfieldd_gordono(l,j,i) +
     &                       gradfieldip(l,j)*dscale(kk)
                        upgradfieldd_gordono(l,j,k) =
     &                       upgradfieldd_gordono(l,j,k) +
     &                       gradfieldkp(l,j)*dscale(kk)
                        do h = 1, 3
                           udhessfield_gordono(h,l,j,i) =
     &                          udhessfield_gordono(h,l,j,i) +
     &                          hessfieldid(h,l,j)
                           udhessfield_gordono(h,l,j,k) =
     &                          udhessfield_gordono(h,l,j,k) +
     &                          hessfieldkd(h,l,j)
                           uphessfield_gordono(h,l,j,i) =
     &                          uphessfield_gordono(h,l,j,i) +
     &                          hessfieldip(h,l,j)
                           uphessfield_gordono(h,l,j,k) =
     &                          uphessfield_gordono(h,l,j,k) +
     &                          hessfieldkp(h,l,j)
                           udhessfieldp_gordono(h,l,j,i) =
     &                          udhessfieldp_gordono(h,l,j,i) +
     &                          hessfieldid(h,l,j)*pscale(kk)
                           udhessfieldp_gordono(h,l,j,k) =
     &                          udhessfieldp_gordono(h,l,j,k) +
     &                          hessfieldkd(h,l,j)*pscale(kk)
                           uphessfieldd_gordono(h,l,j,i) =
     &                          uphessfieldd_gordono(h,l,j,i) +
     &                          hessfieldip(h,l,j)*dscale(kk)
                           uphessfieldd_gordono(h,l,j,k) =
     &                          uphessfieldd_gordono(h,l,j,k) +
     &                          hessfieldkp(h,l,j)*dscale(kk)
                        end do
                     end do
                  end do
                  if (damp_gordonreg) then
                     call dampgordonreg(i,k,rorder,r,scalei,scalek)
                     t2i = t2rr3*scalei(3) + t2rr5*scalei(5)
                     t2k = t2rr3*scalek(3) + t2rr5*scalek(5)
                     t2ik = t2rr3*scaleik(3) + t2rr5*scaleik(5)
                     call cp_ufieldik(i,k,t2i,t2k,t2ik,
     &                    nucfieldid,nucfieldkd,elefieldid,elefieldkd,
     &                    nucfieldip,nucfieldkp,elefieldip,elefieldkp)
                     do j = 1, 3
                        upnucfieldd_gordono(j,i) =
     &                       upnucfieldd_gordono(j,i)+
     &                       nucfieldip(j)*dscale(kk)
                        upnucfieldd_gordono(j,k) =
     &                       upnucfieldd_gordono(j,k)+
     &                       nucfieldkp(j)*dscale(kk)
                     end do
                  end if
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
      udfieldp_thole = udfieldp_tholeo
      upfieldd_thole = upfieldd_tholeo
      udgradfieldp_thole = udgradfieldp_tholeo
      upgradfieldd_thole = upgradfieldd_tholeo
      udhessfieldp_thole = udhessfieldp_tholeo
      uphessfieldd_thole = uphessfieldd_tholeo
c
      udfield_gordon = udfield_gordono
      upfield_gordon = upfield_gordono
      udgradfield_gordon = udgradfield_gordono
      upgradfield_gordon = upgradfield_gordono
      udhessfield_gordon = udhessfield_gordono
      uphessfield_gordon = uphessfield_gordono
c
      udnucfieldp_gordon = udnucfieldp_gordono
      upnucfieldd_gordon = upnucfieldd_gordono
      udfieldp_gordon = udfieldp_gordono
      upfieldd_gordon = upfieldd_gordono
      udgradfieldp_gordon = udgradfieldp_gordono
      upgradfieldd_gordon = upgradfieldd_gordono
      udhessfieldp_gordon = udhessfieldp_gordono
      uphessfieldd_gordon = uphessfieldd_gordono
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
      if ((directdamp.eq."THOLE").or.(mutualdamp.eq."THOLE")) then
         deallocate (udfield_tholeo)
         deallocate (upfield_tholeo)
         deallocate (udgradfield_tholeo)
         deallocate (upgradfield_tholeo)
         deallocate (udhessfield_tholeo)
         deallocate (uphessfield_tholeo)
c
         deallocate (upfieldd_tholeo)
         deallocate (upgradfieldd_tholeo)
         deallocate (uphessfieldd_tholeo)
c
         deallocate (udfieldp_tholeo)
         deallocate (udgradfieldp_tholeo)
         deallocate (udhessfieldp_tholeo)
      end if
      if ((directdamp.eq."GORDON").or.(mutualdamp.eq."GORDON"))then
         deallocate (udfield_gordono)
         deallocate (upfield_gordono)
         deallocate (udgradfield_gordono)
         deallocate (upgradfield_gordono)
         deallocate (udhessfield_gordono)
         deallocate (uphessfield_gordono)
c
         deallocate (upnucfieldd_gordono)
         deallocate (upfieldd_gordono)
         deallocate (upgradfieldd_gordono)
         deallocate (uphessfieldd_gordono)
c
         deallocate (udnucfieldp_gordono)
         deallocate (udfieldp_gordono)
         deallocate (udgradfieldp_gordono)
         deallocate (udhessfieldp_gordono)
      end if
      return
      end
