c
c
c     #############################################################
c     ## COPYRIGHT (C) 2016 by Josh Rackers & Jay William Ponder ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine permfield3  --  electric potential, field,   ##
c     ##  field gradient and field hessian from permanent moments ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "permfield3" computes the electric potential, field,
c     field gradient and field hessian for use in energy and gradient
c     routines as well as for use in computing induced dipoles
c
c     (the 3 stands for up to field hessian)
c
c     assumes multipole components have already been rotated into
c     the global coordinate frame
c
c
      subroutine permfield3
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
c
c     perform dynamic allocation of global arrays
c
      if (.not.allocated(pot)) allocate (pot(npole))
      if (.not.allocated(field)) allocate (field(3,npole))
      if (.not.allocated(gradfield)) allocate (gradfield(3,3,npole))
      if (.not.allocated(hessfield)) allocate (hessfield(3,3,3,npole))
c
c     permanent - permanent exclusion rule
c
      if (.not.allocated(potm)) allocate (potm(npole))
      if (.not.allocated(fieldm)) allocate (fieldm(3,npole))
      if (.not.allocated(gradfieldm)) allocate (gradfieldm(3,3,npole))
      if (.not.allocated(hessfieldm)) allocate (hessfieldm(3,3,3,npole))
c
c     ewald damping
c
      if (.not.allocated(pot_ewald)) allocate (pot_ewald(npole))
      if (.not.allocated(field_ewald)) allocate (field_ewald(3,npole))
      if (.not.allocated(gradfield_ewald)) 
     &     allocate (gradfield_ewald(3,3,npole))
      if (.not.allocated(hessfield_ewald))
     &     allocate (hessfield_ewald(3,3,3,npole))
c
c     thole damping with p and d exclusion rules needed for
c     computing induced dipoles
c
      if (.not.allocated(fieldd_thole)) allocate (fieldd_thole(3,npole))
      if (.not.allocated(fieldp_thole)) allocate (fieldp_thole(3,npole))
      if (.not.allocated(gradfieldd_thole))
     &     allocate (gradfieldd_thole(3,3,npole))
      if (.not.allocated(gradfieldp_thole))
     &     allocate (gradfieldp_thole(3,3,npole))
c
c     gordon charge penetration damping
c
      if (.not.allocated(potm_gordon)) allocate (potm_gordon(npole))
      if (.not.allocated(fieldm_gordon)) 
     &     allocate (fieldm_gordon(3,npole))
      if (.not.allocated(gradfieldm_gordon))
     &     allocate (gradfieldm_gordon(3,3,npole))
      if (.not.allocated(hessfieldm_gordon))
     &     allocate (hessfieldm_gordon(3,3,3,npole))
c
c     gordon charge penetration damping - nuclei
c
      if (.not.allocated(nucpotm_gordon))
     &     allocate (nucpotm_gordon(npole))
      if (.not.allocated(nucfieldm_gordon))
     &     allocate (nucfieldm_gordon(3,npole))
c
c     gordon damping for polarization
c
      if (.not.allocated(fieldd_gordon))
     &     allocate (fieldd_gordon(3,npole))
      if (.not.allocated(fieldp_gordon))
     &     allocate (fieldp_gordon(3,npole))
      if (.not.allocated(gradfieldd_gordon))
     &     allocate (gradfieldd_gordon(3,3,npole))
      if (.not.allocated(gradfieldp_gordon))
     &     allocate (gradfieldp_gordon(3,3,npole))
c
c     gordon damping with nuclear regularization
c
      if (.not.allocated(fieldd_gordonreg))
     &     allocate (fieldd_gordonreg(3,npole))
      if (.not.allocated(fieldp_gordonreg))
     &     allocate (fieldp_gordonreg(3,npole))
c
c     get the electrostatic potential, field and field gradient
c     due to permanent multipoles
c
      if (use_mlist) then
         call permfield3b
      else
         call permfield3a
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine permfield3a  --  potential, field, field    ##
c     ##           gradient and field hessian via double loop    ##
c     ##                                                         ##
c     #############################################################
c
c
c     "permfield3a" computes the direct electrostatic potential, field, 
c     field gradient and field hessian due to permanent multipole 
c     moments via a double loop
c
c
      subroutine permfield3a 
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
      real*8 rr1,rr3,rr5,rr7,rr9,rr11
      real*8 wall,cpu
      real*8 t0,t1(3),t2(3,3),t3(3,3,3),t4(3,3,3,3),t5(3,3,3,3,3)
      real*8 t0i,t1i(3),t2i(3,3),t3i(3,3,3),t4i(3,3,3,3),t5i(3,3,3,3,3)
      real*8 t0k,t1k(3),t2k(3,3),t3k(3,3,3),t4k(3,3,3,3),t5k(3,3,3,3,3)
      real*8 t0ik,t1ik(3),t2ik(3,3),t3ik(3,3,3),t4ik(3,3,3,3)
      real*8 t5ik(3,3,3,3,3)
      real*8 t0rr1,t1rr3(3)
      real*8 t2rr3(3,3),t2rr5(3,3)
      real*8 t3rr5(3,3,3),t3rr7(3,3,3)
      real*8 t4rr5(3,3,3,3),t4rr7(3,3,3,3),t4rr9(3,3,3,3)
      real*8 t5rr7(3,3,3,3,3),t5rr9(3,3,3,3,3),t5rr11(3,3,3,3,3)
      real*8 poti,potk
      real*8 fieldi(3),fieldk(3)
      real*8 gradfieldi(3,3),gradfieldk(3,3)
      real*8 hessfieldi(3,3,3),hessfieldk(3,3,3)
      real*8 nucpoti,nucpotk
      real*8 nucfieldi(3),nucfieldk(3)
      real*8 elepoti,elepotk
      real*8 elefieldi(3),elefieldk(3)
      real*8 elegradfieldi(3,3),elegradfieldk(3,3)
      real*8 elehessfieldi(3,3,3),elehessfieldk(3,3,3)
      real*8, allocatable :: scale(:)
      real*8, allocatable :: scalei(:)
      real*8, allocatable :: scalek(:)
      real*8, allocatable :: scaleik(:)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      logical proceed
      logical usei,usek
      character*6 mode
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         pot(i) = 0.0d0
         potm(i) = 0.0d0
         pot_ewald(i) = 0.0d0
         potm_gordon(i) = 0.0d0
         nucpotm_gordon(i) = 0.0d0
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldm(j,i) = 0.0d0
            field_ewald(j,i) = 0.0d0
            fieldd_thole(j,i) = 0.0d0
            fieldp_thole(j,i) = 0.0d0
            fieldm_gordon(j,i) = 0.0d0
            fieldd_gordon(j,i) = 0.0d0
            fieldp_gordon(j,i) = 0.0d0
            fieldd_gordonreg(j,i) = 0.0d0
            fieldp_gordonreg(j,i) = 0.0d0
            nucfieldm_gordon(j,i) = 0.0d0
            do k = 1, 3
               gradfield(k,j,i) = 0.0d0
               gradfieldm(k,j,i) = 0.0d0
               gradfield_ewald(k,j,i) = 0.0d0
               gradfieldd_thole(k,j,i) = 0.0d0
               gradfieldp_thole(k,j,i) = 0.0d0
               gradfieldm_gordon(k,j,i) = 0.0d0
               gradfieldd_gordon(k,j,i) = 0.0d0
               gradfieldp_gordon(k,j,i) = 0.0d0
               do l = 1, 3
                  hessfield(l,k,j,i) = 0.0d0
                  hessfieldm(l,k,j,i) = 0.0d0
                  hessfield_ewald(l,k,j,i) = 0.0d0
                  hessfieldm_gordon(l,k,j,i) = 0.0d0
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
      allocate (mscale(n))
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set highest order rr and damping terms needed
c     3 = up to field hessian (rr11)
c
      order = 3
      rorder = order*2 + 5
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
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     find the electrostatic potential, field and field gradient 
c     due to permanent multipoles
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
c
c     set m, d and p exclusion rules
c
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
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
               if (use_bounds) call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  rr11 = 9.0d0 * rr9 / r2
c
c     get tmatrix terms separated by powers of rr
c
                  call t0matrixrr1(rr1,t0rr1)
                  call t1matrixrr3(xr,yr,zr,rr3,t1rr3)
                  call t2matrixrr3(xr,yr,zr,rr3,t2rr3)
                  call t2matrixrr5(xr,yr,zr,rr5,t2rr5)
                  call t3matrixrr5(xr,yr,zr,rr5,t3rr5)
                  call t3matrixrr7(xr,yr,zr,rr7,t3rr7)
                  call t4matrixrr5(xr,yr,zr,rr5,t4rr5)
                  call t4matrixrr7(xr,yr,zr,rr7,t4rr7)
                  call t4matrixrr9(xr,yr,zr,rr9,t4rr9)
                  call t5matrixrr7(xr,yr,zr,rr7,t5rr7)
                  call t5matrixrr9(xr,yr,zr,rr9,t5rr9)
                  call t5matrixrr11(xr,yr,zr,rr11,t5rr11)
c
c     call routines that produce potential, field, field gradient
c     for types of damping
c
                  if (damp_none) then
                     t0 = t0rr1
                     t1 = t1rr3
                     t2 = t2rr3 + t2rr5
                     t3 = t3rr5 + t3rr7
                     t4 = t4rr5 + t4rr7 + t4rr9
                     t5 = t5rr7 + t5rr9 + t5rr11
                     call potik(i,k,t0,t1,t2,poti,potk)
                     call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                     call gradfieldik(i,k,t2,t3,t4,
     &                    gradfieldi,gradfieldk)
                     call hessfieldik(i,k,t3,t4,t5,
     &                    hessfieldi,hessfieldk)
                     pot(i) = pot(i) + poti
                     pot(k) = pot(k) + potk
                     potm(i) = potm(i) + poti*mscale(kk)
                     potm(k) = potm(k) + potk*mscale(kk)
                     do j = 1, 3
                        field(j,i) = field(j,i) + fieldi(j)
                        field(j,k) = field(j,k) + fieldk(j)
                        fieldm(j,i) = fieldm(j,i) + fieldi(j)*mscale(kk)
                        fieldm(j,k) = fieldm(j,k) + fieldk(j)*mscale(kk)
                        do l = 1, 3
                           gradfield(l,j,i) = gradfield(l,j,i) +
     &                          gradfieldi(l,j)
                           gradfield(l,j,k) = gradfield(l,j,k) +
     &                          gradfieldk(l,j)
                           gradfieldm(l,j,i) = gradfieldm(l,j,i) +
     &                          gradfieldi(l,j)*mscale(kk)
                           gradfieldm(l,j,k) = gradfieldm(l,j,k) +
     &                          gradfieldk(l,j)*mscale(kk)
                           do h = 1, 3
                              hessfield(h,l,j,i) = hessfield(h,l,j,i) +
     &                             hessfieldi(h,l,j)
                              hessfield(h,l,j,k) = hessfield(h,l,j,k) +
     &                             hessfieldk(h,l,j)
                              hessfieldm(h,l,j,i) = hessfieldm(h,l,j,i)+
     &                             hessfieldi(h,l,j)*mscale(kk)
                              hessfieldm(h,l,j,k) = hessfieldm(h,l,j,k)+
     &                             hessfieldk(h,l,j)*mscale(kk)
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
c     the ewald damping factors already contain their powers of rr
c
                     t0 = t0rr1*scale(1)/rr1
                     t1 = t1rr3*scale(3)/rr3
                     t2 = t2rr3*scale(3)/rr3 + t2rr5*scale(5)/rr5
                     t3 = t3rr5*scale(5)/rr5 + t3rr7*scale(7)/rr7
                     t4 = t4rr5*scale(5)/rr5 + t4rr7*scale(7)/rr7 +
     &                    t4rr9*scale(9)/rr9
                     t5 = t5rr7*scale(7)/rr7 + t5rr9*scale(9)/rr9 +
     &                    t5rr11*scale(11)/rr11
                     call potik(i,k,t0,t1,t2,poti,potk)
                     call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                     call gradfieldik(i,k,t2,t3,t4,
     &                    gradfieldi,gradfieldk)
                     call hessfieldik(i,k,t3,t4,t5,
     &                    hessfieldi,hessfieldk)
                     pot_ewald(i) = pot_ewald(i) + poti
                     pot_ewald(k) = pot_ewald(k) + potk
                     do j = 1, 3
                        field_ewald(j,i) = field_ewald(j,i) + 
     &                       fieldi(j)
                        field_ewald(j,k) = field_ewald(j,k) + 
     &                       fieldk(j)
                        do l = 1, 3
                           gradfield_ewald(l,j,i) = 
     &                          gradfield_ewald(l,j,i) + 
     &                          gradfieldi(l,j)
                           gradfield_ewald(l,j,k) = 
     &                          gradfield_ewald(l,j,k) +
     &                          gradfieldk(l,j)
                           do h = 1, 3
                              hessfield_ewald(h,l,j,i) = 
     &                             hessfield_ewald(h,l,j,i) +
     &                             hessfieldi(h,l,j)
                              hessfield_ewald(h,l,j,k) =
     &                             hessfield_ewald(h,l,j,k) +
     &                             hessfieldk(h,l,j)
                           end do
                        end do
                     end do
                  end if
c
c     thole damping
c
                  if (damp_thole) then
                     call dampthole(i,k,rorder,r,scale)
                     t0 = t0rr1*scale(1)
                     t1 = t1rr3*scale(3)
                     t2 = t2rr3*scale(3) + t2rr5*scale(5)
                     t2 = t2rr3*scale(3) + t2rr5*scale(5)
                     t3 = t3rr5*scale(5) + t3rr7*scale(7)
                     t4 = t4rr5*scale(5) + t4rr7*scale(7) + 
     &                    t4rr9*scale(9)
                     call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                     call gradfieldik(i,k,t2,t3,t4,gradfieldi,
     &                    gradfieldk)
                     do j = 1, 3
                        fieldd_thole(j,i) = fieldd_thole(j,i) +
     &                       fieldi(j)*dscale(kk)
                        fieldd_thole(j,k) = fieldd_thole(j,k) +
     &                       fieldk(j)*dscale(kk)
                        fieldp_thole(j,i) = fieldp_thole(j,i) +
     &                       fieldi(j)*pscale(kk)
                        fieldp_thole(j,k) = fieldp_thole(j,k) +
     &                       fieldk(j)*pscale(kk)
                        do l = 1, 3
                           gradfieldd_thole(l,j,i) =
     &                          gradfieldd_thole(l,j,i) +
     &                          gradfieldi(l,j)*dscale(kk)
                           gradfieldd_thole(l,j,k) =
     &                          gradfieldd_thole(l,j,k) +
     &                          gradfieldk(l,j)*dscale(kk)
                           gradfieldp_thole(l,j,i) =
     &                          gradfieldp_thole(l,j,i) +
     &                          gradfieldi(l,j)*pscale(kk)
                           gradfieldp_thole(l,j,k) =
     &                          gradfieldp_thole(l,j,k) +
     &                          gradfieldk(l,j)*pscale(kk)
                        end do
                     end do
                  end if
c
c     gordon damping - splits charge into nuclear and electronic
c
                  if (damp_gordon) then
                     call dampgordon(i,k,rorder,r,scalei,scalek,scaleik)
                     t0 = t0rr1
                     t1 = t1rr3
                     t2 = t2rr3 + t2rr5
                     t3 = t3rr5 + t3rr7
                     t4 = t4rr5 + t4rr7 + t4rr9
                     t5 = t5rr7 + t5rr9 + t5rr11
                     t0i = t0rr1*scalei(1)
                     t1i = t1rr3*scalei(3)
                     t2i = t2rr3*scalei(3) + t2rr5*scalei(5)
                     t3i = t3rr5*scalei(5) + t3rr7*scalei(7)
                     t4i = t4rr5*scalei(5) + t4rr7*scalei(7) +
     &                     t4rr9*scalei(9)
                     t5i = t5rr7*scalei(7) + t5rr9*scalei(9) + 
     &                     t5rr11*scalei(11)
                     t0k = t0rr1*scalek(1)
                     t1k = t1rr3*scalek(3)
                     t2k = t2rr3*scalek(3) + t2rr5*scalek(5)
                     t3k = t3rr5*scalek(5) + t3rr7*scalek(7)
                     t4k = t4rr5*scalek(5) + t4rr7*scalek(7) +
     &                     t4rr9*scalek(9)
                     t5k = t5rr7*scalek(7) + t5rr9*scalek(9) +
     &                     t5rr11*scalek(11)
                     t0ik = t0rr1*scaleik(1)
                     t1ik = t1rr3*scaleik(3)
                     t2ik = t2rr3*scaleik(3) + t2rr5*scaleik(5)
                     t3ik = t3rr5*scaleik(5) + t3rr7*scaleik(7)
                     t4ik = t4rr5*scaleik(5) + t4rr7*scaleik(7) +
     &                      t4rr9*scaleik(9)
                     t5ik = t5rr7*scaleik(7) + t5rr9*scaleik(9) +
     &                      t5rr11*scaleik(11)
                     call cp_potik(i,k,
     &                  t0,t1,t2,t0i,t1i,t2i,t0k,t1k,t2k,t0ik,t1ik,t2ik,
     &                  nucpoti,nucpotk,elepoti,elepotk)
                     call cp_fieldik(i,k,
     &                  t1,t2,t3,t1i,t2i,t3i,t1k,t2k,t3k,t1ik,t2ik,t3ik,
     &                    nucfieldi,nucfieldk,elefieldi,elefieldk)
                     call cp_gradfieldik(i,k,
     &                  t2,t3,t4,t2i,t3i,t4i,t2k,t3k,t4k,t2ik,t3ik,t4ik,
     &                    elegradfieldi,elegradfieldk)
                     call cp_hessfieldik(i,k,
     &                  t3,t4,t5,t3i,t4i,t5i,t3k,t4k,t5k,t3ik,t4ik,t5ik,
     &                    elehessfieldi,elehessfieldk)
c
c     accumulate fields in global variables
c
                     nucpotm_gordon(i) = nucpotm_gordon(i) + 
     &                    nucpoti*mscale(kk)
                     nucpotm_gordon(k) = nucpotm_gordon(k) + 
     &                    nucpotk*mscale(kk)
                     potm_gordon(i) = potm_gordon(i) +elepoti*mscale(kk)
                     potm_gordon(k) = potm_gordon(k) +elepotk*mscale(kk)
                     do j = 1, 3
                        nucfieldm_gordon(j,i) = nucfieldm_gordon(j,i) + 
     &                       nucfieldi(j)*mscale(kk)
                        nucfieldm_gordon(j,k) = nucfieldm_gordon(j,k) +
     &                       nucfieldk(j)*mscale(kk)
                        fieldm_gordon(j,i) = fieldm_gordon(j,i) +
     &                       elefieldi(j)*mscale(kk)
                        fieldm_gordon(j,k) = fieldm_gordon(j,k) +
     &                       elefieldk(j)*mscale(kk)
                        fieldd_gordon(j,i) = fieldd_gordon(j,i) +
     &                       elefieldi(j)*dscale(kk)
                        fieldd_gordon(j,k) = fieldd_gordon(j,k) +
     &                       elefieldk(j)*dscale(kk)
                        fieldp_gordon(j,i) = fieldp_gordon(j,i) +
     &                       elefieldi(j)*pscale(kk)
                        fieldp_gordon(j,k) = fieldp_gordon(j,k) +
     &                       elefieldk(j)*pscale(kk)
                        do l = 1, 3
                           gradfieldm_gordon(l,j,i) =
     &                          gradfieldm_gordon(l,j,i) +
     &                          elegradfieldi(l,j)*mscale(kk)
                           gradfieldm_gordon(l,j,k) =
     &                          gradfieldm_gordon(l,j,k) +
     &                          elegradfieldk(l,j)*mscale(kk)
                           gradfieldd_gordon(l,j,i) =
     &                          gradfieldd_gordon(l,j,i) +
     &                          elegradfieldi(l,j)*dscale(kk)
                           gradfieldd_gordon(l,j,k) =
     &                          gradfieldd_gordon(l,j,k) +
     &                          elegradfieldk(l,j)*dscale(kk)
                           gradfieldp_gordon(l,j,i) =
     &                          gradfieldp_gordon(l,j,i) +
     &                          elegradfieldi(l,j)*pscale(kk)
                           gradfieldp_gordon(l,j,k) =
     &                          gradfieldp_gordon(l,j,k) +
     &                          elegradfieldk(l,j)*pscale(kk)
                           do h = 1, 3
                              hessfieldm_gordon(h,l,j,i) =
     &                             hessfieldm_gordon(h,l,j,i) +
     &                             elehessfieldi(h,l,j)*mscale(kk)
                              hessfieldm_gordon(h,l,j,k) =
     &                             hessfieldm_gordon(h,l,j,k) +
     &                             elehessfieldk(h,l,j)*mscale(kk)
                           end do
                        end do
                     end do
                  end if
c
c     gordon damping + nuclear regularization
c
                  if (damp_gordonreg) then
                     if (.not. damp_gordon) then
                        call dampgordon(i,k,rorder,r,scalei,scalek,
     &                       scaleik)
                     end if
                     call dampgordonreg(i,k,rorder,r,scalei,scalek)
                     t1 = t1rr3
                     t2 = t2rr3 + t2rr5
                     t3 = t3rr5 + t3rr7
                     t4 = t4rr5 + t4rr7 + t4rr9
                     t1i = t1rr3*scalei(3)
                     t2i = t2rr3*scalei(3) + t2rr5*scalei(5)
                     t3i = t3rr5*scalei(5) + t3rr7*scalei(7)
                     t4i = t4rr5*scalei(5) + t4rr7*scalei(7) +
     &                     t4rr9*scalei(9)
                     t1k = t1rr3*scalek(3)
                     t2k = t2rr3*scalek(3) + t2rr5*scalek(5)
                     t3k = t3rr5*scalek(5) + t3rr7*scalek(7)
                     t4k = t4rr5*scalek(5) + t4rr7*scalek(7) +
     &                     t4rr9*scalek(9)
                     t1ik = t1rr3*scaleik(3)
                     t2ik = t2rr3*scaleik(3) + t2rr5*scaleik(5)
                     t3ik = t3rr5*scaleik(5) + t3rr7*scaleik(7)
                     t4ik = t4rr5*scaleik(5) + t4rr7*scaleik(7) +
     &                      t4rr9*scaleik(9)
                     call cp_fieldik(i,k,
     &                  t1,t2,t3,t1i,t2i,t3i,t1k,t2k,t3k,t1ik,t2ik,t3ik,
     &                    nucfieldi,nucfieldk,elefieldi,elefieldk)
                     do j = 1, 3
                        fieldd_gordonreg(j,i) = fieldd_gordonreg(j,i) +
     &                       elefieldi(j)*dscale(kk)
                        fieldd_gordonreg(j,k) = fieldd_gordonreg(j,k) +
     &                       elefieldk(j)*dscale(kk)
                        fieldp_gordonreg(j,i) = fieldp_gordonreg(j,i) +
     &                       elefieldi(j)*pscale(kk)
                        fieldp_gordonreg(j,k) = fieldp_gordonreg(j,k) +
     &                       elefieldk(j)*pscale(kk)
                     end do
                  end if
c
c     one center gordon damping
c
c                  if (damp_gordonone) then
c                     call dampgordonone(i,k,rorder,r,scale)
c                     t0 = t0rr1*scale(1)
c                     t1 = t1rr3*scale(3)
c                     t2 = t2rr3*scale(3) + t2rr5*scale(5)
c                     t3 = t3rr5*scale(5) + t3rr7*scale(7)
c                     t4 = t4rr5*scale(5) + t4rr7*scale(7) +
c     &                    t4rr9*scale(9)
c                     call potik(i,k,t0,t1,t2,poti,potk)
c                     call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
c                     call gradfieldik(i,k,t2,t3,t4,
c     &                    gradfieldi,gradfieldk)
c                     potm_gordonone(i)=potm_gordonone(i)+poti*mscale(kk)
c                     potm_gordonone(k)=potm_gordonone(k)+potk*mscale(kk)
c                     do j = 1, 3
c                        fieldm_gordonone(j,i) = fieldm_gordonone(j,i) +
c     &                       fieldi(j)*mscale(kk)
c                        fieldm_gordonone(j,k) = fieldm_gordonone(j,k) +
c     &                       fieldk(j)*mscale(kk)
c                        fieldd_gordonone(j,i) = fieldd_gordonone(j,i) +
c     &                       fieldi(j)*dscale(kk)
c                        fieldd_gordonone(j,k) = fieldd_gordonone(j,i) +
c     &                       fieldk(j)*dscale(kk)
c                        fieldp_gordonone(j,i) = fieldp_gordonone(j,i) +
c     &                       fieldi(j)*pscale(kk)
c                        fieldp_gordonone(j,k) = fieldp_gordonone(j,k) +
c     &                       fieldk(j)*pscale(kk)
c                        do l = 1, 3
c                           gradfieldm_gordonone(l,j,i) =
c     &                          gradfieldm_gordonone(l,j,i) +
c     &                          gradfieldi(l,j)*mscale(kk)
c                           gradfieldm_gordonone(l,j,k) =
c     &                          gradfieldm_gordonone(l,j,k) +
c     &                          gradfieldk(l,j)*mscale(kk)
c                           gradfieldp_gordonone(l,j,i) =
c     &                          gradfieldp_gordonone(l,j,i) +
c     &                          gradfieldi(l,j)*pscale(kk)
c                           gradfieldp_gordonone(l,j,k) =
c     &                          gradfieldp_gordonone(l,j,k) +
c     &                          gradfieldk(l,j)*pscale(kk)
c                        end do
c                     end do
c                  end if
c
c     two center gordon damping
c
c                  if (damp_gordontwo) then
c                     call dampgordontwo(i,k,rorder,r,scale)
c                     t0 = t0rr1*scale(1)
c                     t1 = t1rr3*scale(3)
c                     t2 = t2rr3*scale(3) + t2rr5*scale(5)
c                     t3 = t3rr5*scale(5) + t3rr7*scale(7)
c                     t4 = t4rr5*scale(5) + t4rr7*scale(7) +
c     &                    t4rr9*scale(9)
c                     call potik(i,k,t0,t1,t2,poti,potk)
c                     call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
c                     call gradfieldik(i,k,t2,t3,t4,
c     &                    gradfieldi,gradfieldk)
c                     potm_gordontwo(i)=potm_gordontwo(i)+poti*mscale(kk)   
c                     potm_gordontwo(k)=potm_gordontwo(k)+potk*mscale(kk)   
c                     do j = 1, 3
c                        fieldm_gordontwo(j,i) = fieldm_gordontwo(j,i) +
c     &                       fieldi(j)*mscale(kk)
c                        fieldm_gordontwo(j,k) = fieldm_gordontwo(j,k) +
c     &                       fieldk(j)*mscale(kk)
c                        fieldd_gordontwo(j,i) = fieldd_gordontwo(j,i) +
c     &                       fieldi(j)*dscale(kk)
c                        fieldd_gordontwo(j,k) = fieldd_gordontwo(j,i) +
c     &                       fieldk(j)*dscale(kk)
c                        fieldp_gordontwo(j,i) = fieldp_gordontwo(j,i) +
c     &                       fieldi(j)*pscale(kk)
c                        fieldp_gordontwo(j,k) = fieldp_gordontwo(j,k) +
c     &                       fieldk(j)*pscale(kk)
c                        do l = 1, 3
c                           gradfieldm_gordontwo(l,j,i) =
c     &                          gradfieldm_gordontwo(l,j,i) +
c     &                          gradfieldi(l,j)*mscale(kk)
c                           gradfieldm_gordontwo(l,j,k) =
c     &                          gradfieldm_gordontwo(l,j,k) +
c     &                          gradfieldk(l,j)*mscale(kk)
c                           gradfieldp_gordontwo(l,j,i) =
c     &                          gradfieldp_gordontwo(l,j,i) +
c     &                          gradfieldi(l,j)*pscale(kk)
c                           gradfieldp_gordontwo(l,j,k) =
c     &                          gradfieldp_gordontwo(l,j,k) +
c     &                          gradfieldk(l,j)*pscale(kk)
c                        end do
c                     end do
c                  end if
c
c     one center piquemal damping
c
c                  if (damp_piquemalone) then 
c                     call damppiquemalone(i,k,rorder,r,scale)
c                     t0 = t0rr1*scale(1)
c                     t1 = t1rr3*scale(3)
c                     t2 = t2rr3*scale(3) + t2rr5*scale(5)
c                     t3 = t3rr5*scale(5) + t3rr7*scale(7)
c                     t4 = t4rr5*scale(5) + t4rr7*scale(7) +
c     &                    t4rr9*scale(9)
c                     call potik(i,k,t0,t1,t2,poti,potk)
c                     call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
c                     call gradfieldik(i,k,t2,t3,t4,
c     &                    gradfieldi,gradfieldk)
c                     potm_piquemalone(i) = potm_piquemalone(i) + 
c     &                    poti*mscale(kk)
c                     potm_piquemalone(k) = potm_piquemalone(k) + 
c     &                    potk*mscale(kk)
c                     do j = 1, 3
c                        fieldm_piquemalone(j,i)=fieldm_piquemalone(j,i)+
c     &                       fieldi(j)*mscale(kk)
c                        fieldm_piquemalone(j,k)=fieldm_piquemalone(j,k)+
c     &                       fieldk(j)*mscale(kk)
c                        fieldd_piquemalone(j,i)=fieldd_piquemalone(j,i)+
c     &                       fieldi(j)*dscale(kk)
c                        fieldd_piquemalone(j,k)=fieldd_piquemalone(j,i)+
c     &                       fieldk(j)*dscale(kk)
c                        fieldp_piquemalone(j,i)=fieldp_piquemalone(j,i)+
c     &                       fieldi(j)*pscale(kk)
c                        fieldp_piquemalone(j,k)=fieldp_piquemalone(j,k)+
c     &                       fieldk(j)*pscale(kk)
c                        do l = 1, 3
c                           gradfieldm_piquemalone(l,j,i) =
c     &                          gradfieldm_piquemalone(l,j,i) +
c     &                          gradfieldi(l,j)*mscale(kk)
c                           gradfieldm_piquemalone(l,j,k) =
c     &                          gradfieldm_piquemalone(l,j,k) +
c     &                          gradfieldk(l,j)*mscale(kk)
c                           gradfieldp_piquemalone(l,j,i) =
c     &                          gradfieldp_piquemalone(l,j,i) +
c     &                          gradfieldi(l,j)*pscale(kk)
c                           gradfieldp_piquemalone(l,j,k) =
c     &                          gradfieldp_piquemalone(l,j,k) +
c     &                          gradfieldk(l,j)*pscale(kk)
c                        end do
c                     end do
c                  end if
c
c     two center piquemal damping
c
c                  if (damp_piquemaltwo) then
c                     call damppiquemaltwo(i,k,rorder,r,scale)
c                     t0 = t0rr1*scale(1)
c                     t1 = t1rr3*scale(3)
c                     t2 = t2rr3*scale(3) + t2rr5*scale(5)
c                     t3 = t3rr5*scale(5) + t3rr7*scale(7)
c                     t4 = t4rr5*scale(5) + t4rr7*scale(7) +
c     &                    t4rr9*scale(9)
c                     call potik(i,k,t0,t1,t2,poti,potk)
c                     call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
c                     call gradfieldik(i,k,t2,t3,t4,
c     &                    gradfieldi,gradfieldk)
c                     potm_piquemaltwo(i) = potm_piquemaltwo(i) +
c     &                    poti*mscale(kk)
c                     potm_piquemaltwo(k) = potm_piquemaltwo(k) +
c     &                    potk*mscale(kk)
c                     do j = 1, 3
c                        fieldm_piquemaltwo(j,i)=fieldm_piquemaltwo(j,i)+
c     &                       fieldi(j)*mscale(kk)
c                        fieldm_piquemaltwo(j,k)=fieldm_piquemaltwo(j,k)+
c     &                       fieldk(j)*mscale(kk)
c                        fieldd_piquemaltwo(j,i)=fieldd_piquemaltwo(j,i)+
c     &                       fieldi(j)*dscale(kk)
c                        fieldd_piquemaltwo(j,k)=fieldd_piquemaltwo(j,i)+
c     &                       fieldk(j)*dscale(kk)
c                        fieldp_piquemaltwo(j,i)=fieldp_piquemaltwo(j,i)+
c     &                       fieldi(j)*pscale(kk)
c                        fieldp_piquemaltwo(j,k)=fieldp_piquemaltwo(j,k)+
c     &                       fieldk(j)*pscale(kk)
c                        do l = 1, 3
c                           gradfieldm_piquemaltwo(l,j,i) =
c     &                          gradfieldm_piquemaltwo(l,j,i) +
c     &                          gradfieldi(l,j)*mscale(kk)
c                           gradfieldm_piquemaltwo(l,j,k) =
c     &                          gradfieldm_piquemaltwo(l,j,k) +
c     &                          gradfieldk(l,j)*mscale(kk)
c                           gradfieldp_piquemaltwo(l,j,i) =
c     &                          gradfieldp_piquemaltwo(l,j,i) +
c     &                          gradfieldi(l,j)*pscale(kk)
c                           gradfieldp_piquemaltwo(l,j,k) =
c     &                          gradfieldp_piquemaltwo(l,j,k) +
c     &                          gradfieldk(l,j)*pscale(kk)
c                        end do
c                     end do
c                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            dscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            dscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            dscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            dscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
      print *,"USING REPLICA"
c     
c     calculate interaction with other unit cells
c
      do i = 1, npole
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         ii = ipole(i)
c
c     set m, d and p exclusion rules
c
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
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
                     rr11 = 9.0d0 * rr9 / r2
                     call t0matrixrr1(rr1,t0rr1)
                     call t1matrixrr3(xr,yr,zr,rr3,t1rr3)
                     call t2matrixrr3(xr,yr,zr,rr3,t2rr3)
                     call t2matrixrr5(xr,yr,zr,rr5,t2rr5)
                     call t3matrixrr5(xr,yr,zr,rr5,t3rr5)
                     call t3matrixrr7(xr,yr,zr,rr7,t3rr7)
                     call t4matrixrr5(xr,yr,zr,rr5,t4rr5)
                     call t4matrixrr7(xr,yr,zr,rr7,t4rr7)
                     call t4matrixrr9(xr,yr,zr,rr9,t4rr9)
                     call t5matrixrr7(xr,yr,zr,rr7,t5rr7)
                     call t5matrixrr9(xr,yr,zr,rr9,t5rr9)
                     call t5matrixrr11(xr,yr,zr,rr11,t5rr11)
c
c     call routines that produce potential, field, field gradient
c     for types of damping
c
                     if (damp_none) then
                        t0 = t0rr1
                        t1 = t1rr3
                        t2 = t2rr3 + t2rr5
                        t3 = t3rr5 + t3rr7
                        t4 = t4rr5 + t4rr7 + t4rr9
                        t5 = t5rr7 + t5rr9 + t5rr11
                        call potik(i,k,t0,t1,t2,poti,potk)
                        call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                        call gradfieldik(i,k,t2,t3,t4,
     &                       gradfieldi,gradfieldk)
                        call hessfieldik(i,k,t3,t4,t5,
     &                       hessfieldi,hessfieldk)
                        pot(i) = pot(i) + poti
                        pot(k) = pot(k) + potk
                        potm(i) = potm(i) + poti*mscale(kk)
                        potm(k) = potm(k) + potk*mscale(kk)
                        do j = 1, 3
                           field(j,i) = field(j,i) + fieldi(j)
                           field(j,k) = field(j,k) + fieldk(j)
                           fieldm(j,i) = fieldm(j,i) + 
     &                          fieldi(j)*mscale(kk)
                           fieldm(j,k) = fieldm(j,k) + 
     &                          fieldk(j)*mscale(kk)
                           do l = 1, 3
                              gradfield(l,j,i) = gradfield(l,j,i) +
     &                             gradfieldi(l,j)
                              gradfield(l,j,k) = gradfield(l,j,k) +
     &                             gradfieldk(l,j)
                              gradfieldm(l,j,i) = gradfieldm(l,j,i) +
     &                             gradfieldi(l,j)*mscale(kk)
                              gradfieldm(l,j,k) = gradfieldm(l,j,k) +
     &                             gradfieldk(l,j)*mscale(kk)
                              do h = 1, 3
                                 hessfield(h,l,j,i)=hessfield(h,l,j,i) 
     &                                + hessfieldi(h,l,j)
                                 hessfield(h,l,j,k)=hessfield(h,l,j,k) 
     &                                + hessfieldk(h,l,j)
                                 hessfieldm(h,l,j,i)=hessfieldm(h,l,j,i)
     &                                + hessfieldi(h,l,j)*mscale(kk)
                                 hessfieldm(h,l,j,k)=hessfieldm(h,l,j,k)
     &                                + hessfieldk(h,l,j)*mscale(kk)
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
c     the ewald damping factors already contain their powers of rr
c     
                        t0 = t0rr1*scale(1)/rr1
                        t1 = t1rr3*scale(3)/rr3
                        t2 = t2rr3*scale(3)/rr3 + t2rr5*scale(5)/rr5
                        t3 = t3rr5*scale(5)/rr5 + t3rr7*scale(7)/rr7
                        t4 = t4rr5*scale(5)/rr5 + t4rr7*scale(7)/rr7 +
     &                       t4rr9*scale(9)/rr9
                        t5 = t5rr7*scale(7)/rr7 + t5rr9*scale(9)/rr9 +
     &                       t5rr11*scale(11)/rr11
                        call potik(i,k,t0,t1,t2,poti,potk)
                        call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                        call gradfieldik(i,k,t2,t3,t4,
     &                       gradfieldi,gradfieldk)
                        call hessfieldik(i,k,t3,t4,t5,
     &                       hessfieldi,hessfieldk)
                        pot_ewald(i) = pot_ewald(i) + poti
                        pot_ewald(k) = pot_ewald(k) + potk
                        do j = 1, 3
                           field_ewald(j,i) = field_ewald(j,i) + 
     &                          fieldi(j)
                           field_ewald(j,k) = field_ewald(j,k) + 
     &                          fieldk(j)
                           do l = 1, 3
                              gradfield_ewald(l,j,i) = 
     &                             gradfield_ewald(l,j,i) + 
     &                             gradfieldi(l,j)
                              gradfield_ewald(l,j,k) = 
     &                             gradfield_ewald(l,j,k) +
     &                             gradfieldk(l,j)
                              do h = 1, 3
                                 hessfield_ewald(h,l,j,i) = 
     &                                hessfield_ewald(h,l,j,i) +
     &                                hessfieldi(h,l,j)
                                 hessfield_ewald(h,l,j,k) =
     &                                hessfield_ewald(h,l,j,k) +
     &                                hessfieldk(h,l,j)
                              end do
                           end do
                        end do
                     end if
c     
c     thole damping
c     
                     if (damp_thole) then
                        call dampthole(i,k,rorder,r,scale)
                        t0 = t0rr1*scale(1)
                        t1 = t1rr3*scale(3)
                        t2 = t2rr3*scale(3) + t2rr5*scale(5)
                        t2 = t2rr3*scale(3) + t2rr5*scale(5)
                        t3 = t3rr5*scale(5) + t3rr7*scale(7)
                        t4 = t4rr5*scale(5) + t4rr7*scale(7) + 
     &                       t4rr9*scale(9)
                        call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                        call gradfieldik(i,k,t2,t3,t4,gradfieldi,
     &                       gradfieldk)
                        do j = 1, 3
                           fieldd_thole(j,i) = fieldd_thole(j,i) +
     &                          fieldi(j)*dscale(kk)
                           fieldd_thole(j,k) = fieldd_thole(j,k) +
     &                          fieldk(j)*dscale(kk)
                           fieldp_thole(j,i) = fieldp_thole(j,i) +
     &                          fieldi(j)*pscale(kk)
                           fieldp_thole(j,k) = fieldp_thole(j,k) +
     &                          fieldk(j)*pscale(kk)
                           do l = 1, 3
                              gradfieldd_thole(l,j,i) =
     &                             gradfieldd_thole(l,j,i) +
     &                             gradfieldi(l,j)*dscale(kk)
                              gradfieldd_thole(l,j,k) =
     &                             gradfieldd_thole(l,j,k) +
     &                             gradfieldk(l,j)*dscale(kk)
                              gradfieldp_thole(l,j,i) =
     &                             gradfieldp_thole(l,j,i) +
     &                             gradfieldi(l,j)*pscale(kk)
                              gradfieldp_thole(l,j,k) =
     &                             gradfieldp_thole(l,j,k) +
     &                             gradfieldk(l,j)*pscale(kk)
                           end do
                        end do
                     end if
                  end if
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            dscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            dscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            dscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            dscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (dscale)
      deallocate (pscale)
      deallocate (scale)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine permfield3b  --  potential, field, field     ##
c     ##            gradient and field hessian via neighbor list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "permfield3b" computes the direct electrostatic potential, field, 
c     field gradient and field hessian due to permanent multipole 
c     moments using a neighbor list
c
c
      subroutine permfield3b
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
      real*8 rr1,rr3,rr5,rr7,rr9,rr11
      real*8 t0,t1(3),t2(3,3),t3(3,3,3),t4(3,3,3,3),t5(3,3,3,3,3)
      real*8 t0rr1,t1rr3(3)
      real*8 t2rr3(3,3),t2rr5(3,3)
      real*8 t3rr5(3,3,3),t3rr7(3,3,3)
      real*8 t4rr5(3,3,3,3),t4rr7(3,3,3,3),t4rr9(3,3,3,3)
      real*8 t5rr7(3,3,3,3,3),t5rr9(3,3,3,3,3),t5rr11(3,3,3,3,3)
      real*8 poti,potk
      real*8 fieldi(3),fieldk(3)
      real*8 gradfieldi(3,3),gradfieldk(3,3)
      real*8 hessfieldi(3,3,3),hessfieldk(3,3,3)
      real*8, allocatable :: scale(:)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: poto(:)
      real*8, allocatable :: fieldo(:,:)
      real*8, allocatable :: gradfieldo(:,:,:)
      real*8, allocatable :: hessfieldo(:,:,:,:)
      real*8, allocatable :: potmo(:)
      real*8, allocatable :: fieldmo(:,:)
      real*8, allocatable :: gradfieldmo(:,:,:)
      real*8, allocatable :: hessfieldmo(:,:,:,:)
      real*8, allocatable :: pot_ewaldo(:)
      real*8, allocatable :: field_ewaldo(:,:)
      real*8, allocatable :: gradfield_ewaldo(:,:,:)
      real*8, allocatable :: hessfield_ewaldo(:,:,:,:)
      real*8, allocatable :: fieldp_tholeo(:,:)
      real*8, allocatable :: fieldd_tholeo(:,:)
      real*8, allocatable :: gradfieldp_tholeo(:,:,:)
      real*8, allocatable :: gradfieldd_tholeo(:,:,:)
      logical proceed
      logical usei,usek
      character*6 mode
c
c     perform dynamic allocation of some local arrays
c
      allocate (poto(npole))
      allocate (fieldo(3,npole))
      allocate (gradfieldo(3,3,npole))
      allocate (hessfieldo(3,3,3,npole))
c
      allocate (potmo(npole))
      allocate (fieldmo(3,npole))
      allocate (gradfieldmo(3,3,npole))
      allocate (hessfieldmo(3,3,3,npole))
c
      allocate (pot_ewaldo(npole))
      allocate (field_ewaldo(3,npole))
      allocate (gradfield_ewaldo(3,3,npole))
      allocate (hessfield_ewaldo(3,3,3,npole))
c
      allocate (fieldd_tholeo(3,npole))
      allocate (fieldp_tholeo(3,npole))
      allocate (gradfieldd_tholeo(3,3,npole))
      allocate (gradfieldp_tholeo(3,3,npole))
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         pot(i) = 0.0d0
         potm(i) = 0.0d0
         pot_ewald(i) = 0.0d0
c
         poto(i) = 0.0d0
         potmo(i) = 0.0d0
         pot_ewaldo(i) = 0.0d0
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldm(j,i) = 0.0d0
            field_ewald(j,i) = 0.0d0
            fieldd_thole(j,i) = 0.0d0
            fieldp_thole(j,i) = 0.0d0
c
            fieldo(j,i) = 0.0d0
            fieldmo(j,i) = 0.0d0
            field_ewaldo(j,i) = 0.0d0
            fieldd_tholeo(j,i) = 0.0d0
            fieldp_tholeo(j,i) = 0.0d0
            do k = 1, 3
               gradfield(k,j,i) = 0.0d0
               gradfieldm(k,j,i) = 0.0d0
               gradfield_ewald(k,j,i) = 0.0d0
               gradfieldd_thole(k,j,i) = 0.0d0
               gradfieldp_thole(k,j,i) = 0.0d0
c
               gradfieldo(k,j,i) = 0.0d0
               gradfieldmo(k,j,i) = 0.0d0
               gradfield_ewaldo(k,j,i) = 0.0d0
               gradfieldd_tholeo(k,j,i) = 0.0d0
               gradfieldp_tholeo(k,j,i) = 0.0d0
               do l = 1, 3
                  hessfield(l,k,j,i) = 0.0d0
                  hessfieldm(l,k,j,i) = 0.0d0
                  hessfield_ewald(l,k,j,i) = 0.0d0
c
                  hessfieldo(l,k,j,i) = 0.0d0
                  hessfieldmo(l,k,j,i) = 0.0d0
                  hessfield_ewaldo(l,k,j,i) = 0.0d0
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
      rorder = order*3 + 5
      allocate (scale(rorder))
      do i = 1,rorder
          scale(i) = 0.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared)
!$OMP& private(i,j,k,ii,ix,iy,iz,usei,kk,kx,ky,kz,usek,kkk,proceed,
!$OMP& xr,yr,zr,r,r2,rr1,rr3,rr5,rr7,rr9,rr11,fgrp,
!$OMP& poti,potk,fieldi,fieldk,gradfieldi,gradfieldk,
!$OMP& hessfieldi,hessfieldk,
!$OMP& t0rr1,t1rr3,t2rr3,t2rr5,t3rr5,t3rr7,t4rr5,t4rr7,t4rr9,
!$OMP& t5rr7,t5rr9,t5rr11,
!$OMP& t0,t1,t2,t3,t4,t5,scale)
!$OMP& firstprivate(mscale,dscale,pscale)
!$OMP DO reduction(+:poto,fieldo,gradfieldo,hessfieldo,
!$OMP& potmo,fieldmo,gradfieldmo,hessfieldmo,
!$OMP& pot_ewaldo,field_ewaldo,gradfield_ewaldo,hessfield_ewaldo,
!$OMP& fieldd_tholeo,fieldp_tholeo,
!$OMP& gradfieldd_tholeo,gradfieldp_tholeo)
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
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
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
               rr11 = 9.0d0 * rr9 / r2
c
c     get tmatrix terms separated by powers of rr
c
               call t0matrixrr1(rr1,t0rr1)
               call t1matrixrr3(xr,yr,zr,rr3,t1rr3)
               call t2matrixrr3(xr,yr,zr,rr3,t2rr3)
               call t2matrixrr5(xr,yr,zr,rr5,t2rr5)
               call t3matrixrr5(xr,yr,zr,rr5,t3rr5)
               call t3matrixrr7(xr,yr,zr,rr7,t3rr7)
               call t4matrixrr5(xr,yr,zr,rr5,t4rr5)
               call t4matrixrr7(xr,yr,zr,rr7,t4rr7)
               call t4matrixrr9(xr,yr,zr,rr9,t4rr9)
               call t5matrixrr7(xr,yr,zr,rr7,t5rr7)
               call t5matrixrr9(xr,yr,zr,rr9,t5rr9)
               call t5matrixrr11(xr,yr,zr,rr11,t5rr11)
c
c     call routines that produce potential, field, field gradient
c     for types of damping
c
               if (damp_none) then
                  t0 = t0rr1
                  t1 = t1rr3
                  t2 = t2rr3 + t2rr5
                  t3 = t3rr5 + t3rr7
                  t4 = t4rr5 + t4rr7 + t4rr9
                  t5 = t5rr7 + t5rr9 + t5rr11
                  call potik(i,k,t0,t1,t2,poti,potk)
                  call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                  call gradfieldik(i,k,t2,t3,t4,
     &                 gradfieldi,gradfieldk)
                  call hessfieldik(i,k,t3,t4,t5,
     &                 hessfieldi,hessfieldk)
                  poto(i) = poto(i) + poti
                  poto(k) = poto(k) + potk
                  potmo(i) = potmo(i) + poti*mscale(kk)
                  potmo(k) = potmo(k) + potk*mscale(kk)
                  do j = 1, 3
                     fieldo(j,i) = fieldo(j,i) + fieldi(j)
                     fieldo(j,k) = fieldo(j,k) + fieldk(j)
                     fieldmo(j,i) = fieldmo(j,i) + 
     &                    fieldi(j)*mscale(kk)
                     fieldmo(j,k) = fieldmo(j,k) + 
     &                    fieldk(j)*mscale(kk)
                     do l = 1, 3
                        gradfieldo(l,j,i) = gradfieldo(l,j,i) +
     &                       gradfieldi(l,j)
                        gradfieldo(l,j,k) = gradfieldo(l,j,k) +
     &                       gradfieldk(l,j)
                        gradfieldmo(l,j,i) = gradfieldmo(l,j,i) +
     &                       gradfieldi(l,j)*mscale(kk)
                        gradfieldmo(l,j,k) = gradfieldmo(l,j,k) +
     &                       gradfieldk(l,j)*mscale(kk)
                        do h = 1, 3
                           hessfieldo(h,l,j,i) = hessfieldo(h,l,j,i)+
     &                          hessfieldi(h,l,j)
                           hessfieldo(h,l,j,k) = hessfieldo(h,l,j,k)+
     &                          hessfieldk(h,l,j)
                           hessfieldmo(h,l,j,i)=hessfieldmo(h,l,j,i)+
     &                          hessfieldi(h,l,j)*mscale(kk)
                           hessfieldmo(h,l,j,k)=hessfieldmo(h,l,j,k)+
     &                          hessfieldk(h,l,j)*mscale(kk)
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
c     the ewald damping factors already contain their powers of rr
c
                  t0 = t0rr1*scale(1)/rr1
                  t1 = t1rr3*scale(3)/rr3
                  t2 = t2rr3*scale(3)/rr3 + t2rr5*scale(5)/rr5
                  t3 = t3rr5*scale(5)/rr5 + t3rr7*scale(7)/rr7
                  t4 = t4rr5*scale(5)/rr5 + t4rr7*scale(7)/rr7 +
     &                 t4rr9*scale(9)/rr9
                  t5 = t5rr7*scale(7)/rr7 + t5rr9*scale(9)/rr9 +
     &                 t5rr11*scale(11)/rr11
                  call potik(i,k,t0,t1,t2,poti,potk)
                  call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                  call gradfieldik(i,k,t2,t3,t4,
     &                 gradfieldi,gradfieldk)
                  call hessfieldik(i,k,t3,t4,t5,
     &                 hessfieldi,hessfieldk)
                  pot_ewaldo(i) = pot_ewaldo(i) + poti
                  pot_ewaldo(k) = pot_ewaldo(k) + potk
                  do j = 1, 3
                     field_ewaldo(j,i) = field_ewaldo(j,i) + 
     &                    fieldi(j)
                     field_ewaldo(j,k) = field_ewaldo(j,k) + 
     &                    fieldk(j)
                     do l = 1, 3
                        gradfield_ewaldo(l,j,i) = 
     &                       gradfield_ewaldo(l,j,i) + 
     &                       gradfieldi(l,j)
                        gradfield_ewaldo(l,j,k) = 
     &                       gradfield_ewaldo(l,j,k) +
     &                       gradfieldk(l,j)
                        do h = 1, 3
                           hessfield_ewaldo(h,l,j,i) = 
     &                          hessfield_ewaldo(h,l,j,i) +
     &                          hessfieldi(h,l,j)
                           hessfield_ewaldo(h,l,j,k) =
     &                          hessfield_ewaldo(h,l,j,k) +
     &                          hessfieldk(h,l,j)
                        end do
                     end do
                  end do
               end if
c
c     thole damping
c
               if (damp_thole) then
                  call dampthole(i,k,rorder,r,scale)
                  t0 = t0rr1*scale(1)
                  t1 = t1rr3*scale(3)
                  t2 = t2rr3*scale(3) + t2rr5*scale(5)
                  t2 = t2rr3*scale(3) + t2rr5*scale(5)
                  t3 = t3rr5*scale(5) + t3rr7*scale(7)
                  t4 = t4rr5*scale(5) + t4rr7*scale(7) + 
     &                 t4rr9*scale(9)
                  call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                  call gradfieldik(i,k,t2,t3,t4,gradfieldi,
     &                 gradfieldk)
                  do j = 1, 3
                     fieldd_tholeo(j,i) = fieldd_tholeo(j,i) +
     &                    fieldi(j)*dscale(kk)
                     fieldd_tholeo(j,k) = fieldd_tholeo(j,k) +
     &                    fieldk(j)*dscale(kk)
                     fieldp_tholeo(j,i) = fieldp_tholeo(j,i) +
     &                    fieldi(j)*pscale(kk)
                     fieldp_tholeo(j,k) = fieldp_tholeo(j,k) +
     &                    fieldk(j)*pscale(kk)
                     do l = 1, 3
                        gradfieldd_tholeo(l,j,i) =
     &                       gradfieldd_tholeo(l,j,i) +
     &                       gradfieldi(l,j)*dscale(kk)
                        gradfieldd_tholeo(l,j,k) =
     &                       gradfieldd_tholeo(l,j,k) +
     &                       gradfieldk(l,j)*dscale(kk)
                        gradfieldp_tholeo(l,j,i) =
     &                       gradfieldp_tholeo(l,j,i) +
     &                       gradfieldi(l,j)*pscale(kk)
                        gradfieldp_tholeo(l,j,k) =
     &                       gradfieldp_tholeo(l,j,k) +
     &                       gradfieldk(l,j)*pscale(kk)
                     end do
                  end do
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            dscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            dscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            dscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            dscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
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
      pot = poto
      field = fieldo
      gradfield = gradfieldo
      hessfield = hessfieldo
c
      potm = potmo
      fieldm = fieldmo
      gradfieldm = gradfieldmo
      hessfieldm = hessfieldmo
c
      pot_ewald = pot_ewaldo
      field_ewald = field_ewaldo
      gradfield_ewald = gradfield_ewaldo
      hessfield_ewald = hessfield_ewaldo
c
      fieldd_thole = fieldd_tholeo
      fieldp_thole = fieldp_tholeo
      gradfieldd_thole = gradfieldd_tholeo
      gradfieldp_thole = gradfieldp_tholeo
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (dscale)
      deallocate (pscale)
      deallocate (scale)
      deallocate (poto)
      deallocate (fieldo)
      deallocate (gradfieldo)
      deallocate (hessfieldo)
      deallocate (potmo)
      deallocate (fieldmo)
      deallocate (gradfieldmo)
      deallocate (hessfieldmo)
      deallocate (pot_ewaldo)
      deallocate (field_ewaldo)
      deallocate (gradfield_ewaldo)
      deallocate (hessfield_ewaldo)
      deallocate (fieldd_tholeo)
      deallocate (fieldp_tholeo)
      return
      end
