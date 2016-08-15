c
c
c     #############################################################
c     ## COPYRIGHT (C) 2016 by Josh Rackers & Jay William Ponder ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine permfield2  --  electric potential, field,  ##
c     ##  field gradient and higher from permanent moments      ##
c     ##                                                        ##
c     ############################################################
c
c
c     "permfield2" computes the electric potential, field
c     and field gradient for use in energy routines as well as 
c     for use in computing induced dipoles
c
c     (the 2 stands for up to field gradient)
c
c     assumes multipole components have already been rotated into
c     the global coordinate frame
c
c
      subroutine permfield2
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
c
c     permanent - permanent exclusion rule
c
      if (.not.allocated(potm)) allocate (potm(npole))
      if (.not.allocated(fieldm)) allocate (fieldm(3,npole))
      if (.not.allocated(gradfieldm)) allocate (gradfieldm(3,3,npole))
c
c     ewald damping
c
      if (.not.allocated(pot_ewald)) allocate (pot_ewald(npole))
      if (.not.allocated(field_ewald)) allocate (field_ewald(3,npole))
      if (.not.allocated(gradfield_ewald)) 
     &     allocate (gradfield_ewald(3,3,npole))
c
c     thole damping with p and d exclusion rules needed for
c     computing induced dipoles
c
      if (.not.allocated(fieldd_thole)) allocate (fieldd_thole(3,npole))
      if (.not.allocated(fieldp_thole)) allocate (fieldp_thole(3,npole))
c
c     get the electrostatic potential, field and field gradient
c     due to permanent multipoles
c
      if (use_mlist) then
         call permfield2b
      else
         call permfield2a
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine permfield2a  --  potential, field and field  ##
c     ##                              gradient via double loop    ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "permfield2a" computes the direct electrostatic potential, field and 
c     field gradient due to permanent multipole moments via a double loop
c
c
      subroutine permfield2a 
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
      integer i,j,k,l,m
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      integer order,rorder
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr1,rr3,rr5,rr7,rr9,rr11
      real*8 t0,t1(3),t2(3,3),t3(3,3,3),t4(3,3,3,3)
      real*8 t0rr1,t1rr3(3)
      real*8 t2rr3(3,3),t2rr5(3,3)
      real*8 t3rr5(3,3,3),t3rr7(3,3,3)
      real*8 t4rr5(3,3,3,3),t4rr7(3,3,3,3),t4rr9(3,3,3,3)
      real*8 poti,potk
      real*8 fieldi(3),fieldk(3)
      real*8 gradfieldi(3,3),gradfieldk(3,3)
      real*8, allocatable :: scale(:)
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
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldm(j,i) = 0.0d0
            field_ewald(j,i) = 0.0d0
            fieldd_thole(j,i) = 0.0d0
            fieldp_thole(j,i) = 0.0d0
            do k = 1, 3
               gradfield(k,j,i) = 0.0d0
               gradfieldm(k,j,i) = 0.0d0
               gradfield_ewald(k,j,i) = 0.0d0
            end do
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set highest order rr and damping terms needed
c     2 = up to field gradient (rr9)
c
      order = 2
      rorder = order*2 + 5
      allocate (scale(rorder))
      do i = 1, rorder
         scale(i) = 0.0d0
      end do
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
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
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
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
c
c     call routines that produce potential, field, field gradient
c     for types of damping
c
                  print *,"doing damping"
                  if (damp_none) then
                     t0 = t0rr1
                     t1 = t1rr3
                     t2 = t2rr3 + t2rr5
                     t3 = t3rr5 + t3rr7
                     t4 = t4rr5 + t4rr7 + t4rr9
                     call potik(i,k,t0,t1,t2,poti,potk)
                     call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                     call gradfieldik(i,k,t2,t3,t4,
     &                    gradfieldi,gradfieldk)
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
                     call potik(i,k,t0,t1,t2,poti,potk)
                     call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                     call gradfieldik(i,k,t2,t3,t4,
     &                    gradfieldi,gradfieldk)
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
                     do j = 1, 3
                        fieldd_thole(j,i) = fieldd_thole(j,i) +
     &                       fieldi(j)*dscale(kk)
                        fieldd_thole(j,k) = fieldd_thole(j,k) +
     &                       fieldk(j)*dscale(kk)
                        fieldp_thole(j,i) = fieldp_thole(j,i) +
     &                       fieldi(j)*pscale(kk)
                        fieldp_thole(j,k) = fieldp_thole(j,k) +
     &                       fieldk(j)*pscale(kk)
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
c     
c     calculate interaction with other unit cells
c
      do i = 1, npole
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         ii = ipole(i)
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
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
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
                     call t0matrixrr1(rr1,t0rr1)
                     call t1matrixrr3(xr,yr,zr,rr3,t1rr3)
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
                        t0 = t0rr1
                        t1 = t1rr3
                        t2 = t2rr3 + t2rr5
                        t3 = t3rr5 + t3rr7
                        t4 = t4rr5 + t4rr7 + t4rr9
                        call potik(i,k,t0,t1,t2,poti,potk)
                        call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                        call gradfieldik(i,k,t2,t3,t4,
     &                       gradfieldi,gradfieldk)
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
                        call potik(i,k,t0,t1,t2,poti,potk)
                        call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                        call gradfieldik(i,k,t2,t3,t4,
     &                       gradfieldi,gradfieldk)
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
                        do j = 1, 3
                           fieldd_thole(j,i) = fieldd_thole(j,i) +
     &                          fieldi(j)*dscale(kk)
                           fieldd_thole(j,k) = fieldd_thole(j,k) +
     &                          fieldk(j)*dscale(kk)
                           fieldp_thole(j,i) = fieldp_thole(j,i) +
     &                          fieldi(j)*pscale(kk)
                           fieldp_thole(j,k) = fieldp_thole(j,k) +
     &                          fieldk(j)*pscale(kk)
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
c     ##  subroutine permfield2b  --  potential, field and field  ##
c     ##                              gradient via neighbor list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "permfield2b" computes the direct electrostatic potential, field and 
c     field gradient due to permanent multipole moments using a neighbor list
c
c
      subroutine permfield2b
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
      integer i,j,k,l,m
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      integer order,rorder
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr1,rr3,rr5,rr7,rr9,rr11
      real*8 t0,t1(3),t2(3,3),t3(3,3,3),t4(3,3,3,3)
      real*8 t0rr1,t1rr3(3)
      real*8 t2rr3(3,3),t2rr5(3,3)
      real*8 t3rr5(3,3,3),t3rr7(3,3,3)
      real*8 t4rr5(3,3,3,3),t4rr7(3,3,3,3),t4rr9(3,3,3,3)
      real*8 poti,potk
      real*8 fieldi(3),fieldk(3)
      real*8 gradfieldi(3,3),gradfieldk(3,3)
      real*8, allocatable :: scale(:)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: poto(:)
      real*8, allocatable :: fieldo(:,:)
      real*8, allocatable :: gradfieldo(:,:,:)
      real*8, allocatable :: potmo(:)
      real*8, allocatable :: fieldmo(:,:)
      real*8, allocatable :: gradfieldmo(:,:,:)
      real*8, allocatable :: pot_ewaldo(:)
      real*8, allocatable :: field_ewaldo(:,:)
      real*8, allocatable :: gradfield_ewaldo(:,:,:)
      real*8, allocatable :: fieldp_tholeo(:,:)
      real*8, allocatable :: fieldd_tholeo(:,:)
      logical proceed
      logical usei,usek
      character*6 mode
c
c     perform dynamic allocation of some local arrays
c
      allocate (poto(npole))
      allocate (fieldo(3,npole))
      allocate (gradfieldo(3,3,npole))
c
      allocate (potmo(npole))
      allocate (fieldmo(3,npole))
      allocate (gradfieldmo(3,3,npole))
c
      allocate (pot_ewaldo(npole))
      allocate (field_ewaldo(3,npole))
      allocate (gradfield_ewaldo(3,3,npole))
c
      allocate (fieldd_tholeo(3,npole))
      allocate (fieldp_tholeo(3,npole))
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
c
               gradfieldo(k,j,i) = 0.0d0
               gradfieldmo(k,j,i) = 0.0d0
               gradfield_ewaldo(k,j,i) = 0.0d0
            end do
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     set highest order rr and damping terms needed
c     2 = up to field gradient (rr9)
c
      order = 2
      rorder = order*2 + 5
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
!$OMP& xr,yr,zr,r,r2,rr1,rr3,rr5,rr7,rr9,fgrp,
!$OMP& poti,potk,fieldi,fieldk,gradfieldi,gradfieldk,
!$OMP& t0rr1,t1rr3,t2rr3,t2rr5,t3rr5,t3rr7,t4rr5,t4rr7,t4rr9,
!$OMP& t0,t1,t2,t3,t4,scale)
!$OMP& firstprivate(mscale,dscale,pscale)
!$OMP DO reduction(+:poto,fieldo,gradfieldo,
!$OMP& potmo,fieldmo,gradfieldmo,
!$OMP& pot_ewaldo,field_ewaldo,gradfield_ewaldo,
!$OMP& fieldd_tholeo,fieldp_tholeo)
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
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
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
                  call t0matrixrr1(rr1,t0rr1)
                  call t1matrixrr3(xr,yr,zr,rr3,t1rr3)
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
                     t0 = t0rr1
                     t1 = t1rr3
                     t2 = t2rr3 + t2rr5
                     t3 = t3rr5 + t3rr7
                     t4 = t4rr5 + t4rr7 + t4rr9
                     call potik(i,k,t0,t1,t2,poti,potk)
                     call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                     call gradfieldik(i,k,t2,t3,t4,
     &                    gradfieldi,gradfieldk)
                     poto(i) = poto(i) + poti
                     poto(k) = poto(k) + potk
                     potmo(i) = potmo(i) + poti*mscale(kk)
                     potmo(k) = potmo(k) + potk*mscale(kk)
                     do j = 1, 3
                        fieldo(j,i)= fieldo(j,i)+ fieldi(j)
                        fieldo(j,k)= fieldo(j,k)+ fieldk(j)
                        fieldmo(j,i)= fieldmo(j,i)+ fieldi(j)*mscale(kk)
                        fieldmo(j,k)= fieldmo(j,k)+ fieldk(j)*mscale(kk)
                        do l = 1, 3
                           gradfieldo(l,j,i) = gradfieldo(l,j,i) +
     &                          gradfieldi(l,j)
                           gradfieldo(l,j,k) = gradfieldo(l,j,k) +
     &                          gradfieldk(l,j)
                           gradfieldmo(l,j,i) = gradfieldmo(l,j,i) +
     &                          gradfieldi(l,j)*mscale(kk)
                           gradfieldmo(l,j,k) = gradfieldmo(l,j,k) +
     &                          gradfieldk(l,j)*mscale(kk)
                        end do
                     end do
                  end if
c     
c     error function damping for ewald
c     
                  if (damp_ewald) then
                     call dampewald(i,k,rorder,r,r2,scale)
                     t0 = t0rr1*scale(1)/rr1
                     t1 = t1rr3*scale(3)/rr3
                     t2 = t2rr3*scale(3)/rr3 + t2rr5*scale(5)/rr5
                     t3 = t3rr5*scale(5)/rr5 + t3rr7*scale(7)/rr7
                     t4 = t4rr5*scale(5)/rr5 + t4rr7*scale(7)/rr7 +
     &                    t4rr9*scale(9)/rr9
                     call potik(i,k,t0,t1,t2,poti,potk)
                     call fieldik(i,k,t1,t2,t3,fieldi,fieldk)
                     call gradfieldik(i,k,t2,t3,t4,
     &                    gradfieldi,gradfieldk)
                     pot_ewaldo(i) = pot_ewaldo(i) + poti
                     pot_ewaldo(k) = pot_ewaldo(k) + potk
                     do j = 1, 3
                        field_ewaldo(j,i) = field_ewaldo(j,i) + 
     &                       fieldi(j)
                        field_ewaldo(j,k) = field_ewaldo(j,k) + 
     &                       fieldk(j)
                        do l = 1, 3
                           gradfield_ewaldo(l,j,i) = 
     &                          gradfield_ewaldo(l,j,i) + 
     &                          gradfieldi(l,j)
                           gradfield_ewaldo(l,j,k) = 
     &                          gradfield_ewaldo(l,j,k) +
     &                          gradfieldk(l,j)
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
                     do j = 1, 3
                        fieldd_tholeo(j,i) = fieldd_tholeo(j,i) +
     &                       fieldi(j)*dscale(kk)
                        fieldd_tholeo(j,k) = fieldd_tholeo(j,k) +
     &                       fieldk(j)*dscale(kk)
                        fieldp_tholeo(j,i) = fieldp_tholeo(j,i) +
     &                       fieldi(j)*pscale(kk)
                        fieldp_tholeo(j,k) = fieldp_tholeo(j,k) +
     &                       fieldk(j)*pscale(kk)
                     end do
                  end if
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
c
      potm = potmo
      fieldm = fieldmo
      gradfieldm = gradfieldmo
c
      pot_ewald = pot_ewaldo
      field_ewald = field_ewaldo
      gradfield_ewald = gradfield_ewaldo
c
      fieldd_thole = fieldd_tholeo
      fieldp_thole = fieldp_tholeo
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
      deallocate (potmo)
      deallocate (fieldmo)
      deallocate (gradfieldmo)
      deallocate (pot_ewaldo)
      deallocate (field_ewaldo)
      deallocate (gradfield_ewaldo)
      deallocate (fieldd_tholeo)
      deallocate (fieldp_tholeo)
      return
      end
