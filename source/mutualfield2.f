c
c
c     #############################################################
c     ## COPYRIGHT (C) 2016 by Josh Rackers & Jay William Ponder ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mutualfield2  --  electric field gradient     ##
c     ##                               from induce dipoles         ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mutualfield2" computes the electric field gradient 
c     from the induced dipoles
c
c     (the 2 stands for field gradient)
c
c     assumes multipole components have already been rotated into
c     the global coordinate frame
c
c
      subroutine mutualfield2
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
      if (.not.allocated(udgradfield)) allocate (udgradfield(3,3,npole))
      if (.not.allocated(upgradfield)) allocate (upgradfield(3,3,npole))
c
      if (.not.allocated(udgradfield_ewald)) 
     &     allocate (udgradfield_ewald(3,3,npole))
      if (.not.allocated(upgradfield_ewald)) 
     &     allocate (upgradfield_ewald(3,3,npole))
c
      if (.not.allocated(udgradfield_thole)) 
     &     allocate (udgradfield_thole(3,3,npole))
      if (.not.allocated(upgradfield_thole)) 
     &     allocate (upgradfield_thole(3,3,npole))
c
c     get the electrostatic potential, field and field gradient
c     due to the induced dipoles
c
      if (use_mlist) then
         call mutualfield2b
      else
         call mutualfield2a
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mutualfield2a  --  mutual field                ##
c     ##                                gradient via double loop    ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mutualfield2a" computes the mutual 
c     field gradient due to induced dipoles via a double loop
c
c
      subroutine mutualfield2a 
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
      real*8 rr1,rr3,rr5,rr7
      real*8 t3(3,3,3)
      real*8 t3rr5(3,3,3),t3rr7(3,3,3)
      real*8 gradfieldid(3,3),gradfieldkd(3,3)
      real*8 gradfieldip(3,3),gradfieldkp(3,3)
      real*8, allocatable :: scale(:)
      real*8, allocatable :: pscale(:)
      logical proceed
      logical usei,usek
      character*6 mode
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            do k = 1, 3
               udgradfield(k,j,i) = 0.0d0
               upgradfield(k,j,i) = 0.0d0
               udgradfield_ewald(k,j,i) = 0.0d0
               upgradfield_ewald(k,j,i) = 0.0d0
               udgradfield_thole(k,j,i) = 0.0d0
               upgradfield_thole(k,j,i) = 0.0d0
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
      allocate (pscale(n))
c
c     set highest order rr and damping terms needed
c     2 = up to field gradient (rr7 for induced dipoles)
c
      order = 2
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
                  call t3matrixrr5(xr,yr,zr,rr5,t3rr5)
                  call t3matrixrr7(xr,yr,zr,rr7,t3rr7)
c
c     call routines that produce potential, field, field gradient
c     for types of damping
c
                  if (damp_none) then
                     t3 = t3rr5 + t3rr7
                     call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                    gradfieldip,gradfieldkp)
                     do j = 1, 3
                        do l = 1, 3
                           udgradfield(l,j,i) = udgradfield(l,j,i) + 
     &                          gradfieldid(l,j)
                           udgradfield(l,j,k) = udgradfield(l,j,k) + 
     &                          gradfieldkd(l,j)
                           upgradfield(l,j,i) = upgradfield(l,j,i) + 
     &                          gradfieldip(l,j)
                           upgradfield(l,j,k) = upgradfield(l,j,k) + 
     &                          gradfieldkp(l,j)
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
                     t3 = t3rr5*scale(5)/rr5 + t3rr7*scale(7)/rr7
                     call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                    gradfieldip,gradfieldkp)
                     do j = 1, 3
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
                        end do
                     end do
                  end if
c
c     thole damping
c
                  if (damp_thole) then
                     call dampthole(i,k,rorder,r,scale)
                     t3 = t3rr5*scale(5) + t3rr7*scale(7)
                     call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                    gradfieldip,gradfieldkp)
                     do j = 1, 3
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
                        end do
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
c                     call upotik(i,k,t1,potid,potkd,potip,potkp)
c                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
c     &                    fieldkp)
c                     call ugradfieldik(i,k,t3,
c     &                    gradfieldid,gradfieldkd,
c     &                    gradfieldip,gradfieldkp)
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
cc
cc     two center gordon damping
cc
c                  if (damp_gordontwo) then
c                     call dampgordontwo(i,k,rorder,r,scale)
c                     t0 = t0rr1*scale(1)
c                     t1 = t1rr3*scale(3)
c                     t2 = t2rr3*scale(3) + t2rr5*scale(5)
c                     t3 = t3rr5*scale(5) + t3rr7*scale(7)
cc                     t4 = t4rr5*scale(5) + t4rr7*scale(7) +
cc     &                    t4rr9*scale(9)
c                     call upotik(i,k,t1,potid,potkd,potip,potkp)
c                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
c     &                    fieldkp)
c                     call ugradfieldik(i,k,t3,
c     &                    gradfieldid,gradfieldkd,
c     &                    gradfieldip,gradfieldkp)
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
cc                     call damppiquemalone(i,k,rorder,r,scale)
c                     t0 = t0rr1*scale(1)
c                     t1 = t1rr3*scale(3)
c                     t2 = t2rr3*scale(3) + t2rr5*scale(5)
c                     t3 = t3rr5*scale(5) + t3rr7*scale(7)
cc                     t4 = t4rr5*scale(5) + t4rr7*scale(7) +
cc     &                    t4rr9*scale(9)
c                     call upotik(i,k,t1,potid,potkd,potip,potkp)
c                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
c     &                    fieldkp)
c                     call ugradfieldik(i,k,t3,
c     &                    gradfieldid,gradfieldkd,
c     &                    gradfieldip,gradfieldkp)
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
cc                     t4 = t4rr5*scale(5) + t4rr7*scale(7) +
c     &                    t4rr9*scale(9)
c                     call upotik(i,k,t1,potid,potkd,potip,potkp)
c                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
c     &                    fieldkp)
c                     call ugradfieldik(i,k,t3,
c     &                    gradfieldid,gradfieldkd,
c     &                    gradfieldip,gradfieldkp)
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
c                     call t0matrixrr1(rr1,t0rr1)
c                     call t1matrixrr3(xr,yr,zr,rr3,t1rr3)
c                     call t2matrixrr3(xr,yr,zr,rr3,t2rr3)
c                     call t2matrixrr5(xr,yr,zr,rr5,t2rr5)
                     call t3matrixrr5(xr,yr,zr,rr5,t3rr5)
                     call t3matrixrr7(xr,yr,zr,rr7,t3rr7)
c
c     call routines that produce potential, field, field gradient
c     for types of damping
c
c
c     no damping
c
                     if (damp_none) then
                        t3 = t3rr5 + t3rr7
                        call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd
     &                       ,gradfieldip,gradfieldkp)
                        do j = 1, 3
                           do l = 1, 3
                              udgradfield(l,j,i) = udgradfield(l,j,i) + 
     &                             gradfieldid(l,j)
                              udgradfield(l,j,k) = udgradfield(l,j,k) + 
     &                             gradfieldkd(l,j)
                              upgradfield(l,j,i) = upgradfield(l,j,i) + 
     &                             gradfieldip(l,j)
                              upgradfield(l,j,k) = upgradfield(l,j,k) + 
     &                             gradfieldkp(l,j)
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
                        t3 = t3rr5*scale(5)/rr5 + t3rr7*scale(7)/rr7
                        call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd
     &                       ,gradfieldip,gradfieldkp)
                        do j = 1, 3
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
                           end do
                        end do
                     end if
c
c     thole damping
c
                     if (damp_thole) then
                        call dampthole(i,k,rorder,r,scale)
                        t3 = t3rr5*scale(5) + t3rr7*scale(7)
                        call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd
     &                       ,gradfieldip,gradfieldkp)
                        do j = 1, 3
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
c     ##  subroutine mutualfield2b  --  potential, field and field  ##
c     ##                                gradient via neighbor list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mutualfield2b" computes the direct electrostatic potential, field and 
c     field gradient due to induced dipoles using a neighbor list
c
c
      subroutine mutualfield2b
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
      real*8 rr1,rr3,rr5,rr7
      real*8 t3(3,3,3)
      real*8 t3rr5(3,3,3),t3rr7(3,3,3)
      real*8 gradfieldid(3,3),gradfieldkd(3,3)
      real*8 gradfieldip(3,3),gradfieldkp(3,3)
      real*8 test
      real*8, allocatable :: scale(:)
      real*8, allocatable :: udgradfieldo(:,:,:)
      real*8, allocatable :: upgradfieldo(:,:,:)
      real*8, allocatable :: udgradfield_ewaldo(:,:,:)
      real*8, allocatable :: upgradfield_ewaldo(:,:,:)
      real*8, allocatable :: udgradfield_tholeo(:,:,:)
      real*8, allocatable :: upgradfield_tholeo(:,:,:)
      logical proceed
      logical usei,usek
      character*6 mode
c
c     perform dynamic allocation of some local arrays
c
      allocate (udgradfieldo(3,3,npole))
      allocate (upgradfieldo(3,3,npole))
      allocate (udgradfield_ewaldo(3,3,npole))
      allocate (upgradfield_ewaldo(3,3,npole))
      allocate (udgradfield_tholeo(3,3,npole))
      allocate (upgradfield_tholeo(3,3,npole))
c
c     zero out the value of the gradfield at each site
c
      do i = 1, npole
         do j = 1, 3
            do k = 1, 3
               udgradfield(k,j,i) = 0.0d0
               upgradfield(k,j,i) = 0.0d0
               udgradfield_ewald(k,j,i) = 0.0d0
               upgradfield_ewald(k,j,i) = 0.0d0
               udgradfield_thole(k,j,i) = 0.0d0
               upgradfield_thole(k,j,i) = 0.0d0
c
               udgradfieldo(k,j,i) = 0.0d0
               upgradfieldo(k,j,i) = 0.0d0
               udgradfield_ewaldo(k,j,i) = 0.0d0
               upgradfield_ewaldo(k,j,i) = 0.0d0
               udgradfield_tholeo(k,j,i) = 0.0d0
               upgradfield_tholeo(k,j,i) = 0.0d0
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
      rorder = order*2 + 3
      allocate (scale(rorder))
      do i = 1,rorder
          scale(i) = 0.0d0
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared)
!$OMP& private(i,j,k,ii,ix,iy,iz,usei,kk,kx,ky,kz,usek,kkk,proceed,
!$OMP& xr,yr,zr,r,r2,rr1,rr3,rr5,rr7,fgrp,
!$OMP& gradfieldid,gradfieldkd,gradfieldip,gradfieldkp,
!$OMP& t3rr5,t3rr7,
!$OMP& t3,scale)
!$OMP DO reduction(+:udgradfieldo,upgradfieldo,
!$OMP& udgradfield_ewaldo,upgradfield_ewaldo,
!$OMP& udgradfield_tholeo,upgradfield_tholeo)
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
c                  call t0matrixrr1(rr1,t0rr1)
c                  call t1matrixrr3(xr,yr,zr,rr3,t1rr3)
c                  call t2matrixrr3(xr,yr,zr,rr3,t2rr3)
c                  call t2matrixrr5(xr,yr,zr,rr5,t2rr5)
                  call t3matrixrr5(xr,yr,zr,rr5,t3rr5)
                  call t3matrixrr7(xr,yr,zr,rr7,t3rr7)
c
c     call routines that produce potential, field, field gradient
c     for types of damping
c
                  if (damp_none) then
                     t3 = t3rr5 + t3rr7
                     call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                    gradfieldip,gradfieldkp)
                     do j = 1, 3
                        do l = 1, 3
                           udgradfieldo(l,j,i) = udgradfieldo(l,j,i) + 
     &                          gradfieldid(l,j)
                           udgradfieldo(l,j,k) = udgradfieldo(l,j,k) + 
     &                          gradfieldkd(l,j)
                           upgradfieldo(l,j,i) = upgradfieldo(l,j,i) + 
     &                          gradfieldip(l,j)
                           upgradfieldo(l,j,k) = upgradfieldo(l,j,k) + 
     &                          gradfieldkp(l,j)
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
                     t3 = t3rr5*scale(5)/rr5 + t3rr7*scale(7)/rr7
                     call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                    gradfieldip,gradfieldkp)
                     do j = 1, 3
                        do l = 1, 3
                           udgradfield_ewaldo(l,j,i) = 
     &                          udgradfield_ewaldo(l,j,i) + 
     &                          gradfieldid(l,j)
                           udgradfield_ewaldo(l,j,k) = 
     &                          udgradfield_ewaldo(l,j,k) + 
     &                          gradfieldkd(l,j)
                           upgradfield_ewaldo(l,j,i) = 
     &                          upgradfield_ewaldo(l,j,i) +
     &                          gradfieldip(l,j)
                           upgradfield_ewaldo(l,j,k) = 
     &                          upgradfield_ewaldo(l,j,k) +
     &                          gradfieldkp(l,j)
                        end do
                     end do
                  end if
c
c     thole damping
c
                  if (damp_thole) then
                     call dampthole(i,k,rorder,r,scale)
                     t3 = t3rr5*scale(5) + t3rr7*scale(7)
                     call ugradfieldik(i,k,t3,gradfieldid,gradfieldkd,
     &                    gradfieldip,gradfieldkp)
                     do j = 1, 3
                        do l = 1, 3
                           udgradfield_tholeo(l,j,i) = 
     &                          udgradfield_tholeo(l,j,i) +
     &                          gradfieldid(l,j)
                           udgradfield_tholeo(l,j,k) = 
     &                          udgradfield_tholeo(l,j,k) +
     &                          gradfieldkd(l,j)
                           upgradfield_tholeo(l,j,i) = 
     &                          upgradfield_tholeo(l,j,i) +
     &                          gradfieldip(l,j)
                           upgradfield_tholeo(l,j,k) = 
     &                          upgradfield_tholeo(l,j,k) +
     &                          gradfieldkp(l,j)
                        end do
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
      udgradfield = udgradfieldo
      upgradfield = upgradfieldo
      udgradfield_ewald = udgradfield_ewaldo
      upgradfield_ewald = upgradfield_ewaldo
      udgradfield_thole = udgradfield_tholeo
      upgradfield_thole = upgradfield_tholeo
c
c     perform deallocation of some local arrays
c
      deallocate (scale)
      deallocate (udgradfieldo)
      deallocate (upgradfieldo)
      deallocate (udgradfield_ewaldo)
      deallocate (upgradfield_ewaldo)
      deallocate (udgradfield_tholeo)
      deallocate (upgradfield_tholeo)
      return
      end
