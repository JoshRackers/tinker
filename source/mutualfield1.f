c
c
c     #############################################################
c     ## COPYRIGHT (C) 2016 by Josh Rackers & Jay William Ponder ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mutualfield1  --  electric field from induced ##
c     ##                               dipoles                     ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mutualfield1" computes the electric field
c     from the induced dipoles
c
c     (the 1 stands for field)
c
c     assumes multipole components have already been rotated into
c     the global coordinate frame
c
c
      subroutine mutualfield1
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
c
      if (.not.allocated(udfield_ewald)) 
     &     allocate (udfield_ewald(3,npole))
      if (.not.allocated(upfield_ewald)) 
     &     allocate (upfield_ewald(3,npole))
c
      if (.not.allocated(udfield_thole)) 
     &     allocate (udfield_thole(3,npole))
      if (.not.allocated(upfield_thole)) 
     &     allocate (upfield_thole(3,npole))
c
c     get the electrostatic potential, field and field gradient
c     due to the induced dipoles
c
      if (use_mlist) then
         call mutualfield1b
      else
         call mutualfield1a
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mutualfield1a  --  mutual field via a double   ##
c     ##                                loop                        ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mutualfield1a" computes the mutual electrostatic field 
c     due to induced dipoles via a double loop
c
c
      subroutine mutualfield1a 
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
      real*8 t0,t1(3),t2(3,3),t3(3,3,3)
      real*8 t0rr1,t1rr3(3)
      real*8 t2rr3(3,3),t2rr5(3,3)
      real*8 t3rr5(3,3,3),t3rr7(3,3,3)
      real*8 potid,potkd,potip,potkp
      real*8 fieldid(3),fieldkd(3)
      real*8 fieldip(3),fieldkp(3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: scale(:)
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
c     1 = up to field (rr5 for induced dipoles)
c
      order = 1
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
                  call t2matrixrr3(xr,yr,zr,rr3,t2rr3)
                  call t2matrixrr5(xr,yr,zr,rr5,t2rr5)
c
c     call routines that produce potential, field, field gradient
c     for types of damping
c
                  if (damp_none) then
                     t2 = t2rr3 + t2rr5
                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                    fieldkp)
                     do j = 1, 3
                        udfield(j,i) = udfield(j,i) + fieldid(j)
                        udfield(j,k) = udfield(j,k) + fieldkd(j)
                        upfield(j,i) = upfield(j,i) + fieldip(j)
                        upfield(j,k) = upfield(j,k) + fieldkp(j)
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
                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                    fieldkp)
                     do j = 1, 3
                        udfield_ewald(j,i) = udfield_ewald(j,i) + 
     &                       fieldid(j)
                        udfield_ewald(j,k) = udfield_ewald(j,k) + 
     &                       fieldkd(j)
                        upfield_ewald(j,i) = upfield_ewald(j,i) +
     &                       fieldip(j)
                        upfield_ewald(j,k) = upfield_ewald(j,k) +
     &                       fieldkp(j)
                     end do
                  end if
c
c     thole damping
c
                  if (damp_thole) then
                     call dampthole(i,k,rorder,r,scale)
                     t2 = t2rr3*scale(3) + t2rr5*scale(5)
                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                    fieldkp)
                     do j = 1, 3
                        udfield_thole(j,i) = udfield_thole(j,i) +
     &                       fieldid(j)
                        udfield_thole(j,k) = udfield_thole(j,k) +
     &                       fieldkd(j)
                        upfield_thole(j,i) = upfield_thole(j,i) +
     &                       fieldip(j)
                        upfield_thole(j,k) = upfield_thole(j,k) +
     &                       fieldkp(j)
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
                     call t0matrixrr1(rr1,t0rr1)
                     call t1matrixrr3(xr,yr,zr,rr3,t1rr3)
                     call t2matrixrr3(xr,yr,zr,rr3,t2rr3)
                     call t2matrixrr5(xr,yr,zr,rr5,t2rr5)
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
                        t2 = t2rr3 + t2rr5
                        call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                       fieldkp)
                        do j = 1, 3
                           udfield(j,i) = udfield(j,i) + fieldid(j)
                           udfield(j,k) = udfield(j,k) + fieldkd(j)
                           upfield(j,i) = upfield(j,i) + fieldip(j)
                           upfield(j,k) = upfield(j,k) + fieldkp(j)
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
                        call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                       fieldkp)
                        do j = 1, 3
                           udfield_ewald(j,i) = udfield_ewald(j,i) + 
     &                          fieldid(j)
                           udfield_ewald(j,k) = udfield_ewald(j,k) + 
     &                          fieldkd(j)
                           upfield_ewald(j,i) = upfield_ewald(j,i) +
     &                          fieldip(j)
                           upfield_ewald(j,k) = upfield_ewald(j,k) +
     &                          fieldkp(j)
                        end do
                     end if
c
c     thole damping
c
                     if (damp_thole) then
                        call dampthole(i,k,rorder,r,scale)
                        t2 = t2rr3*scale(3) + t2rr5*scale(5)
                        call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                       fieldkp)
                        do j = 1, 3
                           udfield_thole(j,i) = udfield_thole(j,i) +
     &                          fieldid(j)
                           udfield_thole(j,k) = udfield_thole(j,k) +
     &                          fieldkd(j)
                           upfield_thole(j,i) = upfield_thole(j,i) +
     &                          fieldip(j)
                           upfield_thole(j,k) = upfield_thole(j,k) +
     &                          fieldkp(j)
                        end do
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mutualfield1b  --  mutual electric field       ##
c     ##                                via neighbor list           ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mutualfield1b" computes the electric field 
c     due to induced dipoles using a neighbor list
c
c
      subroutine mutualfield1b
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
      real*8 potid,potkd,potip,potkp
      real*8 fieldid(3),fieldkd(3)
      real*8 fieldip(3),fieldkp(3)
      real*8 gradfieldid(3,3),gradfieldkd(3,3)
      real*8 gradfieldip(3,3),gradfieldkp(3,3)
      real*8 test
      real*8, allocatable :: scale(:)
      real*8, allocatable :: udfieldo(:,:)
      real*8, allocatable :: upfieldo(:,:)
      real*8, allocatable :: udfield_ewaldo(:,:)
      real*8, allocatable :: upfield_ewaldo(:,:)
      real*8, allocatable :: udfield_tholeo(:,:)
      real*8, allocatable :: upfield_tholeo(:,:)
      logical proceed
      logical usei,usek
      character*6 mode
c
c     perform dynamic allocation of some local arrays
c
      allocate (udfieldo(3,npole))
      allocate (upfieldo(3,npole))
      allocate (udfield_ewaldo(3,npole))
      allocate (upfield_ewaldo(3,npole))
      allocate (udfield_tholeo(3,npole))
      allocate (upfield_tholeo(3,npole))
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
c
            udfieldo(j,i) = 0.0d0
            upfieldo(j,i) = 0.0d0
            udfield_ewaldo(j,i) = 0.0d0
            upfield_ewaldo(j,i) = 0.0d0
            udfield_tholeo(j,i) = 0.0d0
            upfield_tholeo(j,i) = 0.0d0
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
c      order = 1
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
!$OMP& fieldid,fieldkd,fieldip,fieldkp,
!$OMP& t0rr1,t1rr3,t2rr3,t2rr5,t3rr5,t3rr7,
!$OMP& t0,t1,t2,t3,scale)
!$OMP DO reduction(+:udfieldo,upfieldo,
!$OMP& udfield_ewaldo,upfield_ewaldo,
!$OMP& udfield_tholeo,upfield_tholeo)
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
                  call t2matrixrr3(xr,yr,zr,rr3,t2rr3)
                  call t2matrixrr5(xr,yr,zr,rr5,t2rr5)
c                  call t3matrixrr5(xr,yr,zr,rr5,t3rr5)
c                  call t3matrixrr7(xr,yr,zr,rr7,t3rr7)
c
c     call routines that produce potential, field, field gradient
c     for types of damping
c
                  if (damp_none) then
                     t2 = t2rr3 + t2rr5
                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                    fieldkp)
                     do j = 1, 3
                        udfieldo(j,i) = udfieldo(j,i) + fieldid(j)
                        udfieldo(j,k) = udfieldo(j,k) + fieldkd(j)
                        upfieldo(j,i) = upfieldo(j,i) + fieldip(j)
                        upfieldo(j,k) = upfieldo(j,k) + fieldkp(j)
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
                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                    fieldkp)
                     do j = 1, 3
                        udfield_ewaldo(j,i) = udfield_ewaldo(j,i) + 
     &                       fieldid(j)
                        udfield_ewaldo(j,k) = udfield_ewaldo(j,k) + 
     &                       fieldkd(j)
                        upfield_ewaldo(j,i) = upfield_ewaldo(j,i) +
     &                       fieldip(j)
                        upfield_ewaldo(j,k) = upfield_ewaldo(j,k) +
     &                       fieldkp(j)
                     end do
                  end if
c
c     thole damping
c
                  if (damp_thole) then
                     call dampthole(i,k,rorder,r,scale)
                     t2 = t2rr3*scale(3) + t2rr5*scale(5)
                     call ufieldik(i,k,t2,fieldid,fieldkd,fieldip,
     &                    fieldkp)
                     do j = 1, 3
                        udfield_tholeo(j,i) = udfield_tholeo(j,i) +
     &                       fieldid(j)
                        udfield_tholeo(j,k) = udfield_tholeo(j,k) +
     &                       fieldkd(j)
                        upfield_tholeo(j,i) = upfield_tholeo(j,i) +
     &                       fieldip(j)
                        upfield_tholeo(j,k) = upfield_tholeo(j,k) +
     &                       fieldkp(j)
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
      udfield_ewald = udfield_ewaldo
      upfield_ewald = upfield_ewaldo
      udfield_thole = udfield_tholeo
      upfield_thole = upfield_tholeo
c
c     perform deallocation of some local arrays
c
      deallocate (scale)
      deallocate (udfieldo)
      deallocate (upfieldo)
      deallocate (udfield_ewaldo)
      deallocate (upfield_ewaldo)
      deallocate (udfield_tholeo)
      deallocate (upfield_tholeo)
      return
      end
