c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine exfer3  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "exfer3" calculates any additional user defined potential
c     contribution and also partitions the energy among the atoms
c
c
      subroutine exfer3
      use sizes
      use action
      use analyz
      use atoms
      use energi
c
      use atomid
      use inform
      use math
      use molcul
      use mpole
      use couple
      use chgpen
      use vdw
      use potent
      use chgpot
      use inter
c      use vdwpot
      implicit none
      integer i,k,j,m
      integer ii,kk
      integer rorder
      real*8 r,r2
      real*8 xr,yr,zr
      real*8 e
      real*8 f
      real*8 radi,radk,radik
      real*8 epsi,epsk,epsik
      real*8 radmix,epsik2
      real*8 termi,termk
      logical proceed
c
c
c     zero out energy and partitioning due to extra potential terms
c
c      nex = 0
      exfr = 0.0d0
c      do i = 1, n
c         aex(i) = 0.0d0
c      end do
      f = electric / dielec
      if (.not. use_vdw) then 
         call kvdw
      end if
c
c     i do not want any hydrogen reductions
c
c      allocate (xred(n))
c      allocate (yred(n))
c      allocate (zred(n))
c      do k = 1, nvdw
c         i = ivdw(k)
c         iv = ired(i)
c         rdn = kred(i)
c         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
c         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
c         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
c      end do
c
c     add any user-defined extra potentials and partitioning
c
      do i = 1, nvdw-1
         ii = ivdw(i)
         do k = i+1, nvdw
            kk = ivdw(k)
            proceed = .true.
            if (molcule(kk) .eq. molcule(ii))  proceed = .false.
            if (proceed) then
c
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
c     
               r2 = xr*xr + yr*yr + zr*zr
               r = sqrt(r2)
c
c     get charge transfer parameters
c
               if (use_vdwclass) then
                  epsi = boverlap(vdwclass(ii))
                  epsk = boverlap(vdwclass(kk))
                  radi = overlapr(vdwclass(ii))
                  radk = overlapr(vdwclass(kk))
               else
                  epsi = boverlap(cpclass(ii))
                  epsk = boverlap(cpclass(kk))
                  radi = overlapr(cpclass(ii))
                  radk = overlapr(cpclass(kk))
               end if
c
c     default combining rules
c
c               epsik = epsi * epsk
               epsik = sqrt(epsi * epsk)
               radik = (radi + radk)/2.0d0
c
c     alternative combining rules
c
               if (exmix .eq. 'ARITHMETIC') then
                  epsik = 0.5d0 * (epsi + epsk)
               else if (exmix .eq. 'GEOMETRIC') then
                  epsik = sqrt(epsi) * sqrt(epsk)
               else if (exmix .eq. 'HARMONIC') then
                  epsik = 2.0d0 * (epsi*epsk)/(epsi+epsk)
               else if (exmix .eq. 'HHG') then
                  epsik = 4.0d0 * (epsi*epsk) / 
     &                 (sqrt(epsi) + sqrt(epsk))**2
               else if (exmix .eq. "W-H") then
                  epsik2 = epsi * epsk
                  epsik = sqrt(epsik2)
               else if (exmix .eq."W-H1") then
                  termi = 2.0d0*(radi**3)*(radk**3)
                  termk = (radi**6) + (radk**6)
                  epsik = sqrt(epsi*epsk)*termi/termk
c     bohm ahlrichs (reference?)
               else if (exmix .eq. "BA") then
                  radmix = 2.0d0 * radi * radk / (radi + radk)
                  termi = epsi**(1.0d0/radi)
                  termk = epsk**(1.0d0/radk)
                  epsik = (termi*termk)**(radmix/2.0d0)
               end if
c
c     compute charge transfer energy
c
               if (xfermodel .eq. "EXP") then
                  e = -epsik * 1000.0d0 * exp(-radik*r)
c
c     increment the total exchange energy
c
                  exfr = exfr + e
c
c     increment the total intermolecular energy
c
                  if (molcule(ii) .ne. molcule(kk)) then
                     einter = einter + e
                  end if
               else
c                  print *,"ERROR: SET EXCH-MODEL    ORBITAL2ONLYR"
               end if
            end if
         end do
      end do
      return
      end
