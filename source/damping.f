c
c     #############################################################
c     ## COPYRIGHT (C) 2016 by Josh Rackers & Jay William Ponder ##
c     ##                   All Rights Reserved                   ## 
c     ############################################################# 
c
c     ####################################################################
c     ##                                                                ##
c     ##  damping  --  set of routines to generate damping coefficients ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "damping" generates the damping coefficients that go with 
c     corresponding powers of r for various damping functions
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine dampewald  --  generate ewald damping coefficients ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "dampewald" generates the damping coefficients for the ewald error
c     function damping that go with corresponding powers of r
c
      subroutine dampewald(i,k,rorder,r,r2,scale)
      use sizes
      use ewald
      use math
      implicit none
      integer i,k,j,maxj
      integer rorder
      real*8 r,r2
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 bfac
      real*8 scale(*)
      real*8, allocatable :: bn(:)
c
c     set damping factors to one
c
      do j = 1, rorder
         scale(j) = 1.0d0
      end do
c
c     set max order needed for ewald damping factors
c
      maxj = (rorder - 1)/2
      allocate (bn(0:maxj))
c     
c     compute ewald damping factors
c
      ralpha = aewald * r
      bn(0) = erfc(ralpha) / r
      scale(1) = bn(0)
      alsq2 = 2.0d0 * aewald**2
      alsq2n = 0.0d0
      if (aewald .gt. 0.0d0) alsq2n = 1.0d0 / (sqrtpi*aewald)
      exp2a = exp(-ralpha**2)
      do j = 1, maxj
         bfac = dble(j+j-1)
         alsq2n = alsq2 * alsq2n
         bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
         scale(2*j+1) = bn(j)
      end do
c
c     deallocate local array
c
      deallocate (bn)
      return
      end
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine dampthole  --  generate thole damping coefficients ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "dampthole" generates the damping coefficients for the thole 
c     damping functional form that go with corresponding powers of r
c
      subroutine dampthole(i,k,rorder,r,scale)
      use sizes
      use polar
      implicit none
      integer i,k,j
      integer rorder
      real*8 r
      real*8 damp,expdamp
      real*8 pdi,pti,pdk,ptk
      real*8 pgamma
      real*8 scale(*)
c     
c     set damping factors to one
c
      do j = 1, rorder
         scale(j) = 1.0d0
      end do
c
c     read in damping parameters
c
      pdi = pdamp(i)
      pti = thole(i)
      pdk = pdamp(k)
      ptk = thole(k)
c
c     assign thole damping scale factors
c
      damp = pdi * pdk
      if (damp .ne. 0.0d0) then
         pgamma = min(pti,ptk)
         damp = pgamma * (r/damp)**3
         if (damp .lt. 50.0d0) then
            expdamp = exp(-damp)
            scale(3) = 1.0d0 - expdamp
            scale(5) = 1.0d0 - (1.0d0 + damp)*expdamp
            if (rorder.ge.7) then
            scale(7) = 1.0d0 - (1.0d0 + damp + 0.6d0*damp**2)*expdamp
            end if
            if (rorder.ge.9) then
               scale(9) = 1.0d0 - (1.0d0 + damp + 
     &              (18.0d0/35.0d0)*damp**2 + (9.0d0/35.0d0)*damp**3)*
     &              expdamp
            end if
         end if
      end if
      return
      end
c
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine dampgordon  --  generate gordon damping coefficents  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "dampgordon" generates the damping coefficients for the one-site 
c     gordon damping functional form that go with corresponding powers of r
c
c     one-site scale factors are used for nuclear-electron interactions
c
      subroutine dampgordon(i,k,rorder,r,scalei,scalek,scaleik)
      use sizes
      use chgpen
      implicit none
      integer i,k
      integer rorder
      real*8 r
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8 scalei(*),scalek(*)
      real*8 scaleik(*)
c
c     read in charge penetration damping parameters
c
      alphai = alpha(cpclass(i))
      alphak = alpha(cpclass(k))
c
c     compute common factors for damping
c
      dampi = alphai*r
      dampk = alphak*r
      expdampi = exp(-dampi)
      expdampk = exp(-dampk)
c
c     calculate one-site scale factors
c
      scalei(1) = 1.0d0 - expdampi
      scalek(1) = 1.0d0 - expdampk
      scalei(3) = 1.0d0 - (1.0d0 + dampi)*expdampi
      scalek(3) = 1.0d0 - (1.0d0 +dampk)*expdampk
      scalei(5) = 1.0d0 - (1.0d0 + dampi + 
     &     (1.0d0/3.0d0)*dampi**2)*expdampi
      scalek(5) = 1.0d0 - (1.0d0 +dampk +
     &     (1.0d0/3.0d0)*dampk**2)*expdampk
      if (rorder .ge. 11) then
cccccccccccccccccccccccc need to put in scale(7) for forces
         scalei(7) = 1.0d0 - (1.0d0 + dampi + 0.4d0*dampi**2 +
     &        (1.0d0/15.0d0)*dampi**3)*expdampi
         scalek(7) = 1.0d0 - (1.0d0 + dampk + 0.4d0*dampk**2 +
     &        (1.0d0/15.0d0)*dampk**3)*expdampk
      end if
c
c     calculate two-site scale factors
c
      if (alphai .ne. alphak) then
         termi = alphak**2/(alphak**2 - alphai**2)
         termk = alphai**2/(alphai**2 - alphak**2)
         scaleik(1) = 1.0d0 -termi*expdampi -termk*expdampk
         scaleik(3) = 1.0d0 - termi*(1.0d0 +dampi)*expdampi
     &        - termk*(1.0d0 + dampk)*expdampk
         scaleik(5) = 1.0d0 - termi*(1.0d0 + dampi +
     &        (1.0d0/3.0d0)*dampi**2)*expdampi -
     &        termk*(1.0d0 + dampk +
     &        (1.0d0/3.0d0)*dampk**2)*expdampk
         if (rorder .ge. 7) then
            scaleik(7) = 1.0d0 - termi*(1.0d0 + dampi +
     &           0.4d0*dampi**2 + (1.0d0/15.0d0)*dampi**3)*
     &           expdampi -
     &           termk*(1.0d0 + dampk +
     &           0.4d0*dampk**2 + (1.0d0/15.0d0)*dampk**3)*
     &           expdampk
         end if
         if (rorder .ge. 9) then
            scaleik(9) = 1.0d0 - termi*(1.0d0 + dampi +
     &           (3.0d0/7.0d0)*dampi**2 +
     &           (2.0d0/21.0d0)*dampi**3 +
     &           (1.0d0/105.0d0)*dampi**4)*expdampi -
     &           termk*(1.0d0 + dampk +
     &           (3.0d0/7.0d0)*dampk**2 +
     &           (2.0d0/21.0d0)*dampk**3 +
     &           (1.0d0/105.0d0)*dampk**4)*expdampk
         end if
         if (rorder .ge. 11) then
            scaleik(11) = 1.0d0 - termi*(1.0d0 + dampi +
     &           (4.0d0/9.0d0)*dampi**2 +
     &           (1.0d0/9.0d0)*dampi**3 +
     &           (1.0d0/63.0d0)*dampi**4 + 
     &           (1.0d0/945.0d0)*dampi**5)*expdampi - 
     &           termk*(1.0d0 + dampk +
     &           (4.0d0/9.0d0)*dampk**2+
     &           (1.0d0/9.0d0)*dampk**3+
     &           (1.0d0/63.0d0)*dampk**4 +
     &           (1.0d0/945.0d0)*dampk**5)*expdampk
         end if
      else
         scaleik(1) = 1.0d0 - (1.0d0+0.5d0*dampi)*expdampi
         scaleik(3) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2)
     &        *expdampi
         scaleik(5) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &        + (1.0d0/6.0d0)*dampi**3)*expdampi
         if (rorder .ge. 7) then
            scaleik(7) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &           + (1.0d0/6.0d0)*dampi**3
     &           + (1.0d0/30.0d0)*dampi**4)*expdampi
         end if
         if (rorder .ge. 9) then
            scaleik(9) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &           + (1.0d0/6.0d0)*dampi**3
     &           + (4.0d0/105.0d0)*dampi**4
     &           + (1.0d0/210.0d0)*dampi**5)*expdampi
         end if
         if (rorder .ge. 11) then
            scaleik(11) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2
     &           + (1.0d0/6.0d0)*dampi**3 
     &           + (5.0d0/126.0d0)*dampi**4
     &           + (2.0d0/315.0d0)*dampi**5
     &           + (1.0d0/1890.0d0)*dampi**6)*expdampi
         end if
      end if
      return
      end
c
c
c     #########################################################################
c     ##                                                                     ##
c     ##  subroutine dampgordonreg  --  generate gordon damping coefficents  ##
c     ##                                                                     ##
c     #########################################################################
c
c
c     "dampgordonreg" generates the damping coefficients for 
c     gordon damping functional form times the exp(-ar) nuclear 
c     regularization that go with corresponding powers of r
c
c     one-site scale factors are used for nuclear-electron interactions
c
      subroutine dampgordonreg(i,k,rorder,r,scalei,scalek)
      use sizes
      use chgpen
      implicit none
      integer i,k
      integer rorder
      real*8 r
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8 regulari,regulark
      real*8 rdampi,rdampk
      real*8 rexpdampi,rexpdampk
      real*8 scalei(*),scalek(*)
c
c     read in charge penetration damping parameters
c
      alphai = alpha(cpclass(i))
      alphak = alpha(cpclass(k))
      regulari = regular(cpclass(i))
      regulark = regular(cpclass(k))
c
c     compute common factors for damping
c
      dampi = alphai*r
      dampk = alphak*r
      expdampi = exp(-dampi)
      expdampk = exp(-dampk)
c
c     common factors for regularization
c
      rdampi = regulari*r
      rdampk = regulark*r
      rexpdampi = exp(-rdampi)
      rexpdampk = exp(-rdampk)
c
c     calculate one-site scale factors
c
      scalei(3) = 1.0d0 - (1.0d0 + dampi)*expdampi -
     &     (1.0d0 + rdampk)*rexpdampk + (1.0d0 + dampi +
     &     rdampk)*expdampi*rexpdampk
      scalek(3) = 1.0d0 - (1.0d0 + dampk)*expdampk -
     &     (1.0d0 + rdampi)*rexpdampi + (1.0d0 + dampk +
     &     rdampi)*expdampk*rexpdampi
      return
      end
