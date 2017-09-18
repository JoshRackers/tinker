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
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine dampfunc  --  generate func damping coefficients ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "dampfunc" generates the damping coefficients for the thole 
c     damping functional form that go with corresponding powers of r
c
      subroutine dampfunc(i,k,rorder,r,scale)
      use sizes
      use chgpen
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
      pdk = pdamp(k)
      pti = func10(cpclass(i))
      ptk = func10(cpclass(k))
c
c     assign thole damping scale factors
c
      damp = pdi * pdk
      if (damp .ne. 0.0d0) then
         pgamma = min(pti,ptk)
         damp   = pgamma*(r/(pdi*pdk))**(3.0d0/2.0d0)
         if (damp .lt. 50.0d0) then
            expdamp = exp(-damp)
            scale(3) = 1.0d0 - expdamp
            scale(5) = 1.0d0 - (1.0d0 + (1.0d0/2.0d0)*damp)*expdamp
            if (rorder.ge.7) then
               scale(7) = 1.0d0 - (1.0d0 + (39.0d0/60.0d0)*damp +
     &                   (9.0d0/60.0d0)*damp**2.0d0)*expdamp
            end if
            if (rorder.ge.9) then
               scale(9) = 1.0d0 - (1.0d0 + (609.0d0/840.0d0)*damp + 
     &                 (189.0d0/840.0d0)*damp**2.0d0 + 
     &                 (27.0d0/840.0d0)*damp**3.0d0 )*expdamp
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
      if (use_vdwclass) then
         alphai = alpha(vdwclass(i))
         alphak = alpha(vdwclass(k))
      end if
c
c     compute common factors for damping
c
      dampi = alphai*r
      dampk = alphak*r
      expdampi = exp(-dampi)
      expdampk = exp(-dampk)
c
c     "new" Gordon damping model - model 1
c
      if (gtype .eq. "NEW") then
         scalei(1) = 1.0d0 - (1.0d0 + dampi/2.0d0)*expdampi
         scalek(1) = 1.0d0 - (1.0d0 + dampk/2.0d0)*expdampk
         scalei(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2)*expdampi
         scalek(3) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk**2)*expdampk
         scalei(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &        (1.0d0/6.0d0)*dampi**3)*expdampi
         scalek(5) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk**2 + 
     &        (1.0d0/6.0d0)*dampk**3)*expdampk
         if (rorder .ge. 11) then
            scalei(7) = 100000.0d0
            scalek(7) = 100000.0d0
         end if
         if (abs(alphai - alphak) .ge. 0.0001d0) then
c         if (alphai .ne. alphak) then
            termi = alphak**2/(alphak**2 - alphai**2)
            termk = alphai**2/(alphai**2 - alphak**2)
            scaleik(1) = 1.0d0 - (termi**2)*
     &           (1.0d0 + 2.0d0*termk + 0.5d0*dampi)*expdampi -
     &           (termk**2)*(1.0d0 + 2.0d0*termi + 0.5d0*dampk)*
     &           expdampk
            scaleik(3) = 1.0d0 - (termi**2)*(1.0d0 + dampi + 
     &           0.5d0*dampi**2)*expdampi - 
     &           (termk**2)*(1.0d0 + dampk + 0.5d0*dampk**2)*expdampk
     &           - 2.0d0*(termi**2)*termk*(1.0d0 + dampi)*expdampi
     &           - 2.0d0*(termk**2)*termi*(1.0d0 + dampk)*expdampk
            scaleik(5) = 1.0d0 - (termi**2)*
     &           (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &           (1.0d0/6.0d0)*dampi**3)*expdampi - (termk**2)*
     &           (1.0d0 + dampk + 0.5d0*dampk**2 + 
     &           (1.0d0/6.0d0)*dampk**3)*expdampk - 
     &           2.0d0*(termi**2)*termk*(1.0 + dampi + 
     &           (1.0d0/3.0d0)*dampi**2)*expdampi - 
     &           2.0d0*(termk**2)*termi*(1.0 + dampk +
     &           (1.0d0/3.0d0)*dampk**2)*expdampk
            if (rorder .ge. 7) then
               scaleik(7) = 1.0d0 - (termi**2)*
     &              (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &              (1.0d0/6.0d0)*dampi**3 + (1.0d0/30.0d0)*dampi**4)*
     &              expdampi - 
     &              (termk**2)*(1.0d0 + dampk + 0.5d0*dampk**2 + 
     &              (1.0d0/6.0d0)*dampk**3 + (1.0d0/30.0d0)*dampk**4)*
     &              expdampk - 
     &              2.0d0*(termi**2)*termk*(1.0d0 + dampi +
     &              (2.0d0/5.0d0)*dampi**2 + (1.0d0/15.0d0)*dampi**3)*
     &              expdampi - 
     &              2.0d0*(termk**2)*termi*(1.0d0 + dampk + 
     &              (2.0d0/5.0d0)*dampk**2 + (1.0d0/15.0d0)*dampk**3)*
     &              expdampk
            end if
            if (rorder .ge. 9) then
               scaleik(9) = 1.0d0 - (termi**2)*
     &              (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &              (1.0d0/6.0d0)*dampi**3 + (4.0d0/105.0d0)*dampi**4 +
     &              (1.0d0/210.0d0)*dampi**5)*expdampi - (termk**2)*
     &              (1.0d0 + dampk + 0.5d0*dampk**2 +
     &              (1.0d0/6.0d0)*dampk**3 + (4.0d0/105.0d0)*dampk**4 +
     &              (1.0d0/210.0d0)*dampk**5)*expdampk -
     &              2.0d0*(termi**2)*termk*
     &              (1.0d0 + dampi + (3.0d0/7.0d0)*dampi**2 + 
     &              (2.0d0/21.0d0)*dampi**3 + (1.0d0/105.0d0)*dampi**4)*
     &              expdampi - 
     &              2.0d0*(termk**2)*termi*
     &              (1.0d0 + dampk + (3.0d0/7.0d0)*dampk**2 +
     &              (2.0d0/21.0d0)*dampk**3 + (1.0d0/105.0d0)*dampk**4)*
     &              expdampk
            end if
         else
            scaleik(1) = 1.0d0 - (1.0d0 + (11.0d0/16.0d0)*dampi + 
     &           (3.0d0/16.0d0)*dampi**2 + (1.0d0/48.0d0)*dampi**3)
     &           *expdampi
            scaleik(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &           (7.0d0/48.0d0)*dampi**3 + (1.0d0/48.0d0)*dampi**4)
     &           *expdampi
            scaleik(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &           (1.0d0/6.0d0)*dampi**3 + (1.0d0/24.0d0)*dampi**4 +
     &           (1.0d0/144.0d0)*dampi**5)*expdampi
            if (rorder .ge. 7) then
               scaleik(7) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &              (1.0d0/6.0d0)*dampi**3 + (1.0d0/24.0d0)*dampi**4 + 
     &              (1.0d0/120.0d0)*dampi**5 + (1.0d0/720.0d0)*dampi**6)
     &              *expdampi
            end if
            if (rorder .ge. 9) then
               scaleik(9) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &              (1.0d0/6.0d0)*dampi**3 + (1.0d0/24.0d0)*dampi**4 + 
     &              (1.0d0/120.0d0)*dampi**5 + (1.0d0/720.0d0)*dampi**6+ 
     &              (1.0d0/5040.0d0)*dampi**7)*expdampi
            end if
         end if
      else
c
c     "original" Gordon damping model - model 2
c
c
c     calculate one-site scale factors
c
         scalei(1) = 1.0d0 - expdampi
         scalek(1) = 1.0d0 - expdampk
         scalei(3) = 1.0d0 - (1.0d0 + dampi)*expdampi
         scalek(3) = 1.0d0 - (1.0d0 +dampk)*expdampk
         scalei(5) = 1.0d0 - (1.0d0 + dampi + 
     &        (1.0d0/3.0d0)*dampi**2)*expdampi
         scalek(5) = 1.0d0 - (1.0d0 +dampk +
     &        (1.0d0/3.0d0)*dampk**2)*expdampk
         if (rorder .ge. 11) then
ccccccccccccccccccccccccneed to put in scale(7) for forces
            scalei(7) = 1.0d0 - (1.0d0 + dampi + 0.4d0*dampi**2 +
     &           (1.0d0/15.0d0)*dampi**3)*expdampi
            scalek(7) = 1.0d0 - (1.0d0 + dampk + 0.4d0*dampk**2 +
     &           (1.0d0/15.0d0)*dampk**3)*expdampk
         end if
c
c     calculate two-site scale factors
c
         if (alphai .ne. alphak) then
            termi = alphak**2/(alphak**2 - alphai**2)
            termk = alphai**2/(alphai**2 - alphak**2)
            scaleik(1) = 1.0d0 -termi*expdampi -termk*expdampk
            scaleik(3) = 1.0d0 - termi*(1.0d0 +dampi)*expdampi
     &           - termk*(1.0d0 + dampk)*expdampk
            scaleik(5) = 1.0d0 - termi*(1.0d0 + dampi +
     &           (1.0d0/3.0d0)*dampi**2)*expdampi -
     &           termk*(1.0d0 + dampk +
     &           (1.0d0/3.0d0)*dampk**2)*expdampk
            if (rorder .ge. 7) then
               scaleik(7) = 1.0d0 - termi*(1.0d0 + dampi +
     &              0.4d0*dampi**2 + (1.0d0/15.0d0)*dampi**3)*
     &              expdampi -
     &              termk*(1.0d0 + dampk +
     &              0.4d0*dampk**2 + (1.0d0/15.0d0)*dampk**3)*
     &              expdampk
            end if
            if (rorder .ge. 9) then
               scaleik(9) = 1.0d0 - termi*(1.0d0 + dampi +
     &              (3.0d0/7.0d0)*dampi**2 +
     &              (2.0d0/21.0d0)*dampi**3 +
     &              (1.0d0/105.0d0)*dampi**4)*expdampi -
     &              termk*(1.0d0 + dampk +
     &              (3.0d0/7.0d0)*dampk**2 +
     &              (2.0d0/21.0d0)*dampk**3 +
     &              (1.0d0/105.0d0)*dampk**4)*expdampk
            end if
            if (rorder .ge. 11) then
               scaleik(11) = 1.0d0 - termi*(1.0d0 + dampi +
     &              (4.0d0/9.0d0)*dampi**2 +
     &              (1.0d0/9.0d0)*dampi**3 +
     &              (1.0d0/63.0d0)*dampi**4 + 
     &              (1.0d0/945.0d0)*dampi**5)*expdampi - 
     &              termk*(1.0d0 + dampk +
     &              (4.0d0/9.0d0)*dampk**2+
     &              (1.0d0/9.0d0)*dampk**3+
     &              (1.0d0/63.0d0)*dampk**4 +
     &              (1.0d0/945.0d0)*dampk**5)*expdampk
            end if
         else
            scaleik(1) = 1.0d0 - (1.0d0+0.5d0*dampi)*expdampi
            scaleik(3) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2)
     &           *expdampi
            scaleik(5) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &           + (1.0d0/6.0d0)*dampi**3)*expdampi
            if (rorder .ge. 7) then
               scaleik(7) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &              + (1.0d0/6.0d0)*dampi**3
     &              + (1.0d0/30.0d0)*dampi**4)*expdampi
            end if
            if (rorder .ge. 9) then
               scaleik(9) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &              + (1.0d0/6.0d0)*dampi**3
     &              + (4.0d0/105.0d0)*dampi**4
     &              + (1.0d0/210.0d0)*dampi**5)*expdampi
            end if
            if (rorder .ge. 11) then
               scaleik(11) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2
     &              + (1.0d0/6.0d0)*dampi**3 
     &              + (5.0d0/126.0d0)*dampi**4
     &              + (2.0d0/315.0d0)*dampi**5
     &              + (1.0d0/1890.0d0)*dampi**6)*expdampi
            end if
         end if
      end if
c      print *,"scales",r
c      print *,"1",scalei(1),scalek(1),scaleik(1)
c      print *,"3",scalei(3),scalek(3),scaleik(3)
c      print *,"5",scalei(5),scalek(5),scaleik(5)
c      print *,"7",scaleik(7)
c      print *,"9",scaleik(9)
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
      scalei(1) = (1.0d0 - expdampi)*(1.0d0 - rexpdampi) 
      scalek(1) = (1.0d0 - expdampk)*(1.0d0 - rexpdampk)
      scalei(3) = 1.0d0 - (1.0d0 + dampi)*expdampi -
     &     (1.0d0 + rdampk)*rexpdampk + (1.0d0 + dampi +
     &     rdampk)*expdampi*rexpdampk
      scalek(3) = 1.0d0 - (1.0d0 + dampk)*expdampk -
     &     (1.0d0 + rdampi)*rexpdampi + (1.0d0 + dampk +
     &     rdampi)*expdampk*rexpdampi
      scalei(5) = 1.0d0 - (1.0d0 + dampi +
     &     (1.0d0/3.0d0)*dampi**2)*expdampi - (1.0d0 +
     &     rdampk + (1.0d0/3.0d0)*rdampk**2)*rexpdampk +
     &     (1.0d0 + dampi + rdampk +
     &     (1.0d0/3.0d0)*dampi**2 +
     &     (1.0d0/3.0d0)*rdampk**2 +
     &     (2.0d0/3.0d0)*dampi*rdampk)*expdampi*rexpdampk
      scalek(5) = 1.0d0 - (1.0d0 + dampk +
     &     (1.0d0/3.0d0)*dampk**2)*expdampk - (1.0d0 +
     &     rdampi + (1.0d0/3.0d0)*rdampi**2)*rexpdampi +
     &     (1.0d0 + dampk + rdampi +
     &     (1.0d0/3.0d0)*dampk**2 +
     &     (1.0d0/3.0d0)*rdampi**2 +
     &     (2.0d0/3.0d0)*dampk*rdampi)*expdampk*rexpdampi
      return
      end
c
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine dampgoverlap  --  generate gordon overlap damping coefficents  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "dampgoverlap" generates the damping coefficients for the one-site 
c     gordon damping functional form that go with corresponding powers of r
c
c     one-site scale factors are used for nuclear-electron interactions
c
      subroutine dampgoverlap(i,k,rorder,r,scalei,scalek,scaleik)
      use sizes
      use chgpen
      implicit none
      integer i,k
      integer rorder
      real*8 r,r2
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8 term,pre
      real*8 scalei(*),scalek(*)
      real*8 scaleik(*)
c
c     read in charge penetration damping parameters
c
      alphai = alpha(cpclass(i))
      alphak = alpha(cpclass(k))
      if (use_vdwclass) then
         alphai = alpha(vdwclass(i))
         alphak = alpha(vdwclass(k))
      end if
c
c     compute common factors for damping
c
      dampi = alphai*r
      dampk = alphak*r
      expdampi = exp(-dampi)
      expdampk = exp(-dampk)
      r2 = r**2
c
c     "new" Gordon damping model - model 1
c
      if (gtype .eq. "NEW") then
         scalei(1) = 0.0d0
         scalek(1) = 0.0d0
         scalei(3) = 0.0d0
         scalek(3) = 0.0d0
         scalei(5) = 0.0d0
         scalek(5) = 0.0d0
         scalei(7) = 0.0d0
         scalek(7) = 0.0d0
         if (alphai .ne. alphak) then
            term = alphai**2 - alphak**2
c     pre = 8.0d0*(alphai**3)*(alphak**3)/(r*(term**2))
            pre = 8.0d0*(alphai**3)*(alphak**3)/(term**2) 
            scaleik(1) = pre*(alphai*(r - 4.0d0*alphak/term)
     &           *expdampk + alphak*(r + 4.0d0*alphai/term)
     &           *expdampi)
            scaleik(3) = pre*alphai*alphak*((r2 - 4.0d0/term -
     &           4.0d0*dampk/term)*expdampk + (r2 + 
     &           4.0d0/term + 4.0d0*dampi/term)*expdampi)
            scaleik(5) = pre*alphai*alphak*(((1.0d0/3.0d0)*r2 + 
     &           (1.0d0/3.0d0)*dampk*r2 - 4.0d0/term -
     &           4.0d0*dampk/term - 
     &           4.0d0*(dampk**2)/(3.0d0*term))*expdampk + 
     &           ((1.0d0/3.0d0)*r2 + 
     &           (1.0d0/3.0d0)*dampi*r2 + 4.0d0/term + 
     &           4.0d0*dampi/term + 
     &           4.0d0*(dampi**2)/(3.0d0*term))*expdampi)
            scaleik(7) = pre*alphai*alphak*(((1.0d0/5.0d0)*r2 + 
     &           (1.0d0/5.0d0)*dampk*r2 + 
     &           (1.0d0/15.0d0)*(dampk**2)*r2 - 4.0d0/term
     &           - 4.0d0*dampk/term  
     &           - (8.0d0/5.0d0)*(dampk**2)/term
     &           - (4.0d0/15.0d0)*(dampk**3)/term)
     &           *expdampk
     &           + ((1.0d0/5.0d0)*r2 + (1.0d0/5.0d0)*dampi*r2 +
     &           (1.0d0/15.0d0)*(dampi**2)*r2 + 4.0d0/term +
     &           4.0d0*dampi/term + 
     &           (8.0d0/5.0d0)*(dampi**2)/term +
     &           (4.0d0/15.0d0)*(dampi**3)/term)*expdampi)
            scaleik(9) = pre*alphai*alphak*(((1.0d0/7.0d0)*r2 + 
     &           (1.0d0/7.0d0)*dampk*r2 + 
     &           (2.0d0/35.0d0)*(dampk**2)*r2 + 
     &           (1.0d0/105.0d0)*(dampk**3)*r2 - 
     &           4.0d0/term - 4.0d0*dampk/term - 
     &           (12.0d0/7.0d0)*(dampk**2)/term -
     &           (8.0d0/21.0d0)*(dampk**3)/term - 
     &           (4.0d0/105.0d0)*(dampk**4)/term)*expdampk + 
     &           ((1.0d0/7.0d0)*r2 + (1.0d0/7.0d0)*dampi*r2 + 
     &           (2.0d0/35.0d0)*(dampi**2)*r2 + 
     &           (1.0d0/105.0d0)*(dampi**3)*r2 +
     &           4.0d0/term + 4.0d0*dampi/term +
     &           (12.0d0/7.0d0)*(dampi**2)/term + 
     &           (8.0d0/21.0d0)*(dampi**3)/term +
     &           (4.0d0/105.0d0)*(dampi**4)/term)*expdampi)
         else
            pre = alphai**3
            scaleik(1) = pre*r*(1.0d0 + dampi +
     &           (1.0d0/3.0d0)*(dampi**2))*expdampi
c     oscale1ik = pre*(1.0d0/r)*(1.0d0 + dampi + 
c     &                    (1.0d0/3.0d0)*dampi**2)*expdampi
            scaleik(3) = pre*(1.0d0/3.0d0)*((alphai**2)*r**3 + 
     &           (alphai**3)*r**4)*expdampi
            scaleik(5) = pre*(1.0d0/9.0d0)*(alphai**4)*(r**5)*
     &           expdampi
            scaleik(7) = pre*(1.0d0/45.0d0)*(alphai**5)*(r**6)*
     &           expdampi
            scaleik(9) = pre*(1.0d0/315.0d0)*((alphai**5)*r**6 + 
     &           (alphai**6)*r**7)*expdampi
         end if
c
c     original gordon model (model 2)
c
      else
         scalei(1) = 0.0d0
         scalek(1) = 0.0d0
         scalei(3) = 0.0d0
         scalek(3) = 0.0d0
         scalei(5) = 0.0d0
         scalek(5) = 0.0d0
         scalei(7) = 0.0d0
         scalek(7) = 0.0d0
         if (alphai .ne. alphak) then
            term = 2.0d0*(alphai**2)*(alphak**2)
     &           /(alphai**2 - alphak**2)
            scaleik(1) = term*(expdampk - expdampi)
            scaleik(3) = term*(1.0d0 + dampk)*expdampk -
     &           term*(1.0d0 + dampi)*expdampi
            scaleik(5) = term*(1.0d0 + dampk +
     &           (1.0d0/3.0d0)*dampk**2)*expdampk -
     &           term*(1.0d0 + dampi +
     &           (1.0d0/3.0d0)*dampi**2)*expdampi
            scaleik(7) = term*(1.0d0 + dampk +
     &           (2.0d0/5.0d0)*dampk**2 +
     &           (1.0d0/15.0d0)*dampk**3)*expdampk -
     &           term*(1.0d0 + dampi +
     &           (2.0d0/5.0d0)*dampi**2 +
     &           (1.0d0/15.0d0)*dampi**3)*expdampi
            scaleik(9) = term*(1.0d0 + dampk +
     &           (3.0d0/7.0d0)*dampk**2 +
     &           (2.0d0/21.0d0)*dampk**3 +
     &           (1.0d0/105.0d0)*dampk**4)*expdampk -
     &           term*(1.0d0 + dampi +
     &           (3.0d0/7.0d0)*dampi**2 +
     &           (2.0d0/21.0d0)*dampi**3 +
     &           (1.0d0/105.0d0)*dampi**4)*expdampi
         else
            scaleik(1) = (alphai**3)*r*expdampi
            scaleik(3) = (alphai**4)*r2*expdampi
            scaleik(5) = (alphai**4)*((1.0d0/3.0d0)*r2 +
     &           (1.0d0/3.0d0)*alphai*r**3)*expdampi
            scaleik(7) = (alphai**4)*((1.0d0/5.0d0)*r2 +
     &           (1.0d0/5.0d0)*alphai*r**3 +
     &           (1.0d0/15.0d0)*(alphai**2)*r**4)*expdampi
            scaleik(9) = (alphai**4)*((1.0d0/7.0d0)*r2 +
     &           (1.0d0/7.0d0)*alphai*r**3 +
     &           (2.0d0/35.0d0)*(alphai**2)*r**4 +
     &           (1.0d0/105.0d0)*(alphai**3)*r**5)*expdampi
         end if
      end if
c
c      print *,"scales",r
c      print *,"1",scalei(1),scalek(1),scaleik(1)
c      print *,"3",scalei(3),scalek(3),scaleik(3)
c      print *,"5",scalei(5),scalek(5),scaleik(5)
c      print *,"7",scaleik(7)
c      print *,"9",scaleik(9)
      return
      end
