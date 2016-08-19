c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  chgpen.f  --  specifics of user-defined functional forms  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxele   maximum number of elements from periodic table
c
c     zeta         truhlar penetration parameter
c     alpha        gordon or piquemal penetration parameter
c     beta         piquemal penetration charge overlap parameter
c     regular      nuclear regularization parameter
c     penetration  charge penetration energy functional form
c     num_ele      number of core electrons
c
c
      module chgpen
      use sizes
      implicit none
      real*8 zeta(maxtyp)
      real*8 alpha(maxtyp)
      real*8 alphap(maxtyp)
      real*8 bfactor(maxtyp)
      real*8 regular(maxtyp)
      real*8 cbfactor
      real*8 xpolr(maxtyp)
      real*8 xathl(maxtyp)
      real*8 polfactor(maxtyp)
ccccccccc
c     ignore these
      real*8 malpha
      real*8 dalpha(3)
ccccccccc
      integer cpclass(maxatm)
      character*20 penetration
      character*20 bfactor_mode
      character*20 chgpen_mode
      character*20 num_ele
      character*20 regularize
      save
      end
