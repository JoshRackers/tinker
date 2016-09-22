c
c     #############################################################
c     ## COPYRIGHT (C) 2016 by Josh Rackers & Jay William Ponder ##
c     ##                   All Rights Reserved                   ## 
c     ############################################################# 
c                                     
c     #########################################################  
c     ##                                                     ##             
c     ##  routines below generate orders of the multipole    ## 
c     ##  interaction matrix, "t"                            ##
c     ##                                                     ##
c     #########################################################
c
c
c     "tmatrix" contains routines that generate the undamped
c     and damped multipole interaction matrices for use in 
c     permanent and induced multipole interaction energies 
c     and gradients
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine t0matrixrr1  --  generate undamped t0 multipole  ##
c     ##                           interaction matrix              ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "t0matrixrr1" generates the undamped rank 0 multipole interaction
c     matrix, (1/r)
c
      subroutine t0matrixrr1(rr1,t0rr1)
      use sizes
      implicit none
      real*8 rr1
      real*8 t0rr1
c
      t0rr1 = rr1
      return
      end
c
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine t1matrixrr3  --  generate the undamped t1 multipole    ##
c     ##                           interaction matrix                    ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "t1matrixrr3" generates the rr3 part of the undamped rank 1 
c     multipole interaction matrix, del(1/r)
c
      subroutine t1matrixrr3(xr,yr,zr,rr3,t1rr3)
      use sizes
      implicit none
      real*8 xr,yr,zr
      real*8 rr3
      real*8 t1rr3(3)
c
      t1rr3(1) = -xr*rr3
      t1rr3(2) = -yr*rr3
      t1rr3(3) = -zr*rr3
      return
      end
c
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine t2matrixrr3  --  generate the rr3 part of the undamped ##  
c     ##                           undamped t2 multipole interaction matrix ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "t2matrixrr3" generates the rr3 part of the undamped rank 2 
c     multipole interaction matrix, (del)^2(1/r)
c
      subroutine t2matrixrr3(xr,yr,zr,rr3,t2rr3)
      use sizes
      implicit none
      real*8 xr,yr,zr
      real*8 rr3
      real*8 t2rr3(3,3)
c
c     prefactors already applied to rr
c
      t2rr3(1,1) = -rr3
      t2rr3(2,1) = 0.0d0
      t2rr3(3,1) = 0.0d0
c      t2rr3(2,1) = 0.0d0
      t2rr3(2,2) = -rr3
      t2rr3(3,2) = 0.0d0
c      t2rr3(3,1) = 0.0d0
c      t2rr3(3,2) = 0.0d0
      t2rr3(3,3) = -rr3
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine t2matrixrr5  --  generate the rr5 part of the undamped ##
c     ##                              t2 multipole interaction matrix       ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "t2matrixrr5" generates the rr5 part of the undamped rank 2 
c     multipole interaction matrix, (del)^2(1/r)
c
      subroutine t2matrixrr5(xr,yr,zr,rr5,t2rr5)
      use sizes
      implicit none
      real*8 xr,yr,zr
      real*8 rr5
      real*8 t2rr5(3,3)
c
c     prefactors already applied to rr
c
      t2rr5(1,1) = xr*xr*rr5 
      t2rr5(2,1) = xr*yr*rr5
      t2rr5(3,1) = xr*zr*rr5
c      t2rr5(2,1) = t2rr5(2,1)
      t2rr5(2,2) = yr*yr*rr5
      t2rr5(3,2) = yr*zr*rr5
c      t2rr5(3,1) = t2rr5(3,1)
c      t2rr5(3,2) = t2rr5(3,2)
      t2rr5(3,3) = zr*zr*rr5
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine t3matrixrr5  --  generate the rr5 part of the undamped ##
c     ##                              t3 multipole interaction matrix       ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "t3matrixrr5" generates the rr5 part of the undamped rank 3 
c     multipole interaction matrix, (del)^3(1/r)
c
      subroutine t3matrixrr5(xr,yr,zr,rr5,t3rr5)
      use sizes
      implicit none
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 rr5
      real*8 t3rr5(3,3,3)
      integer i,j,k
c
c     prefactors already applied to rr
c
      xr2 = xr*xr
      yr2 = yr*yr
      zr2 = zr*zr
      t3rr5(1,1,1) = 3.0d0*xr*rr5
      t3rr5(2,1,1) = yr*rr5
      t3rr5(3,1,1) = zr*rr5
c      t3rr5(1,2,1) = t3rr5(2,1,1)
      t3rr5(2,2,1) = xr*rr5
      t3rr5(3,2,1) = 0.0d0
c      t3rr5(1,3,1) = t3rr5(3,1,1)
c      t3rr5(1,3,2) = t3rr5(3,2,1)
      t3rr5(3,3,1) = xr*rr5
c      t3rr5(2,1,1) = t3rr5(2,1,1)
c      t3rr5(2,1,2) = t3rr5(2,2,1)
c      t3rr5(2,1,3) = t3rr5(3,2,1)
c      t3rr5(2,2,1) = t3rr5(2,2,1)
      t3rr5(2,2,2) = 3.0d0*yr*rr5
      t3rr5(3,2,2) = zr*rr5
c      t3rr5(2,3,1) = t3rr5(3,2,1)
c      t3rr5(2,3,2) = t3rr5(3,2,2)
      t3rr5(3,3,2) = yr*rr5
c      t3rr5(3,1,1) = t3rr5(3,1,1)
c      t3rr5(3,1,2) = t3rr5(3,2,1)
c      t3rr5(3,1,3) = t3rr5(3,3,1)
c      t3rr5(3,2,1) = t3rr5(3,2,1)
c      t3rr5(3,2,2) = t3rr5(3,2,2)
c      t3rr5(3,2,3) = t3rr5(3,3,2)
c      t3rr5(3,3,1) = t3rr5(3,3,1)
c      t3rr5(3,3,2) = t3rr5(3,3,2)
      t3rr5(3,3,3) = 3.0d0*zr*rr5
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine t3matrixrr7  --  generate the rr7 part of the undamped ##
c     ##                              t3 multipole interaction matrix       ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "t3matrixrr7" generates the rr7 part of the undamped rank 3 
c     multipole interaction matrix, (del)^3(1/r)
c
      subroutine t3matrixrr7(xr,yr,zr,rr7,t3rr7)
      use sizes
      implicit none
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 rr7
      real*8 t3rr7(3,3,3)
      integer i,j,k
c
c     prefactors already applied to rr
c
      xr2 = xr*xr
      yr2 = yr*yr
      zr2 = zr*zr
      t3rr7(1,1,1) = -xr2*xr*rr7 
      t3rr7(2,1,1) = -xr2*yr*rr7
      t3rr7(3,1,1) = -xr2*zr*rr7
c      t3rr7(1,2,1) = t3rr7(2,1,1)
      t3rr7(2,2,1) = -yr2*xr*rr7
      t3rr7(3,2,1) = -xr*yr*zr*rr7
c      t3rr7(1,3,1) = t3rr7(3,1,1)
c      t3rr7(1,3,2) = t3rr7(3,2,1)
      t3rr7(3,3,1) = -zr2*xr*rr7
c      t3rr7(2,1,1) = t3rr7(2,1,1)
c      t3rr7(2,1,2) = t3rr7(2,2,1)
c      t3rr7(2,1,3) = t3rr7(3,2,1)
c      t3rr7(2,2,1) = t3rr7(2,2,1)
      t3rr7(2,2,2) = -yr2*yr*rr7
      t3rr7(3,2,2) = -yr2*zr*rr7
c      t3rr7(2,3,1) = t3rr7(3,2,1)
c      t3rr7(2,3,2) = t3rr7(3,2,2)
      t3rr7(3,3,2) = -zr2*yr*rr7
c      t3rr7(3,1,1) = t3rr7(3,1,1)
c      t3rr7(3,1,2) = t3rr7(3,2,1)
c      t3rr7(3,1,3) = t3rr7(3,3,1)
c      t3rr7(3,2,1) = t3rr7(3,2,1)
c      t3rr7(3,2,2) = t3rr7(3,2,2)
c      t3rr7(3,2,3) = t3rr7(3,3,2)
c      t3rr7(3,3,1) = t3rr7(3,3,1)
c      t3rr7(3,3,2) = t3rr7(3,3,2)
      t3rr7(3,3,3) = -zr2*zr*rr7
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine t4matrixrr5  --  generate the rr5 part of the undamped ##
c     ##                              t4 multipole interaction matrix       ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "t4matrixrr5" generates the rr5 part of the undamped rank 4 
c     multipole interaction matrix, (del)^4(1/r)
c
      subroutine t4matrixrr5(xr,yr,zr,rr5,t4rr5)
      use sizes
      implicit none
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 xr3,yr3,zr3
      real*8 rr5
      real*8 t4rr5(3,3,3,3)
c
c     prefactors already applied to rr
c
      t4rr5(1,1,1,1) = 3.0d0*rr5
      t4rr5(2,1,1,1) = 0.0d0 
      t4rr5(3,1,1,1) = 0.0d0
c      t4rr5(2,1,1,1) = t4rr5(2,1,1,1)
      t4rr5(2,2,1,1) = rr5
      t4rr5(3,2,1,1) = 0.0d0
c      t4rr5(3,1,1,1) = t4rr5(3,1,1,1)
c      t4rr5(3,2,1,1) = t4rr5(3,2,1,1)
      t4rr5(3,3,1,1) = rr5
c      t4rr5(2,1,1,1) = t4rr5(2,1,1,1)
c      t4rr5(2,2,1,1) = t4rr5(2,2,1,1)
c      t4rr5(3,2,1,1) = t4rr5(3,2,1,1) 
c      t4rr5(2,2,1,1) = t4rr5(2,2,1,1)
      t4rr5(2,2,2,1) = 0.0d0
      t4rr5(3,2,2,1) = 0.0d0
c      t4rr5(3,2,1,1) = t4rr5(3,2,1,1) 
c      t4rr5(3,2,2,1) = t4rr5(3,2,2,1)
      t4rr5(3,3,2,1) = 0.0d0
c      t4rr5(3,1,1,1) = t4rr5(3,1,1,1)
c      t4rr5(3,2,1,1) = t4rr5(3,2,1,1)
c      t4rr5(3,3,1,1) = t4rr5(3,3,1,1)
c      t4rr5(3,2,1,1) = t4rr5(3,2,1,1)
c      t4rr5(3,2,2,1) = t4rr5(3,2,2,1)
c      t4rr5(3,3,2,1) = t4rr5(3,3,2,1)
c      t4rr5(3,3,1,1) = t4rr5(3,3,1,1)
c      t4rr5(3,3,2,1) = t4rr5(3,3,2,1)
      t4rr5(3,3,3,1) = 0.0d0
c      t4rr5(2,1,1,1) = t4rr5(2,1,1,1)
c      t4rr5(2,2,1,1) = t4rr5(2,2,1,1)
c      t4rr5(3,2,1,1) = t4rr5(3,2,1,1)
c      t4rr5(2,2,1,1) = t4rr5(2,2,1,1)
c      t4rr5(2,2,2,1) = t4rr5(2,2,2,1)
c      t4rr5(3,2,2,1) = t4rr5(3,2,2,1)
c      t4rr5(3,2,1,1) = t4rr5(3,2,1,1)
c      t4rr5(3,2,2,1) = t4rr5(3,2,2,1)
c      t4rr5(3,3,2,1) = t4rr5(3,3,2,1)
c      t4rr5(2,2,1,1) = t4rr5(2,2,1,1)
c      t4rr5(2,2,2,1) = t4rr5(2,2,2,1)
c      t4rr5(3,2,2,1) = t4rr5(3,2,2,1)
c      t4rr5(2,2,2,1) = t4rr5(2,2,2,1)
      t4rr5(2,2,2,2) = 3.0d0*rr5
      t4rr5(3,2,2,2) = 0.0d0
c      t4rr5(3,2,2,1) = t4rr5(3,2,2,1)
c      t4rr5(3,2,2,2) = t4rr5(3,2,2,2)
      t4rr5(3,3,2,2) = rr5
c      t4rr5(3,2,1,1) = t4rr5(3,2,1,1)
c      t4rr5(3,2,2,1) = t4rr5(3,2,2,1)
c      t4rr5(3,3,2,1) = t4rr5(3,3,2,1)
c      t4rr5(3,2,2,1) = t4rr5(3,2,2,1)
c      t4rr5(3,2,2,2) = t4rr5(3,2,2,2)
c      t4rr5(3,3,2,2) = t4rr5(3,3,2,2)
c      t4rr5(3,3,2,1) = t4rr5(3,3,2,1)
c      t4rr5(3,3,2,2) = t4rr5(3,3,2,2)
      t4rr5(3,3,3,2) = 0.0d0
c      t4rr5(3,1,1,1) = t4rr5(3,1,1,1)
c      t4rr5(3,2,1,1) = t4rr5(3,2,1,1)
c      t4rr5(3,3,1,1) = t4rr5(3,3,1,1)
c      t4rr5(3,2,1,1) = t4rr5(3,2,1,1)
c      t4rr5(3,2,2,1) = t4rr5(3,2,2,1)
c      t4rr5(3,3,2,1) = t4rr5(3,3,2,1)
c      t4rr5(3,3,1,1) = t4rr5(3,3,1,1)
c      t4rr5(3,3,2,1) = t4rr5(3,3,2,1)
c      t4rr5(3,3,3,1) = t4rr5(3,3,3,1)
c      t4rr5(3,2,1,1) = t4rr5(3,2,1,1)
c      t4rr5(3,2,2,1) = t4rr5(3,2,2,1)
c      t4rr5(3,3,2,1) = t4rr5(3,3,2,1)
c      t4rr5(3,2,2,1) = t4rr5(3,2,2,1)
c      t4rr5(3,2,2,2) = t4rr5(3,2,2,2)
c      t4rr5(3,3,2,2) = t4rr5(3,3,2,2)
c      t4rr5(3,3,2,1) = t4rr5(3,3,2,1)
c      t4rr5(3,3,2,2) = t4rr5(3,3,2,2)
c      t4rr5(3,3,3,2) = t4rr5(3,3,3,2)
c      t4rr5(3,3,1,1) = t4rr5(3,3,1,1)
c      t4rr5(3,3,2,1) = t4rr5(3,3,2,1)
c      t4rr5(3,3,3,1) = t4rr5(3,3,3,1)
c      t4rr5(3,3,2,1) = t4rr5(3,3,2,1)
c      t4rr5(3,3,2,2) = t4rr5(3,3,2,2)
c      t4rr5(3,3,3,2) = t4rr5(3,3,3,2)
c      t4rr5(3,3,3,1) = t4rr5(3,3,3,1)
c      t4rr5(3,3,3,2) = t4rr5(3,3,3,2)
      t4rr5(3,3,3,3) = 3.0d0*rr5
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine t4matrixrr7  --  generate the rr7 part of the undamped ##
c     ##                              t4 multipole interaction matrix       ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "t4matrixrr7" generates the rr7 part of the undamped rank 4 
c     multipole interaction matrix, (del)^4(1/r)
c
      subroutine t4matrixrr7(xr,yr,zr,rr7,t4rr7)
      use sizes
      implicit none
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 xr3,yr3,zr3
      real*8 rr7
      real*8 t4rr7(3,3,3,3)
c
c     prefactors already applied to rr
c
      xr2 = xr*xr
      yr2 = yr*yr
      zr2 = zr*zr
      xr3 = xr2*xr
      yr3 = yr2*yr
      zr3 = zr2*zr
      t4rr7(1,1,1,1) = -6.0d0*xr2*rr7
      t4rr7(2,1,1,1) = -3.0d0*xr*yr*rr7 
      t4rr7(3,1,1,1) = -3.0d0*xr*zr*rr7
c      t4rr7(2,1,1,1) = t4rr7(2,1,1,1)
      t4rr7(2,2,1,1) = -xr2*rr7 - yr2*rr7
      t4rr7(3,2,1,1) = -yr*zr*rr7
c      t4rr7(3,1,1,1) = t4rr7(3,1,1,1)
c      t4rr7(3,2,1,1) = t4rr7(3,2,1,1)
      t4rr7(3,3,1,1) = -xr2*rr7 - zr2*rr7
c      t4rr7(2,1,1,1) = t4rr7(2,1,1,1)
c      t4rr7(2,2,1,1) = t4rr7(2,2,1,1)
c      t4rr7(3,2,1,1) = t4rr7(3,2,1,1) 
c      t4rr7(2,2,1,1) = t4rr7(2,2,1,1)
      t4rr7(2,2,2,1) = -3.0d0*xr*yr*rr7
      t4rr7(3,2,2,1) = -xr*zr*rr7
c      t4rr7(3,2,1,1) = t4rr7(3,2,1,1) 
c      t4rr7(3,2,2,1) = t4rr7(3,2,2,1)
      t4rr7(3,3,2,1) = -xr*yr*rr7
c      t4rr7(3,1,1,1) = t4rr7(3,1,1,1)
c      t4rr7(3,2,1,1) = t4rr7(3,2,1,1)
c      t4rr7(3,3,1,1) = t4rr7(3,3,1,1)
c      t4rr7(3,2,1,1) = t4rr7(3,2,1,1)
c      t4rr7(3,2,2,1) = t4rr7(3,2,2,1)
c      t4rr7(3,3,2,1) = t4rr7(3,3,2,1)
c      t4rr7(3,3,1,1) = t4rr7(3,3,1,1)
c      t4rr7(3,3,2,1) = t4rr7(3,3,2,1)
      t4rr7(3,3,3,1) = -3.0d0*xr*zr*rr7
c      t4rr7(2,1,1,1) = t4rr7(2,1,1,1)
c      t4rr7(2,2,1,1) = t4rr7(2,2,1,1)
c      t4rr7(3,2,1,1) = t4rr7(3,2,1,1)
c      t4rr7(2,2,1,1) = t4rr7(2,2,1,1)
c      t4rr7(2,2,2,1) = t4rr7(2,2,2,1)
c      t4rr7(3,2,2,1) = t4rr7(3,2,2,1)
c      t4rr7(3,2,1,1) = t4rr7(3,2,1,1)
c      t4rr7(3,2,2,1) = t4rr7(3,2,2,1)
c      t4rr7(3,3,2,1) = t4rr7(3,3,2,1)
c      t4rr7(2,2,1,1) = t4rr7(2,2,1,1)
c      t4rr7(2,2,2,1) = t4rr7(2,2,2,1)
c      t4rr7(3,2,2,1) = t4rr7(3,2,2,1)
c      t4rr7(2,2,2,1) = t4rr7(2,2,2,1)
      t4rr7(2,2,2,2) = -6.0d0*yr2*rr7
      t4rr7(3,2,2,2) = -3.0d0*yr*zr*rr7
c      t4rr7(3,2,2,1) = t4rr7(3,2,2,1)
c      t4rr7(3,2,2,2) = t4rr7(3,2,2,2)
      t4rr7(3,3,2,2) = -yr2*rr7 - zr2*rr7
c      t4rr7(3,2,1,1) = t4rr7(3,2,1,1)
c      t4rr7(3,2,2,1) = t4rr7(3,2,2,1)
c      t4rr7(3,3,2,1) = t4rr7(3,3,2,1)
c      t4rr7(3,2,2,1) = t4rr7(3,2,2,1)
c      t4rr7(3,2,2,2) = t4rr7(3,2,2,2)
c      t4rr7(3,3,2,2) = t4rr7(3,3,2,2)
c      t4rr7(3,3,2,1) = t4rr7(3,3,2,1)
c      t4rr7(3,3,2,2) = t4rr7(3,3,2,2)
      t4rr7(3,3,3,2) = -3.0d0*yr*zr*rr7
c      t4rr7(3,1,1,1) = t4rr7(3,1,1,1)
c      t4rr7(3,2,1,1) = t4rr7(3,2,1,1)
c      t4rr7(3,3,1,1) = t4rr7(3,3,1,1)
c      t4rr7(3,2,1,1) = t4rr7(3,2,1,1)
c      t4rr7(3,2,2,1) = t4rr7(3,2,2,1)
c      t4rr7(3,3,2,1) = t4rr7(3,3,2,1)
c      t4rr7(3,3,1,1) = t4rr7(3,3,1,1)
c      t4rr7(3,3,2,1) = t4rr7(3,3,2,1)
c      t4rr7(3,3,3,1) = t4rr7(3,3,3,1)
c      t4rr7(3,2,1,1) = t4rr7(3,2,1,1)
c      t4rr7(3,2,2,1) = t4rr7(3,2,2,1)
c      t4rr7(3,3,2,1) = t4rr7(3,3,2,1)
c      t4rr7(3,2,2,1) = t4rr7(3,2,2,1)
c      t4rr7(3,2,2,2) = t4rr7(3,2,2,2)
c      t4rr7(3,3,2,2) = t4rr7(3,3,2,2)
c      t4rr7(3,3,2,1) = t4rr7(3,3,2,1)
c      t4rr7(3,3,2,2) = t4rr7(3,3,2,2)
c      t4rr7(3,3,3,2) = t4rr7(3,3,3,2)
c      t4rr7(3,3,1,1) = t4rr7(3,3,1,1)
c      t4rr7(3,3,2,1) = t4rr7(3,3,2,1)
c      t4rr7(3,3,3,1) = t4rr7(3,3,3,1)
c      t4rr7(3,3,2,1) = t4rr7(3,3,2,1)
c      t4rr7(3,3,2,2) = t4rr7(3,3,2,2)
c      t4rr7(3,3,3,2) = t4rr7(3,3,3,2)
c      t4rr7(3,3,3,1) = t4rr7(3,3,3,1)
c      t4rr7(3,3,3,2) = t4rr7(3,3,3,2)
      t4rr7(3,3,3,3) = -6.0d0*zr2*rr7
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine t4matrixrr9  --  generate the rr9 part of the undamped ##
c     ##                              t4 multipole interaction matrix       ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "t4matrixrr9" generates the rr9 part of the undamped rank 4 
c     multipole interaction matrix, (del)^4(1/r)
c
      subroutine t4matrixrr9(xr,yr,zr,rr9,t4rr9)
      use sizes
      implicit none
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 xr3,yr3,zr3
      real*8 rr9
      real*8 t4rr9(3,3,3,3)
c
c     prefactors already applied to rr
c
      xr2 = xr*xr
      yr2 = yr*yr
      zr2 = zr*zr
      xr3 = xr2*xr
      yr3 = yr2*yr
      zr3 = zr2*zr
      t4rr9(1,1,1,1) = xr3*xr*rr9
      t4rr9(2,1,1,1) = xr3*yr*rr9 
      t4rr9(3,1,1,1) = xr3*zr*rr9
c      t4rr9(2,1,1,1) = t4rr9(2,1,1,1)
      t4rr9(2,2,1,1) = xr2*yr2*rr9
      t4rr9(3,2,1,1) = xr2*yr*zr*rr9
c      t4rr9(3,1,1,1) = t4rr9(3,1,1,1)
c      t4rr9(3,2,1,1) = t4rr9(3,2,1,1)
      t4rr9(3,3,1,1) = xr2*zr2*rr9
c      t4rr9(2,1,1,1) = t4rr9(2,1,1,1)
c      t4rr9(2,2,1,1) = t4rr9(2,2,1,1)
c      t4rr9(3,2,1,1) = t4rr9(3,2,1,1) 
c      t4rr9(2,2,1,1) = t4rr9(2,2,1,1)
      t4rr9(2,2,2,1) = yr3*xr*rr9
      t4rr9(3,2,2,1) = yr2*xr*zr*rr9
c      t4rr9(3,2,1,1) = t4rr9(3,2,1,1) 
c      t4rr9(3,2,2,1) = t4rr9(3,2,2,1)
      t4rr9(3,3,2,1) = zr2*xr*yr*rr9
c      t4rr9(3,1,1,1) = t4rr9(3,1,1,1)
c      t4rr9(3,2,1,1) = t4rr9(3,2,1,1)
c      t4rr9(3,3,1,1) = t4rr9(3,3,1,1)
c      t4rr9(3,2,1,1) = t4rr9(3,2,1,1)
c      t4rr9(3,2,2,1) = t4rr9(3,2,2,1)
c      t4rr9(3,3,2,1) = t4rr9(3,3,2,1)
c      t4rr9(3,3,1,1) = t4rr9(3,3,1,1)
c      t4rr9(3,3,2,1) = t4rr9(3,3,2,1)
      t4rr9(3,3,3,1) = zr3*xr*rr9
c      t4rr9(2,1,1,1) = t4rr9(2,1,1,1)
c      t4rr9(2,2,1,1) = t4rr9(2,2,1,1)
c      t4rr9(3,2,1,1) = t4rr9(3,2,1,1)
c      t4rr9(2,2,1,1) = t4rr9(2,2,1,1)
c      t4rr9(2,2,2,1) = t4rr9(2,2,2,1)
c      t4rr9(3,2,2,1) = t4rr9(3,2,2,1)
c      t4rr9(3,2,1,1) = t4rr9(3,2,1,1)
c      t4rr9(3,2,2,1) = t4rr9(3,2,2,1)
c      t4rr9(3,3,2,1) = t4rr9(3,3,2,1)
c      t4rr9(2,2,1,1) = t4rr9(2,2,1,1)
c      t4rr9(2,2,2,1) = t4rr9(2,2,2,1)
c      t4rr9(3,2,2,1) = t4rr9(3,2,2,1)
c      t4rr9(2,2,2,1) = t4rr9(2,2,2,1)
      t4rr9(2,2,2,2) = yr3*yr*rr9
      t4rr9(3,2,2,2) = yr3*zr*rr9
c      t4rr9(3,2,2,1) = t4rr9(3,2,2,1)
c      t4rr9(3,2,2,2) = t4rr9(3,2,2,2)
      t4rr9(3,3,2,2) = yr2*zr2*rr9
c      t4rr9(3,2,1,1) = t4rr9(3,2,1,1)
c      t4rr9(3,2,2,1) = t4rr9(3,2,2,1)
c      t4rr9(3,3,2,1) = t4rr9(3,3,2,1)
c      t4rr9(3,2,2,1) = t4rr9(3,2,2,1)
c      t4rr9(3,2,2,2) = t4rr9(3,2,2,2)
c      t4rr9(3,3,2,2) = t4rr9(3,3,2,2)
c      t4rr9(3,3,2,1) = t4rr9(3,3,2,1)
c      t4rr9(3,3,2,2) = t4rr9(3,3,2,2)
      t4rr9(3,3,3,2) = zr3*yr*rr9
c      t4rr9(3,1,1,1) = t4rr9(3,1,1,1)
c      t4rr9(3,2,1,1) = t4rr9(3,2,1,1)
c      t4rr9(3,3,1,1) = t4rr9(3,3,1,1)
c      t4rr9(3,2,1,1) = t4rr9(3,2,1,1)
c      t4rr9(3,2,2,1) = t4rr9(3,2,2,1)
c      t4rr9(3,3,2,1) = t4rr9(3,3,2,1)
c      t4rr9(3,3,1,1) = t4rr9(3,3,1,1)
c      t4rr9(3,3,2,1) = t4rr9(3,3,2,1)
c      t4rr9(3,3,3,1) = t4rr9(3,3,3,1)
c      t4rr9(3,2,1,1) = t4rr9(3,2,1,1)
c      t4rr9(3,2,2,1) = t4rr9(3,2,2,1)
c      t4rr9(3,3,2,1) = t4rr9(3,3,2,1)
c      t4rr9(3,2,2,1) = t4rr9(3,2,2,1)
c      t4rr9(3,2,2,2) = t4rr9(3,2,2,2)
c      t4rr9(3,3,2,2) = t4rr9(3,3,2,2)
c      t4rr9(3,3,2,1) = t4rr9(3,3,2,1)
c      t4rr9(3,3,2,2) = t4rr9(3,3,2,2)
c      t4rr9(3,3,3,2) = t4rr9(3,3,3,2)
c      t4rr9(3,3,1,1) = t4rr9(3,3,1,1)
c      t4rr9(3,3,2,1) = t4rr9(3,3,2,1)
c      t4rr9(3,3,3,1) = t4rr9(3,3,3,1)
c      t4rr9(3,3,2,1) = t4rr9(3,3,2,1)
c      t4rr9(3,3,2,2) = t4rr9(3,3,2,2)
c      t4rr9(3,3,3,2) = t4rr9(3,3,3,2)
c      t4rr9(3,3,3,1) = t4rr9(3,3,3,1)
c      t4rr9(3,3,3,2) = t4rr9(3,3,3,2)
      t4rr9(3,3,3,3) = zr3*zr*rr9
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine t5matrixrr7  --  generate the rr7 part of the undamped ##
c     ##                              t5 multipole interaction matrix       ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "t5matrixrr7" generates the rr7 part of the undamped rank 5 
c     multipole interaction matrix, (del)^5(1/r)
c
      subroutine t5matrixrr7(xr,yr,zr,rr7,t5rr7)
      use sizes
      implicit none
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 xr3,yr3,zr3
      real*8 xr4,yr4,zr4
      real*8 rr7
      real*8 t5rr7(3,3,3,3,3)
c
c     prefactors already applied to rr
c
      xr2 = xr*xr
      yr2 = yr*yr
      zr2 = zr*zr
      xr3 = xr2*xr
      yr3 = yr2*yr
      zr3 = zr2*zr
      xr4 = xr2*xr2
      yr4 = yr2*yr2
      zr4 = zr2*zr2
      t5rr7(1,1,1,1,1) = -15.0d0*xr*rr7
      t5rr7(2,1,1,1,1) = -3.0d0*yr*rr7
      t5rr7(3,1,1,1,1) = -3.0d0*zr*rr7
c      t5rr7(2,1,1,1,1) = t5rr7(2,1,1,1,1)
      t5rr7(2,2,1,1,1) = -3.0d0*xr*rr7
      t5rr7(3,2,1,1,1) = 0.0d0
c      t5rr7(3,1,1,1,1) = t5rr7(3,1,1,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
      t5rr7(3,3,1,1,1) = -3.0d0*xr*rr7
c      t5rr7(2,1,1,1,1) = t5rr7(2,1,1,1,1)
c      t5rr7(2,2,1,1,1) = t5rr7(2,2,1,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(2,2,1,1,1) = t5rr7(2,2,1,1,1)
      t5rr7(2,2,2,1,1) = -3.0d0*yr*rr7
      t5rr7(3,2,2,1,1) = -zr*rr7
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
      t5rr7(3,3,2,1,1) = -yr*rr7
c      t5rr7(3,1,1,1,1) = t5rr7(3,1,1,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,3,1,1,1) = t5rr7(3,3,1,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,1,1,1) = t5rr7(3,3,1,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
      t5rr7(3,3,3,1,1) = -3.0d0*zr*rr7
c      t5rr7(2,1,1,1,1) = t5rr7(2,1,1,1,1)
c      t5rr7(2,2,1,1,1) = t5rr7(2,2,1,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(2,2,1,1,1) = t5rr7(2,2,1,1,1)
c      t5rr7(2,2,2,1,1) = t5rr7(2,2,2,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(2,2,1,1,1) = t5rr7(2,2,1,1,1)
c      t5rr7(2,2,2,1,1) = t5rr7(2,2,2,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(2,2,2,1,1) = t5rr7(2,2,2,1,1)
      t5rr7(2,2,2,2,1) = -3.0d0*xr*rr7
      t5rr7(3,2,2,2,1) = 0.0d0
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
      t5rr7(3,3,2,2,1) = -xr*rr7
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
      t5rr7(3,3,3,2,1) = 0.0d0
c      t5rr7(3,1,1,1,1) = t5rr7(3,1,1,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,3,1,1,1) = t5rr7(3,3,1,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,1,1,1) = t5rr7(3,3,1,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,3,1,1) = t5rr7(3,3,3,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,1,1,1) = t5rr7(3,3,1,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,3,1,1) = t5rr7(3,3,3,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,3,1,1) = t5rr7(3,3,3,1,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
      t5rr7(3,3,3,3,1) = -3.0d0*xr*rr7
c
c      t5rr7(2,1,1,1,1) = t5rr7(2,1,1,1,1)
c      t5rr7(2,2,1,1,1) = t5rr7(2,2,1,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(2,2,1,1,1) = t5rr7(2,2,1,1,1)
c      t5rr7(2,2,2,1,1) = t5rr7(2,2,2,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(2,2,1,1,1) = t5rr7(2,2,1,1,1)
c      t5rr7(2,2,2,1,1) = t5rr7(2,2,2,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(2,2,2,1,1) = t5rr7(2,2,2,1,1)
c      t5rr7(2,2,2,2,1) = t5rr7(2,2,2,2,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(2,2,1,1,1) = t5rr7(2,2,1,1,1)
c      t5rr7(2,2,2,1,1) = t5rr7(2,2,2,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(2,2,2,1,1) = t5rr7(2,2,2,1,1)
c      t5rr7(2,2,2,2,1) = t5rr7(2,2,2,2,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(2,2,2,1,1) = t5rr7(2,2,2,1,1)
c      t5rr7(2,2,2,2,1) = t5rr7(2,2,2,2,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(2,2,2,2,1) = t5rr7(2,2,2,2,1)
      t5rr7(2,2,2,2,2) = -15.0d0*yr*rr7
      t5rr7(3,2,2,2,2) = -3.0d0*zr*rr7
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,2,2,2,2) = t5rr7(3,2,2,2,2)
      t5rr7(3,3,2,2,2) = -3.0d0*yr*rr7
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,2,2,2,2) = t5rr7(3,2,2,2,2)
c      t5rr7(3,3,2,2,2) = t5rr7(3,3,2,2,2)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,2,2) = t5rr7(3,3,2,2,2)
      t5rr7(3,3,3,2,2) = -3.0d0*zr*rr7
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,2,2,2,2) = t5rr7(3,2,2,2,2)
c      t5rr7(3,3,2,2,2) = t5rr7(3,3,2,2,2)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,2,2) = t5rr7(3,3,2,2,2)
c      t5rr7(3,3,3,2,2) = t5rr7(3,3,3,2,2)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,2,2) = t5rr7(3,3,2,2,2)
c      t5rr7(3,3,3,2,2) = t5rr7(3,3,3,2,2)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,3,2,2) = t5rr7(3,3,3,2,2)
      t5rr7(3,3,3,3,2) = -3.0d0*yr*rr7
c
c      t5rr7(3,1,1,1,1) = t5rr7(3,1,1,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,3,1,1,1) = t5rr7(3,3,1,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,1,1,1) = t5rr7(3,3,1,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,3,1,1) = t5rr7(3,3,3,1,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,1,1,1) = t5rr7(3,3,1,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,3,1,1) = t5rr7(3,3,3,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,3,1,1) = t5rr7(3,3,3,1,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,3,3,1) = t5rr7(3,3,3,3,1)
c      t5rr7(3,2,1,1,1) = t5rr7(3,2,1,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,2,2,1,1) = t5rr7(3,2,2,1,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,2,2,2,1) = t5rr7(3,2,2,2,1)
c      t5rr7(3,2,2,2,2) = t5rr7(3,2,2,2,2)
c      t5rr7(3,3,2,2,2) = t5rr7(3,3,2,2,2)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,2,2) = t5rr7(3,3,2,2,2)
c      t5rr7(3,3,3,2,2) = t5rr7(3,3,3,2,2)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,2,2) = t5rr7(3,3,2,2,2)
c      t5rr7(3,3,3,2,2) = t5rr7(3,3,3,2,2)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,3,2,2) = t5rr7(3,3,3,2,2)
c      t5rr7(3,3,3,3,2) = t5rr7(3,3,3,3,2)
c      t5rr7(3,3,1,1,1) = t5rr7(3,3,1,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,3,1,1) = t5rr7(3,3,3,1,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,3,1,1) = t5rr7(3,3,3,1,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,3,3,1) = t5rr7(3,3,3,3,1)
c      t5rr7(3,3,2,1,1) = t5rr7(3,3,2,1,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,2,2,1) = t5rr7(3,3,2,2,1)
c      t5rr7(3,3,2,2,2) = t5rr7(3,3,2,2,2)
c      t5rr7(3,3,3,2,2) = t5rr7(3,3,3,2,2)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,3,2,2) = t5rr7(3,3,3,2,2)
c      t5rr7(3,3,3,3,2) = t5rr7(3,3,3,3,2)
c      t5rr7(3,3,3,1,1) = t5rr7(3,3,3,1,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,3,3,1) = t5rr7(3,3,3,3,1)
c      t5rr7(3,3,3,2,1) = t5rr7(3,3,3,2,1)
c      t5rr7(3,3,3,2,2) = t5rr7(3,3,3,2,2)
c      t5rr7(3,3,3,3,2) = t5rr7(3,3,3,3,2)
c      t5rr7(3,3,3,3,1) = t5rr7(3,3,3,3,1)
c      t5rr7(3,3,3,3,2) = t5rr7(3,3,3,3,2)
      t5rr7(3,3,3,3,3) = -15.0d0*zr*rr7
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine t5matrixrr9  --  generate the rr9 part of the undamped ##
c     ##                              t5 multipole interaction matrix       ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "t5matrixrr9" generates the rr9 part of the undamped rank 5 
c     multipole interaction matrix, (del)^5(1/r)
c
      subroutine t5matrixrr9(xr,yr,zr,rr9,t5rr9)
      use sizes
      implicit none
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 xr3,yr3,zr3
      real*8 xr4,yr4,zr4
      real*8 rr9
      real*8 t5rr9(3,3,3,3,3)
c
c     prefactors already applied to rr
c
      xr2 = xr*xr
      yr2 = yr*yr
      zr2 = zr*zr
      xr3 = xr2*xr
      yr3 = yr2*yr
      zr3 = zr2*zr
      xr4 = xr2*xr2
      yr4 = yr2*yr2
      zr4 = zr2*zr2
      t5rr9(1,1,1,1,1) = 10.0d0*xr3*rr9
      t5rr9(2,1,1,1,1) = 6.0d0*xr2*yr*rr9
      t5rr9(3,1,1,1,1) = 6.0d0*xr2*zr*rr9
c      t5rr9(2,1,1,1,1) = t5rr9(2,1,1,1,1)
      t5rr9(2,2,1,1,1) = xr3*rr9 + 3.0d0*xr*yr2*rr9 
      t5rr9(3,2,1,1,1) = 3.0d0*xr*yr*zr*rr9
c      t5rr9(3,1,1,1,1) = t5rr9(3,1,1,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
      t5rr9(3,3,1,1,1) = xr3*rr9 + 3.0d0*xr*zr2*rr9
c      t5rr9(2,1,1,1,1) = t5rr9(2,1,1,1,1)
c      t5rr9(2,2,1,1,1) = t5rr9(2,2,1,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(2,2,1,1,1) = t5rr9(2,2,1,1,1)
      t5rr9(2,2,2,1,1) = 3.0d0*xr2*yr*rr9 + yr3*rr9
      t5rr9(3,2,2,1,1) = xr2*zr*rr9 + yr2*zr*rr9
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
      t5rr9(3,3,2,1,1) = xr2*yr*rr9 + yr*zr2*rr9
c      t5rr9(3,1,1,1,1) = t5rr9(3,1,1,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,3,1,1,1) = t5rr9(3,3,1,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,1,1,1) = t5rr9(3,3,1,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
      t5rr9(3,3,3,1,1) = 3.0d0*xr2*zr*rr9 + zr3*rr9 
c      t5rr9(2,1,1,1,1) = t5rr9(2,1,1,1,1)
c      t5rr9(2,2,1,1,1) = t5rr9(2,2,1,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(2,2,1,1,1) = t5rr9(2,2,1,1,1)
c      t5rr9(2,2,2,1,1) = t5rr9(2,2,2,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(2,2,1,1,1) = t5rr9(2,2,1,1,1)
c      t5rr9(2,2,2,1,1) = t5rr9(2,2,2,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(2,2,2,1,1) = t5rr9(2,2,2,1,1)
      t5rr9(2,2,2,2,1) = 6.0d0*xr*yr2*rr9
      t5rr9(3,2,2,2,1) = 3.0d0*xr*yr*zr*rr9
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
      t5rr9(3,3,2,2,1) = xr*yr2*rr9 + xr*zr2*rr9
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
      t5rr9(3,3,3,2,1) = 3.0d0*xr*yr*zr*rr9
c      t5rr9(3,1,1,1,1) = t5rr9(3,1,1,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,3,1,1,1) = t5rr9(3,3,1,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,1,1,1) = t5rr9(3,3,1,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,3,1,1) = t5rr9(3,3,3,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,1,1,1) = t5rr9(3,3,1,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,3,1,1) = t5rr9(3,3,3,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,3,1,1) = t5rr9(3,3,3,1,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
      t5rr9(3,3,3,3,1) = 6.0d0*xr*zr2*rr9
c
c      t5rr9(2,1,1,1,1) = t5rr9(2,1,1,1,1)
c      t5rr9(2,2,1,1,1) = t5rr9(2,2,1,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(2,2,1,1,1) = t5rr9(2,2,1,1,1)
c      t5rr9(2,2,2,1,1) = t5rr9(2,2,2,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(2,2,1,1,1) = t5rr9(2,2,1,1,1)
c      t5rr9(2,2,2,1,1) = t5rr9(2,2,2,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(2,2,2,1,1) = t5rr9(2,2,2,1,1)
c      t5rr9(2,2,2,2,1) = t5rr9(2,2,2,2,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(2,2,1,1,1) = t5rr9(2,2,1,1,1)
c      t5rr9(2,2,2,1,1) = t5rr9(2,2,2,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(2,2,2,1,1) = t5rr9(2,2,2,1,1)
c      t5rr9(2,2,2,2,1) = t5rr9(2,2,2,2,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(2,2,2,1,1) = t5rr9(2,2,2,1,1)
c      t5rr9(2,2,2,2,1) = t5rr9(2,2,2,2,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(2,2,2,2,1) = t5rr9(2,2,2,2,1)
      t5rr9(2,2,2,2,2) = 10.0d0*yr3*rr9
      t5rr9(3,2,2,2,2) = 6.0d0*yr2*zr*rr9
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,2,2,2,2) = t5rr9(3,2,2,2,2)
      t5rr9(3,3,2,2,2) = yr3*rr9 + 3.0d0*yr*zr2*rr9
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,2,2,2,2) = t5rr9(3,2,2,2,2)
c      t5rr9(3,3,2,2,2) = t5rr9(3,3,2,2,2)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,2,2) = t5rr9(3,3,2,2,2)
      t5rr9(3,3,3,2,2) = 3.0d0*yr2*zr*rr9 + zr3*rr9
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,2,2,2,2) = t5rr9(3,2,2,2,2)
c      t5rr9(3,3,2,2,2) = t5rr9(3,3,2,2,2)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,2,2) = t5rr9(3,3,2,2,2)
c      t5rr9(3,3,3,2,2) = t5rr9(3,3,3,2,2)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,2,2) = t5rr9(3,3,2,2,2)
c      t5rr9(3,3,3,2,2) = t5rr9(3,3,3,2,2)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,3,2,2) = t5rr9(3,3,3,2,2)
      t5rr9(3,3,3,3,2) = 6.0d0*yr*zr2*rr9
c
c      t5rr9(3,1,1,1,1) = t5rr9(3,1,1,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,3,1,1,1) = t5rr9(3,3,1,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,1,1,1) = t5rr9(3,3,1,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,3,1,1) = t5rr9(3,3,3,1,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,1,1,1) = t5rr9(3,3,1,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,3,1,1) = t5rr9(3,3,3,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,3,1,1) = t5rr9(3,3,3,1,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,3,3,1) = t5rr9(3,3,3,3,1)
c      t5rr9(3,2,1,1,1) = t5rr9(3,2,1,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,2,2,1,1) = t5rr9(3,2,2,1,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,2,2,2,1) = t5rr9(3,2,2,2,1)
c      t5rr9(3,2,2,2,2) = t5rr9(3,2,2,2,2)
c      t5rr9(3,3,2,2,2) = t5rr9(3,3,2,2,2)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,2,2) = t5rr9(3,3,2,2,2)
c      t5rr9(3,3,3,2,2) = t5rr9(3,3,3,2,2)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,2,2) = t5rr9(3,3,2,2,2)
c      t5rr9(3,3,3,2,2) = t5rr9(3,3,3,2,2)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,3,2,2) = t5rr9(3,3,3,2,2)
c      t5rr9(3,3,3,3,2) = t5rr9(3,3,3,3,2)
c      t5rr9(3,3,1,1,1) = t5rr9(3,3,1,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,3,1,1) = t5rr9(3,3,3,1,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,3,1,1) = t5rr9(3,3,3,1,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,3,3,1) = t5rr9(3,3,3,3,1)
c      t5rr9(3,3,2,1,1) = t5rr9(3,3,2,1,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,2,2,1) = t5rr9(3,3,2,2,1)
c      t5rr9(3,3,2,2,2) = t5rr9(3,3,2,2,2)
c      t5rr9(3,3,3,2,2) = t5rr9(3,3,3,2,2)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,3,2,2) = t5rr9(3,3,3,2,2)
c      t5rr9(3,3,3,3,2) = t5rr9(3,3,3,3,2)
c      t5rr9(3,3,3,1,1) = t5rr9(3,3,3,1,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,3,3,1) = t5rr9(3,3,3,3,1)
c      t5rr9(3,3,3,2,1) = t5rr9(3,3,3,2,1)
c      t5rr9(3,3,3,2,2) = t5rr9(3,3,3,2,2)
c      t5rr9(3,3,3,3,2) = t5rr9(3,3,3,3,2)
c      t5rr9(3,3,3,3,1) = t5rr9(3,3,3,3,1)
c      t5rr9(3,3,3,3,2) = t5rr9(3,3,3,3,2)
      t5rr9(3,3,3,3,3) = 10.0d0*zr3*rr9
      return
      end
c
c     ##########################################################################
c     ##                                                                      ##
c     ##  subroutine t5matrixrr11  --  generate the rr11 part of the undamped ##
c     ##                                t5 multipole interaction matrix       ##
c     ##                                                                      ##
c     ##########################################################################
c
c
c     "t5matrixrr11" generates the rr11 part of the undamped rank 5 
c     multipole interaction matrix, (del)^5(1/r)
c
      subroutine t5matrixrr11(xr,yr,zr,rr11,t5rr11)
      use sizes
      implicit none
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 xr3,yr3,zr3
      real*8 xr4,yr4,zr4
      real*8 rr11
      real*8 t5rr11(3,3,3,3,3)
c
c     prefactors already applied to rr
c
      xr2 = xr*xr
      yr2 = yr*yr
      zr2 = zr*zr
      xr3 = xr2*xr
      yr3 = yr2*yr
      zr3 = zr2*zr
      xr4 = xr2*xr2
      yr4 = yr2*yr2
      zr4 = zr2*zr2
      t5rr11(1,1,1,1,1) = -xr4*xr*rr11 
      t5rr11(2,1,1,1,1) = -xr4*yr*rr11 
      t5rr11(3,1,1,1,1) = -xr4*zr*rr11 
c      t5rr11(2,1,1,1,1) = t5rr11(2,1,1,1,1)
      t5rr11(2,2,1,1,1) = -xr3*yr2*rr11 
      t5rr11(3,2,1,1,1) = -xr3*yr*zr*rr11
c      t5rr11(3,1,1,1,1) = t5rr11(3,1,1,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
      t5rr11(3,3,1,1,1) = -xr3*zr2*rr11 
c      t5rr11(2,1,1,1,1) = t5rr11(2,1,1,1,1)
c      t5rr11(2,2,1,1,1) = t5rr11(2,2,1,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(2,2,1,1,1) = t5rr11(2,2,1,1,1)
      t5rr11(2,2,2,1,1) = -xr2*yr3*rr11 
      t5rr11(3,2,2,1,1) = -xr2*yr2*zr*rr11 
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
      t5rr11(3,3,2,1,1) = -xr2*yr*zr2*rr11 
c      t5rr11(3,1,1,1,1) = t5rr11(3,1,1,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,3,1,1,1) = t5rr11(3,3,1,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,1,1,1) = t5rr11(3,3,1,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
      t5rr11(3,3,3,1,1) = -xr2*zr3*rr11 
c      t5rr11(2,1,1,1,1) = t5rr11(2,1,1,1,1)
c      t5rr11(2,2,1,1,1) = t5rr11(2,2,1,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(2,2,1,1,1) = t5rr11(2,2,1,1,1)
c      t5rr11(2,2,2,1,1) = t5rr11(2,2,2,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(2,2,1,1,1) = t5rr11(2,2,1,1,1)
c      t5rr11(2,2,2,1,1) = t5rr11(2,2,2,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(2,2,2,1,1) = t5rr11(2,2,2,1,1)
      t5rr11(2,2,2,2,1) = -xr*yr4*rr11 
      t5rr11(3,2,2,2,1) = -xr*yr3*zr*rr11 
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
      t5rr11(3,3,2,2,1) = -xr*yr2*zr2*rr11 
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
      t5rr11(3,3,3,2,1) = -xr*yr*zr3*rr11 
c      t5rr11(3,1,1,1,1) = t5rr11(3,1,1,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,3,1,1,1) = t5rr11(3,3,1,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,1,1,1) = t5rr11(3,3,1,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,3,1,1) = t5rr11(3,3,3,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,1,1,1) = t5rr11(3,3,1,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,3,1,1) = t5rr11(3,3,3,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,3,1,1) = t5rr11(3,3,3,1,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
      t5rr11(3,3,3,3,1) = -xr*zr4*rr11 
c
c      t5rr11(2,1,1,1,1) = t5rr11(2,1,1,1,1)
c      t5rr11(2,2,1,1,1) = t5rr11(2,2,1,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(2,2,1,1,1) = t5rr11(2,2,1,1,1)
c      t5rr11(2,2,2,1,1) = t5rr11(2,2,2,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(2,2,1,1,1) = t5rr11(2,2,1,1,1)
c      t5rr11(2,2,2,1,1) = t5rr11(2,2,2,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(2,2,2,1,1) = t5rr11(2,2,2,1,1)
c      t5rr11(2,2,2,2,1) = t5rr11(2,2,2,2,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(2,2,1,1,1) = t5rr11(2,2,1,1,1)
c      t5rr11(2,2,2,1,1) = t5rr11(2,2,2,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(2,2,2,1,1) = t5rr11(2,2,2,1,1)
c      t5rr11(2,2,2,2,1) = t5rr11(2,2,2,2,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(2,2,2,1,1) = t5rr11(2,2,2,1,1)
c      t5rr11(2,2,2,2,1) = t5rr11(2,2,2,2,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(2,2,2,2,1) = t5rr11(2,2,2,2,1)
      t5rr11(2,2,2,2,2) = -yr4*yr*rr11 
      t5rr11(3,2,2,2,2) = -yr4*zr*rr11 
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,2,2,2,2) = t5rr11(3,2,2,2,2)
      t5rr11(3,3,2,2,2) = -yr3*zr2*rr11 
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,2,2,2,2) = t5rr11(3,2,2,2,2)
c      t5rr11(3,3,2,2,2) = t5rr11(3,3,2,2,2)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,2,2) = t5rr11(3,3,2,2,2)
      t5rr11(3,3,3,2,2) = -yr2*zr3*rr11 
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,2,2,2,2) = t5rr11(3,2,2,2,2)
c      t5rr11(3,3,2,2,2) = t5rr11(3,3,2,2,2)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,2,2) = t5rr11(3,3,2,2,2)
c      t5rr11(3,3,3,2,2) = t5rr11(3,3,3,2,2)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,2,2) = t5rr11(3,3,2,2,2)
c      t5rr11(3,3,3,2,2) = t5rr11(3,3,3,2,2)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,3,2,2) = t5rr11(3,3,3,2,2)
      t5rr11(3,3,3,3,2) = -yr*zr4*rr11 
c
c      t5rr11(3,1,1,1,1) = t5rr11(3,1,1,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,3,1,1,1) = t5rr11(3,3,1,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,1,1,1) = t5rr11(3,3,1,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,3,1,1) = t5rr11(3,3,3,1,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,1,1,1) = t5rr11(3,3,1,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,3,1,1) = t5rr11(3,3,3,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,3,1,1) = t5rr11(3,3,3,1,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,3,3,1) = t5rr11(3,3,3,3,1)
c      t5rr11(3,2,1,1,1) = t5rr11(3,2,1,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,2,2,1,1) = t5rr11(3,2,2,1,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,2,2,2,1) = t5rr11(3,2,2,2,1)
c      t5rr11(3,2,2,2,2) = t5rr11(3,2,2,2,2)
c      t5rr11(3,3,2,2,2) = t5rr11(3,3,2,2,2)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,2,2) = t5rr11(3,3,2,2,2)
c      t5rr11(3,3,3,2,2) = t5rr11(3,3,3,2,2)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,2,2) = t5rr11(3,3,2,2,2)
c      t5rr11(3,3,3,2,2) = t5rr11(3,3,3,2,2)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,3,2,2) = t5rr11(3,3,3,2,2)
c      t5rr11(3,3,3,3,2) = t5rr11(3,3,3,3,2)
c      t5rr11(3,3,1,1,1) = t5rr11(3,3,1,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,3,1,1) = t5rr11(3,3,3,1,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,3,1,1) = t5rr11(3,3,3,1,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,3,3,1) = t5rr11(3,3,3,3,1)
c      t5rr11(3,3,2,1,1) = t5rr11(3,3,2,1,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,2,2,1) = t5rr11(3,3,2,2,1)
c      t5rr11(3,3,2,2,2) = t5rr11(3,3,2,2,2)
c      t5rr11(3,3,3,2,2) = t5rr11(3,3,3,2,2)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,3,2,2) = t5rr11(3,3,3,2,2)
c      t5rr11(3,3,3,3,2) = t5rr11(3,3,3,3,2)
c      t5rr11(3,3,3,1,1) = t5rr11(3,3,3,1,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,3,3,1) = t5rr11(3,3,3,3,1)
c      t5rr11(3,3,3,2,1) = t5rr11(3,3,3,2,1)
c      t5rr11(3,3,3,2,2) = t5rr11(3,3,3,2,2)
c      t5rr11(3,3,3,3,2) = t5rr11(3,3,3,3,2)
c      t5rr11(3,3,3,3,1) = t5rr11(3,3,3,3,1)
c      t5rr11(3,3,3,3,2) = t5rr11(3,3,3,3,2)
      t5rr11(3,3,3,3,3) = -zr4*zr*rr11 
      return
      end
