c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kextra  --  extra term parameter assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kextra" assigns parameters to any additional user defined
c     potential energy contribution
c
c
      subroutine kextra
      use atoms
      use atomid
      use keys
      use potent
      use polar
      use chgpen
c      use xtrpot
      implicit none
      integer i,k,ii,j,jj,next,m
      real*8 alpha_atom(112)
      real*8 zeta_atom(112)
      real*8 factor
      real*8 polf,dampf
      real*8 pen
      integer c1(29),c2(19),c3(26),c4(3),c5(1)
      integer c6(3),c7(25),c8(6),c9(37)
      integer c10(2),c11(3),c12(3),c13(7)
      integer c14(3),c15(8),c16(3),c17(2),c18(8)
      character*20 keyword
      character*120 record
      character*120 string
c
c     setup atomic multipole parameters if not already done
c
c      if (.not. use_mpole)  call kmpole
c
c     initialize cp parameter arrays
c
      alphaf = 1.0d0
      do k = 1, 18
c         cpclass(k) = 0.0d0
c         alpha(k) = 0.0d0
c         regular(k) = 0.0d0
c         xpolr(k) = 0.0d0
      end do
c
c     check for keyword for charge penetration form
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:12) .eq. 'PENETRATION ') then
c            call getword (record,penetration,next)
            read (string,*,err=10,end=10)  penetration
         end if
         if (keyword(1:15) .eq. 'DIRECT-DAMPING ') then
            read (string,*,err=10,end=10)  directdamp
         end if
         if (keyword(1:15) .eq. 'MUTUAL-DAMPING ') then
            read (string,*,err=10,end=10)  mutualdamp
         end if
         if (keyword(1:15) .eq. 'USE-MUSCALE ') then
            use_muscale = .true.
         end if
         if (keyword(1:8) .eq. 'BFACTOR ') then
            read (string,*,err=10,end=10)  bfactor_mode
         end if
         if (keyword(1:12) .eq. 'CHGPEN-MODE ') then
            read (string,*,err=10,end=10)  chgpen_mode
         end if
         if (keyword(1:9) .eq. 'CORE-ELE ') then
            read (string,*,err=10,end=10)  num_ele
         end if
         if (keyword(1:11) .eq. 'REGULARIZE ') then
            read (string,*,err=10,end=10)  regularize
         end if
 10      continue
      end do
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
c     read in atomic polarizability by class 
         if (keyword(1:6) .eq. 'XPOLR ') then
            ii = 0
            pen = 0.0d0
            string = record(next:120)
            read (string,*,err=20,end=20)  ii,pen
 20         continue
            if (ii .ne. 0)  xpolr(ii) = pen
         end if
c     read in polfactors by class  
         if (keyword(1:10) .eq. 'POLFACTOR ') then
            ii = 0
            pen = 0.0d0
            string = record(next:120)
            read (string,*,err=30,end=30)  ii,pen
 30         continue
            if (ii .ne. 0)  polfactor(ii) = pen
         end if
         if (keyword(1:11) .eq. 'ALPHA-PERM ') then
            ii = 0
            pen = 0.0d0
            string = record(next:120)
            read (string,*,err=40,end=40)  ii,pen
 40         continue
            if (ii .ne. 0)  alpha(ii) = pen
         end if
         if (keyword(1:13) .eq. 'ALPHA-FACTOR ') then
            ii = 0
            pen = 0.0d0
            string = record(next:120)
            read (string,*,err=41,end=41)  pen
 41         continue
            alphaf = pen
         end if
         if (keyword(1:18) .eq. 'FUNC-10 ') then
            ii = 0
            pen = 0.0d0
            string = record(next:120)
            read (string,*,err=42,end=42)  ii,pen
 42         continue
            if (ii .ne. 0)  func10(ii) = pen
         end if
         if (keyword(1:8) .eq. 'REGULAR ') then
            ii = 0
            pen = 0.0d0
            string = record(next:120)
            read (string,*,err=50,end=50)  ii,pen
 50         continue
            if (ii .ne. 0)  regular(ii) = pen
         end if
         if (keyword(1:8) .eq. 'CPCLASS ') then
            ii = 0
            pen = 0.0d0
            string = record(next:120)
            read (string,*,err=60,end=60)  ii,pen
 60         continue
            do k = 1, n
               if (type(k).eq.ii) cpclass(k) = pen
            end do
         end if
         if (keyword(1:11) .eq. 'EX-OVERLAP ') then
            ii = 0
            pen = 0.0d0
            string = record(next:120)
            read (string,*,err=70,end=70)  ii,pen
 70         continue
            if (ii .ne. 0)  overlap(ii) = pen
         end if
      end do
c     
c     set default parameters for gordon charge penetration form
c     
c      if (penetration.eq.'GORDON ') then
c         alpha_atom(1) = 5.00d0
c         alpha_atom(6) = 3.00d0
c         alpha_atom(7) = 2.90d0
c         alpha_atom(8) = 2.80d0
c         alpha_atom(9) = 2.50d0
c         alpha_atom(15) = 2.90d0
c         alpha_atom(16) = 2.80d0
c         alpha_atom(17) = 2.50d0
c         alpha_atom(35) = 2.25d0
cccc
c     hack to test CP-polarization
cccc
         c1 = (/ 404,402,406,408,41,410,411,68,109,
     &        155,415,421,445,449,446,533,529,537,538,545,546,547,
     &        465,468,471,498,501,504,507 /)
c         c1 = (/ 404,402,406,408,41,410,411,68,109,
c     &        155,415,421,445,449,446,533,529,537,538,545,546,547,
c     &        465,468,471,498,501,602,702 /)
         c2 = (/ 2,430,438,107,419,442,39,66,151,514,521,558,
     &        454,458,462,526,539,548,508 /)
         c3 = (/ 218,424,426,428,434,436,515,516,517,522,
     &        523,559,560,561,562,563,564,477,478,479,485,486,487,
     &        493,494,495 /)
         c4 = (/ 451,455,459 /)
         c5 = (/ 448 /)
         c6 = (/ 470,492,503 /)
         c7 = (/ 108,40,67,154,401,413,412,407,409,405,414,
     &        420,506,527,534,535,531,443,447,463,466,469,496,499,502 /)
c         c7 = (/ 108,40,67,154,401,413,412,407,409,405,414,
c     &        420,506,527,534,535,531,443,447,463,466,469,496,499,601 /)
         c8 = (/ 104,152,403,416,524,528 /)
         c9 = (/ 217,423,425,427,431,433,435,439,509,510,512,
     &        518,519,540,541,542,543,549,550,551,552,553,554,556,557,
     &        472,473,474,475,480,481,482,483,488,489,490,491 /)
         c10 = (/ 532,444 /)
         c11 = (/ 467,484,500 /)
         c12 = (/ 106,418,525 /)
         c13 = (/ 422,429,437,511,513,520,555 /)
         c14 = (/ 65,505,536 /)
c         c14 = (/ 65,505,701 /)
         c15 = (/ 1,150,441,38,544,452,456,460 /)
         c16 = (/ 464,476,497 /)
         c17 = (/ 432,440 /)
         c18 = (/ 105,153,417,530,450,453,457,461 /)
         polf = 1.0d0
         dampf = 1.0d0
         do i = 1, n
            ii = type(i)
c     alpha values from Gordon LSQ fit by element
c     H non-polar
            if (any(c1==ii)) then 
c               if (alpha(1).eq.0.0d0) alpha(1) = 4.0d0
               if (alpha(1).eq.0.0d0) alpha(1) = 3.2484d0
c               if (alpha(1).eq.0.0d0) alpha(1) = 3.51814d0
               if (alphap(1).eq.0.0d0) alphap(1) = 2.2233d0
c               if (xpolr(1).eq.0.0d0) xpolr(1) = 0.5588d0
c               if (xpolr(1).eq.0.0d0) xpolr(1) = 0.5588d0 * polf
c               if (polfactor(1).eq.0.0d0) polfactor(1) = 2.5948d0*dampf
c
c     values from lbfgs-POL-Gordon-allpolfactors_0.8.log_step0.001 iteration 171
c
               if (xpolr(1).eq.0.0d0) xpolr(1) = 0.1497d0
               if (polfactor(1).eq.0.0d0) polfactor(1) = 4.0758d0
               cpclass(i) = 1
c     H polar
            else if (any(c2==ii)) then
c               if (alpha(2).eq.0.0d0) alpha(2) = 4.0d0
               if (alpha(2).eq.0.0d0) alpha(2) = 3.2632d0
c               if (alpha(2).eq.0.0d0) alpha(2) = 4.99635d0
               if (alphap(2).eq.0.0d0) alphap(2) = 2.1774d0
c               if (xpolr(2).eq.0.0d0) xpolr(2) = 0.4762d0
c               if (xpolr(2).eq.0.0d0) xpolr(2) = 0.5714d0 * polf
c               if (polfactor(2).eq.0.0d0) polfactor(2) = 2.6123d0*dampf
c               if (xpolr(2).eq.0.0d0) xpolr(2) = 0.0385d0
c               if (xpolr(2).eq.0.0d0) xpolr(2) = 0.0000d0
c               if (polfactor(2).eq.0.0d0) polfactor(2) = 5.6998d0
c     test for PO4H3
c               if (xpolr(2).eq.0.0d0) xpolr(2) = 0.7714d0
c               if (xpolr(2).eq.0.0d0) xpolr(2) = 0.1714d0
c               if (polfactor(2).eq.0.0d0) polfactor(2) = 3000.6123d0
c     test for water
c               if (xpolr(2).eq.0.0d0) xpolr(2) = 0.1723d0 
c               if (polfactor(2).eq.0.0d0) polfactor(2) = 5.5250d0
c     smith+s101
c               if (xpolr(2).eq.0.0d0) xpolr(2) = 0.4741d0
c               if (polfactor(2).eq.0.0d0) polfactor(2) = 3.4631d0
c     fit_polarizability values
               if (xpolr(2).eq.0.0d0) xpolr(2) = 0.4283d0
c               if (polfactor(2).eq.0.0d0) polfactor(2) = 1.5500d0
               if (polfactor(2).eq.0.0d0) polfactor(2) = 5.1d0
               cpclass(i) = 2
c     H aromatic
            else if (any(c3==ii)) then 
c               if (alpha(3).eq.0.0d0) alpha(3) = 4.0d0
               if (alpha(3).eq.0.0d0) alpha(3) = 3.4437d0
c               if (alpha(3).eq.0.0d0) alpha(3) = 3.932770
               if (alphap(3).eq.0.0d0) alphap(3) = 2.3762d0
c               if (xpolr(3).eq.0.0d0) xpolr(3) = 0.6405d0
c               if (xpolr(3).eq.0.0d0) xpolr(3) = 0.6346d0 * polf
c               if (polfactor(3).eq.0.0d0) polfactor(3) = 2.9157d0*dampf
               if (xpolr(3).eq.0.0d0) xpolr(3) = 0.3629d0
               if (polfactor(3).eq.0.0d0) polfactor(3) = 3.6035d0
               cpclass(i) = 3
c     P in Phosphate
            else if (any(c4==ii)) then
c               if (alpha(4).eq.0.0d0) alpha(4) = 3.07d0
               if (alpha(4).eq.0.0d0) alpha(4) = 2.7476d0
c               if (alpha(4).eq.0.0d0) alpha(4) = 3.10969d0
               if (alphap(4).eq.0.0d0) alphap(4) = 1.9386d0
c               if (xpolr(4).eq.0.0d0) xpolr(4) = 0.3656d0
c               if (xpolr(4).eq.0.0d0) xpolr(4) = 0.4135d0 * polf
c               if (polfactor(4).eq.0.0d0) polfactor(4) = 2.1221d0*dampf
               if (xpolr(4).eq.0.0d0) xpolr(4) = 0.5918d0
               if (polfactor(4).eq.0.0d0) polfactor(4) = 2.4756d0
               cpclass(i) = 4
c     S in DMSO
            else if (any(c5==ii)) then 
c               if (alpha(5).eq.0.0d0) alpha(5) = 2.96d0
               if (alpha(5).eq.0.0d0) alpha(5) = 2.6247d0
c               if (alpha(5).eq.0.0d0) alpha(5) = 3.29611d0
               if (alphap(5).eq.0.0d0) alphap(5) = 3.1977d0
c               if (xpolr(5).eq.0.0d0) xpolr(5) = 1.5748d0
c               if (xpolr(5).eq.0.0d0) xpolr(5) = 1.6143d0 * polf
c               if (polfactor(5).eq.0.0d0) polfactor(5) = 2.0974d0*dampf
               if (xpolr(5).eq.0.0d0) xpolr(5) = 1.7085d0
               if (polfactor(5).eq.0.0d0) polfactor(5) = 2.2328d0
               cpclass(i) = 5
c     Br
            else if (any(c6==ii)) then
c               if (alpha(6).eq.0.0d0) alpha(6) = 3.72d0
               if (alpha(6).eq.0.0d0) alpha(6) = 3.6696d0
c               if (alpha(6).eq.0.0d0) alpha(6) = 3.49409d0
               if (alphap(6).eq.0.0d0) alphap(6) = 3.3540d0
c               if (xpolr(6).eq.0.0d0) xpolr(6) = 1.8793d0
c               if (xpolr(6).eq.0.0d0) xpolr(6) = 1.8934d0 * polf
c               if (polfactor(6).eq.0.0d0) polfactor(6) = 2.4807d0*dampf
               if (xpolr(6).eq.0.0d0) xpolr(6) = 1.8507d0
               if (polfactor(6).eq.0.0d0) polfactor(6) = 2.3719d0
               cpclass(i) = 6
c     C Sp3
            else if (any(c7==ii)) then 
c               if (alpha(7).eq.0.0d0) alpha(7) = 3.10d0
               if (alpha(7).eq.0.0d0) alpha(7) = 3.5898d0
c               if (alpha(7).eq.0.0d0) alpha(7) = 3.52327d0
               if (alphap(7).eq.0.0d0) alphap(7) = 1.7416d0
c               if (xpolr(7).eq.0.0d0) xpolr(7) = 0.2668d0
c               if (xpolr(7).eq.0.0d0) xpolr(7) = 0.3036d0 * polf
c               if (polfactor(7).eq.0.0d0) polfactor(7) = 2.4927d0*dampf
               if (xpolr(7).eq.0.0d0) xpolr(7) = 0.7299d0
               if (polfactor(7).eq.0.0d0) polfactor(7) = 2.9134d0
               cpclass(i) = 7
c     C Sp2
            else if (any(c8==ii)) then 
c               if (alpha(8).eq.0.0d0) alpha(8) = 3.10d0
               if (alpha(8).eq.0.0d0) alpha(8) = 3.1286d0
c               if (alpha(8).eq.0.0d0) alpha(8) = 3.20696d0
               if (alphap(8).eq.0.0d0) alphap(8) = 1.7500d0
c               if (xpolr(8).eq.0.0d0) xpolr(8) = 0.2668d0
c               if (xpolr(8).eq.0.0d0) xpolr(8) = 0.2946d0 * polf
c               if (polfactor(8).eq.0.0d0) polfactor(8) = 2.4943d0*dampf
               if (xpolr(8).eq.0.0d0) xpolr(8) = 0.6532d0
               if (polfactor(8).eq.0.0d0) polfactor(8) = 2.6259d0
               cpclass(i) = 8
c     C aromatic
            else if (any(c9==ii)) then 
c               if (alpha(9).eq.0.0d0) alpha(9) = 3.10d0
               if (alpha(9).eq.0.0d0) alpha(9) = 3.2057d0
c               if (alpha(9).eq.0.0d0) alpha(9) = 3.19842d0
               if (alphap(9).eq.0.0d0) alphap(9) = 1.9217d0
c               if (xpolr(9).eq.0.0d0) xpolr(9) = 0.3500d0
c               if (xpolr(9).eq.0.0d0) xpolr(9) = 0.3529d0 * polf
c               if (polfactor(9).eq.0.0d0) polfactor(9) = 2.5030d0*dampf
               if (xpolr(9).eq.0.0d0) xpolr(9) = 0.3064d0
               if (polfactor(9).eq.0.0d0) polfactor(9) = 2.7392d0
               cpclass(i) = 9
c     S in SH
            else if (any(c10==ii)) then
c               if (alpha(10).eq.0.0d0) alpha(10) = 2.96d0
               if (alpha(10).eq.0.0d0) alpha(10) = 3.3112d0
c               if (alpha(10).eq.0.0d0) alpha(10) = 3.15607d0
               if (alphap(10).eq.0.0d0) alphap(10) = 3.7640d0
c               if (xpolr(10).eq.0.0d0) xpolr(10) = 2.6487d0
c               if (xpolr(10).eq.0.0d0) xpolr(10) = 2.6639d0 * polf
c               if (polfactor(10).eq.0.0d0) polfactor(10) =2.4957d0*dampf
               if (xpolr(10).eq.0.0d0) xpolr(10) = 2.7396d0
               if (polfactor(10).eq.0.0d0) polfactor(10) =2.5600d0
               cpclass(i) = 10
c     Cl
            else if (any(c11==ii)) then
c               if (alpha(11).eq.0.0d0) alpha(11) = 3.5d0
               if (alpha(11).eq.0.0d0) alpha(11) = 3.4749d0
c               if (alpha(11).eq.0.0d0) alpha(11) = 3.44190d0
               if (alphap(11).eq.0.0d0) alphap(11) = 3.1989d0
c               if (xpolr(11).eq.0.0d0) xpolr(11) = 1.6454d0
c               if (xpolr(11).eq.0.0d0) xpolr(11) = 1.6608d0 * polf
c               if (polfactor(11).eq.0.0d0) polfactor(11) =2.5015d0*dampf
               if (xpolr(11).eq.0.0d0) xpolr(11) = 1.7132d0
               if (polfactor(11).eq.0.0d0) polfactor(11) = 2.5723d0
               cpclass(i) = 11
c     N Sp2
            else if (any(c12==ii)) then
c               if (alpha(12).eq.0.0d0) alpha(12) = 3.73d0
               if (alpha(12).eq.0.0d0) alpha(12) = 3.7071d0
c               if (alpha(12).eq.0.0d0) alpha(12) = 3.40725d0
               if (alphap(12).eq.0.0d0) alphap(12) = 1.6549d0
c               if (xpolr(12).eq.0.0d0) xpolr(12) = 0.2146d0
c               if (xpolr(12).eq.0.0d0) xpolr(12) = 0.2146d0 * polf
c               if (polfactor(12).eq.0.0d0) polfactor(12) =2.5148d0*dampf
               if (xpolr(12).eq.0.0d0) xpolr(12) = 0.2073d0
               if (polfactor(12).eq.0.0d0) polfactor(12) =2.1638d0
               cpclass(i) = 12
c     N aromatic
            else if (any(c13==ii)) then 
c               if (alpha(13).eq.0.0d0) alpha(13) = 3.73d0
               if (alpha(13).eq.0.0d0) alpha(13) = 3.6358d0
c               if (alpha(13).eq.0.0d0) alpha(13) = 3.66449d0
               if (alphap(13).eq.0.0d0) alphap(13) = 1.6501d0
c               if (xpolr(13).eq.0.0d0) xpolr(13) = 0.2146d0
c               if (xpolr(13).eq.0.0d0) xpolr(13) = 0.2511d0 * polf
c               if (polfactor(13).eq.0.0d0) polfactor(13) =2.5724d0*dampf
               if (xpolr(13).eq.0.0d0) xpolr(13) = 0.3823d0
               if (polfactor(13).eq.0.0d0) polfactor(13) =2.6790d0
               cpclass(i) = 13
c     N Sp3
            else if (any(c14==ii)) then
c               if (alpha(14).eq.0.0d0) alpha(14) = 3.73d0
               if (alpha(14).eq.0.0d0) alpha(14) = 4.0135d0
c               if (alpha(14).eq.0.0d0) alpha(14) = 3.27558d0
               if (alphap(14).eq.0.0d0) alphap(14) = 2.5446d0
c               if (xpolr(14).eq.0.0d0) xpolr(14) = 0.8127d0
c               if (xpolr(14).eq.0.0d0) xpolr(14) = 0.9143d0 * polf
c               if (polfactor(14).eq.0.0d0) polfactor(14) =2.2893d0*dampf
               if (xpolr(14).eq.0.0d0) xpolr(14) = 1.0679d0
               if (polfactor(14).eq.0.0d0) polfactor(14) =1.8082d0
               cpclass(i) = 14
c     O in OH
            else if (any(c15==ii)) then
c               if (alpha(15).eq.0.0d0) alpha(15) = 4.14d0
               if (alpha(15).eq.0.0d0) alpha(15) = 4.1615d0
c               if (alpha(15).eq.0.0d0) alpha(15) = 3.70318d0
               if (alphap(15).eq.0.0d0) alphap(15) = 2.5417d0
c               if (xpolr(15).eq.0.0d0) xpolr(15) = 0.8035d0
c               if (xpolr(15).eq.0.0d0) xpolr(15) = 0.8732d0 * polf
c               if (polfactor(15).eq.0.0d0) polfactor(15) =2.3190d0*dampf
c               if (xpolr(15).eq.0.0d0) xpolr(15) = 0.8612d0
c               if (polfactor(15).eq.0.0d0) polfactor(15) =2.7628d0
c     test for water
c               if (xpolr(15).eq.0.0d0) xpolr(15) = 1.1934d0
c               if (polfactor(15).eq.0.0d0) polfactor(15) =2.9600d0
c     smith+s101
c               if (xpolr(15).eq.0.0d0) xpolr(15) = 1.0775d0
c               if (polfactor(15).eq.0.0d0) polfactor(15) = 2.1994d0
c     fit_polarizability values
               if (xpolr(15).eq.0.0d0) xpolr(15) = 0.7551d0
c               if (polfactor(15).eq.0.0d0) polfactor(15) = 1.6074d0
               if (polfactor(15).eq.0.0d0) polfactor(15) = 2.1d0
               cpclass(i) = 15
c     F
            else if (any(c16==ii)) then
c               if (alpha(16).eq.0.0d0) alpha(16) = 4.49d0
               if (alpha(16).eq.0.0d0) alpha(16) = 4.4675d0
c               if (alpha(16).eq.0.0d0) alpha(16) = 4.42983d0
               if (alphap(16).eq.0.0d0) alphap(16) = 2.2009d0
c               if (xpolr(16).eq.0.0d0) xpolr(16) = 0.4867d0
c               if (xpolr(16).eq.0.0d0) xpolr(16) = 0.5141d0 * polf
c               if (polfactor(16).eq.0.0d0) polfactor(16) =2.5886d0*dampf
               if (xpolr(16).eq.0.0d0) xpolr(16) = 0.5784d0
               if (polfactor(16).eq.0.0d0) polfactor(16) =2.7087d0
               cpclass(i) = 16
c     O in ur
            else if (any(c17==ii)) then
c               if (alpha(17).eq.0.0d0) alpha(17) = 4.14d0
               if (alpha(17).eq.0.0d0) alpha(17) = 4.3778d0
c               if (alpha(17).eq.0.0d0) alpha(17) = 3.76956d0
               if (alphap(17).eq.0.0d0) alphap(17) = 1.5209d0
c               if (xpolr(17).eq.0.0d0) xpolr(17) = 0.1674d0
c               if (xpolr(17).eq.0.0d0) xpolr(17) = 0.1956d0 * polf
c               if (polfactor(17).eq.0.0d0) polfactor(17) =2.6282d0*dampf
               if (xpolr(17).eq.0.0d0) xpolr(17) = 0.9624d0
               if (polfactor(17).eq.0.0d0) polfactor(17) =2.1330d0
               cpclass(i) = 17
c     O in C=O
            else if (any(c18==ii)) then 
c               if (alpha(18).eq.0.0d0) alpha(18) = 4.14d0
               if (alpha(18).eq.0.0d0) alpha(18) = 3.7321d0
c               if (alpha(18).eq.0.0d0) alpha(18) = 3.35067d0
               if (alphap(18).eq.0.0d0) alphap(18) = 1.5427d0
c               if (xpolr(18).eq.0.0d0) xpolr(18) = 0.1842d0
c               if (xpolr(18).eq.0.0d0) xpolr(18) = 0.2210d0 * polf
c               if (polfactor(18).eq.0.0d0) polfactor(18) =2.6154d0*dampf
               if (xpolr(18).eq.0.0d0) xpolr(18) = 0.2871d0
               if (polfactor(18).eq.0.0d0) polfactor(18) =2.5295d0
               cpclass(i) = 18
            end if
c
c     this allows polarizability and thole damping to vary by cpclass
c            xpolr(cpclass(i)) = polarity(i)
c            xathl(cpclass(i)) = thole(i)
c            polarity(i) = xpolr(cpclass(i))
c            thole(i) = xathl(cpclass(i))
c            alphap(cpclass(i)) = alpha(cpclass(i))
c            alphap(cpclass(i)) = polfactor(cpclass(i)) *
c     &           xpolr(cpclass(i))**(1.0d0/3.0d0)
c            if (xpolr(cpclass(i)).eq.0.0d0) then
c               xpolr(cpclass(i)) = polarity(i)
c            end if
ccccc            polarity(i) = xpolr(cpclass(i))
c            print *,"atom",i,"cpclass",cpclass(i)
c            print *,"xpolr",xpolr(cpclass(i))
c            print *,"polfactor",polfactor(cpclass(i))
c            print *,"alphap",alphap(cpclass(i))
            if (regular(cpclass(i)).eq.0.0d0) regular(cpclass(i))= 3.0d0
         end do
c      end if
      return
      end
