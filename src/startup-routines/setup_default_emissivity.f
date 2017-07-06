      subroutine eirene_setup_default_emissivity

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comsig
      use eirmod_comusr

      implicit none

      TYPE(TCONTRIB) :: CNT
      
      integer :: i, no_compo, iat, iml, ipl, nat, npl, nml

      no_lines = 6
!pbh3+      no_compo = 6  
      no_compo = 5
      MOD_ADDV = 0

      ALLOCATE (EMIS_LINES(NO_LINES))
      EMIS_LINES%LINE_NAME = REPEAT(' ',80)
      EMIS_LINES%NO_COMPO = 0
 

************************************************
* BALMER ALPHA
************************************************

      EMIS_LINES(1)%LINE_NAME = 'BA_ALPHA'
      EMIS_LINES(1)%NO_COMPO = NO_COMPO
      EMIS_LINES(1)%L1 = 3
      EMIS_LINES(1)%L2 = 2
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(1)%FAC = 4.410E7
      EMIS_LINES(1)%ENERGY = 1.8889_DP
      EMIS_LINES(1)%IADV_TOTAL = NADVI + NO_COMPO+1 
      
      ALLOCATE (EMIS_LINES(1)%COMPO(NO_COMPO))

C  CONTRIBUTION LINEAR IN H   -ATOM      DENSITY
C  H(n=3)/H(n=1)
  
      EMIS_LINES(1)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(1)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(1)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(1)%COMPO(1)%NO_CONTRIB = NAT  
      CNT%ISP          = 1
      CNT%ITP          = 1 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.5a   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO
      
C  CONTRIBUTION LINEAR IN H+  -ION       DENSITY
C  H(n=3)/H+
  
      EMIS_LINES(1)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(1)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(1)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(1)%COMPO(2)%NO_CONTRIB = NPL  
      CNT%ISP          = 1
      CNT%ITP          = 4 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.8a   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2  -MOLEC.    DENSITY
C  H(n=3)/H2(g)
  
      EMIS_LINES(1)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(1)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(1)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(1)%COMPO(3)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.5a   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(3)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2+ -MOLEC.ION DENSITY
C  H(n=3)/H2+(g)
  
      EMIS_LINES(1)%COMPO(4)%COMPO_NAME = 
     .     'DIATOMIC NEUTRAL HYDR. MOL ION'
      EMIS_LINES(1)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(1)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(1)%COMPO(4)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.14a   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL  '
      CNT%RAT_H2       = 'H.12' 
      CNT%RAT_REACTION = '2.0c     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(4)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H-  -NEG. ION  DENSITY
C  H(n=3)/H-
  
      EMIS_LINES(1)%COMPO(5)%COMPO_NAME = 
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(1)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(1)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(1)%COMPO(5)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '7.2a     '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL   '
      CNT%RAT_H2       = 'H.11' 
      CNT%RAT_REACTION = '7.0a     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(5)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H3+ -MOL. ION  DENSITY
C  H(n=3)/H3+
  
!pbh3+      EMIS_LINES(1)%COMPO(6)%COMPO_NAME = 
!pbh3+     .     'TRIATOMIC HYDR. ION'
!pbh3+      EMIS_LINES(1)%COMPO(6)%IADV = NADVI + 6

!pbh3+      NML = COUNT(NCHARM == 2)
!pbh3+      ALLOCATE (EMIS_LINES(1)%COMPO(6)%CONTRIB(NML))
!pbh3+      EMIS_LINES(1)%COMPO(6)%NO_CONTRIB = NML  
!pbh3+      CNT%ISP          = 1
!pbh3+      CNT%ITP          = 2
!pbh3+      CNT%IRATIO       = 1
!pbh3+      CNT%FNAME        = 'AMJUEL  '
!pbh3+      CNT%H2           = 'H.12'
!pbh3+      CNT%REACTION     = '2.2.15a  '
!pbh3+      CNT%CR           = 'OT ' 
!pbh3+      CNT%FRATIO       = 'AMJUEL  '
!pbh3+      CNT%RAT_H2       = 'H.11' 
!pbh3+      CNT%RAT_REACTION = '4.0a     '
!pbh3+      CNT%RAT_CR       = 'OT '
!pbh3+      CNT%IRC          = 0
!pbh3+      CNT%IRC_RAT      = 0

!pbh3+      IML = 0
!pbh3+      DO I = 1, NMOLI
!pbh3+        IF (NCHARM(I) == 2) THEN
!pbh3+          IML = IML + 1
!pbh3+          CNT%ISP = I
!pbh3+          EMIS_LINES(1)%COMPO(6)%CONTRIB(IPL) = CNT
!pbh3+        END IF
!pbh3+      END DO


************************************************
* BALMER BETA
************************************************
      
      EMIS_LINES(2)%LINE_NAME = 'BA_BETA'
      EMIS_LINES(2)%NO_COMPO = NO_COMPO
      EMIS_LINES(2)%L1 = 4
      EMIS_LINES(2)%L2 = 2
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(2)%FAC = 8.419E6
      EMIS_LINES(2)%ENERGY = 2.5500_DP
      EMIS_LINES(2)%IADV_TOTAL = NADVI + NO_COMPO+1 
      
      ALLOCATE (EMIS_LINES(2)%COMPO(NO_COMPO))

C  CONTRIBUTION LINEAR IN H   -ATOM      DENSITY
C  H(n=4)/H(n=1)
  
      EMIS_LINES(2)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(2)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(2)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(2)%COMPO(1)%NO_CONTRIB = NAT  
      CNT%ISP          = 1
      CNT%ITP          = 1 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.5c   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO
      
C  CONTRIBUTION LINEAR IN H+  -ION       DENSITY
C  H(n=4)/H+
  
      EMIS_LINES(2)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(2)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(2)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(2)%COMPO(2)%NO_CONTRIB = NPL  
      CNT%ISP          = 1
      CNT%ITP          = 4 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.8c   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2  -MOLEC.    DENSITY
C  H(n=4)/H2(g)
  
      EMIS_LINES(2)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(2)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(2)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(2)%COMPO(3)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.5c   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(3)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2+ -MOLEC.ION DENSITY
C  H(n=4)/H2+(g)
  
      EMIS_LINES(2)%COMPO(4)%COMPO_NAME = 
     .     'DIATOMIC NEUTRAL HYDR. MOL ION'
      EMIS_LINES(2)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(2)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(2)%COMPO(4)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.14c   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL  '
      CNT%RAT_H2       = 'H.12' 
      CNT%RAT_REACTION = '2.0c     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(4)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H-  -NEG. ION  DENSITY
C  H(n=4)/H-
  
      EMIS_LINES(2)%COMPO(5)%COMPO_NAME = 
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(2)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(2)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(2)%COMPO(5)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '7.2c      '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL  '
      CNT%RAT_H2       = 'H.11' 
      CNT%RAT_REACTION = '7.0a     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(5)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H3+ -MOL. ION  DENSITY
C  H(n=4)/H3+
  
!pbh3+      EMIS_LINES(2)%COMPO(6)%COMPO_NAME = 
!pbh3+     .     'TRIATOMIC HYDR. ION'
!pbh3+      EMIS_LINES(2)%COMPO(6)%IADV = NADVI + 6

!pbh3+      NML = COUNT(NCHARM == 2)
!pbh3+      ALLOCATE (EMIS_LINES(2)%COMPO(6)%CONTRIB(NML))
!pbh3+      EMIS_LINES(2)%COMPO(6)%NO_CONTRIB = NML  
!pbh3+      CNT%ISP          = 1
!pbh3+      CNT%ITP          = 2
!pbh3+      CNT%IRATIO       = 1
!pbh3+      CNT%FNAME        = 'AMJUEL  '
!pbh3+      CNT%H2           = 'H.12'
!pbh3+      CNT%REACTION     = '2.2.15c  '
!pbh3+      CNT%CR           = 'OT ' 
!pbh3+      CNT%FRATIO       = 'AMJUEL  '
!pbh3+      CNT%RAT_H2       = 'H.11' 
!pbh3+      CNT%RAT_REACTION = '4.0a     '
!pbh3+      CNT%RAT_CR       = 'OT '
!pbh3+      CNT%IRC          = 0
!pbh3+      CNT%IRC_RAT      = 0

!pbh3+      IML = 0
!pbh3+      DO I = 1, NMOLI
!pbh3+        IF (NCHARM(I) == 2) THEN
!pbh3+          IML = IML + 1
!pbh3+          CNT%ISP = I
!pbh3+          EMIS_LINES(2)%COMPO(6)%CONTRIB(IPL) = CNT
!pbh3+        END IF
!pbh3+      END DO


************************************************
* BALMER GAMMA
************************************************
      
      EMIS_LINES(3)%LINE_NAME = 'BA_GAMMA'
      EMIS_LINES(3)%NO_COMPO = NO_COMPO
      EMIS_LINES(3)%L1 = 5
      EMIS_LINES(3)%L2 = 2
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(3)%FAC = 2.530E6
      EMIS_LINES(3)%ENERGY = 2.8560_DP
      EMIS_LINES(3)%IADV_TOTAL = NADVI + NO_COMPO+1
      
      ALLOCATE (EMIS_LINES(3)%COMPO(NO_COMPO))

C  CONTRIBUTION LINEAR IN H   -ATOM      DENSITY
C  H(n=5)/H(n=1)
  
      EMIS_LINES(3)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(3)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(3)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(3)%COMPO(1)%NO_CONTRIB = NAT  
      CNT%ISP          = 1
      CNT%ITP          = 1 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.5d   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO
      
C  CONTRIBUTION LINEAR IN H+  -ION       DENSITY
C  H(n=5)/H+
  
      EMIS_LINES(3)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(3)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(3)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(3)%COMPO(2)%NO_CONTRIB = NPL  
      CNT%ISP          = 1
      CNT%ITP          = 4 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.8d   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2  -MOLEC.    DENSITY
C  H(n=5)/H2(g)
  
      EMIS_LINES(3)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(3)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(3)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(3)%COMPO(3)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.5d   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(3)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2+ -MOLEC.ION DENSITY
C  H(n=5)/H2+(g)
  
      EMIS_LINES(3)%COMPO(4)%COMPO_NAME = 
     .     'DIATOMIC NEUTRAL HYDR. MOL ION'
      EMIS_LINES(3)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(3)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(3)%COMPO(4)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.14d   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL  '
      CNT%RAT_H2       = 'H.12' 
      CNT%RAT_REACTION = '2.0c     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(4)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H-  -NEG. ION  DENSITY
C  H(n=5)/H-
  
      EMIS_LINES(3)%COMPO(5)%COMPO_NAME = 
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(3)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(3)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(3)%COMPO(5)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '7.2d      '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL  '
      CNT%RAT_H2       = 'H.11' 
      CNT%RAT_REACTION = '7.0a     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(5)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H3+ -MOL. ION  DENSITY
C  H(n=5)/H3+
  
!pbh3+      EMIS_LINES(3)%COMPO(6)%COMPO_NAME = 
!pbh3+     .     'TRIATOMIC HYDR. ION'
!pbh3+      EMIS_LINES(3)%COMPO(6)%IADV = NADVI + 6

!pbh3+      NML = COUNT(NCHARM == 2)
!pbh3+      ALLOCATE (EMIS_LINES(3)%COMPO(6)%CONTRIB(NML))
!pbh3+      EMIS_LINES(3)%COMPO(6)%NO_CONTRIB = NML  
!pbh3+      CNT%ISP          = 1
!pbh3+      CNT%ITP          = 2
!pbh3+      CNT%IRATIO       = 1
!pbh3+      CNT%FNAME        = 'AMJUEL  '
!pbh3+      CNT%H2           = 'H.12'
!pbh3+      CNT%REACTION     = '2.2.15d  '
!pbh3+      CNT%CR           = 'OT ' 
!pbh3+      CNT%FRATIO       = 'AMJUEL  '
!pbh3+      CNT%RAT_H2       = 'H.11' 
!pbh3+      CNT%RAT_REACTION = '4.0a     '
!pbh3+      CNT%RAT_CR       = 'OT '
!pbh3+      CNT%IRC          = 0
!pbh3+      CNT%IRC_RAT      = 0

!pbh3+      IML = 0
!pbh3+      DO I = 1, NMOLI
!pbh3+        IF (NCHARM(I) == 2) THEN
!pbh3+          IML = IML + 1
!pbh3+          CNT%ISP = I
!pbh3+          EMIS_LINES(3)%COMPO(6)%CONTRIB(IPL) = CNT
!pbh3+        END IF
!pbh3+      END DO
      


************************************************
* BALMER DELTA
************************************************
     
      EMIS_LINES(4)%LINE_NAME = 'BA_DELTA'
      EMIS_LINES(4)%NO_COMPO = NO_COMPO
      EMIS_LINES(4)%L1 = 6
      EMIS_LINES(4)%L2 = 2
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(4)%FAC = 9.732E5
      EMIS_LINES(4)%ENERGY = 3.0222_DP
      EMIS_LINES(4)%IADV_TOTAL = NADVI + NO_COMPO+1
      
      ALLOCATE (EMIS_LINES(4)%COMPO(NO_COMPO))

C  CONTRIBUTION LINEAR IN H   -ATOM      DENSITY
C  H(n=6)/H(n=1)
  
      EMIS_LINES(4)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(4)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(4)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(4)%COMPO(1)%NO_CONTRIB = NAT  
      CNT%ISP          = 1
      CNT%ITP          = 1 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.5e   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO
      
C  CONTRIBUTION LINEAR IN H+  -ION       DENSITY
C  H(n=6)/H+
  
      EMIS_LINES(4)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(4)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(4)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(4)%COMPO(2)%NO_CONTRIB = NPL  
      CNT%ISP          = 1
      CNT%ITP          = 4 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.8e   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2  -MOLEC.    DENSITY
C  H(n=6)/H2(g)
  
      EMIS_LINES(4)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(4)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(4)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(4)%COMPO(3)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.5e   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(3)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2+ -MOLEC.ION DENSITY
C  H(n=6)/H2+(g)
  
      EMIS_LINES(4)%COMPO(4)%COMPO_NAME = 
     .     'DIATOMIC NEUTRAL HYDR. MOL ION'
      EMIS_LINES(4)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(4)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(4)%COMPO(4)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.14e   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL  '
      CNT%RAT_H2       = 'H.12' 
      CNT%RAT_REACTION = '2.0c     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(4)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H-  -NEG. ION  DENSITY
C  H(n=6)/H-
  
      EMIS_LINES(4)%COMPO(5)%COMPO_NAME = 
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(4)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(4)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(4)%COMPO(5)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '7.2e      '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL  '
      CNT%RAT_H2       = 'H.11' 
      CNT%RAT_REACTION = '7.0a     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(5)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H3+ -MOL. ION  DENSITY
C  H(n=6)/H3+
  
!pbh3+      EMIS_LINES(4)%COMPO(6)%COMPO_NAME = 
!pbh3+     .     'TRIATOMIC HYDR. ION'
!pbh3+      EMIS_LINES(4)%COMPO(6)%IADV = NADVI + 6

!pbh3+      NML = COUNT(NCHARM == 2)
!pbh3+      ALLOCATE (EMIS_LINES(4)%COMPO(6)%CONTRIB(NML))
!pbh3+      EMIS_LINES(4)%COMPO(6)%NO_CONTRIB = NML  
!pbh3+      CNT%ISP          = 1
!pbh3+      CNT%ITP          = 2
!pbh3+      CNT%IRATIO       = 1
!pbh3+      CNT%FNAME        = 'AMJUEL  '
!pbh3+      CNT%H2           = 'H.12'
!pbh3+      CNT%REACTION     = '2.2.15e  '
!pbh3+      CNT%CR           = 'OT ' 
!pbh3+      CNT%FRATIO       = 'AMJUEL  '
!pbh3+      CNT%RAT_H2       = 'H.11' 
!pbh3+      CNT%RAT_REACTION = '4.0a     '
!pbh3+      CNT%RAT_CR       = 'OT '
!pbh3+      CNT%IRC          = 0
!pbh3+      CNT%IRC_RAT      = 0

!pbh3+      IML = 0
!pbh3+      DO I = 1, NMOLI
!pbh3+        IF (NCHARM(I) == 2) THEN
!pbh3+          IML = IML + 1
!pbh3+          CNT%ISP = I
!pbh3+          EMIS_LINES(4)%COMPO(6)%CONTRIB(IPL) = CNT
!pbh3+        END IF
!pbh3+      END DO     


************************************************
* LYMAN ALPHA
************************************************
      
      EMIS_LINES(5)%LINE_NAME = 'LY_ALPHA'
      EMIS_LINES(5)%NO_COMPO = NO_COMPO
      EMIS_LINES(5)%L1 = 2
      EMIS_LINES(5)%L2 = 1
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(5)%FAC = 4.699E8
      EMIS_LINES(5)%ENERGY = 10.2375_DP
      EMIS_LINES(5)%IADV_TOTAL = NADVI + NO_COMPO+1
      
      ALLOCATE (EMIS_LINES(5)%COMPO(NO_COMPO))

C  CONTRIBUTION LINEAR IN H   -ATOM      DENSITY
C  H(n=2)/H(n=1)
  
      EMIS_LINES(5)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(5)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(5)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(5)%COMPO(1)%NO_CONTRIB = NAT  
      CNT%ISP          = 1
      CNT%ITP          = 1 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.5b   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO
      
C  CONTRIBUTION LINEAR IN H+  -ION       DENSITY
C  H(n=2)/H+
  
      EMIS_LINES(5)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(5)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(5)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(5)%COMPO(2)%NO_CONTRIB = NPL  
      CNT%ISP          = 1
      CNT%ITP          = 4 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.8b   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2  -MOLEC.    DENSITY
C  H(n=2)/H2(g)
  
      EMIS_LINES(5)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(5)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(5)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(5)%COMPO(3)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.5b   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(3)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2+ -MOLEC.ION DENSITY
C  H(n=2)/H2+(g)
  
      EMIS_LINES(5)%COMPO(4)%COMPO_NAME = 
     .     'DIATOMIC NEUTRAL HYDR. MOL ION'
      EMIS_LINES(5)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(5)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(5)%COMPO(4)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.14b   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL  '
      CNT%RAT_H2       = 'H.12' 
      CNT%RAT_REACTION = '2.0c     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(4)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H-  -NEG. ION  DENSITY
C  H(n=2)/H-
  
      EMIS_LINES(5)%COMPO(5)%COMPO_NAME = 
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(5)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(5)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(5)%COMPO(5)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '7.2b      '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL  '
      CNT%RAT_H2       = 'H.11' 
      CNT%RAT_REACTION = '7.0a     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(5)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H3+ -MOL. ION  DENSITY
C  H(n=2)/H3+
  
!pbh3+      EMIS_LINES(5)%COMPO(6)%COMPO_NAME = 
!pbh3+     .     'TRIATOMIC HYDR. ION'
!pbh3+      EMIS_LINES(5)%COMPO(6)%IADV = NADVI + 6

!pbh3+      NML = COUNT(NCHARM == 2)
!pbh3+      ALLOCATE (EMIS_LINES(5)%COMPO(6)%CONTRIB(NML))
!pbh3+      EMIS_LINES(5)%COMPO(6)%NO_CONTRIB = NML  
!pbh3+      CNT%ISP          = 1
!pbh3+      CNT%ITP          = 2
!pbh3+      CNT%IRATIO       = 1
!pbh3+      CNT%FNAME        = 'AMJUEL  '
!pbh3+      CNT%H2           = 'H.12'
!pbh3+      CNT%REACTION     = '2.2.15b  '
!pbh3+      CNT%CR           = 'OT ' 
!pbh3+      CNT%FRATIO       = 'AMJUEL  '
!pbh3+      CNT%RAT_H2       = 'H.11' 
!pbh3+      CNT%RAT_REACTION = '4.0a     '
!pbh3+      CNT%RAT_CR       = 'OT '
!pbh3+      CNT%IRC          = 0
!pbh3+      CNT%IRC_RAT      = 0

!pbh3+      IML = 0
!pbh3+      DO I = 1, NMOLI
!pbh3+        IF (NCHARM(I) == 2) THEN
!pbh3+          IML = IML + 1
!pbh3+          CNT%ISP = I
!pbh3+          EMIS_LINES(5)%COMPO(6)%CONTRIB(IPL) = CNT
!pbh3+        END IF
!pbh3+      END DO


************************************************
* LYMAN BETA
************************************************
      
      EMIS_LINES(6)%LINE_NAME = 'LY_BETA'
      EMIS_LINES(6)%NO_COMPO = NO_COMPO
      EMIS_LINES(6)%L1 = 3
      EMIS_LINES(6)%L2 = 1
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(6)%FAC = 5.575E7
      EMIS_LINES(6)%ENERGY = 12.089_DP
      EMIS_LINES(6)%IADV_TOTAL = NADVI + NO_COMPO+1
      
      ALLOCATE (EMIS_LINES(6)%COMPO(NO_COMPO))

C  CONTRIBUTION LINEAR IN H   -ATOM      DENSITY
C  H(n=3)/H(n=1)
  
      EMIS_LINES(6)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(6)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(6)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(6)%COMPO(1)%NO_CONTRIB = NAT  
      CNT%ISP          = 1
      CNT%ITP          = 1 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.5a   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO
      
C  CONTRIBUTION LINEAR IN H+  -ION       DENSITY
C  H(n=3)/H+
  
      EMIS_LINES(6)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(6)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(6)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(6)%COMPO(2)%NO_CONTRIB = NPL  
      CNT%ISP          = 1
      CNT%ITP          = 4 
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.1.8a   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2  -MOLEC.    DENSITY
C  H(n=3)/H2(g)
  
      EMIS_LINES(6)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(6)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(6)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(6)%COMPO(3)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.5a   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = ''
      CNT%RAT_H2       = '' 
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(3)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2+ -MOLEC.ION DENSITY
C  H(n=3)/H2+(g)
  
      EMIS_LINES(6)%COMPO(4)%COMPO_NAME = 
     .     'DIATOMIC NEUTRAL HYDR. MOL ION'
      EMIS_LINES(6)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(6)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(6)%COMPO(4)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '2.2.14a   '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL  '
      CNT%RAT_H2       = 'H.12' 
      CNT%RAT_REACTION = '2.0c     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(4)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H-  -NEG. ION  DENSITY
C  H(n=3)/H-
  
      EMIS_LINES(6)%COMPO(5)%COMPO_NAME = 
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(6)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(6)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(6)%COMPO(5)%NO_CONTRIB = NML  
      CNT%ISP          = 1
      CNT%ITP          = 2
      CNT%IRATIO       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H2           = 'H.12'
      CNT%REACTION     = '7.2a      '
      CNT%CR           = 'OT ' 
      CNT%FRATIO       = 'AMJUEL  '
      CNT%RAT_H2       = 'H.11' 
      CNT%RAT_REACTION = '7.0a     '
      CNT%RAT_CR       = 'OT '
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(5)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H3+ -MOL. ION  DENSITY
C  H(n=2)/H3+
  
!pbh3+      EMIS_LINES(6)%COMPO(6)%COMPO_NAME = 
!pbh3+     .     'TRIATOMIC HYDR. ION'
!pbh3+      EMIS_LINES(6)%COMPO(6)%IADV = NADVI + 6

!pbh3+      NML = COUNT(NCHARM == 2)
!pbh3+      ALLOCATE (EMIS_LINES(6)%COMPO(6)%CONTRIB(NML))
!pbh3+      EMIS_LINES(6)%COMPO(6)%NO_CONTRIB = NML  
!pbh3+      CNT%ISP          = 1
!pbh3+      CNT%ITP          = 2
!pbh3+      CNT%IRATIO       = 1
!pbh3+      CNT%FNAME        = 'AMJUEL  '
!pbh3+      CNT%H2           = 'H.12'
!pbh3+      CNT%REACTION     = '2.2.15a  '
!pbh3+      CNT%CR           = 'OT ' 
!pbh3+      CNT%FRATIO       = 'AMJUEL  '
!pbh3+      CNT%RAT_H2       = 'H.11' 
!pbh3+      CNT%RAT_REACTION = '4.0a     '
!pbh3+      CNT%RAT_CR       = 'OT '
!pbh3+      CNT%IRC          = 0
!pbh3+      CNT%IRC_RAT      = 0

!pbh3+      IML = 0
!pbh3+      DO I = 1, NMOLI
!pbh3+        IF (NCHARM(I) == 2) THEN
!pbh3+          IML = IML + 1
!pbh3+          CNT%ISP = I
!pbh3+          EMIS_LINES(6)%COMPO(6)%CONTRIB(IPL) = CNT
!pbh3+        END IF
!pbh3+      END DO

      
      end subroutine eirene_setup_default_emissivity
