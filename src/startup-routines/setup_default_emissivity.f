cdr  Oct. 18: further comments: 

cdr           This routine implicity makes some assumptions regarding the
cdr           species in input block 4a,b,c,d:
cdr           H2  (type=2)
cdr           H+  (type=4)
cdr           H   (type=1)
cdr   to be checked: contributions? multiple isotopes ?


      subroutine eirene_setup_default_emissivity

cdr originally programmed by PB 2017
cdr Routine is called from subr. INPUT.f, input block 12.
cdr april 18:  the calculation of volumetric line emissitivies
cdr            and their storing on additional tallies ADDV
cdr            has been generalized,
cdr            replacing the former 6 routines:
cdr            ba_alpha.f, ba_beta.f, ba_gamma.f, ba_delta.f,
cdr            ly_alpha.f, ly_beta.f

cdr  this present routine:
cdr  Try to reproduce the old version of these 6 routines,
cdr  by using the new structures EMIS_LINES%....
cdr
cdr  number of lines       6     (BA_AL, BA_BET, ..., LY_BET)
cdr  number of components: 6     (COUPLING TO H, H+,H2,H2+,H-,H3+)
cdr  number of contributions:  detected from input file,
cdr                            as in old ba... ly... routines
cdr                           (there sum over contributions only
cdr                            on ADDV tallies),
cdr  hard-coded here: use pop. coeffs from amjuel H.12, and
cdr                   use ratios for short living radicals (H2+, H3+, H-)
cdr                   from amjuel H.11 and H.12
cdr  hard coded: 
cdr              H2+, H- and H3+ QSS states, because of hard coded density ratios.
cdr              H2+ must be produced from both EI and IC processes, because
cdr              hard coded ratio H.12 2.0c is used here.
cdr
cdr  tbd:  make consistent notation "component vs. contribution":  DONE !
cdr  
cdr  The ADDV tallies are filled later,
cdr  in calls to emissivity.f from MCARLO, per stratum. 

cdr  Apparently we need at least one chord and nchtal=2, (even if unused)
cdr  to fill the ADDV arrays, because of a hidden link between emissivity options
cdr  and line-of-sight options.
cdr

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comsig
      use eirmod_comusr

      implicit none

      TYPE(TCONTRIB) :: CNT

      integer :: i, NUM_compo, iat, iml, ipl, nat, npl, nml
      real(dp) :: ry = 13.605

      NUM_lines    = 6
      NUM_compo    = 6
c     NUM_contrib  = inferred from input file, species specification block 4.
      MOD_ADDV = 0

! NOT USED IN DEFAULT MODEL
      CNT%IZ = 0
      CNT%IZ_RAT = 0
      CNT%ELEMENT = '  '
      CNT%RAT_ELEMENT = '  '

      ALLOCATE (EMIS_LINES(NUM_LINES))
      EMIS_LINES%LINE_NAME = REPEAT(' ',80)
      EMIS_LINES%NUM_COMPO = 0


************************************************
* BALMER ALPHA, LINE NO. 1, all 6 components
************************************************

      EMIS_LINES(1)%LINE_NAME = 'BA_ALPHA'
      EMIS_LINES(1)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(1)%EINSTEIN = 4.410E7
c  transition energy
      EMIS_LINES(1)%TRANS_EN = RY *
     .                        (1._dp/(2._DP*2._DP)-1._DP/(3._DP*3._DP))
C  identifier of Line:
      EMIS_LINES(1)%ENERGY = 1.8889_DP
      EMIS_LINES(1)%POP_ESC = 1.0_DP
      EMIS_LINES(1)%IROW_ESC = 0
      EMIS_LINES(1)%ICOL_ESC = 0
      EMIS_LINES(1)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(1)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=3)/H(n=1)

      EMIS_LINES(1)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(1)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(1)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(1)%COMPO(1)%NUM_CONTRIB = NAT
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.5a   '
      CNT%CR           = 'OT '

cdr all reaction data are the same for all contributions.
cdr only CNT%ISP (species index) may differ for different contributions.
      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  COMPONENT 2: LINEAR IN H+/D+/T+ ION DENSITY
C  ALL BULK ION (ITYP=4) CONTRIBUTIONS WITH
C                        NUCLEAR CHARGE NUMBER=1 AND CHARGE STATE NUMBER=1
C  H(n=3)/H+

      EMIS_LINES(1)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(1)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(1)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(1)%COMPO(2)%NUM_CONTRIB = NPL
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 4
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.8a   '
      CNT%CR           = 'OT '

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  COMPONENT 3: LINEAR IN "H2" MOLEC. DENSITY
C  ALL MOLECULE (ITYP=2) CONTRIBUTIONS WITH
C                        NUCLEAR CHARGE NUMBER=2
C  H(n=3)/H2(g)

      EMIS_LINES(1)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(1)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(1)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(1)%COMPO(3)%NUM_CONTRIB = NML
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 2
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.2.5a   '
      CNT%CR           = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 4: LINEAR IN "H2+" MOLEC. ION DENSITY
C  H(n=3)/H2+(g)

      EMIS_LINES(1)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(1)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(1)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(1)%COMPO(4)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.14a  '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.12'
      CNT%RAT_REACTION(1) = '2.0c     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 5: LINEAR IN H- NEG. ION DENSITY
C  H(n=3)/H-

      EMIS_LINES(1)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(1)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(1)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(1)%COMPO(5)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '7.2a     '
      CNT%CR              = 'OT '

      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '7.0a     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 6: LINEAR IN H3+ MOL. ION DENSITY
C  H(n=3)/H3+

      EMIS_LINES(1)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(1)%COMPO(6)%IADV = NADVI + 6

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(1)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(1)%COMPO(6)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 2
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.15a  '
      CNT%CR              = 'OT '

      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '4.0a     '
      CNT%RAT_CR(1)       = 'OT '

      CNT%ISP(2)          = 1
      CNT%ITP(2)          = 2
      CNT%ISP(3)          = 1
      CNT%ITP(3)          = 5
      CNT%FRATIO(2)       = 'AMJUEL  '
      CNT%RAT_H123(2)     = 'H.12'
      CNT%RAT_REACTION(2) = '2.0c     '
      CNT%RAT_CR(2)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO


************************************************
* BALMER BETA, LINE NO. 2
************************************************

      EMIS_LINES(2)%LINE_NAME = 'BA_BETA'
      EMIS_LINES(2)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(2)%EINSTEIN = 8.419E6
      EMIS_LINES(2)%TRANS_EN = RY *
     .                        (1._dp/(2._DP*2._DP)-1._DP/(4._DP*4._DP))
      EMIS_LINES(2)%ENERGY = 2.5500_DP
      EMIS_LINES(2)%POP_ESC = 1.0_DP
      EMIS_LINES(2)%IROW_ESC = 0
      EMIS_LINES(2)%ICOL_ESC = 0
      EMIS_LINES(2)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(2)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=4)/H(n=1)

      EMIS_LINES(2)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(2)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(2)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(2)%COMPO(1)%NUM_CONTRIB = NAT
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.5c   '
      CNT%CR           = 'OT '

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  COMPONENT 2: LINEAR IN H+ ION DENSITY
C  H(n=4)/H+

      EMIS_LINES(2)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(2)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(2)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(2)%COMPO(2)%NUM_CONTRIB = NPL
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 4
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.8c   '
      CNT%CR           = 'OT '

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  COMPONENT 3: LINEAR IN H2 MOLEC. DENSITY
C  H(n=4)/H2(g)

      EMIS_LINES(2)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(2)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(2)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(2)%COMPO(3)%NUM_CONTRIB = NML
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 2
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.2.5c   '
      CNT%CR           = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 4: LINEAR IN H2+ MOLEC. ION DENSITY
C  H(n=4)/H2+(g)

      EMIS_LINES(2)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(2)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(2)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(2)%COMPO(4)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.14c  '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.12'
      CNT%RAT_REACTION(1) = '2.0c     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 5: LINEAR IN H- NEG. ION DENSITY
C  H(n=4)/H-

      EMIS_LINES(2)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(2)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(2)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(2)%COMPO(5)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '7.2c     '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '7.0a     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 6: LINEAR IN H3+ MOL. ION DENSITY
C  H(n=4)/H3+

      EMIS_LINES(2)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(2)%COMPO(6)%IADV = NADVI + 6

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(2)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(2)%COMPO(6)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 2
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.15c  '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '4.0a     '
      CNT%RAT_CR(1)       = 'OT '
      CNT%ISP(2)          = 1
      CNT%ITP(2)          = 2
      CNT%ISP(3)          = 1
      CNT%ITP(3)          = 5
      CNT%FRATIO(2)       = 'AMJUEL  '
      CNT%RAT_H123(2)     = 'H.12'
      CNT%RAT_REACTION(2) = '2.0c     '
      CNT%RAT_CR(2)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
         IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO


************************************************
* BALMER GAMMA
************************************************

      EMIS_LINES(3)%LINE_NAME = 'BA_GAMMA'
      EMIS_LINES(3)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(3)%EINSTEIN = 2.530E6
      EMIS_LINES(3)%TRANS_EN = RY *
     .                        (1._dp/(2._DP*2._DP)-1._DP/(5._DP*5._DP))
      EMIS_LINES(3)%ENERGY = 2.8560_DP
      EMIS_LINES(3)%POP_ESC = 1.0_DP
      EMIS_LINES(3)%IROW_ESC = 0
      EMIS_LINES(3)%ICOL_ESC = 0
      EMIS_LINES(3)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(3)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=5)/H(n=1)

      EMIS_LINES(3)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(3)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(3)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(3)%COMPO(1)%NUM_CONTRIB = NAT
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 1
      CNT%IRATIO       = 0
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.5d   '
      CNT%CR           = 'OT '

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  COMPONENT 2: LINEAR IN H+ ION DENSITY
C  H(n=5)/H+

      EMIS_LINES(3)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(3)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(3)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(3)%COMPO(2)%NUM_CONTRIB = NPL
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 4
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.8d   '
      CNT%CR           = 'OT '

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  COMPONENT 3: LINEAR IN H2  -MOLEC.    DENSITY
C  H(n=5)/H2(g)

      EMIS_LINES(3)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(3)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(3)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(3)%COMPO(3)%NUM_CONTRIB = NML
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 2
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.2.5d   '
      CNT%CR           = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2+ MOLEC. ION DENSITY
C  H(n=5)/H2+(g)

      EMIS_LINES(3)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(3)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(3)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(3)%COMPO(4)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.14d  '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.12'
      CNT%RAT_REACTION(1) = '2.0c     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 5: LINEAR IN H- NEG. ION DENSITY
C  H(n=5)/H-

      EMIS_LINES(3)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(3)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(3)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(3)%COMPO(5)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '7.2d     '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '7.0a     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 6: LINEAR IN H3+ MOL. ION DENSITY
C  H(n=5)/H3+

      EMIS_LINES(3)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(3)%COMPO(6)%IADV = NADVI + 6

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(3)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(3)%COMPO(6)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 2
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.15d  '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '4.0a     '
      CNT%RAT_CR(1)       = 'OT '
      CNT%ISP(2)          = 1
      CNT%ITP(2)          = 2
      CNT%ISP(3)          = 1
      CNT%ITP(3)          = 5
      CNT%FRATIO(2)       = 'AMJUEL  '
      CNT%RAT_H123(2)     = 'H.12'
      CNT%RAT_REACTION(2) = '2.0c     '
      CNT%RAT_CR(2)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO



************************************************
* BALMER DELTA
************************************************

      EMIS_LINES(4)%LINE_NAME = 'BA_DELTA'
      EMIS_LINES(4)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(4)%EINSTEIN = 9.732E5
      EMIS_LINES(4)%TRANS_EN = RY *
     .                        (1._dp/(2._DP*2._DP)-1._DP/(6._DP*6._DP))
      EMIS_LINES(4)%ENERGY = 3.0222_DP
      EMIS_LINES(4)%POP_ESC = 1.0_DP
      EMIS_LINES(4)%IROW_ESC = 0
      EMIS_LINES(4)%ICOL_ESC = 0
      EMIS_LINES(4)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(4)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=6)/H(n=1)

      EMIS_LINES(4)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(4)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(4)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(4)%COMPO(1)%NUM_CONTRIB = NAT
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.5e   '
      CNT%CR           = 'OT '

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  COMPONENT 2: LINEAR IN H+ ION DENSITY
C  H(n=6)/H+

      EMIS_LINES(4)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(4)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(4)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(4)%COMPO(2)%NUM_CONTRIB = NPL
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 4
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.8e   '
      CNT%CR           = 'OT '

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  COMPONENT 3: LINEAR IN H2 MOLEC. DENSITY
C  H(n=6)/H2(g)

      EMIS_LINES(4)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(4)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(4)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(4)%COMPO(3)%NUM_CONTRIB = NML
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 2
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.2.5e   '
      CNT%CR           = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2+ MOLEC. ION DENSITY
C  H(n=6)/H2+(g)

      EMIS_LINES(4)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(4)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(4)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(4)%COMPO(4)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.14e  '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.12'
      CNT%RAT_REACTION(1) = '2.0c     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 5: LINEAR IN H- NEG. ION DENSITY
C  H(n=6)/H-

      EMIS_LINES(4)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(4)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(4)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(4)%COMPO(5)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '7.2e     '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '7.0a     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H3+ MOL. ION DENSITY
C  H(n=6)/H3+

      EMIS_LINES(4)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(4)%COMPO(6)%IADV = NADVI + 6

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(4)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(4)%COMPO(6)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 2
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.15e  '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '4.0a     '
      CNT%RAT_CR(1)       = 'OT '
      CNT%ISP(2)          = 1
      CNT%ITP(2)          = 2
      CNT%ISP(3)          = 1
      CNT%ITP(3)          = 5
      CNT%FRATIO(2)       = 'AMJUEL  '
      CNT%RAT_H123(2)     = 'H.12'
      CNT%RAT_REACTION(2) = '2.0c     '
      CNT%RAT_CR(2)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO


************************************************
* LYMAN ALPHA
************************************************

      EMIS_LINES(5)%LINE_NAME = 'LY_ALPHA'
      EMIS_LINES(5)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(5)%EINSTEIN = 4.699E8
      EMIS_LINES(5)%TRANS_EN = RY *
     .                        (1._dp/(1._DP*1._DP)-1._DP/(2._DP*2._DP))
      EMIS_LINES(5)%ENERGY = 10.2375_DP
      EMIS_LINES(5)%POP_ESC = 1.0_DP
      EMIS_LINES(5)%IROW_ESC = 0
      EMIS_LINES(5)%ICOL_ESC = 0
      EMIS_LINES(5)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(5)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=2)/H(n=1)

      EMIS_LINES(5)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(5)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(5)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(5)%COMPO(1)%NUM_CONTRIB = NAT
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.5b   '
      CNT%CR           = 'OT '

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H+ ION DENSITY
C  H(n=2)/H+

      EMIS_LINES(5)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(5)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(5)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(5)%COMPO(2)%NUM_CONTRIB = NPL
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 4
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.8b   '
      CNT%CR           = 'OT '

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2 MOLEC. DENSITY
C  H(n=2)/H2(g)

      EMIS_LINES(5)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(5)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(5)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(5)%COMPO(3)%NUM_CONTRIB = NML
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 2
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.2.5b   '
      CNT%CR           = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2+ MOLEC. ION DENSITY
C  H(n=2)/H2+(g)

      EMIS_LINES(5)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(5)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(5)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(5)%COMPO(4)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.14b  '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.12'
      CNT%RAT_REACTION(1) = '2.0c     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H- NEG. ION DENSITY
C  H(n=2)/H-

      EMIS_LINES(5)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(5)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(5)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(5)%COMPO(5)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '7.2b     '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '7.0a     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H3+ MOL. ION DENSITY
C  H(n=2)/H3+

      EMIS_LINES(5)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(5)%COMPO(6)%IADV = NADVI + 6

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(5)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(5)%COMPO(6)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 2
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.15b  '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '4.0a     '
      CNT%RAT_CR(1)       = 'OT '
      CNT%ISP(2)          = 1
      CNT%ITP(2)          = 2
      CNT%ISP(3)          = 1
      CNT%ITP(3)          = 5
      CNT%FRATIO(2)       = 'AMJUEL  '
      CNT%RAT_H123(2)     = 'H.12'
      CNT%RAT_REACTION(2) = '2.0c     '
      CNT%RAT_CR(2)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO


************************************************
* LYMAN BETA
************************************************

      EMIS_LINES(6)%LINE_NAME = 'LY_BETA'
      EMIS_LINES(6)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(6)%EINSTEIN = 5.575E7
      EMIS_LINES(6)%TRANS_EN = RY *
     .                        (1._dp/(1._DP*1._DP)-1._DP/(3._DP*3._DP))
      EMIS_LINES(6)%ENERGY = 12.089_DP
      EMIS_LINES(6)%POP_ESC  = 1.0_DP
      EMIS_LINES(6)%IROW_ESC = 0
      EMIS_LINES(6)%ICOL_ESC = 0
      EMIS_LINES(6)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(6)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=3)/H(n=1)

      EMIS_LINES(6)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(6)%COMPO(1)%IADV = NADVI + 1

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(6)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(6)%COMPO(1)%NUM_CONTRIB = NAT
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 1
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.5a   '
      CNT%CR           = 'OT '

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H+ ION DENSITY
C  H(n=3)/H+

      EMIS_LINES(6)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(6)%COMPO(2)%IADV = NADVI + 2

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(6)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(6)%COMPO(2)%NUM_CONTRIB = NPL
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 4
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.1.8a   '
      CNT%CR           = 'OT '

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2 MOLEC. DENSITY
C  H(n=3)/H2(g)

      EMIS_LINES(6)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(6)%COMPO(3)%IADV = NADVI + 3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(6)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(6)%COMPO(3)%NUM_CONTRIB = NML
      CNT%ISP          = -1
      CNT%ITP          = -1
      CNT%FRATIO       = ''
      CNT%RAT_H123     = ''
      CNT%RAT_REACTION = ''
      CNT%RAT_CR       = ''
      CNT%IRC          = 0
      CNT%IRC_RAT      = 0

      CNT%IRATIO       = 0
      CNT%ISP(1)       = 1
      CNT%ITP(1)       = 2
      CNT%FNAME        = 'AMJUEL  '
      CNT%H123         = 'H.12'
      CNT%REACTION     = '2.2.5a   '
      CNT%CR           = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H2+ MOLEC. ION DENSITY
C  H(n=3)/H2+(g)

      EMIS_LINES(6)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(6)%COMPO(4)%IADV = NADVI + 4

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(6)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(6)%COMPO(4)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.14a  '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.12'
      CNT%RAT_REACTION(1) = '2.0c     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H- NEG. ION DENSITY
C  H(n=3)/H-

      EMIS_LINES(6)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(6)%COMPO(5)%IADV = NADVI + 5

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(6)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(6)%COMPO(5)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 1
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '7.2a     '
      CNT%CR              = 'OT '
      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '7.0a     '
      CNT%RAT_CR(1)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  CONTRIBUTION LINEAR IN H3+ MOL. ION DENSITY
C  H(n=2)/H3+

      EMIS_LINES(6)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(6)%COMPO(6)%IADV = NADVI + 6

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(6)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(6)%COMPO(6)%NUM_CONTRIB = NML
      CNT%ISP             = -1
      CNT%ITP             = -1
      CNT%FRATIO          = ''
      CNT%RAT_H123        = ''
      CNT%RAT_REACTION    = ''
      CNT%RAT_CR          = ''
      CNT%IRC             = 0
      CNT%IRC_RAT         = 0

      CNT%IRATIO          = 2
      CNT%ISP(1)          = 1
      CNT%ITP(1)          = 2
      CNT%FNAME           = 'AMJUEL  '
      CNT%H123            = 'H.12'
      CNT%REACTION        = '2.2.15a  '
      CNT%CR              = 'OT '

      CNT%FRATIO(1)       = 'AMJUEL  '
      CNT%RAT_H123(1)     = 'H.11'
      CNT%RAT_REACTION(1) = '4.0a     '
      CNT%RAT_CR(1)       = 'OT '

      CNT%ISP(2)          = 1
      CNT%ITP(2)          = 2
      CNT%ISP(3)          = 1
      CNT%ITP(3)          = 5
      CNT%FRATIO(2)       = 'AMJUEL  '
      CNT%RAT_H123(2)     = 'H.12'
      CNT%RAT_REACTION(2) = '2.0c     '
      CNT%RAT_CR(2)       = 'OT '

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO


      end subroutine eirene_setup_default_emissivity
