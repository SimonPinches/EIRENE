cdr Nov. 16: MODULE FOR ALL ATOMIC/MOLECULAR/PHOTONIC DATA STRUCTURES.
cdr
cdr  MXCOLLS --> MSTOR0

      MODULE EIRMOD_COMXS
 
!  jan-05: natprc_2,..... introduced
!  07.12.05: bugfix: IFTFLG is now available for default reactions too
!                    via dimensioning IFTFLG(-10:NREAC)
!  30.08.06: data structure for reaction data redefined
!  12.10.06: modcol revised
!  19.12.06: test functions added which allow to test if a rate-coefficient
!            is defined via ADAS database
!  02.03.07: remove ESCD2* arrays
!  02.03.07: fourth secondary group specifier introduced
!
!  02.03.07: IMESS added in photon line reaction data in order to allow for
!            a complete printout of input data in case of HYDKIN default
!            database option
!
cdr sometime between 2004 and 2007 the atomic data structure was revised.
cdr 
cdr  now it is on REACDAT.  Commenting, cleanup started: jan 2016.
!
!  24.03.15: number of default reactions increased from 10 to 11, REACDAT(-11)...
cdr23.04.15: only text, comments.... continued: Nov. 15, still not complete
cdr  JAN  16:  additional species index for eplds-->eplei, eplpi
!pb  APR  16:  ipplds -> ipplei, pplds -> pplei
!pb  APR  16:  ipatds -> ipatei, patds -> patei, eatds -> eatei
!pb  APR  16:  ipmlds -> ipmlei, pmlds -> pmlei, emlds -> emlei
!pb  APR  16:  ipiods -> ipioei, piods -> pioei, eiods -> eioei
!pb  APR  16:  pelds  -> pelei,  eelds -> eelei
!pb  MAY  16:  tabds1 -> tabei1
!pb  MAY  16:  nrds   -> nrei
!pb  JUL  16:  ehvds1 -> ehvei1
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
 
      IMPLICIT NONE
 
      PRIVATE
 
      PUBLIC :: EIRENE_ALLOC_COMXS, EIRENE_DEALLOC_COMXS,
     .          EIRENE_INIT_CMDTA,  EIRENE_WRITE_CMDTA,   
     .          EIRENE_READ_CMDTA,
     .          EIRENE_WRITE_CMAMF, EIRENE_READ_CMAMF, 
     .          EIRENE_CMDTA_XDR, EIRENE_CMAMF_XDR,
cdr
     .          LINE_DATA, EIRENE_GET_REACTION, POLY_DATA,
     .          REACTION_DATA, EIRENE_SET_REACTION_DATA, ADAS_DATA,
     .          FIT_FORMS, EIRENE_IS_RTC_ADAS, EIRENE_IS_RTCEW_ADAS, 
     .          EIRENE_IS_RTCMW_ADAS,
     .          REACTION_INPUT_LINE, HYDKIN_DATA
 
      TYPE LINE_DATA
        REAL(DP) :: E0, E1, AIK, G1, G2, C2, C3, C4, C6, B12, B21
        REAL(DP) :: C6A(12)
        INTEGER :: IGND, IRCART, IPROFILETYPE, IFREMD, NRJPRT, IMESS
        INTEGER :: IPLSC6(12)
        CHARACTER(2) :: KENN(12)
        CHARACTER(50) :: REACNAME
      END TYPE LINE_DATA
 
      TYPE ADAS_DATA
        INTEGER :: NDENS, NTEMP
        REAL(DP), POINTER :: DENS(:), TEMP(:), FIT(:,:)
        REAL(DP), POINTER :: DDE(:), DTE(:)
      END TYPE ADAS_DATA
 
      TYPE POLY_DATA
        INTEGER :: IFEXMN, IFEXMX
        REAL(DP), POINTER :: DBLPOL(:,:)
        REAL(DP) :: RCMN, RCMX, FPARM(6)
      END TYPE POLY_DATA
 
      TYPE HYDKIN_DATA
        INTEGER :: NTEMPS
        REAL(DP), POINTER :: TEMPS(:), RATES(:), RATIO(:)
        CHARACTER(50) :: REAC_STRING, REACNAME
        CHARACTER(100) :: RPRT
      END TYPE HYDKIN_DATA
 
      TYPE FIT_FORMS
        INTEGER :: IFIT
        TYPE(POLY_DATA),   POINTER :: POLY
        TYPE(ADAS_DATA),   POINTER :: ADAS
        TYPE(LINE_DATA),   POINTER :: LINE
        TYPE(HYDKIN_DATA), POINTER :: HYD
      END TYPE FIT_FORMS
 
      TYPE REACTION_DATA
        TYPE(FIT_FORMS), POINTER :: POT, CRS, RTC, RTCMW, RTCEW,
     T                              OTH, PHR
        LOGICAL :: LPOT, LCRS, LRTC, LRTCMW, LRTCEW, LOTH, LPHR
        INTEGER :: NOSEC
      END TYPE REACTION_DATA
 
      TYPE REACTION_INPUT_LINE
        INTEGER :: NO, MT, MP, IZ, JFEXMN, JFEXMX, NCONST
        REAL(DP) :: RMN, RMX, DPP, FP(6), CONST(9)
        CHARACTER(8) :: FILE
        CHARACTER(50) :: REAC_STRING
        CHARACTER(4) :: H_SELECT
        CHARACTER(3) :: REACTYP
        CHARACTER(2) :: ELEMENT
      END TYPE REACTION_INPUT_LINE
 
      TYPE(LINE_DATA), POINTER, PUBLIC, SAVE :: REACTION
      INTEGER, PUBLIC, SAVE :: IDREAC, IRLINES
 
      TYPE(REACTION_DATA), ALLOCATABLE, PUBLIC, SAVE :: REACDAT(:)
 
      TYPE(REACTION_INPUT_LINE), ALLOCATABLE, PUBLIC, SAVE ::
     .                           REACLINES(:)
 
 
      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        XSTOR(:,:), XSTORV(:)

cdr  local (on the flight) atomic-moleculer reaction data 
      REAL(DP), PUBLIC, POINTER, SAVE ::
c  reaction rates, by reaction
     R SIGVCX(:),   SIGVPI(:),   SIGVEI(:),   SIGVEL(:),
c  energy exchange rates, by reaction
     R ESIGCX(:,:), ESIGPI(:,:), ESIGEI(:,:), ESIGEL(:,:),
c  momentum exchange rates, by reaction
     R VSIGCX(:),   VSIGPI(:),   VSIGEL(:),
c  totals
     R SIGCXT,      SIGPIT,      SIGEIT,      SIGELT,      SIGTOT,
     R SIGBGK,
c  invers mean free path
     R ZMFPI
 
      REAL(DP), PUBLIC, SAVE :: ZMFPTHI, TDGTEMX
 
csw added OTHER (OT) reactions
      REAL(DP), PUBLIC, POINTER, SAVE :: SIGVOT(:),   ESIGOT(:,:),
     R SIGOTT
 
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R TABEI1(:,:),   TABRC1(:,:),
     R TABPI3(:,:,:), TABCX3(:,:,:), TABEL3(:,:,:),
     R FDLMPI(:),     FDLMCX(:),     FDLMEL(:),
     R ADDPI(:,:),    ADDCX(:,:),    ADDEL(:,:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R FACRRC(:,:), FACRPI(:,:), FACREL(:,:), FACREI(:,:), FACRCX(:,:) 

c  secondaries, species distribution, for EI and PI processes 
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R PELEI(:),  PATEI(:,:), PMLEI(:,:), PIOEI(:,:), PPLEI(:,:),
     R PELPI(:),  PATPI(:,:), PMLPI(:,:), PIOPI(:,:), PPLPI(:,:),
c  ...and cummulated distributions thereof, for species sampling
     R P2ND(:,:), P2NP(:,:),  P2NDS(:),   P2NPI(:)
 
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R EELEI1(:,:),   EELRC1(:,:),   EELPI1(:,:), !  missing: eelot1,  el and cx processes have no secondary electrons
     R EHVEI1(:,:),   EHVPI3(:,:,:),
     R EPLPI3(:,:,:), EPLCX3(:,:,:), EPLEL3(:,:,:), EPLOT3(:,:,:)
 
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R EATEI(:,:,:), EMLEI(:,:,:), EIOEI(:,:,:), EPLEI(:,:,:),
     R EATPI(:,:,:), EMLPI(:,:,:), EIOPI(:,:,:), EPLPI(:,:,:)
 
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I MODCOL(:,:,:),
     I IESTCX(:,:), IESTEL(:,:), IESTPI(:,:), IESTEI(:,:),
     I NAEII(:),    NMEII(:),    NIEII(:),
     I NACXI(:),    NMCXI(:),    NICXI(:),
     I NAELI(:),    NMELI(:),    NIELI(:),
     I NAPII(:),    NMPII(:),    NIPII(:),
     I NAEIIM(:),   NMEIIM(:),   NIEIIM(:),
     I NACXIM(:),   NMCXIM(:),   NICXIM(:),
     I NAELIM(:),   NMELIM(:),   NIELIM(:),
     I NAPIIM(:),   NMPIIM(:),   NIPIIM(:),
     I NPRCI(:),    NPRCIM(:),
     I NPBGKA(:),   NPBGKM(:),   NPBGKI(:), NPBGKP(:,:)
 
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NATPRC(:),  NMLPRC(:), NIOPRC(:), NPLPRC(:), NPHPRC(:),
     I NATPRC_2(:),  NMLPRC_2(:), NIOPRC_2(:), NPLPRC_2(:), NPHPRC_2(:),
     I N1STX(:,:), N2NDX(:,:)
 
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NSEACX(:,:,:), NSEMCX(:,:,:), NSEICX(:,:,:),
     I NSEAEL(:,:,:), NSEMEL(:,:,:), NSEIEL(:,:,:),
     I NSEPRC(:,:)
 
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NREACX(:),NREAPI(:),NREAEL(:),
     I NREAEI(:),JEREAEI(:),NREARC(:),JEREARC(:),
     I NELREI(:),JELREI(:),NREAHV(:),NELREL(:),
     I NELRRC(:),JELRRC(:),NELRPI(:),JELRPI(:),NELRCX(:),
     I NELROT(:),NREAOT(:),NREACT(:),NRHVPI(:),
     I IPATEI(:,:),IPMLEI(:,:),
     I IPIOEI(:,:),IPPLEI(:,:),
     I IPATPI(:,:),IPMLPI(:,:),
     I IPIOPI(:,:),IPPLPI(:,:),
     I LGACX(:,:,:),LGMCX(:,:,:),
     I LGICX(:,:,:),
     I LGAEI(:,:),    LGMEI(:,:),
     I LGIEI(:,:),
     I LGAEL(:,:,:),LGMEL(:,:,:),
     I LGIEL(:,:,:),
     I LGPRC(:,:),
     I LGAPI(:,:,:),LGMPI(:,:,:),
     I LGIPI(:,:,:)
 
      INTEGER, PUBLIC, SAVE ::
     I NRPII, NREII, NRCXI, NRELI, NRRCI, NRBGI
 
      INTEGER, PUBLIC, SAVE ::
     I NSTOR1, NSTOR,  NSTORV, NTAB, NDAT, NMDTA, MMDTA, NAMF, MAMF,
     I MSTOR0, MSTOR1, MSTOR2
 
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R DELPOT(:),   FACREA(:,:),
     R FREACA(:,:), FREACM(:,:), FREACI(:,:), FREACP(:,:), FREACPH(:,:),
     R FLDLMA(:,:), FLDLMM(:,:), FLDLMI(:,:), FLDLMP(:,:), FLDLMPH(:,:),
     R EELECA(:,:), EELECM(:,:), EELECI(:,:), EELECP(:,:), EELECPH(:,:),
     R EBULKA(:,:), EBULKM(:,:), EBULKI(:,:), EBULKP(:,:), EBULKPH(:,:),
     R ESCD1A(:,:), ESCD1M(:,:), ESCD1I(:,:), ESCD1P(:,:), ESCD1PH(:,:)
 
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I ISWR(:),     MODCLF(:),   MASSP(:),    MASST(:),
     I IFTFLG(:,:),
     I NRCP(:),     NRCA(:),     NRCM(:),     NRCI(:),     NRCPH(:),
     I IREACA(:,:), IREACM(:,:), IREACI(:,:), IREACP(:,:), IREACPH(:,:),
     I IBULKA(:,:), IBULKM(:,:), IBULKI(:,:), IBULKP(:,:), IBULKPH(:,:),
     I ISCD1A(:,:), ISCD1M(:,:), ISCD1I(:,:), ISCD1P(:,:), ISCD1PH(:,:),
     I ISCD2A(:,:), ISCD2M(:,:), ISCD2I(:,:), ISCD2P(:,:), ISCD2PH(:,:),
     I ISCD3A(:,:), ISCD3M(:,:), ISCD3I(:,:), ISCD3P(:,:), ISCD3PH(:,:),
     I ISCD4A(:,:), ISCD4M(:,:), ISCD4I(:,:), ISCD4P(:,:), ISCD4PH(:,:),
     I ISCDEA(:,:), ISCDEM(:,:), ISCDEI(:,:), ISCDEP(:,:), ISCDEPH(:,:),
     I IESTMA(:,:), IESTMM(:,:), IESTMI(:,:), IESTMPH(:,:),
     I IBGKA (:,:), IBGKM (:,:), IBGKI (:,:), IBGKPH (:,:)
 
      INTEGER, PUBLIC, SAVE ::
     I NREACI

      INTEGER, PUBLIC, SAVE :: MAXSPC(0:4)
 
      CHARACTER(50), PUBLIC, ALLOCATABLE, SAVE :: REAC_NAME(:)
 
      CONTAINS
 
 
      SUBROUTINE EIRENE_ALLOC_COMXS (ICAL)
CDR
C  AUTOMATTED ALLOCATION OF STORAGE FOR A&M DATA STRUCTURES AND ARRAYS. 
C  CALLED FROM:  ALLOCATE_MODULES.F
C  ICAL=1: ... 
C  ICAL=2: ...
 
      INTEGER, INTENT(IN) :: ICAL
      INTEGER, PARAMETER :: IL = SELECTED_INT_KIND(15)
      INTEGER(IL) :: MEM

      IF (ICAL == 1) THEN
 
        IF (ALLOCATED(XSTORV)) RETURN
 
        NSTORV = 8
C
        NAMF=9*11*11+11+2*(2*11+6*11)+
     P       NREAC*(9*11+18+ 6*NPHOT+ 6*NATM+ 6*NMOL+ 6*NION+ 6*NPLS)
        MAMF=NREAC*(      8+ 8*NPHOT+ 8*NATM+ 8*NMOL+ 8*NION+ 6*NPLS)+
     P       1*NATM+ 1*NMOL+ 1*NION+ 1*NPLS+ 1 +
     P       (12+NREAC)*10
 
        ALLOCATE (XSTORV(NSTORV))
 
        SIGCXT  => XSTORV(1)
        SIGPIT  => XSTORV(2)
        SIGEIT  => XSTORV(3)
        SIGELT  => XSTORV(4)
        SIGOTT  => XSTORV(5)
        SIGTOT  => XSTORV(6)
        SIGBGK  => XSTORV(7)
        ZMFPI   => XSTORV(8)
 
        ALLOCATE (NAEII(NATM))
        ALLOCATE (NMEII(NMOL))
        ALLOCATE (NIEII(NION))
        ALLOCATE (NACXI(NATM))
        ALLOCATE (NMCXI(NMOL))
        ALLOCATE (NICXI(NION))
        ALLOCATE (NAELI(NATM))
        ALLOCATE (NMELI(NMOL))
        ALLOCATE (NIELI(NION))
        ALLOCATE (NAPII(NATM))
        ALLOCATE (NMPII(NMOL))
        ALLOCATE (NIPII(NION))
        ALLOCATE (NPRCI(NPLS))
        ALLOCATE (NAEIIM(NATM))
        ALLOCATE (NMEIIM(NMOL))
        ALLOCATE (NIEIIM(NION))
        ALLOCATE (NACXIM(NATM))
        ALLOCATE (NMCXIM(NMOL))
        ALLOCATE (NICXIM(NION))
        ALLOCATE (NAELIM(NATM))
        ALLOCATE (NMELIM(NMOL))
        ALLOCATE (NIELIM(NION))
        ALLOCATE (NAPIIM(NATM))
        ALLOCATE (NMPIIM(NMOL))
        ALLOCATE (NIPIIM(NION))
        ALLOCATE (NPRCIM(NPLS))
        ALLOCATE (NPBGKA(NATM))
        ALLOCATE (NPBGKM(NMOL))
        ALLOCATE (NPBGKI(NION))
        ALLOCATE (NPBGKP(NPLS,2))
 
        ALLOCATE (NSEACX(NATM,NPLS,5))
        ALLOCATE (NSEMCX(NMOL,NPLS,5))
        ALLOCATE (NSEICX(NION,NPLS,5))
        ALLOCATE (NSEAEL(NATM,NPLS,5))
        ALLOCATE (NSEMEL(NMOL,NPLS,5))
        ALLOCATE (NSEIEL(NION,NPLS,5))
 
        ALLOCATE (DELPOT(NREAC))
        ALLOCATE (FACREA(-11:NREAC,2))
        ALLOCATE (FREACA(NATM,NREAC))
        ALLOCATE (FREACM(NMOL,NREAC))
        ALLOCATE (FREACI(NION,NREAC))
        ALLOCATE (FREACP(NPLS,NREAC))
        ALLOCATE (FREACPH(NPHOT,NREAC))
        ALLOCATE (FLDLMA(NATM,NREAC))
        ALLOCATE (FLDLMM(NMOL,NREAC))
        ALLOCATE (FLDLMI(NION,NREAC))
        ALLOCATE (FLDLMP(NPLS,NREAC))
        ALLOCATE (FLDLMPH(NPHOT,NREAC))
        ALLOCATE (EELECA(NATM,NREAC))
        ALLOCATE (EELECM(NMOL,NREAC))
        ALLOCATE (EELECI(NION,NREAC))
        ALLOCATE (EELECP(NPLS,NREAC))
        ALLOCATE (EELECPH(NPHOT,NREAC))
        ALLOCATE (EBULKA(NATM,NREAC))
        ALLOCATE (EBULKM(NMOL,NREAC))
        ALLOCATE (EBULKI(NION,NREAC))
        ALLOCATE (EBULKP(NPLS,NREAC))
        ALLOCATE (EBULKPH(NPHOT,NREAC))
        ALLOCATE (ESCD1A(NATM,NREAC))
        ALLOCATE (ESCD1M(NMOL,NREAC))
        ALLOCATE (ESCD1I(NION,NREAC))
        ALLOCATE (ESCD1P(NPLS,NREAC))
        ALLOCATE (ESCD1PH(NPHOT,NREAC))
 
        ALLOCATE (ISWR(NREAC))
        ALLOCATE (MODCLF(NREAC))
        ALLOCATE (MASSP(NREAC))
        ALLOCATE (MASST(NREAC))
        ALLOCATE (IFTFLG(-11:NREAC,0:5))
        ALLOCATE (NRCP(NPLS))
        ALLOCATE (NRCA(NATM))
        ALLOCATE (NRCM(NMOL))
        ALLOCATE (NRCI(NION))
        ALLOCATE (NRCPH(NPHOT))

        ALLOCATE (IREACA(NATM,NREAC))
        ALLOCATE (IREACM(NMOL,NREAC))
        ALLOCATE (IREACI(NION,NREAC))
        ALLOCATE (IREACP(NPLS,NREAC))
        ALLOCATE (IREACPH(NPHOT,NREAC))

        ALLOCATE (IBULKA(NATM,NREAC))
        ALLOCATE (IBULKM(NMOL,NREAC))
        ALLOCATE (IBULKI(NION,NREAC))
        ALLOCATE (IBULKP(NPLS,NREAC))
        ALLOCATE (IBULKPH(NPHOT,NREAC))

        ALLOCATE (ISCD1A(NATM,NREAC))
        ALLOCATE (ISCD1M(NMOL,NREAC))
        ALLOCATE (ISCD1I(NION,NREAC))
        ALLOCATE (ISCD1P(NPLS,NREAC))
        ALLOCATE (ISCD1PH(NPHOT,NREAC))

        ALLOCATE (ISCD2A(NATM,NREAC))
        ALLOCATE (ISCD2M(NMOL,NREAC))
        ALLOCATE (ISCD2I(NION,NREAC))
        ALLOCATE (ISCD2P(NPLS,NREAC))
        ALLOCATE (ISCD2PH(NPHOT,NREAC))

        ALLOCATE (ISCD3A(NATM,NREAC))
        ALLOCATE (ISCD3M(NMOL,NREAC))
        ALLOCATE (ISCD3I(NION,NREAC))
        ALLOCATE (ISCD3P(NPLS,NREAC))
        ALLOCATE (ISCD3PH(NPHOT,NREAC))

        ALLOCATE (ISCD4A(NATM,NREAC))
        ALLOCATE (ISCD4M(NMOL,NREAC))
        ALLOCATE (ISCD4I(NION,NREAC))
        ALLOCATE (ISCD4P(NPLS,NREAC))
        ALLOCATE (ISCD4PH(NPHOT,NREAC))

        ALLOCATE (ISCDEA(NATM,NREAC))
        ALLOCATE (ISCDEM(NMOL,NREAC))
        ALLOCATE (ISCDEI(NION,NREAC))
        ALLOCATE (ISCDEP(NPLS,NREAC))
        ALLOCATE (ISCDEPH(NPHOT,NREAC))

        ALLOCATE (IESTMA(NATM,NREAC))
        ALLOCATE (IESTMM(NMOL,NREAC))
        ALLOCATE (IESTMI(NION,NREAC))
        ALLOCATE (IESTMPH(NPHOT,NREAC))

        ALLOCATE (IBGKA (NATM,NREAC))
        ALLOCATE (IBGKM (NMOL,NREAC))
        ALLOCATE (IBGKI (NION,NREAC))
        ALLOCATE (IBGKPH(NPHOT,NREAC))
 
        ALLOCATE (REAC_NAME(NREAC))
cdr  -11 ... -1   : internal default atomic-molecular data
cdr    1 ... NREAC: atomic/molecular data read from external data files, input block 4
        ALLOCATE (REACDAT(-11:NREAC))
        ALLOCATE (REACLINES(NREAC_LINES))

        MEM = (NSTORV+NAMF)*8_IL + (MAMF+
     .                      9_IL*(NATM+NMOL+NION)+4_IL*NPLS+
     .                      10_IL*NPLS*(NATM+NMOL+NION))*4_IL +
     .                      NREAC*LEN(REAC_NAME(1))
 
        WRITE (55+IFOFF,'(A,T25,I15)')
     .        ' COMXS(1) ', MEM
 
 
      ELSE IF (ICAL == 2) THEN
 
        IF (ALLOCATED(XSTOR)) RETURN
C  DIMENSION OF FULL REACTION SPECIFIC ARRAYS: CFLAG, MODCOL,.... 
        MSTOR0 = MAX(NRPI, NREI, NRCX, NREL, NREC, NROT)
C  FIRST DIMENSION OF XSTOR ARRAY
        MSTOR1 = MAX(NRCX, NRPI, NREI, NREL, NROT)
C  SECOND DIMENSION OF XSTOR ARRAY
        MSTOR2 = 24
 
        NSTOR1 = NREL+NRCX+NRPI+NREI
        NSTOR  = NSTOR1+
     .           2*(NREL+NRCX+NRPI)+5*NREI+
     .           NREL+NRCX+NRPI
C
        NTAB=NSTORDR*(NREI+NREC)+
     P       NSTORDR*NSTORDT*(NRCX+NREL+NRPI)+
     P       (NPLS+1)*(NRPI+NRCX+NREL)+
     P       2*(NREC+NRPI+NREL+NREI+NRCX)
C
        NDAT=NSTORDR*(2*NREI+NREC+NRPI+
     P       NSTORDT*(NRCX+NREL+2*NRPI+NROT))+
     P      (NREI+NRPI)*
     P      (NATMP+NMOLP+NIONP+NPLSP+1)+
     P      (NRPI+NREI)*(NSPZP+1)+
     P       NRPI*2*(NATMP+NMOLP+NIONP+1)+
     P       NREI*2*(NATMP+NMOLP+NIONP+1)
C
        NMDTA=NTAB+NDAT
C
        MMDTA=7*5*MSTOR0+3*(NRCX+NREL+NRPI+NREI)+6+
     P        5*NREC+
     P        6*NRCX+
     P        10*NREC+
     P        2*NRCX+4*NRPI+2*NREL+5*NREI+4*NREC+NREAC+2*NROT+
     P        (NREI+NRPI)*
     P        (NATMP+NMOLP+NIONP+NPLSP)+
C  LG... ARRAYS
     P          (NATMP+NMOLP+NIONP      )*(NREI+1)+
     P        2*(NATMP+NMOLP+NIONP      )*(NRCX+1)+
     P        2*(NATMP+NMOLP+NIONP      )*(NREL+1)+
     P          (                  NPLSP)*(NREC+1)+
     P        2*(NATMP+NMOLP+NIONP      )*(NRPI+1)
 
 

 
        ALLOCATE (XSTOR(MSTOR1,MSTOR2))
 
        SIGVCX => XSTOR(:,1)
        SIGVPI => XSTOR(:,2)
        SIGVEI => XSTOR(:,3)
        SIGVEL => XSTOR(:,4)
        SIGVOT => XSTOR(:,22)
 
        ESIGCX => XSTOR(:,5:6)
        ESIGPI => XSTOR(:,7:11)
        ESIGEI => XSTOR(:,12:16)
        ESIGEL => XSTOR(:,17:18)
        ESIGOT => XSTOR(:,23:24)
 
        VSIGCX => XSTOR(:,19)
        VSIGPI => XSTOR(:,20)
        VSIGEL => XSTOR(:,21)
cdr     vsigei  : fehlt noch
cdr     vsigot  : fehlt noch
 
 
        ALLOCATE (TABEI1(NREI,NSTORDR))
        ALLOCATE (TABRC1(NREC,NSTORDR))
        ALLOCATE (TABPI3(NRPI,NSTORDR,NSTORDT))
        ALLOCATE (TABCX3(NRCX,NSTORDR,NSTORDT))
        ALLOCATE (TABEL3(NREL,NSTORDR,NSTORDT))
        ALLOCATE (FDLMPI(NRPI))
        ALLOCATE (FDLMCX(NRCX))
        ALLOCATE (FDLMEL(NREL))

c  factors for scaling reaction rates to other target (ipls) masses
c  (only for heavy particle impact reaction) 
        ALLOCATE (ADDPI(NRPI,NPLS))
        ALLOCATE (ADDCX(NRCX,NPLS))
        ALLOCATE (ADDEL(NREL,NPLS))

c  factors for scaling reaction processes
        ALLOCATE (FACRRC(NREC,2)) 
        ALLOCATE (FACRPI(NRPI,2)) 
        ALLOCATE (FACREL(NREL,2)) 
        ALLOCATE (FACREI(NREI,2))
        ALLOCATE (FACRCX(NRCX,2)) 
 
c  secondaries, EI processes 
        ALLOCATE (PELEI(NREI))
        ALLOCATE (PATEI(NREI,0:NATM))
        ALLOCATE (PMLEI(NREI,0:NMOL))
        ALLOCATE (PIOEI(NREI,0:NION))
        ALLOCATE (PPLEI(NREI,0:NPLS))
        ALLOCATE (P2ND(NREI,0:NSPZ))
        ALLOCATE (P2NDS(NREI))
c  secondaries, PI processes
        ALLOCATE (PELPI(NRPI))
        ALLOCATE (PATPI(NRPI,0:NATM))
        ALLOCATE (PMLPI(NRPI,0:NMOL))
        ALLOCATE (PIOPI(NRPI,0:NION))
        ALLOCATE (PPLPI(NRPI,0:NPLS))       
        ALLOCATE (P2NP(NRPI,0:NSPZ))        
        ALLOCATE (P2NPI(NRPI))
 
        ALLOCATE (EELEI1(NREI,NSTORDR))
        ALLOCATE (EHVEI1(NREI,NSTORDR))


        ALLOCATE (EELRC1(NREC,NSTORDR))

        ALLOCATE (EELPI1(NRPI,NSTORDR))
        ALLOCATE (EHVPI3(NRPI,NSTORDR,NSTORDT))
        ALLOCATE (EPLPI3(NRPI,NSTORDR,NSTORDT))

        ALLOCATE (EPLCX3(NRCX,NSTORDR,NSTORDT))
        ALLOCATE (EPLEL3(NREL,NSTORDR,NSTORDT))
 

        ALLOCATE (EPLOT3(NROT,NSTORDR,NSTORDT))
 
        ALLOCATE (EATPI(NRPI,0:NATM,2))
        ALLOCATE (EMLPI(NRPI,0:NMOL,2))
        ALLOCATE (EIOPI(NRPI,0:NION,2))
        ALLOCATE (EPLPI(NRPI,0:NPLS,2))

        ALLOCATE (EATEI(NREI,0:NATM,2))
        ALLOCATE (EMLEI(NREI,0:NMOL,2))
        ALLOCATE (EIOEI(NREI,0:NION,2))
        ALLOCATE (EPLEI(NREI,0:NPLS,2))
 
        ALLOCATE (MODCOL(7,0:4,MSTOR0))

c   flags for collision or tracklength estimators, 
c   for particle (1), momentum (2) and energy (3) source rates, resp.
        ALLOCATE (IESTCX(NRCX,3))
        ALLOCATE (IESTEL(NREL,3))
        ALLOCATE (IESTPI(NRPI,3))
        ALLOCATE (IESTEI(NREI,3))
 
        ALLOCATE (NATPRC(NREC))
        ALLOCATE (NMLPRC(NREC))
        ALLOCATE (NIOPRC(NREC))
        ALLOCATE (NPLPRC(NREC))
        ALLOCATE (NPHPRC(NREC))
        ALLOCATE (NATPRC_2(NREC))
        ALLOCATE (NMLPRC_2(NREC))
        ALLOCATE (NIOPRC_2(NREC))
        ALLOCATE (NPLPRC_2(NREC))
        ALLOCATE (NPHPRC_2(NREC))

        ALLOCATE (N1STX(NRCX,3))
        ALLOCATE (N2NDX(NRCX,3))
 
        ALLOCATE (NSEPRC(NREC,5))
 
        ALLOCATE (NREACX(NRCX))
        ALLOCATE (NREAPI(NRPI))
        ALLOCATE (NREAEL(NREL))
        ALLOCATE (NREAEI(NREI))
        ALLOCATE (JEREAEI(NREI))
        ALLOCATE (NREARC(NREC))
        ALLOCATE (JEREARC(NREC))
        ALLOCATE (NELREI(NREI))
        ALLOCATE (JELREI(NREI))
        ALLOCATE (NREAHV(NREI))
        ALLOCATE (NELREL(NREL))
        ALLOCATE (NELRRC(NREC))
        ALLOCATE (JELRRC(NREC))
        ALLOCATE (NELRPI(NRPI))
        ALLOCATE (JELRPI(NRPI))
        ALLOCATE (NELRCX(NRCX))
        ALLOCATE (NELROT(NROT))
        ALLOCATE (NREAOT(NROT))
        ALLOCATE (NREACT(NREAC))
        ALLOCATE (NRHVPI(NRPI))
c  again: some arrays for species distribution of secondaries
c         derived from P..DS and P..PI, above. 
c         for speeding up scoring in update, collide 
        ALLOCATE (IPATEI(NREI,0:NATM))
        ALLOCATE (IPMLEI(NREI,0:NMOL))
        ALLOCATE (IPIOEI(NREI,0:NION))
        ALLOCATE (IPPLEI(NREI,0:NPLS))
        ALLOCATE (IPATPI(NRPI,0:NATM))
        ALLOCATE (IPMLPI(NRPI,0:NMOL))
        ALLOCATE (IPIOPI(NRPI,0:NION))
        ALLOCATE (IPPLPI(NRPI,0:NPLS))
c
        ALLOCATE (LGACX(0:NATM,0:NRCX,0:1))
        ALLOCATE (LGMCX(0:NMOL,0:NRCX,0:1))
        ALLOCATE (LGICX(0:NION,0:NRCX,0:1))
        ALLOCATE (LGAEI(0:NATM,0:NREI))
        ALLOCATE (LGMEI(0:NMOL,0:NREI))
        ALLOCATE (LGIEI(0:NION,0:NREI))
        ALLOCATE (LGAEL(0:NATM,0:NREL,0:1))
        ALLOCATE (LGMEL(0:NMOL,0:NREL,0:1))
        ALLOCATE (LGIEL(0:NION,0:NREL,0:1))
        ALLOCATE (LGPRC(0:NPLS,0:NREC))
        ALLOCATE (LGAPI(0:NATM,0:NRPI,0:1))
        ALLOCATE (LGMPI(0:NMOL,0:NRPI,0:1))
        ALLOCATE (LGIPI(0:NION,0:NRPI,0:1))

        MEM = (MSTOR1*MSTOR2+NMDTA)*8_IL +
     .                      MMDTA*4_IL
 
        WRITE (55+IFOFF,'(A,T25,I15)')
     .        ' COMXS(2) ', MEM
 
      END IF
 
      CALL EIRENE_INIT_CMDTA (ICAL)
 
      RETURN
      END SUBROUTINE EIRENE_ALLOC_COMXS
 
 
      SUBROUTINE EIRENE_DEALLOC_COMXS
 
      IF (.NOT.ALLOCATED(XSTOR)) RETURN
 
      DEALLOCATE (XSTOR)
      DEALLOCATE (XSTORV)
 
 
      DEALLOCATE (TABEI1)
      DEALLOCATE (TABRC1)
      DEALLOCATE (TABPI3)
      DEALLOCATE (TABCX3)
      DEALLOCATE (TABEL3)
      DEALLOCATE (FDLMPI)
      DEALLOCATE (FDLMCX)
      DEALLOCATE (FDLMEL)
      DEALLOCATE (ADDPI)
      DEALLOCATE (ADDCX)
      DEALLOCATE (ADDEL)

      DEALLOCATE (FACRRC) 
      DEALLOCATE (FACRPI) 
      DEALLOCATE (FACREL) 
      DEALLOCATE (FACREI)
      DEALLOCATE (FACRCX) 
 
      DEALLOCATE (PELEI)
      DEALLOCATE (PATEI)
      DEALLOCATE (PMLEI)
      DEALLOCATE (PIOEI)
      DEALLOCATE (PPLEI)
      DEALLOCATE (PELPI)
      DEALLOCATE (PATPI)
      DEALLOCATE (PMLPI)
      DEALLOCATE (PIOPI)
      DEALLOCATE (PPLPI)
      DEALLOCATE (P2ND)
      DEALLOCATE (P2NP)
      DEALLOCATE (P2NDS)
      DEALLOCATE (P2NPI)
 
      DEALLOCATE (EELEI1)
      DEALLOCATE (EHVEI1)
      DEALLOCATE (EELRC1)
      DEALLOCATE (EELPI1)
      DEALLOCATE (EHVPI3)
      DEALLOCATE (EPLPI3)
      DEALLOCATE (EPLCX3)
      DEALLOCATE (EPLEL3)
      DEALLOCATE (EPLOT3)
 
      DEALLOCATE (EATPI)
      DEALLOCATE (EMLPI)
      DEALLOCATE (EIOPI)
      DEALLOCATE (EPLPI)

      DEALLOCATE (EATEI)
      DEALLOCATE (EMLEI)
      DEALLOCATE (EIOEI)
      DEALLOCATE (EPLEI)
 
      DEALLOCATE (MODCOL)
      DEALLOCATE (IESTCX)
      DEALLOCATE (IESTEL)
      DEALLOCATE (IESTPI)
      DEALLOCATE (IESTEI)

      DEALLOCATE (NAEII)
      DEALLOCATE (NMEII)
      DEALLOCATE (NIEII)
      DEALLOCATE (NACXI)
      DEALLOCATE (NMCXI)
      DEALLOCATE (NICXI)
      DEALLOCATE (NAELI)
      DEALLOCATE (NMELI)
      DEALLOCATE (NIELI)
      DEALLOCATE (NAPII)
      DEALLOCATE (NMPII)
      DEALLOCATE (NIPII)
      DEALLOCATE (NPRCI)
      DEALLOCATE (NAEIIM)
      DEALLOCATE (NMEIIM)
      DEALLOCATE (NIEIIM)
      DEALLOCATE (NACXIM)
      DEALLOCATE (NMCXIM)
      DEALLOCATE (NICXIM)
      DEALLOCATE (NAELIM)
      DEALLOCATE (NMELIM)
      DEALLOCATE (NIELIM)
      DEALLOCATE (NAPIIM)
      DEALLOCATE (NMPIIM)
      DEALLOCATE (NIPIIM)
      DEALLOCATE (NPRCIM)
      DEALLOCATE (NPBGKA)
      DEALLOCATE (NPBGKM)
      DEALLOCATE (NPBGKI)
      DEALLOCATE (NPBGKP)
 
      DEALLOCATE (NATPRC_2)
      DEALLOCATE (NMLPRC_2)
      DEALLOCATE (NIOPRC_2)
      DEALLOCATE (NPLPRC_2)
      DEALLOCATE (NPHPRC_2)
      DEALLOCATE (NATPRC)
      DEALLOCATE (NMLPRC)
      DEALLOCATE (NIOPRC)
      DEALLOCATE (NPLPRC)
      DEALLOCATE (NPHPRC)
      DEALLOCATE (N1STX)
      DEALLOCATE (N2NDX)
 
      DEALLOCATE (NSEACX)
      DEALLOCATE (NSEMCX)
      DEALLOCATE (NSEICX)
      DEALLOCATE (NSEAEL)
      DEALLOCATE (NSEMEL)
      DEALLOCATE (NSEIEL)
      DEALLOCATE (NSEPRC)
 
      DEALLOCATE (NREACX)
      DEALLOCATE (NREAPI)
      DEALLOCATE (NREAEL)
      DEALLOCATE (NREAEI)
      DEALLOCATE (JEREAEI)
      DEALLOCATE (NREARC)
      DEALLOCATE (JEREARC)
      DEALLOCATE (NELREI)
      DEALLOCATE (JELREI)
      DEALLOCATE (NREAHV)
      DEALLOCATE (NELREL)
      DEALLOCATE (NELRRC)
      DEALLOCATE (JELRRC)
      DEALLOCATE (NELRPI)
      DEALLOCATE (JELRPI)
      DEALLOCATE (NELRCX)
      DEALLOCATE (NELROT)
      DEALLOCATE (NREAOT)
      DEALLOCATE (NREACT)
      DEALLOCATE (NRHVPI)
      DEALLOCATE (IPATEI)
      DEALLOCATE (IPMLEI)
      DEALLOCATE (IPIOEI)
      DEALLOCATE (IPPLEI)
      DEALLOCATE (IPATPI)
      DEALLOCATE (IPMLPI)
      DEALLOCATE (IPIOPI)
      DEALLOCATE (IPPLPI)
      DEALLOCATE (LGACX)
      DEALLOCATE (LGMCX)
      DEALLOCATE (LGICX)
      DEALLOCATE (LGAEI)
      DEALLOCATE (LGMEI)
      DEALLOCATE (LGIEI)
      DEALLOCATE (LGAEL)
      DEALLOCATE (LGMEL)
      DEALLOCATE (LGIEL)
      DEALLOCATE (LGPRC)
      DEALLOCATE (LGAPI)
      DEALLOCATE (LGMPI)
      DEALLOCATE (LGIPI)
 
      DEALLOCATE (DELPOT)
      DEALLOCATE (FACREA)
      DEALLOCATE (FREACA)
      DEALLOCATE (FREACM)
      DEALLOCATE (FREACI)
      DEALLOCATE (FREACP)
      DEALLOCATE (FREACPH)
      DEALLOCATE (FLDLMA)
      DEALLOCATE (FLDLMM)
      DEALLOCATE (FLDLMI)
      DEALLOCATE (FLDLMP)
      DEALLOCATE (FLDLMPH)
      DEALLOCATE (EELECA)
      DEALLOCATE (EELECM)
      DEALLOCATE (EELECI)
      DEALLOCATE (EELECP)
      DEALLOCATE (EELECPH)
      DEALLOCATE (EBULKA)
      DEALLOCATE (EBULKM)
      DEALLOCATE (EBULKI)
      DEALLOCATE (EBULKP)
      DEALLOCATE (EBULKPH)
      DEALLOCATE (ESCD1A)
      DEALLOCATE (ESCD1M)
      DEALLOCATE (ESCD1I)
      DEALLOCATE (ESCD1P)
      DEALLOCATE (ESCD1PH)
 
      DEALLOCATE (ISWR)
      DEALLOCATE (MODCLF)
      DEALLOCATE (MASSP)
      DEALLOCATE (MASST)
      DEALLOCATE (IFTFLG)
      DEALLOCATE (NRCP)
      DEALLOCATE (NRCA)
      DEALLOCATE (NRCM)
      DEALLOCATE (NRCI)
      DEALLOCATE (NRCPH)
      DEALLOCATE (IREACA)
      DEALLOCATE (IREACM)
      DEALLOCATE (IREACI)
      DEALLOCATE (IREACP)
      DEALLOCATE (IREACPH)
      DEALLOCATE (IBULKA)
      DEALLOCATE (IBULKM)
      DEALLOCATE (IBULKI)
      DEALLOCATE (IBULKP)
      DEALLOCATE (IBULKPH)
      DEALLOCATE (ISCD1A)
      DEALLOCATE (ISCD1M)
      DEALLOCATE (ISCD1I)
      DEALLOCATE (ISCD1P)
      DEALLOCATE (ISCD1PH)
      DEALLOCATE (ISCD2A)
      DEALLOCATE (ISCD2M)
      DEALLOCATE (ISCD2I)
      DEALLOCATE (ISCD2P)
      DEALLOCATE (ISCD2PH)
      DEALLOCATE (ISCD3A)
      DEALLOCATE (ISCD3M)
      DEALLOCATE (ISCD3I)
      DEALLOCATE (ISCD3P)
      DEALLOCATE (ISCD3PH)
      DEALLOCATE (ISCD4A)
      DEALLOCATE (ISCD4M)
      DEALLOCATE (ISCD4I)
      DEALLOCATE (ISCD4P)
      DEALLOCATE (ISCD4PH)
      DEALLOCATE (ISCDEA)
      DEALLOCATE (ISCDEM)
      DEALLOCATE (ISCDEI)
      DEALLOCATE (ISCDEP)
      DEALLOCATE (ISCDEPH)
      DEALLOCATE (IESTMA)
      DEALLOCATE (IESTMM)
      DEALLOCATE (IESTMI)
      DEALLOCATE (IESTMPH)
      DEALLOCATE (IBGKA )
      DEALLOCATE (IBGKM )
      DEALLOCATE (IBGKI )
      DEALLOCATE (IBGKPH)
      DEALLOCATE (REAC_NAME)
 
!pb      DEALLOCATE (REACDAT)
      CALL EIRENE_FREE_REACDAT
      DEALLOCATE (REACLINES)
 
      RETURN
      END SUBROUTINE EIRENE_DEALLOC_COMXS
 
 
      SUBROUTINE EIRENE_INIT_CMDTA (ICAL)
cdr  initialize (nullify) A&M data
cdr  ical=1:  ??
cdr  ical=2:  ??
 
      INTEGER, INTENT(IN) :: ICAL
      INTEGER :: IREAC, IL
 
      IF (ICAL == 1) THEN
        NAEII   = 0
        NMEII   = 0
        NIEII   = 0
        NACXI   = 0
        NMCXI   = 0
        NICXI   = 0
        NAELI   = 0
        NMELI   = 0
        NIELI   = 0
        NAPII   = 0
        NMPII   = 0
        NIPII   = 0
        NPRCI   = 0
        NAEIIM  = 0
        NMEIIM  = 0
        NIEIIM  = 0
        NACXIM  = 0
        NMCXIM  = 0
        NICXIM  = 0
        NAELIM  = 0
        NMELIM  = 0
        NIELIM  = 0
        NAPIIM  = 0
        NMPIIM  = 0
        NIPIIM  = 0
        NPRCIM  = 0
        NPBGKA  = 0
        NPBGKM  = 0
        NPBGKI  = 0
        NPBGKP  = 0
 
        NSEACX  = 0
        NSEMCX  = 0
        NSEICX  = 0
        NSEAEL  = 0
        NSEMEL  = 0
        NSEIEL  = 0
 
        DELPOT  = 0._DP
        FACREA(:,1) = 1._DP
        FACREA(:,2) = 0._DP
        FREACA  = 0._DP
        FREACM  = 0._DP
        FREACI  = 0._DP
        FREACP  = 0._DP
        FREACPH = 0._DP
        FLDLMA  = 0._DP
        FLDLMM  = 0._DP
        FLDLMI  = 0._DP
        FLDLMP  = 0._DP
        FLDLMPH = 0._DP
        EELECA  = 0._DP
        EELECM  = 0._DP
        EELECI  = 0._DP
        EELECP  = 0._DP
        EELECPH = 0._DP
        EBULKA  = 0._DP
        EBULKM  = 0._DP
        EBULKI  = 0._DP
        EBULKP  = 0._DP
        EBULKPH = 0._DP
        ESCD1A  = 0._DP
        ESCD1M  = 0._DP
        ESCD1I  = 0._DP
        ESCD1P  = 0._DP
        ESCD1PH = 0._DP
 
        ZMFPTHI = 0._DP
        TDGTEMX = 1.E-30_DP
 
        ISWR    = 0
        MODCLF  = 0
        MASSP   = 0
        MASST   = 0
        IFTFLG  = 0
        NRCP    = 0
        NRCA    = 0
        NRCM    = 0
        NRCI    = 0
        NRCPH   = 0
        IREACA  = 0
        IREACM  = 0
        IREACI  = 0
        IREACP  = 0
        IREACPH = 0
        IBULKA  = 0
        IBULKM  = 0
        IBULKI  = 0
        IBULKP  = 0
        IBULKPH = 0
        ISCD1A  = 0
        ISCD1M  = 0
        ISCD1I  = 0
        ISCD1P  = 0
        ISCD1PH = 0
        ISCD2A  = 0
        ISCD2M  = 0
        ISCD2I  = 0
        ISCD2P  = 0
        ISCD2PH = 0
        ISCD3A  = 0
        ISCD3M  = 0
        ISCD3I  = 0
        ISCD3P  = 0
        ISCD3PH = 0
        ISCD4A  = 0
        ISCD4M  = 0
        ISCD4I  = 0
        ISCD4P  = 0
        ISCD4PH = 0
        ISCDEA  = 0
        ISCDEM  = 0
        ISCDEI  = 0
        ISCDEP  = 0
        ISCDEPH = 0
        IESTMA  = 0
        IESTMM  = 0
        IESTMI  = 0
        IESTMPH = 0
        IBGKA   = 0
        IBGKM   = 0
        IBGKI   = 0
        IBGKPH  = 0
 
        REAC_NAME = REPEAT(' ',LEN(REAC_NAME(1)))
 
        idreac  = 0
 
        XSTORV  = 0._DP
 
        DO IREAC= -11, NREAC
          REACDAT(IREAC)%LPOT   = .FALSE.
          REACDAT(IREAC)%LCRS   = .FALSE.
          REACDAT(IREAC)%LRTC   = .FALSE.
          REACDAT(IREAC)%LRTCMW = .FALSE.
          REACDAT(IREAC)%LRTCEW = .FALSE.
          REACDAT(IREAC)%LOTH   = .FALSE.
          REACDAT(IREAC)%LPHR   = .FALSE.
          REACDAT(IREAC)%NOSEC  = 0
          NULLIFY(REACDAT(IREAC)%POT)
          NULLIFY(REACDAT(IREAC)%CRS)
          NULLIFY(REACDAT(IREAC)%RTC)
          NULLIFY(REACDAT(IREAC)%RTCMW)
          NULLIFY(REACDAT(IREAC)%RTCEW)
          NULLIFY(REACDAT(IREAC)%OTH)
          NULLIFY(REACDAT(IREAC)%PHR)
        END DO
 
        DO IL = 1, NREAC_LINES
          REACLINES(IL)%NO = 0
          REACLINES(IL)%MT = 0
          REACLINES(IL)%MP = 0
          REACLINES(IL)%IZ = 0
          REACLINES(IL)%JFEXMN = 0
          REACLINES(IL)%JFEXMX = 0
          REACLINES(IL)%NCONST = 0
          REACLINES(IL)%RMN = 0._DP
          REACLINES(IL)%RMX = 0._DP
          REACLINES(IL)%DPP = 0._DP
          REACLINES(IL)%FP = 0._DP
          REACLINES(IL)%CONST = 0._DP
          REACLINES(IL)%FILE = REPEAT(' ',8)
          REACLINES(IL)%REAC_STRING = REPEAT(' ',50)
          REACLINES(IL)%H_SELECT = REPEAT(' ',4)
          REACLINES(IL)%REACTYP = REPEAT(' ',3)
          REACLINES(IL)%ELEMENT = REPEAT(' ',2)
        END DO
 
        IRLINES = 0
 
      ELSE IF (ICAL == 2) THEN
 
        XSTOR  = 0._DP
 
        TABEI1  = 0._DP
        TABRC1  = 0._DP
        TABPI3  = 0._DP
        TABCX3  = 0._DP
        TABEL3  = 0._DP
        FDLMPI  = 0._DP
        FDLMCX  = 0._DP
        FDLMEL  = 0._DP
        ADDPI   = 0._DP
        ADDCX   = 0._DP
        ADDEL   = 0._DP

        FACRRC(:,1) = 1._DP
        FACRRC(:,2) = 0._DP
        FACRPI(:,1) = 1._DP 
        FACRPI(:,2) = 0._DP 
        FACREL(:,1) = 1._DP 
        FACREL(:,2) = 0._DP 
        FACREI(:,1) = 1._DP
        FACREI(:,2) = 0._DP
        FACRCX(:,1) = 1._DP
        FACRCX(:,2) = 0._DP
 
        PELEI   = 0._DP
        PATEI   = 0._DP
        PMLEI   = 0._DP
        PIOEI   = 0._DP
        PPLEI   = 0._DP
        PELPI   = 0._DP
        PATPI   = 0._DP
        PMLPI   = 0._DP
        PIOPI   = 0._DP
        PPLPI   = 0._DP
        P2ND    = 0._DP
        P2NP    = 0._DP
        P2NDS   = 0._DP
        P2NPI   = 0._DP
 
        EELEI1  = 0._DP
        EHVEI1  = 0._DP
        EELRC1  = 0._DP
        EELPI1  = 0._DP
        EHVPI3  = 0._DP
        EPLPI3  = 0._DP
        EPLCX3  = 0._DP
        EPLEL3  = 0._DP
        EPLOT3  = 0._DP
 
        EATPI   = 0._DP
        EMLPI   = 0._DP
        EIOPI   = 0._DP
        EPLPI   = 0._DP
        EATEI   = 0._DP
        EMLEI   = 0._DP
        EIOEI   = 0._DP
        EPLEI   = 0._DP
 
        MODCOL  = 0
        IESTCX  = 0
        IESTEL  = 0
        IESTPI  = 0
        IESTEI  = 0
 
        NATPRC  = 0
        NMLPRC  = 0
        NIOPRC  = 0
        NPLPRC  = 0
        NPHPRC  = 0
        NATPRC_2= 0
        NMLPRC_2= 0
        NIOPRC_2= 0
        NPLPRC_2= 0
        NPHPRC_2= 0
        N1STX   = 0
        N2NDX   = 0
 
        NRPII   = 0
        NREII   = 0
        NRCXI   = 0
        NRELI   = 0
        NRRCI   = 0
        NRBGI   = 0
        NSEPRC  = 0
 
        NREACX  = 0
        NREAPI  = 0
        NREAEL  = 0
        NREAEI  = 0
        JEREAEI = 0
        NREARC  = 0
        JEREARC = 0
        NELREI  = 0
        JELREI  = 0
        NREAHV  = 0
        NELREL  = 0
        NELRRC  = 0
        JELRRC  = 0
        NELRPI  = 0
        JELRPI  = 0
        NELRCX  = 0
        NELROT  = 0
        NREAOT  = 0
        NREACT  = 0
        NRHVPI  = 0
        IPATEI  = 0
        IPMLEI  = 0
        IPIOEI  = 0
        IPPLEI  = 0
        IPATPI  = 0
        IPMLPI  = 0
        IPIOPI  = 0
        IPPLPI  = 0
        LGACX   = 0
        LGMCX   = 0
        LGICX   = 0
        LGAEI   = 0
        LGMEI   = 0
        LGIEI   = 0
        LGAEL   = 0
        LGMEL   = 0
        LGIEL   = 0
        LGPRC   = 0
        LGAPI   = 0
        LGMPI   = 0
        LGIPI   = 0
 
      END IF
 
      RETURN
      END SUBROUTINE EIRENE_INIT_CMDTA
 
 
      SUBROUTINE EIRENE_WRITE_CMDTA
cdr  read and write A&M data onto fort 13., controlled by NFILEL option (input block 1)
 
      WRITE (13+IFOFF)
     . TABEI1 ,TABRC1 ,TABPI3 ,TABCX3 ,TABEL3 ,
     . FDLMPI ,FDLMCX ,FDLMEL ,
     . ADDPI  ,ADDCX  ,ADDEL  ,
     . FACRRC ,FACRPI ,FACREL ,FACREI ,FACRCX ,
 
     . PELEI  ,PATEI  ,PMLEI  ,PIOEI  ,PPLEI  ,
     . PELPI  ,PATPI  ,PMLPI  ,PIOPI  ,PPLPI  ,
     . P2ND   ,P2NP   ,P2NDS  ,P2NPI  ,
 
     . EELEI1 ,EELRC1 ,EELPI1 ,
     . EHVEI1 ,EHVPI3 ,
     . EPLPI3 ,EPLCX3 ,EPLEL3 ,EPLOT3 ,
 
     . EATPI  ,EMLPI  ,EIOPI  ,EPLPI  ,
     . EATEI  ,EMLEI  ,EIOEI  ,EPLEI
 
      WRITE (13+IFOFF)
     . MODCOL ,IESTCX ,IESTEL ,IESTPI ,IESTEI ,
     . NAEII  ,NMEII  ,NIEII  ,NACXI  ,NMCXI  ,NICXI  ,
     . NAELI  ,NMELI  ,NIELI  ,NAPII  ,NMPII  ,NIPII  ,NPRCI  ,
     . NAEIIM ,NMEIIM ,NIEIIM ,NACXIM ,NMCXIM ,NICXIM ,
     . NAELIM ,NMELIM ,NIELIM ,NAPIIM ,NMPIIM ,NIPIIM ,NPRCIM ,
     . NPBGKA ,NPBGKM ,NPBGKI ,NPBGKP ,
     . NATPRC ,NMLPRC ,NIOPRC ,NPLPRC ,NPHPRC ,
     . NATPRC_2 ,NMLPRC_2 ,NIOPRC_2 ,NPLPRC_2 ,NPHPRC_2 ,
     . N1STX  ,N2NDX  ,
 
     . NRPII  ,NREII  ,NRCXI  ,NRELI  ,NRRCI  ,NRBGI  ,
 
     . NSEACX ,NSEMCX ,NSEICX ,NSEAEL ,NSEMEL ,NSEIEL ,NSEPRC ,
     . NREACX ,NREAPI ,NREAEL ,NREAEI ,JEREAEI,NREARC ,JEREARC,
     . NELREI ,JELREI ,NREAHV ,NELREL ,NELRRC ,JELRRC ,NELRPI ,JELRPI ,
     . NELRCX ,NELROT ,NREAOT ,NREACT ,NRHVPI ,
     . IPATEI ,IPMLEI ,IPIOEI ,IPPLEI ,IPATPI ,IPMLPI ,IPIOPI ,IPPLPI ,
     . LGACX  ,LGMCX  ,LGICX  ,LGAEI  ,LGMEI  ,LGIEI  ,
     . LGAEL  ,LGMEL  ,LGIEL  ,LGPRC  ,LGAPI  ,LGMPI  ,LGIPI
 
      RETURN
      END SUBROUTINE EIRENE_WRITE_CMDTA
 
 
      SUBROUTINE EIRENE_READ_CMDTA
 
      READ (13+IFOFF)
     . TABEI1 ,TABRC1 ,TABPI3 ,TABCX3 ,TABEL3 ,
     . FDLMPI ,FDLMCX ,FDLMEL ,
     . ADDPI  ,ADDCX  ,ADDEL  ,
     . FACRRC ,FACRPI ,FACREL ,FACREI ,FACRCX ,
 
     . PELEI  ,PATEI  ,PMLEI  ,PIOEI  ,PPLEI  ,
     . PELPI  ,PATPI  ,PMLPI  ,PIOPI  ,PPLPI  ,
     . P2ND   ,P2NP   ,P2NDS  ,P2NPI  ,
 
     . EELEI1 ,EELRC1 ,EELPI1 ,
     . EHVEI1 ,EHVPI3 ,
     . EPLPI3 ,EPLCX3 ,EPLEL3 ,EPLOT3 ,
 
     . EATPI  ,EMLPI  ,EIOPI  ,EPLPI  ,
     . EATEI  ,EMLEI  ,EIOEI  ,EPLEI
 
      READ (13+IFOFF)
     . MODCOL ,IESTCX ,IESTEL ,IESTPI ,IESTEI ,
     . NAEII  ,NMEII  ,NIEII  ,NACXI  ,NMCXI  ,NICXI  ,
     . NAELI  ,NMELI  ,NIELI  ,NAPII  ,NMPII  ,NIPII  ,NPRCI  ,
     . NAEIIM ,NMEIIM ,NIEIIM ,NACXIM ,NMCXIM ,NICXIM ,
     . NAELIM ,NMELIM ,NIELIM ,NAPIIM ,NMPIIM ,NIPIIM ,NPRCIM ,
     . NPBGKA ,NPBGKM ,NPBGKI ,NPBGKP ,
     . NATPRC ,NMLPRC ,NIOPRC ,NPLPRC ,NPHPRC ,
     . NATPRC_2 ,NMLPRC_2 ,NIOPRC_2 ,NPLPRC_2 ,NPHPRC_2 ,
     . N1STX  ,N2NDX  ,
 
     . NRPII  ,NREII  ,NRCXI  ,NRELI  ,NRRCI  ,NRBGI  ,
 
     . NSEACX ,NSEMCX ,NSEICX ,NSEAEL ,NSEMEL ,NSEIEL ,NSEPRC ,
     . NREACX ,NREAPI ,NREAEL ,NREAEI ,JEREAEI,NREARC ,JEREARC,
     . NELREI ,JELREI ,NREAHV ,NELREL ,NELRRC ,JELRRC ,NELRPI ,JELRPI ,
     . NELRCX ,NELROT ,NREAOT ,NREACT ,NRHVPI ,
     . IPATEI ,IPMLEI ,IPIOEI ,IPPLEI ,IPATPI ,IPMLPI ,IPIOPI ,IPPLPI ,
     . LGACX  ,LGMCX  ,LGICX  ,LGAEI  ,LGMEI  ,LGIEI  ,
     . LGAEL  ,LGMEL  ,LGIEL  ,LGPRC  ,LGAPI  ,LGMPI  ,LGIPI
 
      RETURN
      END SUBROUTINE EIRENE_READ_CMDTA
 
 
      SUBROUTINE EIRENE_CMDTA_XDR (IUN)
 
      INTEGER, INTENT(IN) :: IUN
      INTEGER :: IHELP(1)
c
      CALL FXDRDBL (IUN,TABEI1,NREI*NSTORDR)
      CALL FXDRDBL (IUN,TABRC1,NREC*NSTORDR)
      CALL FXDRDBL (IUN,TABPI3,NRPI*NSTORDR*NSTORDT)
      CALL FXDRDBL (IUN,TABCX3,NRCX*NSTORDR*NSTORDT)
      CALL FXDRDBL (IUN,TABEL3,NREL*NSTORDR*NSTORDT)
      CALL FXDRDBL (IUN,FDLMPI,NRPI)
      CALL FXDRDBL (IUN,FDLMCX,NRCX)
      CALL FXDRDBL (IUN,FDLMEL,NREL)
      CALL FXDRDBL (IUN,ADDPI,NRPI*NPLS)
      CALL FXDRDBL (IUN,ADDCX,NRCX*NPLS)
      CALL FXDRDBL (IUN,ADDEL,NREL*NPLS)
      CALL FXDRDBL (IUN,FACRRC,NREC*2)
      CALL FXDRDBL (IUN,FACRPI,NRPI*2)
      CALL FXDRDBL (IUN,FACREL,NREL*2)
      CALL FXDRDBL (IUN,FACREI,NREI*2)
      CALL FXDRDBL (IUN,FACRCX,NRCX*2)
 
      CALL FXDRDBL (IUN,PELEI,NREI)
      CALL FXDRDBL (IUN,PATEI,NREI*(NATM+1))
      CALL FXDRDBL (IUN,PMLEI,NREI*(NMOL+1))
      CALL FXDRDBL (IUN,PIOEI,NREI*(NION+1))
      CALL FXDRDBL (IUN,PPLEI,NREI*(NPLS+1))
      CALL FXDRDBL (IUN,PELPI,NRPI)
      CALL FXDRDBL (IUN,PATPI,NRPI*(NATM+1))
      CALL FXDRDBL (IUN,PMLPI,NRPI*(NMOL+1))
      CALL FXDRDBL (IUN,PIOPI,NRPI*(NION+1))
      CALL FXDRDBL (IUN,PPLPI,NRPI*(NPLS+1))
      CALL FXDRDBL (IUN,P2ND,NREI*(NSPZ+1))
      CALL FXDRDBL (IUN,P2NP,NRPI*(NSPZ+1))
      CALL FXDRDBL (IUN,P2NDS,NREI)
      CALL FXDRDBL (IUN,P2NPI,NRPI)
 
      CALL FXDRDBL (IUN,EELEI1,NREI*NSTORDR)
      CALL FXDRDBL (IUN,EHVEI1,NREI*NSTORDR)
      CALL FXDRDBL (IUN,EELRC1,NREC*NSTORDR)
      CALL FXDRDBL (IUN,EELPI1,NRPI*NSTORDR)
      CALL FXDRDBL (IUN,EHVPI3,NRPI*NSTORDR*NSTORDT)
      CALL FXDRDBL (IUN,EPLPI3,NRPI*NSTORDR*NSTORDT)
      CALL FXDRDBL (IUN,EPLCX3,NRCX*NSTORDR*NSTORDT)
      CALL FXDRDBL (IUN,EPLEL3,NREL*NSTORDR*NSTORDT)
      CALL FXDRDBL (IUN,EPLOT3,NROT*NSTORDR*NSTORDT)

      CALL FXDRDBL (IUN,EATPI,NRPI*(NATM+1)*2)
      CALL FXDRDBL (IUN,EMLPI,NRPI*(NMOL+1)*2)
      CALL FXDRDBL (IUN,EIOPI,NRPI*(NION+1)*2)
      CALL FXDRDBL (IUN,EPLPI,NRPI*(NPLS+1)*2)

      CALL FXDRDBL (IUN,EATEI,NREI*(NATM+1)*2)
      CALL FXDRDBL (IUN,EMLEI,NREI*(NMOL+1)*2)
      CALL FXDRDBL (IUN,EIOEI,NREI*(NION+1)*2)
      CALL FXDRDBL (IUN,EPLEI,NREI*(NPLS+1)*2)
 
c
      CALL FXDRINT (IUN,MODCOL ,7*5*MSTOR0)
      CALL FXDRINT (IUN,IESTCX ,3*NRCX)
      CALL FXDRINT (IUN,IESTEL ,3*NREL)
      CALL FXDRINT (IUN,IESTPI ,3*NRPI)
      CALL FXDRINT (IUN,IESTEI ,3*NREI)
      CALL FXDRINT (IUN,NAEII  ,NATM)
      CALL FXDRINT (IUN,NMEII  ,NMOL)
      CALL FXDRINT (IUN,NIEII  ,NION)
      CALL FXDRINT (IUN,NACXI  ,NATM)
      CALL FXDRINT (IUN,NMCXI  ,NMOL)
      CALL FXDRINT (IUN,NICXI  ,NION)
      CALL FXDRINT (IUN,NAELI  ,NATM)
      CALL FXDRINT (IUN,NMELI  ,NMOL)
      CALL FXDRINT (IUN,NIELI  ,NION)
      CALL FXDRINT (IUN,NAPII  ,NATM)
      CALL FXDRINT (IUN,NMPII  ,NMOL)
      CALL FXDRINT (IUN,NIPII  ,NION)
      CALL FXDRINT (IUN,NPRCI  ,NPLS)
      CALL FXDRINT (IUN,NAEIIM ,NATM)
      CALL FXDRINT (IUN,NMEIIM ,NMOL)
      CALL FXDRINT (IUN,NIEIIM ,NION)
      CALL FXDRINT (IUN,NACXIM ,NATM)
      CALL FXDRINT (IUN,NMCXIM ,NMOL)
      CALL FXDRINT (IUN,NICXIM ,NION)
      CALL FXDRINT (IUN,NAELIM ,NATM)
      CALL FXDRINT (IUN,NMELIM ,NMOL)
      CALL FXDRINT (IUN,NIELIM ,NION)
      CALL FXDRINT (IUN,NAPIIM ,NATM)
      CALL FXDRINT (IUN,NMPIIM ,NMOL)
      CALL FXDRINT (IUN,NIPIIM ,NION)
      CALL FXDRINT (IUN,NPRCIM ,NPLS)
      CALL FXDRINT (IUN,NPBGKA ,NATM)
      CALL FXDRINT (IUN,NPBGKM ,NMOL)
      CALL FXDRINT (IUN,NPBGKI ,NION)
      CALL FXDRINT (IUN,NPBGKP ,2*NPLS)
      CALL FXDRINT (IUN,NATPRC ,NREC)
      CALL FXDRINT (IUN,NMLPRC ,NREC)
      CALL FXDRINT (IUN,NIOPRC ,NREC)
      CALL FXDRINT (IUN,NPLPRC ,NREC)
      CALL FXDRINT (IUN,NPHPRC ,NREC)
      CALL FXDRINT (IUN,NATPRC_2 ,NREC)
      CALL FXDRINT (IUN,NMLPRC_2 ,NREC)
      CALL FXDRINT (IUN,NIOPRC_2 ,NREC)
      CALL FXDRINT (IUN,NPLPRC_2 ,NREC)
      CALL FXDRINT (IUN,NPHPRC_2 ,NREC)
      CALL FXDRINT (IUN,N1STX  ,3*NRCX)
      CALL FXDRINT (IUN,N2NDX  ,3*NRCX)
 
      IHELP(1) = NRPII
      CALL FXDRINT (IUN, IHELP ,1)
      IHELP(1) = NREII
      CALL FXDRINT (IUN, IHELP ,1)
      IHELP(1) = NRCXI
      CALL FXDRINT (IUN, IHELP ,1)
      IHELP(1) = NRELI
      CALL FXDRINT (IUN, IHELP ,1)
      IHELP(1) = NRRCI
      CALL FXDRINT (IUN, IHELP ,1)
      IHELP(1) = NRBGI
      CALL FXDRINT (IUN, IHELP ,1)
 
      CALL FXDRINT (IUN,NSEACX ,5*NATM*NPLS)
      CALL FXDRINT (IUN,NSEMCX ,5*NMOL*NPLS)
      CALL FXDRINT (IUN,NSEICX ,5*NION*NPLS)
      CALL FXDRINT (IUN,NSEAEL ,5*NATM*NPLS)
      CALL FXDRINT (IUN,NSEMEL ,5*NMOL*NPLS)
      CALL FXDRINT (IUN,NSEIEL ,5*NION*NPLS)
      CALL FXDRINT (IUN,NSEPRC ,5*NREC)
      CALL FXDRINT (IUN,NREACX ,NRCX)
      CALL FXDRINT (IUN,NREAPI ,NRPI)
      CALL FXDRINT (IUN,NREAEL ,NREL)
      CALL FXDRINT (IUN,NREAEI ,NREI)
      CALL FXDRINT (IUN,JEREAEI,NREI)
      CALL FXDRINT (IUN,NREARC ,NREC)
      CALL FXDRINT (IUN,JEREARC,NREC)
      CALL FXDRINT (IUN,NELREI ,NREI)
      CALL FXDRINT (IUN,JELREI ,NREI)
      CALL FXDRINT (IUN,NREAHV ,NREI)
      CALL FXDRINT (IUN,NELREL ,NREL)
      CALL FXDRINT (IUN,NELRRC ,NREC)
      CALL FXDRINT (IUN,JELRRC ,NREC)
      CALL FXDRINT (IUN,NELRPI ,NRPI)
      CALL FXDRINT (IUN,JELRPI ,NRPI)
      CALL FXDRINT (IUN,NELRCX ,NRCX)
      CALL FXDRINT (IUN,NELROT ,NROT)
      CALL FXDRINT (IUN,NREAOT ,NROT)
      CALL FXDRINT (IUN,NREACT ,NREAC)
      CALL FXDRINT (IUN,NRHVPI ,NRPI)
      CALL FXDRINT (IUN,IPATEI ,NREI*(NATM+1))
      CALL FXDRINT (IUN,IPMLEI ,NREI*(NMOL+1))
      CALL FXDRINT (IUN,IPIOEI ,NREI*(NION+1))
      CALL FXDRINT (IUN,IPPLEI ,NREI*(NPLS+1))
      CALL FXDRINT (IUN,IPATPI ,NRPI*(NATM+1))
      CALL FXDRINT (IUN,IPMLPI ,NRPI*(NMOL+1))
      CALL FXDRINT (IUN,IPIOPI ,NRPI*(NION+1))
      CALL FXDRINT (IUN,IPPLPI ,NRPI*(NPLS+1))
 
      CALL FXDRINT (IUN,LGACX  ,2*(NATM+1)*(NRCX+1))
      CALL FXDRINT (IUN,LGMCX  ,2*(NMOL+1)*(NRCX+1))
      CALL FXDRINT (IUN,LGICX  ,2*(NION+1)*(NRCX+1))
      CALL FXDRINT (IUN,LGAEI  ,(NATM+1)*(NREI+1))
      CALL FXDRINT (IUN,LGMEI  ,(NMOL+1)*(NREI+1))
      CALL FXDRINT (IUN,LGIEI  ,(NION+1)*(NREI+1))
      CALL FXDRINT (IUN,LGAEL  ,2*(NATM+1)*(NREL+1))
      CALL FXDRINT (IUN,LGMEL  ,2*(NMOL+1)*(NREL+1))
      CALL FXDRINT (IUN,LGIEL  ,2*(NION+1)*(NREL+1))
      CALL FXDRINT (IUN,LGPRC  ,(NPLS+1)*(NREC+1))
      CALL FXDRINT (IUN,LGAPI  ,2*(NATM+1)*(NRPI+1))
      CALL FXDRINT (IUN,LGMPI  ,2*(NMOL+1)*(NRPI+1))
      CALL FXDRINT (IUN,LGIPI  ,2*(NION+1)*(NRPI+1))
 
      RETURN
      END SUBROUTINE EIRENE_CMDTA_XDR
 
 
      SUBROUTINE EIRENE_WRITE_CMAMF
 
      INTEGER :: IR
 
      WRITE (13+IFOFF)
     . DELPOT, FACREA,
     . FREACA, FREACM, FREACI, FREACP, FREACPH,
     . FLDLMA, FLDLMM, FLDLMI, FLDLMP, FLDLMPH,
     . EELECA, EELECM, EELECI, EELECP, EELECPH,
     . EBULKA, EBULKM, EBULKI, EBULKP, EBULKPH,
     . ESCD1A, ESCD1M, ESCD1I, ESCD1P, ESCD1PH,
 
     . NREACI, ISWR,   MODCLF, MASSP,  MASST,  IFTFLG,
     . NRCP,   NRCA,   NRCM,   NRCI, NRCPH,
     . IREACA, IREACM, IREACI, IREACP, IREACPH,
     . IBULKA, IBULKM, IBULKI, IBULKP, IBULKPH,
     . ISCD1A, ISCD1M, ISCD1I, ISCD1P, ISCD1PH,
     . ISCD2A, ISCD2M, ISCD2I, ISCD2P, ISCD2PH,
     . ISCD3A, ISCD3M, ISCD3I, ISCD3P, ISCD3PH,
     . ISCD4A, ISCD4M, ISCD4I, ISCD4P, ISCD4PH,
     . ISCDEA, ISCDEM, ISCDEI, ISCDEP, ISCDEPH,
     . IESTMA, IESTMM, IESTMI, IESTMPH, IBGKA , IBGKM , IBGKI, IBGKPH
 
      DO IR=1,NREACI
        WRITE (13+IFOFF) 
     .             REACDAT(IR)%LPOT,REACDAT(IR)%LCRS,
     .             REACDAT(IR)%LRTC,
     .             REACDAT(IR)%LRTCMW,REACDAT(IR)%LRTCEW,
     .             REACDAT(IR)%LOTH,REACDAT(IR)%LPHR,
     .             REACDAT(IR)%NOSEC
        IF (REACDAT(IR)%LPOT)   CALL
     .  EIRENE_WRITE_FIT_FORM(REACDAT(IR)%POT)
        IF (REACDAT(IR)%LCRS)   CALL
     .  EIRENE_WRITE_FIT_FORM(REACDAT(IR)%CRS)
        IF (REACDAT(IR)%LRTC)   CALL
     .  EIRENE_WRITE_FIT_FORM(REACDAT(IR)%RTC)
        IF (REACDAT(IR)%LRTCMW) CALL
     .  EIRENE_WRITE_FIT_FORM(REACDAT(IR)%RTCMW)
        IF (REACDAT(IR)%LRTCEW) CALL
     .  EIRENE_WRITE_FIT_FORM(REACDAT(IR)%RTCEW)
        IF (REACDAT(IR)%LOTH)   CALL
     .  EIRENE_WRITE_FIT_FORM(REACDAT(IR)%OTH)
        IF (REACDAT(IR)%LPHR)   CALL
     .  EIRENE_WRITE_FIT_FORM(REACDAT(IR)%PHR)
      END DO
 
      RETURN
 
      CONTAINS
 
        SUBROUTINE EIRENE_WRITE_FIT_FORM (RP)
        TYPE(FIT_FORMS),POINTER :: RP
 
        WRITE (13+IFOFF) RP%IFIT
 
        IF (RP%IFIT < 0) THEN
! DATA FOR PHOTONIC LINE SHAPE, BROADENING.  SPECIAL FORMAT
          WRITE (13+IFOFF) RP%LINE%E0, RP%LINE%E1, RP%LINE%AIK,
     .               RP%LINE%G1, RP%LINE%G2, RP%LINE%C2,
     .               RP%LINE%C3, RP%LINE%C4, RP%LINE%C6,
     .               RP%LINE%B12, RP%LINE%B21, RP%LINE%C6A
          WRITE (13+IFOFF) RP%LINE%IGND, RP%LINE%IRCART, 
     .               RP%LINE%IPROFILETYPE,
     .               RP%LINE%IFREMD, RP%LINE%NRJPRT, RP%LINE%IPLSC6,
     .               RP%LINE%IMESS
          WRITE (13+IFOFF) RP%LINE%REACNAME, RP%LINE%KENN
 
        ELSE IF (RP%IFIT <= 2) THEN
! DATA FOR FIT EXPRESSIONS  (E.G. POLYNOMIAL, IN CASE OF HYDHEL DATABASE)
          WRITE (13+IFOFF) RP%POLY%IFEXMN,RP%POLY%IFEXMX,
     .               UBOUND(RP%POLY%DBLPOL)
          WRITE (13+IFOFF) RP%POLY%DBLPOL,
     .               RP%POLY%RCMN,RP%POLY%RCMX,RP%POLY%FPARM
 
        ELSE
! DATA FROM FORMATTED TABLES (E.G. ADAS)
        END IF
 
        END SUBROUTINE EIRENE_WRITE_FIT_FORM
 
      END SUBROUTINE EIRENE_WRITE_CMAMF
 
 
      SUBROUTINE EIRENE_READ_CMAMF
 
      INTEGER :: IR
 
      READ (13+IFOFF)
     . DELPOT, FACREA,
     . FREACA, FREACM, FREACI, FREACP, FREACPH,
     . FLDLMA, FLDLMM, FLDLMI, FLDLMP, FLDLMPH,
     . EELECA, EELECM, EELECI, EELECP, EELECPH,
     . EBULKA, EBULKM, EBULKI, EBULKP, EBULKPH,
     . ESCD1A, ESCD1M, ESCD1I, ESCD1P, ESCD1PH,
 
     . NREACI, ISWR,   MODCLF, MASSP,  MASST,  IFTFLG,
     . NRCP,   NRCA,   NRCM,   NRCI, NRCPH,
     . IREACA, IREACM, IREACI, IREACP, IREACPH,
     . IBULKA, IBULKM, IBULKI, IBULKP, IBULKPH,
     . ISCD1A, ISCD1M, ISCD1I, ISCD1P, ISCD1PH,
     . ISCD2A, ISCD2M, ISCD2I, ISCD2P, ISCD2PH,
     . ISCD3A, ISCD3M, ISCD3I, ISCD3P, ISCD3PH,
     . ISCD4A, ISCD4M, ISCD4I, ISCD4P, ISCD4PH,
     . ISCDEA, ISCDEM, ISCDEI, ISCDEP, ISCDEPH,
     . IESTMA, IESTMM, IESTMI, IESTMPH, IBGKA , IBGKM , IBGKI, IBGKPH
 
      DO IR=1,NREACI
        READ (13+IFOFF) 
     .            REACDAT(IR)%LPOT,REACDAT(IR)%LCRS,
     .            REACDAT(IR)%LRTC,
     .            REACDAT(IR)%LRTCMW,REACDAT(IR)%LRTCEW,
     .            REACDAT(IR)%LOTH,REACDAT(IR)%LPHR,
     .            REACDAT(IR)%NOSEC
 
        IF (REACDAT(IR)%LPOT) CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%POT)
        IF (REACDAT(IR)%LCRS) CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%CRS)
        IF (REACDAT(IR)%LRTC) CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%RTC)
        IF (REACDAT(IR)%LRTCMW) 
     .    CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%RTCMW)
        IF (REACDAT(IR)%LRTCEW) 
     .    CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%RTCEW)
        IF (REACDAT(IR)%LOTH) CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%OTH)
        IF (REACDAT(IR)%LPHR) CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%PHR)
      END DO
 
 
      RETURN
      CONTAINS
 
        SUBROUTINE EIRENE_READ_FIT_FORM (RP)
        TYPE(FIT_FORMS),POINTER :: RP
        INTEGER :: ND, ND2
 
        READ (13+IFOFF) RP%IFIT
 
        IF (RP%IFIT < 0) THEN
! DATA FOR PHOTONIC LINE
          IF (.NOT.ASSOCIATED(RP%LINE)) ALLOCATE (RP%LINE)
          READ (13+IFOFF) RP%LINE%E0, RP%LINE%E1, RP%LINE%AIK,
     .              RP%LINE%G1, RP%LINE%G2, RP%LINE%C2,
     .              RP%LINE%C3, RP%LINE%C4, RP%LINE%C6,
     .              RP%LINE%B12, RP%LINE%B21, RP%LINE%C6A
          READ (13+IFOFF) RP%LINE%IGND, RP%LINE%IRCART, 
     .              RP%LINE%IPROFILETYPE,
     .              RP%LINE%IFREMD, RP%LINE%NRJPRT, RP%LINE%IPLSC6,
     .              RP%LINE%IMESS
          READ (13+IFOFF) RP%LINE%REACNAME, RP%LINE%KENN
 
        ELSE IF (RP%IFIT <= 2) THEN
! DATA FOR FIT EXPRESSIONS  (E.G. POLYNOMIAL, IN CASE OF HYDHEL DATABASE)
          IF (.NOT.ASSOCIATED(RP%POLY)) ALLOCATE (RP%POLY)
          IF (ASSOCIATED(RP%POLY%DBLPOL)) DEALLOCATE (RP%POLY%DBLPOL)
 
          READ (13+IFOFF) RP%POLY%IFEXMN,RP%POLY%IFEXMX,ND,ND2
          ALLOCATE (RP%POLY%DBLPOL(ND,ND2))
          READ (13+IFOFF) RP%POLY%DBLPOL,
     .              RP%POLY%RCMN,RP%POLY%RCMX,RP%POLY%FPARM
 
        ELSE
! DATA FROM FORMATTED TABLES (E.G. ADAS)
        END IF
 
        END SUBROUTINE EIRENE_READ_FIT_FORM
 
      END SUBROUTINE EIRENE_READ_CMAMF
 
 
      SUBROUTINE EIRENE_CMAMF_XDR (IUN,IFLG)
 
      INTEGER, INTENT(IN) :: IUN,IFLG
      INTEGER :: IR, IHELP(1)
      LOGICAL :: LHELP(7)
 
      CALL FXDRDBL (IUN,DELPOT,NREAC)
      CALL FXDRDBL (IUN,FACREA,(NREAC+11)*2)
      CALL FXDRDBL (IUN,FREACA,NATM*NREAC)
      CALL FXDRDBL (IUN,FREACM,NMOL*NREAC)
      CALL FXDRDBL (IUN,FREACI,NION*NREAC)
      CALL FXDRDBL (IUN,FREACP,NPLS*NREAC)
      CALL FXDRDBL (IUN,FREACPH,NPHOT*NREAC)
      CALL FXDRDBL (IUN,FLDLMA,NATM*NREAC)
      CALL FXDRDBL (IUN,FLDLMM,NMOL*NREAC)
      CALL FXDRDBL (IUN,FLDLMI,NION*NREAC)
      CALL FXDRDBL (IUN,FLDLMP,NPLS*NREAC)
      CALL FXDRDBL (IUN,FLDLMPH,NPHOT*NREAC)
      CALL FXDRDBL (IUN,EELECA,NATM*NREAC)
      CALL FXDRDBL (IUN,EELECM,NMOL*NREAC)
      CALL FXDRDBL (IUN,EELECI,NION*NREAC)
      CALL FXDRDBL (IUN,EELECP,NPLS*NREAC)
      CALL FXDRDBL (IUN,EELECPH,NPHOT*NREAC)
      CALL FXDRDBL (IUN,EBULKA,NATM*NREAC)
      CALL FXDRDBL (IUN,EBULKM,NMOL*NREAC)
      CALL FXDRDBL (IUN,EBULKI,NION*NREAC)
      CALL FXDRDBL (IUN,EBULKP,NPLS*NREAC)
      CALL FXDRDBL (IUN,EBULKPH,NPHOT*NREAC)
      CALL FXDRDBL (IUN,ESCD1A,NATM*NREAC)
      CALL FXDRDBL (IUN,ESCD1M,NMOL*NREAC)
      CALL FXDRDBL (IUN,ESCD1I,NION*NREAC)
      CALL FXDRDBL (IUN,ESCD1P,NPLS*NREAC)
      CALL FXDRDBL (IUN,ESCD1PH,NPHOT*NREAC)
 
      IHELP(1) = NREACI
      CALL FXDRINT (IUN,IHELP,1)
      CALL FXDRINT (IUN,ISWR,NREAC)
      CALL FXDRINT (IUN,MODCLF,NREAC)
      CALL FXDRINT (IUN,MASSP,NREAC)
      CALL FXDRINT (IUN,MASST,NREAC)
      CALL FXDRINT (IUN,IFTFLG,6*NREAC)
      CALL FXDRINT (IUN,NRCP,NPLS)
      CALL FXDRINT (IUN,NRCA,NATM)
      CALL FXDRINT (IUN,NRCM,NMOL)
      CALL FXDRINT (IUN,NRCI,NION)
      CALL FXDRINT (IUN,NRCPH,NPHOT)
      CALL FXDRINT (IUN,IREACA,NATM*NREAC)
      CALL FXDRINT (IUN,IREACM,NMOL*NREAC)
      CALL FXDRINT (IUN,IREACI,NION*NREAC)
      CALL FXDRINT (IUN,IREACP,NPLS*NREAC)
      CALL FXDRINT (IUN,IREACPH,NPHOT*NREAC)
      CALL FXDRINT (IUN,IBULKA,NATM*NREAC)
      CALL FXDRINT (IUN,IBULKM,NMOL*NREAC)
      CALL FXDRINT (IUN,IBULKI,NION*NREAC)
      CALL FXDRINT (IUN,IBULKP,NPLS*NREAC)
      CALL FXDRINT (IUN,IBULKPH,NPHOT*NREAC)
      CALL FXDRINT (IUN,ISCD1A,NATM*NREAC)
      CALL FXDRINT (IUN,ISCD1M,NMOL*NREAC)
      CALL FXDRINT (IUN,ISCD1I,NION*NREAC)
      CALL FXDRINT (IUN,ISCD1P,NPLS*NREAC)
      CALL FXDRINT (IUN,ISCD1PH,NPHOT*NREAC)
      CALL FXDRINT (IUN,ISCD2A,NATM*NREAC)
      CALL FXDRINT (IUN,ISCD2M,NMOL*NREAC)
      CALL FXDRINT (IUN,ISCD2I,NION*NREAC)
      CALL FXDRINT (IUN,ISCD2P,NPLS*NREAC)
      CALL FXDRINT (IUN,ISCD2PH,NPHOT*NREAC)
      CALL FXDRINT (IUN,ISCD3A,NATM*NREAC)
      CALL FXDRINT (IUN,ISCD3M,NMOL*NREAC)
      CALL FXDRINT (IUN,ISCD3I,NION*NREAC)
      CALL FXDRINT (IUN,ISCD3P,NPLS*NREAC)
      CALL FXDRINT (IUN,ISCD3PH,NPHOT*NREAC)
      CALL FXDRINT (IUN,ISCD4A,NATM*NREAC)
      CALL FXDRINT (IUN,ISCD4M,NMOL*NREAC)
      CALL FXDRINT (IUN,ISCD4I,NION*NREAC)
      CALL FXDRINT (IUN,ISCD4P,NPLS*NREAC)
      CALL FXDRINT (IUN,ISCD4PH,NPHOT*NREAC)
      CALL FXDRINT (IUN,ISCDEA,NATM*NREAC)
      CALL FXDRINT (IUN,ISCDEM,NMOL*NREAC)
      CALL FXDRINT (IUN,ISCDEI,NION*NREAC)
      CALL FXDRINT (IUN,ISCDEP,NPLS*NREAC)
      CALL FXDRINT (IUN,ISCDEPH,NPHOT*NREAC)
      CALL FXDRINT (IUN,IESTMA,NATM*NREAC)
      CALL FXDRINT (IUN,IESTMM,NMOL*NREAC)
      CALL FXDRINT (IUN,IESTMI,NION*NREAC)
      CALL FXDRINT (IUN,IESTMPH,NPHOT*NREAC)
      CALL FXDRINT (IUN,IBGKA ,NATM*NREAC)
      CALL FXDRINT (IUN,IBGKM ,NMOL*NREAC)
      CALL FXDRINT (IUN,IBGKI, NION*NREAC)
      CALL FXDRINT (IUN,IBGKPH,NPHOT*NREAC)
 
      IF (IFLG == 0) THEN
        DO IR=1,NREACI
          LHELP = (/ REACDAT(IR)%LPOT, REACDAT(IR)%LCRS,
     .               REACDAT(IR)%LRTC, REACDAT(IR)%LRTCMW,
     .               REACDAT(IR)%LRTCEW, REACDAT(IR)%LOTH,
     .               REACDAT(IR)%LPHR /)
          CALL FXDRLOG(IUN,LHELP,7)
          IHELP(1) = REACDAT(IR)%NOSEC
          CALL FXDRINT (IUN,IHELP,1)
          IF (REACDAT(IR)%LPOT) CALL
     .  EIRENE_WXDR_FIT_FORM(REACDAT(IR)%POT)
          IF (REACDAT(IR)%LCRS) CALL
     .  EIRENE_WXDR_FIT_FORM(REACDAT(IR)%CRS)
          IF (REACDAT(IR)%LRTC) CALL
     .  EIRENE_WXDR_FIT_FORM(REACDAT(IR)%RTC)
          IF (REACDAT(IR)%LRTCMW) CALL
     .  EIRENE_WXDR_FIT_FORM(REACDAT(IR)%RTCMW)
          IF (REACDAT(IR)%LRTCEW) CALL
     .  EIRENE_WXDR_FIT_FORM(REACDAT(IR)%RTCEW)
          IF (REACDAT(IR)%LOTH) CALL
     .  EIRENE_WXDR_FIT_FORM(REACDAT(IR)%OTH)
          IF (REACDAT(IR)%LPHR) CALL
     .  EIRENE_WXDR_FIT_FORM(REACDAT(IR)%PHR)
        END DO
      ELSE
        LHELP = .FALSE.
        DO IR=1,NREACI
          CALL FXDRLOG(IUN,LHELP,7)
          REACDAT(IR)%LPOT = LHELP(1)
          REACDAT(IR)%LCRS = LHELP(2)
          REACDAT(IR)%LRTC = LHELP(3)
          REACDAT(IR)%LRTCMW = LHELP(4)
          REACDAT(IR)%LRTCEW = LHELP(5)
          REACDAT(IR)%LOTH = LHELP(6)
          REACDAT(IR)%LPHR = LHELP(7)
          CALL FXDRINT (IUN,IHELP,1)
          REACDAT(IR)%NOSEC = IHELP(1)
          IF (REACDAT(IR)%LPOT) CALL
     .  EIRENE_RXDR_FIT_FORM(REACDAT(IR)%POT)
          IF (REACDAT(IR)%LCRS) CALL
     .  EIRENE_RXDR_FIT_FORM(REACDAT(IR)%CRS)
          IF (REACDAT(IR)%LRTC) CALL
     .  EIRENE_RXDR_FIT_FORM(REACDAT(IR)%RTC)
          IF (REACDAT(IR)%LRTCMW) CALL
     .  EIRENE_RXDR_FIT_FORM(REACDAT(IR)%RTCMW)
          IF (REACDAT(IR)%LRTCEW) CALL
     .  EIRENE_RXDR_FIT_FORM(REACDAT(IR)%RTCEW)
          IF (REACDAT(IR)%LOTH) CALL
     .  EIRENE_RXDR_FIT_FORM(REACDAT(IR)%OTH)
          IF (REACDAT(IR)%LPHR) CALL
     .  EIRENE_RXDR_FIT_FORM(REACDAT(IR)%PHR)
        END DO
      END IF
 
      RETURN
 
      CONTAINS
 
        SUBROUTINE EIRENE_WXDR_FIT_FORM(RP)
        TYPE(FIT_FORMS),POINTER :: RP
        INTEGER :: IHELP(18), ND, ND2
        REAL(DP) :: RHELP(23)
        CHARACTER(24) :: CKENN
 
        IHELP(1) = RP%IFIT
        CALL FXDRINT (IUN,IHELP,1)
 
        IF (RP%IFIT < 0) THEN
! DATA FOR PHOTONIC LINE
          RHELP(1:11) = (/ RP%LINE%E0, RP%LINE%E1, RP%LINE%AIK,
     .               RP%LINE%G1, RP%LINE%G2, RP%LINE%C2,
     .               RP%LINE%C3, RP%LINE%C4, RP%LINE%C6,
     .               RP%LINE%B12, RP%LINE%B21 /)
          RHELP(12:23) = RP%LINE%C6A(1:12)
 
          IHELP(1:6) = (/ RP%LINE%IGND, RP%LINE%IRCART,
     .                    RP%LINE%IPROFILETYPE,
     .                    RP%LINE%IFREMD, RP%LINE%NRJPRT,
     .                    RP%LINE%IMESS /)
          IHELP(7:18) = RP%LINE%IPLSC6
          CKENN = TRANSFER(RP%LINE%KENN,CKENN)
 
          CALL FXDRDBL (IUN,RHELP,23)
          CALL FXDRINT (IUN,IHELP,18)
          CALL FXDRCHR (IUN,RP%LINE%REACNAME)
          CALL FXDRCHR (IUN,CKENN)
 
        ELSE IF (RP%IFIT <= 2) THEN
! DATA FOR FIT EXPRESSIONS  (E.G. POLYNOMIAL, IN CASE OF HYDHEL DATABASE)
          ND = UBOUND(RP%POLY%DBLPOL,1)
          ND2 = UBOUND(RP%POLY%DBLPOL,2)
          IHELP(1:4) = (/ RP%POLY%IFEXMN, RP%POLY%IFEXMX, ND, ND2 /)
 
          RHELP(1) = RP%POLY%RCMN
          RHELP(2) = RP%POLY%RCMX
          RHELP(3:8) = RP%POLY%FPARM
 
          CALL FXDRINT (IUN,IHELP,4)
          CALL FXDRDBL (IUN,RP%POLY%DBLPOL,ND*ND2)
          CALL FXDRDBL (IUN,RHELP,8)
 
        ELSE
! DATA FROM FORMATTED TABLES (E.G. ADAS)
        END IF
 
        END SUBROUTINE EIRENE_WXDR_FIT_FORM
 
 
        SUBROUTINE EIRENE_RXDR_FIT_FORM(RP)
        TYPE(FIT_FORMS),POINTER :: RP
        INTEGER :: IHELP(18), ND, ND2
        REAL(DP) :: RHELP(23)
        CHARACTER(24) :: CKENN
 
        CALL FXDRINT (IUN,IHELP,1)
        RP%IFIT = IHELP(1)
 
        IF (RP%IFIT < 0) THEN
! DATA FOR PHOTONIC LINE
          IF (.NOT.ASSOCIATED(RP%LINE)) ALLOCATE (RP%LINE)
          CALL FXDRDBL (IUN,RHELP,23)
          CALL FXDRINT (IUN,IHELP,18)
          CALL FXDRCHR (IUN,RP%LINE%REACNAME)
          CALL FXDRCHR (IUN,CKENN)
 
          RP%LINE%KENN = TRANSFER(CKENN,RP%LINE%KENN)
 
          RP%LINE%E0 = RHELP(1)
          RP%LINE%E1 = RHELP(2)
          RP%LINE%AIK = RHELP(3)
          RP%LINE%G1 = RHELP(4)
          RP%LINE%G2 = RHELP(5)
          RP%LINE%C2 = RHELP(6)
          RP%LINE%C3 = RHELP(7)
          RP%LINE%C4 = RHELP(8)
          RP%LINE%C6 = RHELP(9)
          RP%LINE%B12 = RHELP(10)
          RP%LINE%B21 = RHELP(11)
          RP%LINE%C6A(1:12) = RHELP(12:23)
 
          RP%LINE%IGND  = IHELP(1)
          RP%LINE%IRCART = IHELP(2)
          RP%LINE%IPROFILETYPE = IHELP(3)
          RP%LINE%IFREMD = IHELP(4)
          RP%LINE%NRJPRT = IHELP(5)
          RP%LINE%IMESS = IHELP(6)
          RP%LINE%IPLSC6(1:12) = IHELP(7:18)
 
 
        ELSE IF (RP%IFIT <= 2) THEN
! DATA FOR POLYNOMIAL FIT
 
          CALL FXDRINT (IUN,IHELP,4)
          RP%POLY%IFEXMN = IHELP(1)
          RP%POLY%IFEXMX = IHELP(2)
          ND = IHELP(3)
          ND2 = IHELP(4)
 
          IF (.NOT.ASSOCIATED(RP%POLY)) ALLOCATE (RP%POLY)
          IF (ASSOCIATED(RP%POLY%DBLPOL)) DEALLOCATE (RP%POLY%DBLPOL)
          ALLOCATE (RP%POLY%DBLPOL(ND,ND2))
          CALL FXDRDBL (IUN,RP%POLY%DBLPOL,ND*ND2)
 
          CALL FXDRDBL (IUN,RHELP,8)
          RP%POLY%RCMN = RHELP(1)
          RP%POLY%RCMX = RHELP(2)
          RP%POLY%FPARM = RHELP(3:8)
 
        ELSE
! DATA FROM FORMATTED TABLES (E.G. ADAS)
        END IF
 
        END SUBROUTINE EIRENE_RXDR_FIT_FORM
 
      END SUBROUTINE EIRENE_CMAMF_XDR
 
 
      SUBROUTINE EIRENE_GET_REACTION (IR)
 
      INTEGER, INTENT(IN) :: IR
 
      REACTION => REACDAT(IR)%PHR%LINE
 
      IDREAC = IR
 
      RETURN
      END SUBROUTINE EIRENE_GET_REACTION
 
 
      SUBROUTINE EIRENE_SET_REACTION_DATA
     .           (IR,ISW,IFTFLG,RDATA,IUNOUT,LTEST,
     .            RCMIN, RCMAX, FP, JFEXMN, JFEXMX)
 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IR, ISW, IFTFLG, IUNOUT
      INTEGER, OPTIONAL, INTENT(IN) :: JFEXMN, JFEXMX
      REAL(DP), INTENT(IN) :: RDATA(9,*)
      REAL(DP), OPTIONAL, INTENT(IN) :: RCMIN, RCMAX, FP(6)
      LOGICAL, INTENT(IN) :: LTEST
      INTEGER :: NDIM, NDIM2, I, J, IFIT
      REAL(DP) :: CTEST
      TYPE(POLY_DATA), POINTER :: REA
      INTEGER, SAVE :: ISW2D(7) = (/ 3, 4, 6, 7, 9, 10, 12 /)
 
      NDIM = 9
      IF (MOD(IFTFLG,100) == 10) NDIM = 1
 
      NDIM2 = 1
      IF (COUNT(ISW2D == ISW) > 0) NDIM2=9
 
      IF (LTEST) THEN
!pb        DO J = 1, NDIM2
!pb          CTEST = SUM(ABS(RDATA(1:NDIM,J)))
!pb          IF (CTEST.LE.1.E-30_DP) THEN
!pb            WRITE (iunout,*)
!pb     .            'ERROR IN SUBROUTINE EIRENE_SET_REACTION_DATA:',
!pb     .            ' ZERO FIT COEFFICIENTS'
!pb            WRITE (iunout,*) 'J,IR = ',J,IR,'  EXIT CALLED!'
!pb            CALL EIRENE_EXIT_OWN(1)
!pb          END IF
!pb        END DO
        CTEST = SUM(ABS(RDATA(1:NDIM,1:NDIM2)))
        IF (CTEST.LE.1.E-30_DP) THEN
           WRITE (iunout,*)
     .            'ERROR IN SUBROUTINE EIRENE_SET_REACTION_DATA:',
     .            ' ZERO FIT COEFFICIENTS'
           WRITE (iunout,*) 'IR = ',IR,'  EXIT CALLED!'
           CALL EIRENE_EXIT_OWN(1)
        END IF
      END IF
 
      ALLOCATE (REA)
      ALLOCATE (REA%DBLPOL(1:NDIM,1:NDIM2))
 
      REA%RCMN = -HUGE(1._DP)
      REA%RCMX =  HUGE(1._DP)
      REA%FPARM = 0._DP
      REA%IFEXMN = 0
      REA%IFEXMX = 0
 
      REA%DBLPOL(1:NDIM,1:NDIM2) = RDATA(1:NDIM,1:NDIM2)
      IFIT = 1
      IF (NDIM2 == 9) IFIT = 2
 
      IF (PRESENT(RCMIN)) REA%RCMN = RCMIN
      IF (PRESENT(RCMAX)) REA%RCMX = RCMAX
      IF (PRESENT(FP)) REA%FPARM = FP
      IF (PRESENT(JFEXMN)) REA%IFEXMN = JFEXMN
      IF (PRESENT(JFEXMX)) REA%IFEXMX = JFEXMX
 
      SELECT CASE(ISW)
 
      CASE (0)
        IF (REACDAT(IR)%LPOT) THEN
          WRITE (IUNOUT,*) ' POTENTIAL ALREADY SPECIFIED FOR REACTION',
     .                       IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        ALLOCATE (REACDAT(IR)%POT)
        NULLIFY(REACDAT(IR)%POT%ADAS)
        NULLIFY(REACDAT(IR)%POT%LINE)
        NULLIFY(REACDAT(IR)%POT%HYD)
        REACDAT(IR)%POT%POLY => REA
        REACDAT(IR)%LPOT = .TRUE.
        REACDAT(IR)%POT%IFIT = IFIT
 
      CASE (1)
        IF (REACDAT(IR)%LCRS) THEN
          WRITE (IUNOUT,*) ' CROSS SECTION ALREADY SPECIFIED',
     .                     ' FOR REACTION', IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        ALLOCATE (REACDAT(IR)%CRS)
        NULLIFY(REACDAT(IR)%CRS%ADAS)
        NULLIFY(REACDAT(IR)%CRS%LINE)
        NULLIFY(REACDAT(IR)%CRS%HYD)
        REACDAT(IR)%CRS%POLY => REA
        REACDAT(IR)%LCRS = .TRUE.
        REACDAT(IR)%CRS%IFIT = IFIT
 
      CASE (2:4)
        IF (REACDAT(IR)%LRTC) THEN
          WRITE (IUNOUT,*) ' RATE COEFFICIENT ALREADY SPECIFIED',
     .                     ' FOR REACTION', IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        ALLOCATE (REACDAT(IR)%RTC)
        NULLIFY(REACDAT(IR)%RTC%ADAS)
        NULLIFY(REACDAT(IR)%RTC%LINE)
        NULLIFY(REACDAT(IR)%RTC%HYD)
        REACDAT(IR)%RTC%POLY => REA
        REACDAT(IR)%LRTC = .TRUE.
        REACDAT(IR)%RTC%IFIT = IFIT
 
      CASE (5:7)
        IF (REACDAT(IR)%LRTCMW) THEN
          WRITE (IUNOUT,*) ' MOMEMTUM WEIGHTED RATE COEFFICIENT',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        ALLOCATE (REACDAT(IR)%RTCMW)
        NULLIFY(REACDAT(IR)%RTCMW%ADAS)
        NULLIFY(REACDAT(IR)%RTCMW%LINE)
        NULLIFY(REACDAT(IR)%RTCMW%HYD)
        REACDAT(IR)%RTCMW%POLY => REA
        REACDAT(IR)%LRTCMW = .TRUE.
        REACDAT(IR)%RTCMW%IFIT = IFIT
 
      CASE (8:10)
        IF (REACDAT(IR)%LRTCEW) THEN
          WRITE (IUNOUT,*) ' ENERGY WEIGHTED RATE COEFFICIENT',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        ALLOCATE (REACDAT(IR)%RTCEW)
        NULLIFY(REACDAT(IR)%RTCEW%ADAS)
        NULLIFY(REACDAT(IR)%RTCEW%LINE)
        NULLIFY(REACDAT(IR)%RTCEW%HYD)
        REACDAT(IR)%RTCEW%POLY => REA
        REACDAT(IR)%LRTCEW = .TRUE.
        REACDAT(IR)%RTCEW%IFIT = IFIT
 
      CASE (11:12)
        IF (REACDAT(IR)%LOTH) THEN
          WRITE (IUNOUT,*) ' OTHER POLYNOMIAL FIT COEFFICIENTS',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        ALLOCATE (REACDAT(IR)%OTH)
        NULLIFY(REACDAT(IR)%OTH%ADAS)
        NULLIFY(REACDAT(IR)%OTH%LINE)
        NULLIFY(REACDAT(IR)%OTH%HYD)
        REACDAT(IR)%OTH%POLY => REA
        REACDAT(IR)%LOTH = .TRUE.
        REACDAT(IR)%OTH%IFIT = IFIT
 
      CASE DEFAULT
        WRITE (IUNOUT,*) ' WRONG REACTION TYPE SPCIFIED '
        WRITE (IUNOUT,*) ' REACTION NO. ', IR
        WRITE (IUNOUT,*) ' REACTION TYPE H.', ISW
        CALL EIRENE_EXIT_OWN(1)
      END SELECT
 
      RETURN
      END SUBROUTINE EIRENE_SET_REACTION_DATA
 
 
      FUNCTION EIRENE_IS_RTC_ADAS (IREAC) RESULT(RES)
 
      INTEGER, INTENT(IN) :: IREAC
      LOGICAL :: RES
 
      RES = .FALSE.
      IF (REACDAT(IREAC)%LRTC) THEN
        RES = REACDAT(IREAC)%RTC%IFIT == 3
      END IF
 
      END FUNCTION EIRENE_IS_RTC_ADAS
 
 
      FUNCTION EIRENE_IS_RTCEW_ADAS (IREAC) RESULT(RES)
 
      INTEGER, INTENT(IN) :: IREAC
      LOGICAL :: RES
 
      RES = .FALSE.
      IF (REACDAT(IREAC)%LRTCEW) THEN
        RES = REACDAT(IREAC)%RTCEW%IFIT == 3
      END IF
 
      END FUNCTION EIRENE_IS_RTCEW_ADAS
 
 
      FUNCTION EIRENE_IS_RTCMW_ADAS (IREAC) RESULT(RES)
 
      INTEGER, INTENT(IN) :: IREAC
      LOGICAL :: RES
 
      RES = .FALSE.
      IF (REACDAT(IREAC)%LRTCMW) THEN
        RES = REACDAT(IREAC)%RTCMW%IFIT == 3
      END IF
 
      END FUNCTION EIRENE_IS_RTCMW_ADAS
  

      SUBROUTINE EIRENE_FREE_REACDAT

      type(fit_forms), pointer :: rea
      integer :: ir

      DO IR = -11, NREAC

        IF (REACDAT(IR)%LPOT) THEN
           rea => REACDAT(IR)%POT
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF

        IF (REACDAT(IR)%LCRS) THEN
           rea => REACDAT(IR)%CRS
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF

        IF (REACDAT(IR)%LRTC) THEN
           rea => REACDAT(IR)%RTC
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF

        IF (REACDAT(IR)%LRTCMW) THEN
           rea => REACDAT(IR)%RTCMW
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF

        IF (REACDAT(IR)%LRTCEW) THEN
           rea => REACDAT(IR)%RTCEW
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF

        IF (REACDAT(IR)%LOTH) THEN
           rea => REACDAT(IR)%OTH
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF

        IF (REACDAT(IR)%LPHR) THEN
           rea => REACDAT(IR)%PHR
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF

      END DO

      deallocate (reacdat)

      END SUBROUTINE EIRENE_FREE_REACDAT



      SUBROUTINE EIRENE_FREE_FIT_FORM(RP)

      TYPE(FIT_FORMS),POINTER :: RP

      IF (ASSOCIATED(RP%POLY)) THEN
         DEALLOCATE (RP%POLY%DBLPOL)
         DEALLOCATE (RP%POLY)
      END IF

      IF (ASSOCIATED(RP%ADAS)) THEN
         DEALLOCATE (RP%ADAS%DENS)
         DEALLOCATE (RP%ADAS%TEMP)
         DEALLOCATE (RP%ADAS%FIT)
         DEALLOCATE (RP%ADAS%DDE)
         DEALLOCATE (RP%ADAS%DTE)
         DEALLOCATE (RP%ADAS)
      END IF

      IF (ASSOCIATED(RP%LINE)) THEN
         DEALLOCATE (RP%LINE)
      END IF

      IF (ASSOCIATED(RP%HYD)) THEN
         DEALLOCATE (RP%HYD%TEMPS)
         DEALLOCATE (RP%HYD%RATES)
         DEALLOCATE (RP%HYD%RATIO)
         DEALLOCATE (RP%HYD)
      END IF

      END SUBROUTINE EIRENE_FREE_FIT_FORM


      END MODULE EIRMOD_COMXS
