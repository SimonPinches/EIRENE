
!pb  18.12.06:  NPARTC and NPARTT reduced because of cancelation of XNUE
!    20.06.07:  NUM_PARM = maximum number of parameters introduced
cdr  input tallies   ntali, increased from 21 to 22 (electr. potential)
cdr  surface tallies ntals, increased from 59 to 79 (more sputter tallies)
cpb  surface tallies ntals, increased from 79 to 84 (even more sputter tallies)

cdr  naming conventions for variance tallies also for spectra tallies
cdr  spcint --> spcs
cdr 21.09.15:  NPARTT REDUCED FROM 12 TO 11 (XGENER NOT ON CENSUS)
cdr  Dec. 15:  species-resolved energy tallies for pl (bulk ion) energy balance.
!pb  May  16:  nrds -> nrei
cdr  May  17: eliminate NCOP, NCOPI, only use NCPV, NCPVI
cdr           tbd: similar: eliminate NBGK, NBGKI,  only use  NBGV, NBGVI
cdr  July 17: remove NTALW  (was same as NTALS), NAIN added to N1MX
cpb  Dec. 17: remove type SPECT_ARRAY, not needed in Fortran 2003
cdr   dec.17: add nspztotw, at same place as formerly NTALW was.
cdr           fully corresponds to vol tally parameter nspztot,
cdr           but is for surface tally pointers
cpb  input tallies   ntali, increased from 22 to 24 (BVIN, PARMOM)
cdr  jan.18:  added: NUM_LINES, NADV_ADD
cdr  nov.18:  notational cleanup: separate OT from PH (photonic) processes
cpb  dez.18:  ntalg, increased from 24 to 30, free slots for future use
c
      MODULE EIRMOD_PARMMOD
c
c handling of storage for dynamically allocated arrays
c contains:
c    set_parmmod(ical),  ical=1,2,3
c    collect_parm
c    distrib_parm

      USE EIRMOD_PRECISION

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_SET_PARMMOD, EIRENE_COLLECT_PARM,
     P          EIRENE_DISTRIB_PARM,
     P          EIRENE_SPECTRUM,
     P          ASSIGNMENT(=)

      INTEGER, PUBLIC, PARAMETER ::
     P         NUM_PARM=200,
     P         NPARTC=12, NPARTT=11,
     P         MPARTC=14, MPARTT=10
C> Unit number for memory usage output file
      INTEGER, PUBLIC, SAVE :: IUNMEM = 55
C> Unit number for RAPS output file
      INTEGER, PUBLIC, SAVE :: IUNRAPS = 65
C> Unit number for RAPS vector field output file
      INTEGER, PUBLIC, SAVE :: IUNRAPSVEC = 60
      integer, public, save :: IFOFF = 0
C> Indicates whether output files 'output.*' should be appended or overwritten
      LOGICAL, PUBLIC, SAVE :: LOUTAPP = .FALSE.

      INTEGER, PUBLIC, SAVE ::
     I N1ST,   N2ND,   N3RD,   NADD,   NTOR,
     I NRTAL,  NLIM,   NSTS,
     I NPLG,   NPPART, NKNOT,  NTRI,   NTETRA, NCOORD,
     I NOPTIM, NOPTM1, NCORNER

      INTEGER, PUBLIC, SAVE ::
     I NSTRA,  NSRFS,  NSTEP

      INTEGER, PUBLIC, SAVE ::
     I NATM,   NMOL,   NION,   NPLS,   NPHOT,  NADV,   NADS,
     I NCLV,   NSNV,   NALV,   NALS,   NAIN,   NCPV,   NBGK,
     I NPLSTI, NPLSV,
     I NADSPC, NBACK_SPEC ,NADSPC_S, NADSPC_C, NADSPC_D,NADSPC_CD,
     I NADV_ADD

      INTEGER, PUBLIC, SAVE ::
     I NSD,    NSDW,   NCV

      INTEGER, PUBLIC, SAVE ::
     I NREAC,  NREC,   NREI,   NRCX,   NREL,   NRPI,   NRPH

      INTEGER, PUBLIC, SAVE ::
     I NHD1,   NHD2,   NHD3,   NHD4,   NHD5,   NHD6

      INTEGER, PUBLIC, SAVE ::
     I NCHOR,  NCHEN, NUM_LINES

      INTEGER, PUBLIC, SAVE ::
     I NDX,    NDY,    NFL,    NDXP,   NDYP,   NPTRGT

      INTEGER, PUBLIC, SAVE ::
     I NPRNL

      INTEGER, PUBLIC, SAVE ::
     I NTRJ


      INTEGER, PUBLIC, SAVE ::
     I NGEOM_USR, NCOUP_INPUT, NSMSTRA, NSTORAM, NGSTAL

      INTEGER, PUBLIC, SAVE ::
     I NRAD,   NSWIT,   N1F,    N2F,    N3F,    NGITT,  NGITTP,
     I NRADS,  N2NDPLGS, N1STS,  N2NDS,  NTRIS,  NKNOTS, NRTALS

      INTEGER, PUBLIC, SAVE ::
     I NGTSFT, NLIMPS, NLMPGS

      INTEGER, PUBLIC, SAVE ::
     I NBGV,   NBMAX,   NPTAL,  NCPV_STAT, NSCOP

      INTEGER, PUBLIC, SAVE ::
     I NSTRAP

      INTEGER, PUBLIC, SAVE ::
     I NIONP,  NATMP,  NMOLP,
     I NPLSP,  NPHOTP, NADVP,  NADSP,
     I NCLVP,  NALVP,  NALSP,
     I NSNVP,  NCPVP,  NBGVP,
     I NTALI,  NTALG,  NTALN,  NTALO,  NTALV,
     I NTALA,  NTALC,  NTALT,
     I NTALM,  NTALB,  NTALR,
     I NTALS,  NTLSA,  NTLSR,  NSPZTOTW,
     I N1MX,   N2MX,   NSPZ,   NSPZP, NSPZMC, NCOLMC, NSPZTOT

      INTEGER, PUBLIC, SAVE ::
     I NVOLTL, NVLTLP,
     I NSRFTL, NSFTLP,
     I NINPTL

      INTEGER, PUBLIC, SAVE ::
     I NH0,    NH1,    NH2,    NH3

      INTEGER, PUBLIC, SAVE ::
     I NHSTOR, NSTORDT, NSTORDR

      INTEGER, PUBLIC, SAVE ::
     I NPLT,   NVLPR,  NSRPR

      INTEGER, PUBLIC, SAVE :: INT_PARM(NUM_PARM)



      PRIVATE :: EIRENE_SPEC_TO_SPEC
      TYPE EIRENE_SPECTRUM
        REAL(DP) :: SPCMIN, SPCMAX, SPCDEL, SPCDELI, ESP_MIN,
     .              ESP_MAX, ESP_00, SPC_XPLT, SPC_YPLT, SPC_SAME,
     .              SPCVX, SPCVY, SPCVZ
        REAL(DP) :: SPCS, SGMS, STVS, GGS
        INTEGER :: NSPC, ISPCTYP, ISPCSRF, IPRTYP, IPRSP, IMETSP,
     .             ISRFCLL, IDIREC
        LOGICAL :: LOG
        REAL(DP), DIMENSION(:), POINTER :: SPC, SDV, SGM, STV, GG
      END TYPE EIRENE_SPECTRUM


      INTERFACE ASSIGNMENT(=)  ! DEFINE ASSIGNMENT
        MODULE PROCEDURE EIRENE_SPEC_TO_SPEC
      END INTERFACE



      CONTAINS


      SUBROUTINE EIRENE_SET_PARMMOD(ICAL)
C  ical=1:  called directly from eirene.f after "find_param.f", before input.f
C  ical=2:  called from inside "setamd.f",  prepare allocatable storage for comxs, comsou, czt1
C  ical=3:  called from inside "input.f", prepare allocatable storage for cgeom, comusr

      INTEGER, INTENT(IN) :: ICAL


      IF (ICAL == 1) THEN
C.......................................................................
C  CALLED AFTER FIND_PARAM.F, AND BEFORE INPUT.F
C.......................................................................
C
C  GEOMETRY
C
        NRAD=MAX(N1ST*N2ND*N3RD,NTRI*N3RD,NTETRA)+NADD+1
        IF (NRTAL==0) NRTAL=NRAD
        IF (NOPTIM < 0) NOPTIM = NRAD

C  NSWIT: ELIMINATE SOME ARRAYS IN CASE OF LEVGEO=10 OPTION
C                     (GEOMETRY ARRAYS OUTSIDE EIRENE-CODE)
C  ngeom_usr=1:  use eirene grid tallies RSURF, PSURF,...  (default)
c  ngeom_usr=0:  eirene grid tallies are eliminated, no storage, (e.g. in case of levgeo=10)
        NSWIT=1-NGEOM_USR

C  IDENTIFY: WHICH GRIDS ARE THERE? N1F=0 OR N1F=1, IF N1ST=1, OR IF N1ST GT 1, RESP.
        N1F=1-1/N1ST
        N2F=1-1/N2ND
        N3F=1-1/N3RD
C
c  ngitt is the largest possible 2d dimension (grid on a coordinate surface or line)
        IF (NGITT <= 1) NGITT=N1ST*N2ND*N3F+N1ST*N3RD*N2F+N2ND*N3RD*N1F
        NGITTP=NGITT+1
C
C STORAGE FOR GRIDS. SWITCH OFF GRID STORAGE IN CASE OF LEVGEO=10 (external geometry package)
        NRADS=NSWIT*    NRAD+                  (1-NSWIT)*1
        NRTALS=NSWIT*   NRTAL+                 (1-NSWIT)*1
        N2NDPLGS=NSWIT*(N2ND*N2F+NPLG*(1-N2F))+(1-NSWIT)*1
        N1STS=NSWIT*    N1ST+                  (1-NSWIT)*1
        N2NDS=NSWIT*    N2ND+                  (1-NSWIT)*1
        NTRIS=NSWIT*    NTRI+                  (1-NSWIT)*1
        NKNOTS=NSWIT*   NKNOT+                 (1-NSWIT)*1
C
C TALLIES
C
        NLIMPS=NLIM+NSTS
C
C  GENERATION LIMIT TALLIES
C
        NBMAX=10
        NPTAL=30

C  PRIMARY SOURCE
        NSTRAP=NSTRA+1

C  SPECIES AND TALLIES  NTALV: TOTAL NUMBER OF VOLUME OUTPUT TALLIES
C                           NTALA: INDEX OF THE ADDITIONAL TRACKLENGTH
C                                  ESTIMATED TALLY
C                           NTALC: INDEX OF THE ADDITIONAL COLLISION
C                                  ESTIMATED TALLY
C                           NTALT: INDEX OF THE TIME DEP. TALLY
C                                  (SNAPSHOT ESTIMATOR)
C                           NTALM: INDEX OF THE TALLIES FOR COUPLING,
C                                  (E.G. MOMENTUM SOURCES)
C                           NTALB: INDEX OF THE BGK TALLY
C                           NTALR: INDEX OF THE ALGEBRAIC TALLY

C                       NTALS: TOTAL NUMBER OF SURFACE OUTPUT TALLIES
C                           NTLSA: INDEX OF THE ADDITIONAL TALLY
C                                  (TRACKLENGTH AND COLLISION ESTIMATORS
C                                   ARE IDENTICAL FOR SURFACE AVERAGES)
C                           NTLSR: INDEX OF THE ALGEBRAIC TALLY

C                       NTALI: TOTAL NUMBER OF INPUT TALLIES
C                           NTALG: NUMBER OF INPUT TALLIES, EXCLUDING THE OPT. GRADIENT TALLIES
C                           NTALN: INDEX OF THE ADDITIONAL INPUT TALLIES
C                           NTALO: INDEX OF THE CELL VOLUME TALLY

        NIONP=NION+1
        NATMP=NATM+1
        NMOLP=NMOL+1
        NPLSP=NPLS+1
        NPHOTP=NPHOT+1
        NADVP=NADV+1
        NADSP=NADS+1
        NCLVP=NCLV+1
        NALVP=NALV+1
        NALSP=NALS+1
        NSNVP=NSNV+1

        NTALG=30        ! number of VOLUME INPUT TALLIES,  WITHOUT COUNTING GRADIENT TALLIES:
c                         INCREASED IN 2014 FROM 21 TO 22
c                         INCREASED IN 2018 FROM 22 TO 24
c                         INCREASED IN 2018 FROM 24 TO 30
        NTALI=NTALG*4   ! total number of VOLUME INPUT TALLIES
c                         INCLUDE ALL POSSIBLE GRADIENT TALLIES d(TL)/dX, d(TL)/dY,  d(TL)/dZ...

c  additional volume-averaged INPUT tallies
        NTALN=12  ! (ADIN: ADDITIONAL INPUT TALLIES)
        NTALO=14  ! (CELL VOLUME)

        NTALV=100  ! total number of VOLUME-AVERAGED OUTPUT TALLIES
c  additional volume-averaged output tallies
        NTALA=57
        NTALC=58
        NTALT=59
        NTALM=60
        NTALB=61
        NTALR=62

! SURFACE-AVERAGED OUTPUT TALLIES: INCREASED IN 2014 FROM 59 TO 84 (MORE SPUTTER TALLIES)
        NTALS=84
c  additional surface-averaged output tallies
        NTLSA=NTALS-2
        NTLSR=NTALS-1

C  MAX SPECIES INDEX IN SURFACE-AVERAGED OUTPUT TALLIES
        N2MX=MAX(NPHOT,NATM,NMOL,NION,NPLS,NADS,NALS)

        NSPZ=NPHOT+NATM+NMOL+NION+NPLS  ! TOTAL NUMBER OF MC SPECIES PLUS BULK
        NSPZP=NSPZ+1
        NSPZMC=NPHOT+NATM+NMOL+NION     ! TOTAL NUMBER OF MC SPECIES


C  TOTAL NUMBER OF SURFACE-AVERAGED TALLIES
C  SET IN SETPRM ACCORDING TO THE LIVING TALLIES SPECIFIED IN LIVTALS
        NSFTLP=17*NATMP+17*NMOLP+17*NIONP+17*NPHOTP+7*NPLSP+6+
     P        1*NADSP+1*NALSP+1*NSPZP

C  SURFACE REFLECTION DATA
        NHD1=12
        NHD2=7
        NHD3=5
        NHD4=5
        NHD5=5

C  ATOMIC DATA STORAGE. CURRENTLY ONLY TWO OPTIONS
C  NSTORAM=0     : --> NHSTOR=0 --> NSTORDT=1,       NSTORDR=1
C  NSTORAM=9     : --> NHSTOR=1 --> NSTORDT=NSTORAM, NSTORDR=NRAD
        NHSTOR=1-1/(NSTORAM+1)
        NSTORDT=NHSTOR*NSTORAM+(1-NHSTOR)*1
        NSTORDR=NHSTOR*NRAD+   (1-NHSTOR)*1

! NUMBER OF TRAJECTORIES THAT CAN BE STORED
        NTRJ = 1

C
      ELSE IF (ICAL == 2) THEN

c  set some derived storage parameters
        NBGV=NBGK*3
        NCPVP=NCPV+1
        NBGVP=NBGV+1
        NCOLMC=NPLS+NREI+NREC

C  N1MX: storage parameter for species text for output tallies, and scltal in mcarlo.f

        N1MX=    NSPZ+NADV+NALV+NCLV+NCPV+NBGV+NSNV+NAIN

!pb     N1MX=MAX(NPHOT,NATM,NMOL,NION,NPLS,NADV,NALV,NCLV,NCPV,NBGV,
!    .           NSNV,NAIN)
cdr  same MEANING as n1mx?.  Check: why not n1mx=max(....)

C  NSPZTOT: storage parameter for LMETSP(NSPZTOT) array, for standard deviation estimators
        NSPZTOT = NSPZ+NADV+NALV+NCLV+NCPV+NBGV+NSNV

C  NSPZTOTW: storage parameter for LMETSPW(NSPZTOTW) array, for standard deviation estimators
        NSPZTOTW= NSPZ+NADS+NALS

C  TOTAL NUMBER OF VOLUME-AVERAGED OUTPUT TALLIES
C  SET IN SETPRM ACCORDING TO LIVING TALLIES SPECIFIED IN LIVTALV

        NVLTLP=6*NATMP+6*NMOLP+6*NIONP+6*NPHOTP+1*NADVP+1*NCLVP+
     P         1*NSNVP+1*NCPVP+1*NALVP+1*NBGVP+
     P         4*NPLSP+28+3*(NATMP+NMOLP+NIONP+NPHOTP)+
C
     P         NATMP+NMOLP+NIONP+NPHOTP+NPLSP+5*NPLSP+
     P         3*(NATMP+NMOLP+NIONP+NPHOTP)+4*NPLSP

!pb arrays in module CSDVI_COP no longer needed
!pb     NCPV_STAT=(NCPV+NPLS+2)*NSWIT+1
        NCPV_STAT=1
        NSCOP=NCPV_STAT*NRTALS

      ELSE IF (ICAL == 3) THEN

C  SPATIALLY RESOLVED SURFACE TALLIES?

        NGTSFT=NGSTAL*NGITT
        NLMPGS=NLIM+NSTS+NGTSFT*NSTS

      END IF

      RETURN
      END SUBROUTINE EIRENE_SET_PARMMOD


      SUBROUTINE EIRENE_COLLECT_PARM

      INT_PARM = 0

      INT_PARM(  1) = N1ST
      INT_PARM(  2) = N2ND
      INT_PARM(  3) = N3RD
      INT_PARM(  4) = NADD
      INT_PARM(  5) = NTOR
      INT_PARM(  6) = NRTAL
      INT_PARM(  7) = NLIM
      INT_PARM(  8) = NSTS
      INT_PARM(  9) = NPLG
      INT_PARM( 10) = NPPART
      INT_PARM( 11) = NKNOT
      INT_PARM( 12) = NTRI
      INT_PARM( 13) = NTETRA
      INT_PARM( 14) = NCOORD
      INT_PARM( 15) = NOPTIM
      INT_PARM( 16) = NOPTM1

      INT_PARM( 17) = NSTRA
      INT_PARM( 18) = NSRFS
      INT_PARM( 19) = NSTEP

c  leading dimensions of output tally arrays
      INT_PARM( 20) = NATM
      INT_PARM( 21) = NMOL
      INT_PARM( 22) = NION
      INT_PARM( 23) = NPLS
      INT_PARM( 24) = NPHOT
      INT_PARM( 25) = NADV
      INT_PARM( 26) = NADS
      INT_PARM( 27) = NCLV
      INT_PARM( 28) = NSNV
      INT_PARM( 29) = NALV
      INT_PARM( 30) = NALS
      INT_PARM( 31) = NAIN  !  this is an input tally!
      INT_PARM( 32) = NCPV
      INT_PARM( 33) = NBGK  !dr  tbd: more logical: put NBGV here, eliminate NBGK

c  variances, covariances
      INT_PARM( 34) = NSD
      INT_PARM( 35) = NSDW
      INT_PARM( 36) = NCV
c  collision processes
      INT_PARM( 37) = NREAC
      INT_PARM( 38) = NREC
      INT_PARM( 39) = NREI
      INT_PARM( 40) = NRCX
      INT_PARM( 41) = NREL
      INT_PARM( 42) = NRPI
      INT_PARM(134) = NRPH
c  surface reflection model
      INT_PARM( 43) = NHD1
      INT_PARM( 44) = NHD2
      INT_PARM( 45) = NHD3
      INT_PARM( 46) = NHD4
      INT_PARM( 47) = NHD5
      INT_PARM( 48) = NHD6
c  lines of sight integrals (post-processing)
      INT_PARM( 49) = NCHOR
      INT_PARM( 50) = NCHEN

c  parameters for 2d cfd- code coupling, 2d polygonal grid, no. of fluids, target sources
      INT_PARM( 51) = NDX
      INT_PARM( 52) = NDY
      INT_PARM( 53) = NFL
      INT_PARM( 54) = NDXP
      INT_PARM( 55) = NDYP
      INT_PARM( 56) = NPTRGT

      INT_PARM( 57) = NPRNL

c  storage reduction parameters
      INT_PARM( 58) = NGEOM_USR
      INT_PARM( 59) = NCOUP_INPUT
      INT_PARM( 60) = NSMSTRA
      INT_PARM( 61) = NSTORAM
      INT_PARM( 62) = NGSTAL

C     INT_PARM( 63) =        !free, not in use.

      INT_PARM( 64) = NRAD
      INT_PARM( 65) = NSWIT
      INT_PARM( 66) = N1F
      INT_PARM( 67) = N2F
      INT_PARM( 68) = N3F
      INT_PARM( 69) = NGITT
      INT_PARM( 70) = NGITTP

c  STORAGE FOR GEOMETRY (GRID) ARRAYS, ELIMINATED IN CASE LEVGEO=10 (EXTERNAL GEOMETRY)
      INT_PARM( 71) = NRADS
      INT_PARM( 72) = N2NDPLGS
      INT_PARM( 73) = N1STS
      INT_PARM( 74) = N2NDS
      INT_PARM( 75) = NTRIS
      INT_PARM( 76) = NKNOTS
      INT_PARM( 77) = NRTALS
c
      INT_PARM( 78) = NGTSFT
      INT_PARM( 79) = NLIMPS
      INT_PARM( 80) = NLMPGS

C     INT_PARM( 81) =        !dr free, not in use.
      INT_PARM( 82) = NBGV   !dr either nbgk or nbgv should be made redundant
      INT_PARM( 83) = NBMAX
      INT_PARM( 84) = NPTAL
      INT_PARM( 85) = NCPV_STAT
      INT_PARM( 86) = NSCOP

      INT_PARM( 87) = NSTRAP

      INT_PARM( 88) = NIONP
      INT_PARM( 89) = NATMP
      INT_PARM( 90) = NMOLP
      INT_PARM( 91) = NPLSP
      INT_PARM( 92) = NPHOTP
      INT_PARM( 93) = NADVP
      INT_PARM( 94) = NADSP
      INT_PARM( 95) = NCLVP
      INT_PARM( 96) = NALVP
      INT_PARM( 97) = NALSP
      INT_PARM( 98) = NSNVP
      INT_PARM( 99) = NCPVP
      INT_PARM(100) = NBGVP

      INT_PARM(101) = NTALI
      INT_PARM(102) = NTALN
      INT_PARM(103) = NTALO
      INT_PARM(104) = NTALV
      INT_PARM(105) = NTALA
      INT_PARM(106) = NTALC
      INT_PARM(107) = NTALT
      INT_PARM(108) = NTALM
      INT_PARM(109) = NTALB
      INT_PARM(110) = NTALR

      INT_PARM(111) = NTALS
      INT_PARM(112) = NTLSA
      INT_PARM(113) = NTLSR
C     INT_PARM(114) = NTALW   !    OUT, WAS SAME AS NTALS
      INT_PARM(114) = NSPZTOTW

      INT_PARM(115) = N1MX
      INT_PARM(116) = N2MX
      INT_PARM(117) = NSPZ
      INT_PARM(118) = NSPZP
      INT_PARM(119) = NSPZMC
      INT_PARM(120) = NCOLMC
      INT_PARM(121) = NSPZTOT


      INT_PARM(122) = NVOLTL
      INT_PARM(123) = NVLTLP
      INT_PARM(124) = NSRFTL
      INT_PARM(125) = NSFTLP


      INT_PARM(126) = NH0
      INT_PARM(127) = NH1
      INT_PARM(128) = NH2
      INT_PARM(129) = NH3

      INT_PARM(130) = NHSTOR
      INT_PARM(131) = NSTORDT
      INT_PARM(132) = NSTORDR

      INT_PARM(133) = NPLT

      INT_PARM(135) = NADSPC

      INT_PARM(136) = NPLSTI
      INT_PARM(137) = NPLSV
      INT_PARM(138) = NTRJ
      INT_PARM(139) = NBACK_SPEC

cdr   INT_PARM(140) = not in use

      INT_PARM(141) = NCORNER
      INT_PARM(142) = NVLPR
      INT_PARM(143) = NSRPR

      INT_PARM(144) = NADSPC_S
      INT_PARM(145) = NADSPC_C
      INT_PARM(146) = NADSPC_D
      INT_PARM(147) = NADSPC_CD

      INT_PARM(148) = NUM_LINES
      INT_PARM(149) = NADV_ADD

      INT_PARM(150) = NINPTL
      INT_PARM(151) = NTALG   
 
      RETURN
      END SUBROUTINE EIRENE_COLLECT_PARM


      SUBROUTINE EIRENE_DISTRIB_PARM

      N1ST        = INT_PARM(  1)
      N2ND        = INT_PARM(  2)
      N3RD        = INT_PARM(  3)
      NADD        = INT_PARM(  4)
      NTOR        = INT_PARM(  5)
      NRTAL       = INT_PARM(  6)
      NLIM        = INT_PARM(  7)
      NSTS        = INT_PARM(  8)
      NPLG        = INT_PARM(  9)
      NPPART      = INT_PARM( 10)
      NKNOT       = INT_PARM( 11)
      NTRI        = INT_PARM( 12)
      NTETRA      = INT_PARM( 13)
      NCOORD      = INT_PARM( 14)
      NOPTIM      = INT_PARM( 15)
      NOPTM1      = INT_PARM( 16)

      NSTRA       = INT_PARM( 17)
      NSRFS       = INT_PARM( 18)
      NSTEP       = INT_PARM( 19)

c  species indices (1st dimension) of output tallies
      NATM        = INT_PARM( 20)
      NMOL        = INT_PARM( 21)
      NION        = INT_PARM( 22)
      NPLS        = INT_PARM( 23)
      NPHOT       = INT_PARM( 24)
      NADV        = INT_PARM( 25)
      NADS        = INT_PARM( 26)
      NCLV        = INT_PARM( 27)
      NSNV        = INT_PARM( 28)
      NALV        = INT_PARM( 29)
      NALS        = INT_PARM( 30)
      NAIN        = INT_PARM( 31)   ! actually: an input tally !
      NCPV        = INT_PARM( 32)
      NBGK        = INT_PARM( 33)

      NSD         = INT_PARM( 34)
      NSDW        = INT_PARM( 35)
      NCV         = INT_PARM( 36)

      NREAC       = INT_PARM( 37)
      NREC        = INT_PARM( 38)
      NREI        = INT_PARM( 39)
      NRCX        = INT_PARM( 40)
      NREL        = INT_PARM( 41)
      NRPI        = INT_PARM( 42)
      NRPH        = INT_PARM(134)

      NHD1        = INT_PARM( 43)
      NHD2        = INT_PARM( 44)
      NHD3        = INT_PARM( 45)
      NHD4        = INT_PARM( 46)
      NHD5        = INT_PARM( 47)
      NHD6        = INT_PARM( 48)

      NCHOR       = INT_PARM( 49)
      NCHEN       = INT_PARM( 50)

      NDX         = INT_PARM( 51)
      NDY         = INT_PARM( 52)
      NFL         = INT_PARM( 53)
      NDXP        = INT_PARM( 54)
      NDYP        = INT_PARM( 55)
      NPTRGT      = INT_PARM( 56)

      NPRNL       = INT_PARM( 57)


      NGEOM_USR   = INT_PARM( 58)
      NCOUP_INPUT = INT_PARM( 59)
      NSMSTRA     = INT_PARM( 60)
      NSTORAM     = INT_PARM( 61)
      NGSTAL      = INT_PARM( 62)
C                 = INT_PARM( 63)

      NRAD        = INT_PARM( 64)
      NSWIT       = INT_PARM( 65)
      N1F         = INT_PARM( 66)
      N2F         = INT_PARM( 67)
      N3F         = INT_PARM( 68)
      NGITT       = INT_PARM( 69)
      NGITTP      = INT_PARM( 70)
c
      NRADS       = INT_PARM( 71)
      N2NDPLGS    = INT_PARM( 72)
      N1STS       = INT_PARM( 73)
      N2NDS       = INT_PARM( 74)
      NTRIS       = INT_PARM( 75)
      NKNOTS      = INT_PARM( 76)
      NRTALS      = INT_PARM( 77)

      NGTSFT      = INT_PARM( 78)
      NLIMPS      = INT_PARM( 79)
      NLMPGS      = INT_PARM( 80)

C     NCPV        = INT_PARM( 81)  !dr  out, NCOP eliminted, only NCPV retained.
      NBGV        = INT_PARM( 82)
      NBMAX       = INT_PARM( 83)
      NPTAL       = INT_PARM( 84)
      NCPV_STAT   = INT_PARM( 85)
      NSCOP       = INT_PARM( 86)

      NSTRAP      = INT_PARM( 87)

      NIONP       = INT_PARM( 88)
      NATMP       = INT_PARM( 89)
      NMOLP       = INT_PARM( 90)
      NPLSP       = INT_PARM( 91)
      NPHOTP      = INT_PARM( 92)
      NADVP       = INT_PARM( 93)
      NADSP       = INT_PARM( 94)
      NCLVP       = INT_PARM( 95)
      NALVP       = INT_PARM( 96)
      NALSP       = INT_PARM( 97)
      NSNVP       = INT_PARM( 98)
      NCPVP       = INT_PARM( 99)
      NBGVP       = INT_PARM(100)

      NTALI       = INT_PARM(101)
      NTALN       = INT_PARM(102)
      NTALO       = INT_PARM(103)
      NTALV       = INT_PARM(104)
      NTALA       = INT_PARM(105)
      NTALC       = INT_PARM(106)
      NTALT       = INT_PARM(107)
      NTALM       = INT_PARM(108)
      NTALB       = INT_PARM(109)
      NTALR       = INT_PARM(110)
C
      NTALS       = INT_PARM(111)
      NTLSA       = INT_PARM(112)
      NTLSR       = INT_PARM(113)
c     NTALW       = INT_PARM(114)  !dr out, was same as ntals
      NSPZTOTW    = INT_PARM(114)


      N1MX        = INT_PARM(115)
      N2MX        = INT_PARM(116)
      NSPZ        = INT_PARM(117)
      NSPZP       = INT_PARM(118)
      NSPZMC      = INT_PARM(119)
      NCOLMC      = INT_PARM(120)
      NSPZTOT     = INT_PARM(121)

      NVOLTL      = INT_PARM(122)
      NVLTLP      = INT_PARM(123)
      NSRFTL      = INT_PARM(124)
      NSFTLP      = INT_PARM(125)


      NH0         = INT_PARM(126)
      NH1         = INT_PARM(127)
      NH2         = INT_PARM(128)
      NH3         = INT_PARM(129)

      NHSTOR      = INT_PARM(130)
      NSTORDT     = INT_PARM(131)
      NSTORDR     = INT_PARM(132)

      NPLT        = INT_PARM(133)

      NADSPC      = INT_PARM(135)

      NPLSTI      = INT_PARM(136)
      NPLSV       = INT_PARM(137)
      NTRJ        = INT_PARM(138)
      NBACK_SPEC  = INT_PARM(139)
cdr   not in use  = INT_PARM(140)

      NCORNER     = INT_PARM(141)
      NVLPR       = INT_PARM(142)
      NSRPR       = INT_PARM(143)

      NADSPC_S    = INT_PARM(144)
      NADSPC_C    = INT_PARM(145)
      NADSPC_D    = INT_PARM(146)
      NADSPC_CD   = INT_PARM(147)

      NUM_LINES   = INT_PARM(148)
      NADV_ADD    = INT_PARM(149)

      NINPTL      = INT_PARM(150)
      NTALG       = INT_PARM(151)

      RETURN
      END SUBROUTINE EIRENE_DISTRIB_PARM


      SUBROUTINE EIRENE_SPEC_TO_SPEC (SPECA, SPECB)

      TYPE(EIRENE_SPECTRUM), INTENT(OUT) :: SPECA
      TYPE(EIRENE_SPECTRUM), INTENT(IN) :: SPECB

      SPECA%SPCMIN  = SPECB%SPCMIN
      SPECA%SPCMAX  = SPECB%SPCMAX
      SPECA%SPCDEL  = SPECB%SPCDEL
      SPECA%SPCDELI = SPECB%SPCDELI
      SPECA%ESP_MIN = SPECB%ESP_MIN
      SPECA%ESP_MAX = SPECB%ESP_MAX
      SPECA%ESP_00  = SPECB%ESP_00
      SPECA%SPC_XPLT= SPECB%SPC_XPLT
      SPECA%SPC_YPLT= SPECB%SPC_YPLT
      SPECA%SPC_SAME= SPECB%SPC_SAME
      SPECA%SPCVX   = SPECB%SPCVX
      SPECA%SPCVY   = SPECB%SPCVY
      SPECA%SPCVZ   = SPECB%SPCVZ

      SPECA%SPCS    = SPECB%SPCS
      SPECA%SGMS    = SPECB%SGMS
      SPECA%STVS    = SPECB%STVS
      SPECA%GGS     = SPECB%GGS

      SPECA%NSPC    = SPECB%NSPC
      SPECA%ISPCTYP = SPECB%ISPCTYP
      SPECA%ISPCSRF = SPECB%ISPCSRF
      SPECA%IPRTYP  = SPECB%IPRTYP
      SPECA%IPRSP   = SPECB%IPRSP
      SPECA%IMETSP  = SPECB%IMETSP
      SPECA%ISRFCLL = SPECB%ISRFCLL
      SPECA%IDIREC  = SPECB%IDIREC
      SPECA%LOG     = SPECB%LOG

      if (associated(speca%spc)) then
        if (size(speca%spc) < specb%nspc+2) deallocate(speca%spc)
      end if
      if (.not.associated(speca%spc))
     .  allocate(speca%spc(0:specb%nspc+1))
      SPECA%SPC     = SPECB%SPC

      IF (ASSOCIATED(SPECB%SDV)) THEN
        if (associated(speca%sdv)) then
          if (size(speca%sdv) < specb%nspc+2) then
            deallocate(speca%sdv)
            deallocate(speca%sgm)
            deallocate(speca%stv)
            deallocate(speca%gg)
          end if
        end if
        if (.not.associated(speca%sdv)) then
          allocate(speca%sdv(0:specb%nspc+1))
          allocate(speca%sgm(0:specb%nspc+1))
          allocate(speca%stv(0:specb%nspc+1))
          allocate(speca%gg(0:specb%nspc+1))
        end if
        SPECA%SDV     = SPECB%SDV
        SPECA%SGM     = SPECB%SGM
        SPECA%STV     = SPECB%STV
        SPECA%GG      = SPECB%GG
      END IF
      END SUBROUTINE EIRENE_SPEC_TO_SPEC



      END MODULE EIRMOD_PARMMOD
