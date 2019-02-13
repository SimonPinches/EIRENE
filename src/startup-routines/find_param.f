Cdr  Purpose: find storage parameters NPARM for a number of allocatable
c             arrays.
c             Set default for NPARM             (for example NATM=0, no atomic species in this run)
c             Read input file, to find NPARMI,  (for example NATMI)
c             Then set NPARM=MAX(NPARM,NPARMI)  (for example NATM=MAX(NATM,NATMI))
c
c    Later these parameters NPARMI may be modified, but NPARM >= NPARMI must be assured always.
C
!pb  11.12.06:  allow letters 'f' or 't' in case name of fem or tetrahedron
!pb             calculation
!pb  27.12.06:  bug fix: increase NSTS in case of time-dependent mode
!pb  15.01.07:  additional line in input block 4 defining HYDKIN model
!pb  02.03.07:  NUMSEC=4 introduced
!pb  20.03.07:  include input block written by HYDKIN model
!pb  22.03.07:  input for NLFEM and NLTET corrected.
!dr  16.01.14:  default NOPTIM changed from 1 to NRAD (automatically), some printout rearranged
!cd  29.10.14:  reading external file for block 4&5: allow comment lines at the beginning of file
!               (same in find_param.f)
cdr  2.2.15:    nflr renamed to nfr (number of TRIM A_on_B files), now same name as in input.f
cdr             because nflr (common CREF) is later used in RDTRIM and REFDAT with a slightly other meaning.
cdr  Jan 2016:  storage for second dimension only if nlpol=true
cdr             to be tested: storage for nplg, if nlpol=false?
cdr             storage for third dimension only if nltor=true
cdr             to be tested:  storage for nltra, if nltor=false?
cdr             to be done: check for comment lines *... synchronized with input.f?
!pb  June 16:   default for NPLSTI changed from 1 to NPLS
!pb  MAY  16:   nrds -> nrei
cdr  March 17:  NPTRGT printed. May have been changed in call to if0parm, block 14.
CDR  May 2017:  try to fix NSTRAI, NSRFSI, consistent with input.f
cdr             same thing: NCPVI, NCPV  (and eliminate old parameters NCOP, NCOPI)
cdr  July 17 :  lmulti, lmulvi:  automatic options for multiple ion temperatures,
cdr                              multiple ion velocities in case of BGK non-lin. collisions
cdr  July 17 :  initialize 2D CFD code coupling parameters NDX,....
c               move NRAD=... after call to if0prm, because of 3D CFD (emc3) coupling
cdr  Jun 18  : various corrections, comments in new (generalized) block 12 options.
cdr            nadv=nadv+10: now out, is contained in more general storage settings.
C
      SUBROUTINE EIRENE_FIND_PARAM
C
C   SET DIMENSIONS (PARAMETERS) FOR ALLOCATABLE ARRAYS (TALLIES, GRIDS, ETC..) AND SET SOME FURTHER DEFAULT VALUES
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMSOU, ONLY: NSTRAI
      USE EIRMOD_COMPRT, ONLY: IUNIN, IUNOUT

      IMPLICIT NONE

      INTEGER :: INDGRD(3), INDPRO(12), IDUM(12)
      INTEGER, ALLOCATABLE :: INDSRC(:), IEIGEN(:)
      INTEGER :: NFR, ISOR, NSRFSI, NRADD,
     .           NREACI, NSTSI, NLIMI, NVOLPL, NSP, ICO,
     .           NSURPR, NVOLPR, NPRNLI, NCHORI,
     .           NCHENI, NSIGI_BGK, NSIGSI, ID, NSIGVI,
     .           NSIGI_COP, NR1ST, NRSEP, NTIME0,
     .           NP1, NP2, NRKNOT, NRPLG, NPPLG,
     .           NITER0, NTPER, NTTRA, NCOOR, NTET,
     .           NT3RD, NTSEP, NTRII, NP2ND, NPPER, NPSEP, NPPLA,
     .           NSIGCI, IREAD, NCOPII, NCOPIE, NREAC_ADD,
     .           NRC, NRE, NLINES, LL, NB1, NB2, NB3, NS1,
     .           NS2, NS3, INM1, INM2, INM3, INMDL, IEND, ITOK, IER,
     .           N_REAC, N_SPEC, N_ATOMS, N_MOL, N_IONS, N_TESTIONS,
     .           N_BULKIONS, NB4, NS4, INM4, IUNIN_SAVE, I1, NPRMUL,
     .           IATM, IMOL, IION, IPHOT, IPLS,
     .           ISTRA, ISPZ,
     .           NUMSEC, IC, NINITL_READ,
     .           NCHTAL, MOD_ADDV, NUM_COMPO,
     .           NUM_CONTRIB, ISP, ITP, IRATIO,
     .           I, J, K,
     .           ILINE, JCOMP, KCONTR, IREAC_ADD,
     .           NREAC_LINES
      REAL(DP) :: SORIND, SORLIM, DUMM1, ROA, ZAA, ZZA, ZGA, YAA, YYA,
     .            ZIA, YP, XP, YIA, YGA, EMIN1, EMAX1
      LOGICAL :: NLSCL, NLTEST, NLANA, NLDRFT, NLCRR, NLERG, NLIDENT,
     .           NLONE, NLMOVIE, LINCL45, NLCASCAD, NLDFST,
     .           NLOLDRAN, NLOCTREE, NLWRMSH
      LOGICAL :: NLSLB, NLCRC,  NLELL, NLTRI,  NLPLG, NLFEM, NLTET,
     .           NLGEN
      LOGICAL :: NLRAD, NLPOL,  NLTOR, NLADD,  NLMLT, NLTRIM
      LOGICAL :: NLTRA, NLTRT, NLTRZ
      LOGICAL :: PLTL2D, PLTL3D, LRPSCUT, LHYDDEF, LADAPT
      LOGICAL :: LDEFSTOR
      LOGICAL :: NLEMIS
      LOGICAL :: LMULTI, LMULVI   ! multiple ion temperatures (per species) multiple ion velocities (per species)
      CHARACTER(420) :: CASENAME, FILENAME, ULINE
      character(420) :: ZEILE, FILE45
      CHARACTER(12) :: HYDKIN_DEFAULT, CHR, CADAPT
      CHARACTER(4) :: CLAB
      CHARACTER(2) :: COR, CREP
      CHARACTER(15), ALLOCATABLE :: HYDSPEC(:), BULK_NAME(:),
     .                              PART_NAME(:)
      CHARACTER(15) :: BNAME
      CHARACTER(1000) :: HLINE
      CHARACTER(8) :: FNAME, FRATIO
C
C  SET DEFAULT VALUES FOR STORAGE PARAMETERS
C
C  GEOMETRY
      N1ST=1
      N2ND=1
      N3RD=1
!pb   NADD=1
      NADD=0
      NTOR=1
      NRTAL=0
      NLIM=1
      NSTS=1
      NPLG=1
      NPPART=1
      NKNOT=1
      NTRI=1
      NTRII=0
      NTETRA=1
      NCOORD=1
C  PRIMARY SOURCE
      NSTRA=1
      NSRFS=1
      NSTEP=1
      NGITT=0
      NPTRGT=1
C  TEST PARTICLE SPECIES AND TALLIES
      NATM=0
      NMOL=0
      NION=0
      NPLS=1
      NPHOT=0
      NADV=0
      NADS=0
      NCLV=0
      NSNV=0
      NALV=0
      NALS=0
      NAIN=1
      NCPV=0
      NBGK=0
      NADSPC=0
C  BACKGROUND PARTICLE SPECIES AND TALLIES
      LMULTI = .FALSE.
      LMULVI = .FALSE.
      NPLSTI=1
      NPLSV=1
C  STATISTICS
      NSD=1
      NSDW=1
      NCV=1
C  ATOMIC DATA
      NREAC=1
      NREAC_ADD=0
      NREC=1
      NREI=1
      NRCX=1
      NREL=1
      NRPI=1
C  LINE-OF-SIGHT DIAGNOSTICS
      NCHOR=0
      NCHEN=0
      NPRNL=0
      NVLPR=0
      NSRPR=0
C  2D CFD INTERFACE

cdr STORAGE FOR COUPLING TO 2D CFD CODES: COMMON 2D GRID SIZE, BACKGROUND FLUIDS, TARGET RECYCLING
C   CORRECT PARAMETERS  MUST BE SET IN INTERFACING TO CFD CODE, MODULE.  INTERFACES/COUPLE_../IF0PRM.f
c   ONLY INITIALIZE HERE  (e.g. for stand alone runs or code coupling to other than 2D CFD codes)
      NDX=0
      NDY=0
      NFL=0
      NPTRGT=0
      NDXP=1
      NDYP=1


C  OPTIMIZATION OF GEOMETRICAL CALCULATIONS: STORAGE FOR IGJUM3(NCELL,NSURF)
C  ALSO AFFECTS ARRAYS NLIMI(NCELL), NLIME(NCELL) OPTIMIZATION OF CALLS TO TIMEA.F
C  NOPTIM=1   IGJUM3 AND NLIMI, NLIME ARRAYS ARE REMOVED, NO OPTIMIZATION
C  ELSE:  STORAGE IS PROVIDED, CH3 OPTIONS CAN BE USED,  IGJUM3(NOPTIM,NSURF), ETC.
      NOPTIM=1  ! DEFAULT WILL BE AUTOMATICALLY SET TO NRAD, BELOW, LDEFSTOR

C  BIT ARITHMETIC FOR (LARGE) IGJUM.. ARRAYS: ONLY VALUES 0 OR 1 ARE ON THESE ARRAYS
C  NOPTM1=1   USE REGULAR INTEGER ARITHMETIC (8 BIT PER INTEGER)
C  NOPTM1= ???   DO WHAT ??  DEFAULT  ?  LDEFSTOR ?
      NOPTM1=1

C  USER-DEFINED GEOMETRY (LEVGEO=10)?
C  NGEOM_USR = 1  ==> STORAGE PROVIDED FOR USER-DEFINED GEOMETRY OPTION
C  NGEOM_USR = 0  ==> NO STORAGE FOR USER-DEFINED GEOMETRY
      NGEOM_USR=0   ! LDEFSTOR:  SHOULD BE MADE DEPENDENT ON WHETHER LEVGEO=10 OR NOT

C  Input from coupling routine
C  NCOUP_INPUT = 0  ==> NO PLASMA INPUT FROM COUPLING ROUTINE
C  NCOUP_INPUT = 1  ==> PLASMA INPUT FROM COUPLING ROUTINE
      NCOUP_INPUT = 1

C  Sum over Strata
C  NSMSTRA = 0  ==> SUM OVER STRATA IS NOT PERFORMED
C  NSMSTRA = 1  ==> SUM OVER STRATA IS PERFORMED
      NSMSTRA=1
C
C  CALCULATION OR STORAGE OF ATOMIC DATA
C  NSTORAM=0,9.    =0: MINIMUM STORAGE, MAXIMUM CALCULATION
C                  =9: MAXIMUM STORAGE, MINIMUM CALCULATION
      NSTORAM=9
C
C  SPATIALLY RESOLVED SURFACE TALLIES?
C  NGSTAL = 0  ==> NO SPACE FOR SPATIALLY RESOLVED SURFACE TALLIES
C  NGSTAL = 1  ==> SPATIALLY RESOLVED SURFACE TALLIES ARE COMPUTED
      NGSTAL=0

C  NUMBER OF BACKGROUND SPECTRA
      NBACK_SPEC=0


c  NEXT:  BROWSE INPUT FILE AND IDENTIFY THE REAL STORAGE NEEDS.
c   e.g. NPARMI, then set the storage (for allocatable arrays): NPARM = MAX(NPARM,NPARMI)
c   in most cases then: NPARM=NPARMI


C
C  UNIT NUMBER FOR INPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15
C
      IF (IUNIN.EQ.5.OR.IUNIN.EQ.8.OR.IUNIN.EQ.10.OR.
     .    IUNIN.EQ.11.OR.IUNIN.EQ.12.OR.IUNIN.EQ.13.OR.
     .    IUNIN.EQ.14.OR.IUNIN.EQ.15) THEN
        WRITE (IUNOUT,*) 'INVALID INPUT STREAM IUNIN: ',IUNIN
        WRITE (IUNOUT,*) 'ERROR EXIT FROM FIND_PARAM.F      '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

      REWIND IUNIN
C
      CALL EIRENE_LEER(3)

      WRITE (IUNOUT,*) 'PRINTOUT FROM EIRENE PRE-PROCESSING:'
      WRITE (IUNOUT,*) 'BROWSE INPUT FOR STORAGE NEEDS (FIND_PARAM.F)'
      WRITE (IUNOUT,*) 'IUNIN = ', IUNIN
      CALL EIRENE_LEER(1)
C
C  read and write header
      READ (IUNIN,'(A72)') ZEILE
      WRITE (IUNOUT,'(A72)') ZEILE
      CALL EIRENE_LEER(1)

c  skip further comments in header
      DO WHILE (ZEILE(1:1).EQ.'*')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      READ (ZEILE,6666) NPRLL,NMODE,NTCPU,NFILE,NITER0,NITER,
     .                  NTIME0,NTIME

      READ (IUNIN,'(A72)') ZEILE
      LDEFSTOR = .FALSE.   ! INDICATES: NO STORAGE OPTIMIZATION INPUT CARD
      IF ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') + INDEX(ZEILE,'T') +
     .     INDEX(ZEILE,'t')) == 0) THEN
        LDEFSTOR = .TRUE.  ! INDICATES: STORAGE OPTIMIZATION INPUT CARD IS READ
C   READ OPTIONAL INPUT CARD FOR STORAGE HANDLING.
C   OTHERWISE: USE DEFAULTS DEFINED ABOVE.
        READ (ZEILE,6666) NOPTIM,NOPTM1,NGEOM_USR,NCOUP_INPUT,
     .                    NSMSTRA,NSTORAM,NGSTAL,NRTAL,NREAC_ADD
        READ (IUNIN,'(A72)') ZEILE
      ENDIF
C
      call fix_logical_input(zeile,14)
      READ (ZEILE,6665) NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
     .                  NLERG,NLIDENT,NLONE,NLMOVIE,NLDFST,
     .                  NLOLDRAN,NLCASCAD,NLOCTREE,NLWRMSH

C  NSTORAM IS REDEFINED, FINALLY EITHER =0  (A&M STORAGE SAVE MODE)
C                                    OR =9  (FULL A&M STORAGE MODE, =DEFAULT)
      NSTORAM = MIN(NSTORAM,9)
      IF (NSTORAM < 9) NSTORAM = 0
      NOPTM1 = MAX(NOPTM1,1)

      WRITE (iunout,*) '*** 1. DATA FOR OPERATING MODE'

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C
C  READ DATA FOR STANDARD MESH, 200---299


      WRITE (iunout,*) '*** 2. DATA FOR VOXEL GRID GENERATION '

C
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) (INDGRD(J),J=1,3)
C
C INPUT SUB-BLOCK 2A
C
C  RADIAL MESH
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLRAD
      IREAD=0
      IF (NLRAD) THEN
C
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,8)
        READ (ZEILE,6665) NLSLB,NLCRC,NLELL,NLTRI,NLPLG,
     .                    NLFEM,NLTET,NLGEN
        READ (IUNIN,6666) NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT,NCOOR
        N1ST = MAX(N1ST,NR1ST)
        IF (NLPLG) NPLG = MAX(NPLG,NRPLG)
        IF (NLPLG) NPPART = MAX(NPPART,NPPLG)
        NKNOT = MAX(NKNOT,NRKNOT)
        NCOORD = MAX(NCOORD,NCOOR)
        IF (INDGRD(1).LE.5) THEN
          IF (NLSLB.OR.NLCRC.OR.NLELL.OR.NLTRI) THEN
            READ (IUNIN,*)
            IF (NLELL.OR.NLTRI) THEN
              READ (IUNIN,*)
              READ (IUNIN,*)
              IF (NLTRI) THEN
                READ (IUNIN,*)
              ENDIF
            ENDIF
          ENDIF
          IF (NLPLG) THEN
            READ (IUNIN,*)
            READ (IUNIN,6666) (NP1,NP2,K=1,NPPLG)
            DO 212 I=1,NR1ST
              READ (IUNIN,6664) (XP,YP,    J=1,NRPLG)
  212       CONTINUE
          ENDIF

          IF (NLFEM .OR. NLTET) THEN
            READ (IUNIN,'(A72)') ZEILE
            CLAB = ZEILE(1:4)
            CALL EIRENE_UPPERCASE(CLAB)
            IF (INDEX(ZEILE,'CASE')==0) THEN
              WRITE (IUNOUT,*) ' ERROR IN GEOMETRY SPECIFICATION '
              WRITE (IUNOUT,*)
     .          ' TRIANGLE OR TETRAHEDRON GRID SWITCHED ON '
              WRITE (IUNOUT,*) ' BUT NO CASENAME SPECIFIED '
              CALL EIRENE_EXIT_OWN(1)
            END IF

            READ (ZEILE(6:),'(A66)') CASENAME
            CASENAME=ADJUSTL(CASENAME)
            LL=LEN_TRIM(CASENAME)
          END IF

          IF (NLFEM) THEN
            NTRII=NR1ST
            NTRI=MAX(NTRI,NR1ST)

            FILENAME=CASENAME(1:LL) // '.npco_char'
            OPEN (UNIT=30+ifoff,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .            FORM='FORMATTED')
            ZEILE='*   '
            DO WHILE (ZEILE(1:1) == '*')
              READ (30+ifoff,'(A)') ZEILE
            END DO

            READ (ZEILE,*) NRKNOT
            CLOSE (UNIT=30+ifoff)
            NKNOT=NRKNOT

            FILENAME=CASENAME(1:LL) // '.elemente'
            OPEN (UNIT=30+ifoff,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .            FORM='FORMATTED')

            ZEILE='*   '
            DO WHILE (ZEILE(1:1) == '*')
              READ (30+ifoff,'(A)') ZEILE
            END DO

            READ (ZEILE,*) NTRII
            CLOSE (UNIT=30+ifoff)
            NTRI=NTRII+1
            NR1ST=NTRI

            FILENAME=CASENAME(1:LL) // '.neighbors'
            OPEN (UNIT=30+ifoff,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .            FORM='FORMATTED')

            ZEILE='*   '
            DO WHILE (ZEILE(1:1) == '*')
               READ (30+ifoff,'(A)') ZEILE
            END DO

            NGITT=1
            DO I=1,NTRII
              READ (30+ifoff,*) ID, NB1, NS1, INM1,
     .                        NB2, NS2, INM2,
     .                        NB3, NS3, INM3
              IF (INM1 /= 0) NGITT = NGITT + 1
              IF (INM2 /= 0) NGITT = NGITT + 1
              IF (INM3 /= 0) NGITT = NGITT + 1
            END DO

            CLOSE (UNIT=30+ifoff)

          ENDIF

          IF (NLTET) THEN
            NTET=NR1ST
            NTETRA=MAX(NTETRA,NR1ST)

            FILENAME=CASENAME(1:LL) // '.npco_char'
            OPEN (UNIT=30+ifoff,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .            FORM='FORMATTED')
            ZEILE='*   '
            DO WHILE (ZEILE(1:1) == '*')
              READ (30+ifoff,'(A)') ZEILE
            END DO

            READ (ZEILE,*) NCOOR
            CLOSE (UNIT=30+ifoff)
            NCOORD=NCOOR

            FILENAME=CASENAME(1:LL) // '.elemente'
            OPEN (UNIT=30+ifoff,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .            FORM='FORMATTED')

            ZEILE='*   '
            DO WHILE (ZEILE(1:1) == '*')
              READ (30+ifoff,'(A)') ZEILE
            END DO

            READ (ZEILE,*) NTET
            CLOSE (UNIT=30+ifoff)
            NTETRA=NTET+1
            NR1ST=NTETRA

            FILENAME=CASENAME(1:LL) // '.neighbors'
            OPEN (UNIT=30+ifoff,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .            FORM='FORMATTED')

            ZEILE='*   '
            DO WHILE (ZEILE(1:1) == '*')
               READ (30+ifoff,'(A100)') ZEILE
            END DO

            NGITT=1
            DO I=1,NTET
              READ (30+ifoff,*) ID, NB1, NS1, INM1,
     .                        NB2, NS2, INM2,
     .                        NB3, NS3, INM3,
     .                        NB4, NS4, INM4
              IF (INM1 /= 0) NGITT = NGITT + 1
              IF (INM2 /= 0) NGITT = NGITT + 1
              IF (INM3 /= 0) NGITT = NGITT + 1
              IF (INM4 /= 0) NGITT = NGITT + 1
            END DO

            CLOSE (UNIT=30+ifoff)

          ENDIF
        ELSEIF (INDGRD(1).EQ.6) THEN
C  IS THERE ONE MORE LINE, OR IS NLPOL THE NEXT VARIABLE
          READ (IUNIN,'(A72)') ZEILE
          DO WHILE ((ZEILE(1:1) .NE. '*') .AND.
     .         ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') +
     .         INDEX(ZEILE,'T') + INDEX(ZEILE,'t')) == 0))
             READ (IUNIN,'(A72)') ZEILE
          END DO
          IREAD=1
        ENDIF
      ENDIF
      NTRIS=NTRII
      NKNOTS=NKNOT
      N1ST = MAX(N1ST,NR1ST)
C
C  POLOIDAL MESH
C
C INPUT SUB-BLOCK 2B
C
      IF (IREAD == 0) READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLPOL
C
      READ (IUNIN,*)
      READ (IUNIN,6666) NP2ND,NPSEP,NPPLA,NPPER
      IF (INDGRD(2).LE.5) THEN
        READ (IUNIN,6664) YIA,YGA,YAA,YYA
      ENDIF
cdr  storage for 2nd coordinate grid only, if NLPOL=T
      IF (NLPOL) N2ND = MAX(N2ND,NP2ND)
C
C  TOROIDAL MESH
C
C INPUT SUB-BLOCK 2C
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLTOR
C
      READ (IUNIN,'(A72)') ZEILE
      call fix_logical_input(zeile,3)
      READ (ZEILE,6665) NLTRZ,NLTRA,NLTRT
      READ (IUNIN,6666) NT3RD,NTSEP,NTTRA,NTPER
      IF (INDGRD(3).LE.5) THEN
        READ (IUNIN,6664) ZIA,ZGA,ZAA,ZZA,ROA
      ELSEIF (INDGRD(3).EQ.6) THEN
      ENDIF
      IF (NLTOR.AND.NLTRA) NTTRA=NT3RD
cdr  storage for 3nd coordinate grid only, if NLTOR=T
      IF (NLTOR) N3RD = MAX(N3RD,NT3RD)
CDR STORAGE FOR TOROIDAL EFFECTS, EVEN IF
CDR NO TOROIDAL RESOLUTION IS USED
      NTOR = MAX(NTOR,NTTRA)
C
C  MESH MULTIPLICATION
C
C INPUT SUB-BLOCK 2D
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLMLT
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE ((ZEILE(1:1) .NE. '*') .AND.
     .          ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') +
     .            INDEX(ZEILE,'T') + INDEX(ZEILE,'t')) == 0))
         READ (IUNIN,'(A72)') ZEILE
      END DO
      IF (ZEILE(1:1) .EQ. '*') READ (IUNIN,'(A72)') ZEILE
C
C  ADDITIONAL CELLS OUTSIDE STANDARD MESH
C
      READ (ZEILE,6665) NLADD
C
      IF (NLADD) THEN
        READ (IUNIN,6666) NRADD
        NADD = MAX(NADD,NRADD)
      ENDIF

C  FIND START OF NEXT INPUT BLOCK: 3A

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C
C  READING FOR INPUT BLOCK 2 DONE
C
C
      WRITE (iunout,*) '*** 3. DATA FOR BREP SURFACES'
C  READ DATA FOR NON-DEFAULT SURFACE MODELS ON STANDARD SURFACES
C  300--349
C
      WRITE (iunout,*) '*** 3A. DATA FOR NON-DEFAULT STANDARD SURFACES'
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ(ZEILE,6666) NSTSI
      NSTS = MAX(NSTS,NSTSI)

C  FIND START OF NEXT INPUT BLOCK: 3B

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR ADDITIONAL SURFACES 350--399
C
      WRITE (iunout,*) '*** 3B. DATA FOR ADDITIONAL SURFACES           '
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) NLIMI
      NLIM = MAX(NLIM,NLIMI)

C  FIND START OF NEXT INPUT BLOCK: 4

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR SPECIES SPECIFICATION AND ATOMIC PHYSICS MODULE
C  400--499
C
      WRITE (iunout,*) '*** 4. DATA FOR SPECIES SPECIFICATION AND'
      WRITE (iunout,*) '       ATOMIC PHYSICS MODULE'

! CHECK FOR INCLUDE LINE
      READ (IUNIN,'(A420)') ZEILE
      ULINE = ZEILE
      IREAD=1
      CALL EIRENE_UPPERCASE(ULINE)
      I1 = INDEX(ULINE,'INCLUDE')
      LINCL45 = .FALSE.
      IF (I1 > 0) THEN
c   read block 4 and block 5 from external include file, stream fort.2
        IREAD = 0
        CALL EIRENE_READ_TOKEN(ZEILE(I1+7:),' ',FILE45,ITOK,IER,.FALSE.)
        LINCL45 = .TRUE.
        IUNIN_SAVE = IUNIN
        IUNIN = 2+ifoff
        OPEN (IUNIN,FILE=FILE45,FORM='FORMATTED',ACCESS='SEQUENTIAL')
c  read comment lines on external A&M data file FILE45, stream fort.2
  401   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 401
        IREAD=1  ! now ZEILE contains the first non-comment line from fort.2
        GOTO 402
      END IF
C
C
      IF (IREAD == 0) READ (IUNIN,*)
      WRITE (iunout,*)
     .  '       ATOMIC REACTION CARDS, NREACI DATA FIELDS'
      READ (IUNIN,'(A420)') ZEILE

  402 CALL EIRENE_UPPERCASE(ZEILE)
      IEND=INDEX(ZEILE,'DEFAULT')
      LHYDDEF =.FALSE.
cdr ............................................
      IF (IEND > 0) THEN
        LHYDDEF=.TRUE.
cdr  here error exit: LHYDDEF: unfinished option, proprietary version only...
        WRITE (IUNOUT,*) 'INVALID OPTION LHYDDEF IN INPUT BLOCK 4 '
        WRITE (IUNOUT,*) 'USE LHYDDEF ONLY IN PROPRIETARY VERSIONS'
        WRITE (IUNOUT,*) 'ERROR EXIT FROM FIND_PARAM.F      '
        CALL EIRENE_EXIT_OWN(1)

        CALL EIRENE_READ_TOKEN
     .       (ZEILE(IEND+7:),' ',HYDKIN_DEFAULT,ITOK,IER,.FALSE.)

        IEND = IEND + 7 + ITOK
        CALL EIRENE_READ_TOKEN(ZEILE(IEND+1:),' ',CADAPT,ITOK,IER,
     .                  .FALSE.)
        READ (IUNIN,'(A72)') ZEILE
      END IF
cdr ....................................
      READ (ZEILE,*) NREACI
      NREAC = MAX(NREAC,NREACI)+NREAC_ADD
C
cdr  count the number of reaction cards read here.
      NREAC_LINES=0
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        NREAC_LINES=NREAC_LINES+1
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
      ALLOCATE (PART_NAME(500))
      PART_NAME=REPEAT(' ',15)

cdr  start reading species specification block 4a,4b,4c,4d
      WRITE (iunout,*)
     .  '*** 4A. NEUTRAL ATOMS SPECIES CARDS, NATMI SPECIES'
      READ (IUNIN,*) NATMI
      NATM = MAX(NATM,NATMI)

      ISPZ = 0
      DO IATM=1,NATMI
        READ (IUNIN,'(A72)') ZEILE
        ISPZ = ISPZ + 1
        PART_NAME(ISPZ)(1:8) = ZEILE(4:11)
cdr     READ (ZEILE(12:17),'(2I3)') NMASSA(IATM),NCHARA(IATM)
        READ (ZEILE(30:35),'(2I3)') NUMSEC,NRC
        DO K=1,NRC
cdr  read 2 cards per reaction assigned to IATM.  I.e.:  NRC*NATMI*2 cards
cpb......................................
cdr:  try to identify if there are so-called NONLINEAR BGK collisions, input flag IBGK:
cdr:  to be generalized: there may be other reactions, which require multiple Ti, Vi profiles
          READ (IUNIN,'(12I6)') IDUM(1:12)
          IF (NUMSEC < 3) THEN
            LMULTI = LMULTI .OR. (IDUM(7) /= 0)
            LMULVI = LMULVI .OR. (IDUM(7) /= 0)
          ELSEIF (NUMSEC == 3) THEN
            LMULTI = LMULTI .OR. (IDUM(8) /= 0)
            LMULVI = LMULVI .OR. (IDUM(8) /= 0)
          ELSEIF (NUMSEC == 4) THEN
            LMULTI = LMULTI .OR. (IDUM(9) /= 0)
            LMULVI = LMULVI .OR. (IDUM(9) /= 0)
          END IF
cpb.......................................
          READ (IUNIN,*)
        END DO
      END DO
C
C  READ NEUTRAL MOLECULES SPECIES CARDS
C
      WRITE (iunout,*)
     .  '*** 4B. NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,*) NMOLI
      NMOL = MAX(NMOL,NMOLI)

      DO IMOL=1,NMOLI
        READ (IUNIN,'(A72)') ZEILE
        ISPZ = ISPZ + 1
        PART_NAME(ISPZ)(1:8) = ZEILE(4:11)
cdr     READ (ZEILE(12:14),'(I3)') NMASSM(IMOL)
        READ (ZEILE(30:35),'(2I3)') NUMSEC,NRC
        DO K=1,NRC
cpb......................................
cdr:  try to identify if there are so-called BGK collisions, input flag IBGK:
cdr:  to be generalized: there may be other reactions, which require multiple Ti profiles
          READ (IUNIN,'(12I6)') IDUM(1:12)
          IF (NUMSEC < 3) THEN
            LMULTI = LMULTI .OR. (IDUM(7) /= 0)
            LMULVI = LMULVI .OR. (IDUM(7) /= 0)
          ELSEIF (NUMSEC == 3) THEN
            LMULTI = LMULTI .OR. (IDUM(8) /= 0)
            LMULVI = LMULVI .OR. (IDUM(8) /= 0)
          ELSEIF (NUMSEC == 4) THEN
            LMULTI = LMULTI .OR. (IDUM(9) /= 0)
            LMULVI = LMULVI .OR. (IDUM(9) /= 0)
          END IF
cpb.......................................
          READ (IUNIN,*)
        END DO
      END DO
C
C  READ TEST PARTICLE IONS SPECIES CARDS
C
      WRITE (iunout,*)
     . '*** 4C. TEST IONS SPECIES CARDS, NIONI SPECIES'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,*) NIONI
      NION = MAX(NION,NIONI)

      DO IION=1,NIONI
        READ (IUNIN,'(A72)') ZEILE
        ISPZ = ISPZ + 1
        PART_NAME(ISPZ)(1:8) = ZEILE(4:11)
cdr     READ (ZEILE(12:14),'(I3)') NMASSI(IION)
        READ (ZEILE(30:35),'(2I3)') NUMSEC,NRC
        DO K=1,NRC
cpb......................................
cdr:  try to identify if there are so-called BGK collisions, input flag IBGK::
cdr:  to be generalized: there may be other reactions, which require multiple Ti profiles
          READ (IUNIN,'(12I6)') IDUM(1:12)
          IF (NUMSEC < 3) THEN
            LMULTI = LMULTI .OR. (IDUM(7) /= 0)
            LMULVI = LMULVI .OR. (IDUM(7) /= 0)
          ELSEIF (NUMSEC == 3) THEN
            LMULTI = LMULTI .OR. (IDUM(8) /= 0)
            LMULVI = LMULVI .OR. (IDUM(8) /= 0)
          ELSEIF (NUMSEC == 4) THEN
            LMULTI = LMULTI .OR. (IDUM(9) /= 0)
            LMULVI = LMULVI .OR. (IDUM(9) /= 0)
          END IF
cpb.......................................
          READ (IUNIN,*)
        END DO
      END DO

C  FIND START OF NEXT INPUT BLOCK: 4D
      READ (IUNIN,'(A72)') ZEILE
      NPHOTI=0
      IF (ZEILE(1:3) == '***') GOTO 500
      WRITE (iunout,*) '*4D.   PHOTONS SPECIES CARDS, NPHOTI SPECIES '
      READ (IUNIN,*) NPHOTI
      NPHOT = MAX(NPHOT,NPHOTI)

      DO IPHOT=1,NPHOTI
        READ (IUNIN,'(A72)') ZEILE
        ISPZ = ISPZ + 1
        PART_NAME(ISPZ)(1:8) = ZEILE(4:11)
cdr
        READ (ZEILE(30:35),'(2I3)') NUMSEC,NRC
        DO K=1,NRC
cpb......................................
cdr:  try to identify if there are so-called BGK collisions, input flag IBGK:
cdr:  to be generalized: there may be other reactions, which require multiple Ti profiles
          READ (IUNIN,'(12I6)') IDUM(1:12)
          IF (NUMSEC < 3) THEN
            LMULTI = LMULTI .OR. (IDUM(7) /= 0)
            LMULVI = LMULVI .OR. (IDUM(7) /= 0)
          ELSEIF (NUMSEC == 3) THEN
            LMULTI = LMULTI .OR. (IDUM(8) /= 0)
            LMULVI = LMULVI .OR. (IDUM(8) /= 0)
          ELSEIF (NUMSEC == 4) THEN
            LMULTI = LMULTI .OR. (IDUM(9) /= 0)
            LMULVI = LMULVI .OR. (IDUM(9) /= 0)
          END IF
cpb.......................................
          READ (IUNIN,*)
        END DO
      END DO
C
C  READ DATA FOR PLASMA BACKGROUND, 500--599
C
  500 WRITE (iunout,*) '*** 5. DATA FOR PLASMA BACKGROUND'
C
C  READ BULK IONS SPECIES CARDS
C
      WRITE (iunout,*) '*** 5A. BULK ION SPECIES CARDS, NPLSI SPECIES'

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ(ZEILE,6666) NPLSI
      NPLS = MAX(NPLS,NPLSI)

      ALLOCATE (BULK_NAME(NPLSI))
      BULK_NAME = REPEAT(' ',15)
      ICO = 0
      DO IPLS=1,NPLSI
        READ (IUNIN,'(A72)') ZEILE
        BULK_NAME(IPLS)(1:8) = ZEILE(4:11)
        ISPZ = ISPZ + 1
        PART_NAME(ISPZ)(1:8) = ZEILE(4:11)
        ULINE = ZEILE
        CALL EIRENE_UPPERCASE(ULINE)
        INMDL=INDEX(ULINE,'FORT')+
     .        INDEX(ULINE,'SAHA')+
     .        INDEX(ULINE,'CORONA')+
     .        INDEX(ULINE,'BOLTZMANN')+
     .        INDEX(ULINE,'COLRAD')+
     .        INDEX(ULINE,'CONSTANT')
        IF (INMDL > 0) ICO = ICO + 1
        READ (ZEILE(33:35),'(I3)') NRC
        DO K=1,NRC
          READ (IUNIN,*)
          READ (IUNIN,*)
        END DO

        IF (INMDL > 0) THEN
          NRE=0
          IF (VERIFY(ULINE(INMDL+11:),' ') > 0)
     .      READ(ZEILE(INMDL+11:),*) NRE
          NRE=MAX(NRE,1)
          NLINES=1
          IF (INDEX(ULINE(INMDL:),'COLRAD') > 0) NLINES=NRE
          DO K=1,NLINES
             READ (IUNIN,*)
          END DO
        END IF
      END DO

cdr  this next line is probably not needed.
cdr  ico > 0 indicates: at least one bulk species has a special
cdr  background data model,  fort.., saha, corona, ...etc...
      IF (ICO > 0) NREAC=NREAC+1

cdr ..................................................................
      IF (LHYDDEF) THEN

c   unfinished option, only for proprietary version of code.
cdr here should come an error exit: unfinished option, proprietary version only...

        HYDKIN_DEFAULT=ADJUSTL(HYDKIN_DEFAULT)
        LL=LEN_TRIM(HYDKIN_DEFAULT)
        FILENAME=HYDKIN_DEFAULT(1:LL) // '.reactions'
        OPEN (UNIT=27+ifoff,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .        FORM='FORMATTED')
        DO
          READ (27+ifoff,'(A80)') ZEILE
          CALL EIRENE_UPPERCASE(ZEILE)
          IF (INDEX(ZEILE,'N_REAC') /= 0) EXIT
        END DO
        READ (ZEILE,*) CHR,n_reac
        READ (27+ifoff,*) CHR,n_spec
        READ (27+ifoff,*) CHR,n_atoms
        READ (27+ifoff,*) CHR,n_ions
        READ (27+ifoff,*) CHR,n_mol

        ALLOCATE (HYDSPEC(N_SPEC))
        ALLOCATE (IEIGEN(N_SPEC))

        READ (27+ifoff,*) HYDSPEC(1:N_SPEC)
        READ (27+ifoff,'(A1000)') HLINE
        READ (HLINE(52:),*) IEIGEN(1:N_SPEC)

        CLOSE (UNIT=27+ifoff)

        cor = '  '
        crep = '  '
        ladapt = len_trim(cadapt) > 0
        if (ladapt) then
          ic = index(cadapt,'-->')
          if (ic == 0) then
            ladapt = .false.
          else
            cor = cadapt(1:ic-1)
            crep = cadapt(ic+3:ic+4)
          end if
        end if

!pb        N_BULKIONS = COUNT((SCAN(HYDSPEC(1:N_SPEC),'+') > 0) .AND.
!pb     .                     (IEIGEN(1:N_SPEC) == 0))
!pb        N_TESTIONS = N_IONS - N_BULKIONS

!PB        N_BULKIONS = COUNT(IEIGEN(1:N_SPEC) == 0)

        N_BULKIONS = 0
        ILOOP: DO I = 1, N_SPEC
          IF (IEIGEN(I) == 0) THEN
            BNAME = REPEAT(' ',15)
            call EIRENE_remove_char (hydspec(i),BNAME,'_^')
            if (ladapt)
     .        call EIRENE_replace_string(BNAME,cor,crep,iunout)
            DO IPLS = 1, ISPZ
              IF (PART_NAME(IPLS)(1:8) == BNAME(1:8)) CYCLE ILOOP
            END DO
! SPECIES NAME IS UNKNOWN IN INPUTFILE ==> INCREASE NO. OF BULKS
            N_BULKIONS = N_BULKIONS + 1
          END IF
        END DO ILOOP

        N_TESTIONS = COUNT((SCAN(HYDSPEC(1:N_SPEC),'+-') > 0) .AND.
     .                     (IEIGEN(1:N_SPEC) /= 0))

        NATM = NATM + N_ATOMS
        NMOL = NMOL + N_MOL
        NION = NION + N_TESTIONS
        NPLS = NPLS + N_BULKIONS
        NREAC = NREAC + N_REAC
        NREAC_LINES = NREAC_LINES + N_REAC

        DEALLOCATE (HYDSPEC)
        DEALLOCATE (IEIGEN)

      END IF  ! LHYDDEF

      DEALLOCATE (BULK_NAME)
      DEALLOCATE (PART_NAME)
cdr ..................................................................

      READ (IUNIN,'(A72)') ZEILE
      WRITE (IUNOUT,*) '*** 5B. PLASMA BACKGROUND DATA'
      DO WHILE (ZEILE(1:1) == '*')
         READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) (INDPRO(J),J=1,12)

cdr to be done: synchronisation of options for Ti and Vi.
cdr these next 2 lines for Ti(ipls) have historically been just opposite to V_IN options.
!pb   NPLSTI = 1
!pb   IF ((INDPRO(2) < 0) .OR. (MOD(INDPRO(2),100) > 9)) NPLSTI=NPLS

      NPLSTI = NPLS  !dr now same as for V_IN.  Good
cdr   IF (MOD(ABS(INDPRO(2)),100) > 9) NPLSTI = 1  this should be here, to synchronize with Vi
      IF (INDPRO(2)<0) NPLSTI=1  !dr  different still from V_IN logic. Bad

      IF ((NPLS > 1) .AND. (NPLSTI == 1)) THEN
        WRITE (IUNOUT,*) 'WARNING FROM FIND_PARAM'
        WRITE (IUNOUT,*) 'TIIN PROVIDED FOR ONE SPECIES ONLY'
        WRITE (IUNOUT,*) 'DUE TO INDPRO(2) < 0'  !dr  or:  > 10 ???
        IF (LMULTI) THEN
          WRITE (IUNOUT,*) 'DIMENSION OF TIIN OVERWRITTEN'
          WRITE (IUNOUT,*) 'BECAUSE BGK REACTIONS ARE PRESENT'
          NPLSTI = NPLS
        END IF
        WRITE (IUNOUT,*) ' NPLSTI = ',NPLSTI
      END IF

cdr these next 2 lines for V_IN(ipls)
      NPLSV = NPLS
      IF (MOD(ABS(INDPRO(4)),100) > 9) NPLSV = 1

      IF ((NPLS > 1) .AND. (NPLSV == 1)) THEN
        WRITE (IUNOUT,*) 'WARNING FROM FIND_PARAM'
        WRITE (IUNOUT,*) 'V_IN PROVIDED FOR ONE SPECIES ONLY'
        WRITE (IUNOUT,*) 'DUE TO INDPRO(4) > 10'  !dr above, for Ti, we say:  < 0
        IF (LMULVI) THEN
          WRITE (IUNOUT,*) 'DIMENSION OF V_IN ARRAYS OVERWRITTEN'
          WRITE (IUNOUT,*) 'BECAUSE BGK REACTIONS ARE  PRESENT'
          NPLSV = NPLS
        END IF
        WRITE (IUNOUT,*) ' NPLSV = ',NPLSV
      END IF

C  FIND START OF NEXT INPUT BLOCK: 6

      IF (LINCL45) THEN
        CLOSE (IUNIN)
        IUNIN = IUNIN_SAVE
        LINCL45 =.FALSE.
      END IF

      DO
        READ (IUNIN,'(A72)') ZEILE
!pb     IF ((ZEILE(1:3) == '***') .AND.
!    ,      (INDEX(ZEILE,'6.') > 0)) EXIT
!pb  ,      (INDEX(ZEILE,'6') > 0)) EXIT
        IF (ZEILE(1:5) == '*** 6') EXIT
      END DO
C
C  READ  DATA FOR REFLECTION MODEL  600--699
C
      WRITE (iunout,*) '*** 6. GENERAL DATA FOR REFLECTION MODEL'
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLTRIM
      NFR=0
      IF (NLTRIM) THEN
        READ (IUNIN,'(A72)') ZEILE
        IF (INDEX(ZEILE,'PATH')+INDEX(ZEILE,'path').NE.0) THEN
C  PATH SPECIFICATION FOR DATABASE FOUND
          READ (IUNIN,'(A72)') ZEILE
          DO WHILE ((INDEX(ZEILE,'ON')+INDEX(ZEILE,'on')) > 0)
            NFR=NFR+1
            READ (IUNIN,'(A72)') ZEILE
          END DO
        ENDIF
      ENDIF

C  NHD6 FOR STORAGE ON ALLOCATABLE ARRAYS: number of TRIM target-projectile combinations
      IF (NFR > 0) THEN
        NHD6 = NFR
      ELSE
        NHD6 = 12
      END IF

C  FIND START OF NEXT INPUT BLOCK: 7

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR PRIMARY SOURCE  700--799
C
      WRITE (iunout,*) '*** 7. DATA FOR PRIMARY SOURCES, NSTRAI STRATA'
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) NSTRAI
      NSTRA = MAX(NSTRA,NSTRAI)

CDR  TRY TO SET NSTEP, THE NUMBER OF STEP FUNCTIONS FOR SOURCE SAMPLING
cdr  set nstep = smallest stratum number, which receives primary source data from external code.
cdr  this must be highly case specfic. To be reconsidered !!
      ALLOCATE (INDSRC(NSTRAI))
      READ (IUNIN,6666) (INDSRC(J),J=1,NSTRAI)
      IF (ANY(INDSRC == 6)) THEN
        DO J=NSTRAI,1,-1
          IF (INDSRC(J) == 6) THEN
            NSTEP = J
            EXIT
          END IF
        END DO
      END IF

      READ (IUNIN,*)
      DO ISTRA=1,NSTRAI
        IF (INDSRC(ISTRA) == 6) CYCLE
C * 7ABCD...: STRATUM NAME
        READ (IUNIN,'(A72)') ZEILE
        WRITE (IUNOUT,'(A1,A72)') ' ',ZEILE
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,6666) NSRFSI
        NSRFS = MAX(NSRFS,NSRFSI)

        DO I=1,NSRFSI
          READ (IUNIN,*)
          READ (IUNIN,6664) DUMM1, SORLIM, SORIND
          ISOR = INT(SORLIM)
          DO WHILE (ISOR > 0)
cdr  here NSTEP is set to the largest step function number specified on SORIND
            ID = MOD(ISOR,10)
            IF ((ID == 4).OR.(ID==5)) NSTEP = MAX(NSTEP,INT(SORIND))
            ISOR = ISOR / 10
          END DO
          READ (IUNIN,*)
          READ (IUNIN,*)
        END DO
        READ (IUNIN,*)
        READ (IUNIN,*)
C
      END DO
      DEALLOCATE (INDSRC)


      READ (IUNIN,'(A72)') ZEILE
C
C     READ ADDITIONAL DATA FOR SOME SPECIFIC ZONES
C
      WRITE (iunout,*) '*** 8. ADDITIONAL DATA FOR SPECIFIC ZONES'

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO

C
C  READ DATA FOR STATISTICS AND NONANALOG MODEL, 900--999
C
      WRITE (iunout,*) '*** 9. DATA FOR STATISTIC AND NONANALOG MODEL'

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  DATA FOR STANDARD DEVIATION
      WRITE (iunout,*) '       CARDS FOR STANDARD DEVIATION'
      READ (IUNIN,6666) NSIGVI,NSIGSI,NSIGCI,NSIGI_BGK,NSIGI_COP
      NSD = MAX(NSD,NSIGVI)
      NSDW = MAX(NSDW,NSIGSI)
      NCV = MAX(NCV,NSIGCI)

C  FIND START OF NEXT INPUT BLOCK: 10

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C   READ DATA FOR ADDITIONAL AND SURFACE-AVERAGED TALLIES
      WRITE (iunout,*)
     . '*** 10. DATA FOR ADDITIONAL TALLIES, COLLISION'
      WRITE (iunout,*)
     . '        ESTIMATORS AND ALGEBRAIC EXPRESSIONS'
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) NADVI,NCLVI,NALVI,NADSI,NALSI,NADSPC

      NADV = MAX(NADV,NADVI)
      NCLV = MAX(NCLV,NCLVI)
      NALV = MAX(NALV,NALVI)
      NADS = MAX(NADS,NADSI)
      NALS = MAX(NALS,NALSI)
C
      WRITE (iunout,*) '*** 10A. DATA FOR ADDITIONAL TALLIES'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10B. DATA FOR COLLISION ESTIMATORS'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10C. DATA FOR ALGEBRAIC EXPRESSIONS'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10D. DATA FOR ADDITIONAL SURFACE TALLIES'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10E. DATA FOR ALGEBRAIC SURFACE TALLIES'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10F. DATA FOR SPECTRA'
      READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:3) .NE. '***') THEN
        READ (IUNIN,'(A72)') ZEILE
        DO WHILE (ZEILE(1:3) .NE. '***')
          READ (IUNIN,'(A72)') ZEILE
        END DO
      END IF
C
C   READ DATA FOR NUMERICAL AND GRAPHICAL OUTPUT 1100--1199
C
      WRITE (iunout,*)
     .  '*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT'
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C   READ TRCSRC (60 LOGICALS PER LINE)
      do j=0, NSTRAI, 60
        READ (IUNIN,*)
      end do

      READ (IUNIN,6666) NVOLPR
      NVLPR=NVOLPR
C  ERGODIC OPTION NEEDS PRINTOUT OF VOLUME, AND ONE, TWO OR THREE FURTHER TALLIES AT LEAST
      IF (NLERG) NVLPR=MAX(4,NVLPR)
      DO J=1,NVOLPR
        READ (IUNIN,*)
      END DO
C
      READ (IUNIN,6666) NSURPR
      NSRPR=NSURPR
C  ERGODIC OPTION NEEDS PRINTOUT AT LEAST FROM TIME-HORIZON
      IF (NLERG) NSRPR=MAX(1,NSRPR)
      DO J=1,NSURPR
        READ (IUNIN,*)
      END DO

C  SKIP READING ALSO POSSIBLE LINES FOR DELIBERATE DE-ACTIVATION OR RE-ACTIVATION OF TALLIES
c     to be written:  allow for comment lines here

C  SEARCH START OF PLOTTING INPUT: NEXT LINE WITH F OR T
      READ (IUNIN,'(A72)') ZEILE
      CALL EIRENE_UPPERCASE(ZEILE)
      DO WHILE ((SCAN(ZEILE,'FT') == 0) .OR. (ZEILE(1:1) == '*'))
        READ (IUNIN,'(A72)') ZEILE
        CALL EIRENE_UPPERCASE(ZEILE)
      END DO
C
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C  FIRST CARD FOR (LOGICAL) PLOTTING FLAGS NOW ON 'zeile'
      READ (ZEILE(16:16),'(L1)') LRPSCUT  ! needed below for further storage considerations

c  there are 13 geometry plotting input cards following, before plotting of volume tallies
      DO J=1,13
        READ (IUNIN,*)
      END DO
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  DATA FOR PLOTS OF VOLUME-AVERAGED TALLIES
      READ (ZEILE,6666) NVOLPL
C
C
      NPLT = 1
      IF (NVOLPL > 0) THEN
C   READ PLTSRC (60 LOGICALS PER LINE)
        DO J=0, NSTRAI, 60
          READ (IUNIN,*)
        END DO

cdr wrong place for this card here
        IF (LRPSCUT) READ (IUNIN,*) !dr if the "raps-cut option flags" would be
                                    !dr read only below (3d plots and nlraps) then
                                    !dr this exception would not be needed at all.

        DO J=1,NVOLPL
          READ (IUNIN,'(A72)') ZEILE
          DO WHILE (ZEILE(1:1) .EQ. '*')
            READ (IUNIN,'(A72)') ZEILE
          END DO
          READ (ZEILE,6666) NSP
          NPLT = MAX(NPLT, NSP)
          READ (IUNIN,'(A72)') ZEILE
          call fix_logical_input(zeile,2)
          READ (ZEILE,6665) PLTL2D,PLTL3D
          READ (IUNIN,*)
          IF (PLTL2D) THEN
            READ (IUNIN,*)
            DO I=1,NSP
              READ (IUNIN,*)
            END DO
          ENDIF
          IF (PLTL3D) THEN
            READ (IUNIN,*)
            READ (IUNIN,*)
            DO I=1,NSP
              READ (IUNIN,*)
            END DO
            READ (IUNIN,*)
          ENDIF
        END DO
      END IF
C
C  SKIP INPUT LINES, UNTIL INPUT BLOCK 12 STARTS
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR DIAGNOSTIC MODULE  1200--1299
C
      WRITE (iunout,*) '*** 12. DATA FOR DIAGNOSTIC MODULE'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      IREAD=1

c  optional input cards: 'DEFINE_LINES'

c  read up to NUM_LINES transitions (volumetric line emissions),
c  Each LINE may consist of NUM_CONTRIB
c  for different parent (donor) state components.
c  Identify the block of LINES and COMPONENTS available in this run
C  by an extra input card containing 'DEFINE_LINES'
      ULINE=ZEILE
      CALL EIRENE_UPPERCASE(ULINE)
      NADV_ADD = 0
      NLEMIS = .FALSE.
      NUM_LINES = 0

      IF (INDEX(ULINE,'DEFINE_LINES') > 0) THEN
CDR AT LEAST ONE (OR MORE) VOLUMETRIC LINE EMISSIVITY TALLY DEFINED IN INPUT BLOCK 12
cdr as additional output tally ADDV(...).
        IREAC_ADD = 0
        NLEMIS = .TRUE.
c  read number of lines, and the flag MOD_ADDV for storage mode on ADDV tallies
        READ (IUNIN,6666) NUM_LINES, MOD_ADDV
        DO ILINE=1, NUM_LINES
          READ (IUNIN,'(A80)') ZEILE
          DO WHILE (ZEILE(1:1) == '*')
            READ (IUNIN,'(A80)') ZEILE
          END DO
          READ (IUNIN,6666) NUM_COMPO  ! components of line ILINE
          READ (IUNIN,*)
          IF (MOD_ADDV == 0) THEN
cdr  minimal storage, but each time when a new lines comes,
cdr  the emissivity profiles on ADDV must be re-calculated
            NADV_ADD = MAX(NADV_ADD, (NUM_COMPO + 1))
          ELSE
cdr  all possible emissivity profiles are kept on ADDV tallies.
            NADV_ADD =     NADV_ADD +(NUM_COMPO + 1)
          END IF
          DO JCOMP=1, NUM_COMPO
            READ (IUNIN,*)
            READ (IUNIN,*) NUM_CONTRIB     ! contributions to component JCOMP for line ILINE
            IREAC_ADD = IREAC_ADD + NUM_CONTRIB
cdr  specify all required contributions explicitly.
cdr  In the old default with was automatically detected
cdr     from mass and charge states/numbers of hydrogenic particles.
cdr     And only one set of emission data for all contributions was used,
cdr     plus one or two population ratios.
cdr     Now we provide storage for one additional AM data set for each contribution,
cdr     plus one or two population ratios.
            DO KCONTR = 1, NUM_CONTRIB
              READ (IUNIN,'(3I6,1X,A6)') ISP, ITP, IRATIO, FNAME
cdr skip one more input line in case of TAB2D or ADAS input
              IF (INDEX(FNAME,'ADAS')  .NE. 0 .OR.
     .            INDEX(FNAME,'TAB2D') .NE. 0) READ (IUNIN,*)
cdr do we require a QSS population ratio for this contribution?
              IF (IRATIO > 0) THEN
                IREAC_ADD = IREAC_ADD + 1
                READ (IUNIN,'(18X,1X,A6)') FRATIO
cdr skip one more input line in case of TAB2D or ADAS input
                IF (INDEX(FRATIO,'ADAS')  .NE. 0  .OR.
     .              INDEX(FRATIO,'TAB2D') .NE. 0)  READ (IUNIN,*)
cdr do we require a second QSS population ratio for this contribution?
                IF (IRATIO == 2) THEN
                  IREAC_ADD = IREAC_ADD + 1
                  READ (IUNIN,*)
                  READ (IUNIN,'(18X,1X,A6)') FRATIO
                  IF (INDEX(FRATIO,'ADAS')  .NE. 0 .OR.
     .                INDEX(FRATIO,'TAB2D') .NE. 0)  READ (IUNIN,*)
                END IF
              END IF  !  IRATIO
            END DO    !  NUM_CONTRIB   (POSSIBLE D, H, T CONTRIBUTE TO GROUND STATE EMISSIVITY)
          END DO      !  NUM_COMPO     (E.G.  GROUND STATE
        END DO        !  NUM_LINES     (E.G. BA-ALPHA)

c  STORAGE FOR ADDITIONAL TALLIES NADV_ADD, AND REACTIONS IREAC_ADD (LINE EMISSIVITIES)
        NADV = NADV + NADV_ADD
        NREAC = NREAC + IREAC_ADD

        READ (IUNIN,'(A72)') ZEILE
      END IF

      READ (ZEILE,6666) NCHORI,NCHENI
      NCHOR = MAX(NCHOR,NCHORI)
      NCHEN = MAX(NCHEN,NCHENI)

cdr this next condition for old default: better also check for nchtal=2 ??
      NLEMIS = NLEMIS .OR. (NCHOR > 0)
      IF (NLEMIS.AND.(NUM_LINES == 0)) THEN
! USE OLD HYDROGENIC DEFAULT LINES FOR EMISSIVITY
        MOD_ADDV = 0
        NADV=NADV +7
        NUM_LINES = 6
        NUM_COMPO = 6
! USE MAXIMUM POSSIBLE NUMBER OF CONTRIBUTIONS, AS NCHAR AND NCHRG ARE NOT YET AVAILABLE
        NUM_CONTRIB =
     .   NATMI + NMOLI + 2*NMOLI + 2*NMOLI + 3*NMOLI + NPLSI
cdr  ?? perhaps: in old default only one line possible at a time?
cdr  ?? but why then: num_lines=6 rather than num_lines=1 ?
        NREAC = NREAC + NUM_CONTRIB*NUM_COMPO  !dr: this must be way too large

      END IF

C  SKIP READING REST OF THIS BLOCK
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR TIME-DEPENDENT AND NONLINEAR MODE  1300--1399
C
      IREAD=0
      WRITE (iunout,*)
     .  '*** 13. DATA FOR ITERATIVE AND TIME DEP. OPTION'

C
      READ (IUNIN,6666) NPRNLI, NINITL_READ, NPRMUL
      IF (NPRMUL > 1) NPRNLI = NPRNLI * NPRMUL
      NPRNL = MAX(NPRNL,NPRNLI)

      if ((NTIME.GE.1.AND.NPRNLI > 0).OR.NLERG) THEN
        NSTSI=NSTSI+1
        NSTRAI=NSTRAI+1
      ENDIF
      NSTS = MAX(NSTS,NSTSI)
      NSTRA = MAX(NSTRA,NSTRAI)

C   SNAPSHOT TALLIES AND CENSUS ARRAY
      NSNVI=0
      IF (NPRNLI > 0) THEN
        READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1).NE.'*') THEN
          READ (IUNIN,*)
C   READ DATA FOR SNAPSHOT TALLIES
          READ (IUNIN,*)
          WRITE (iunout,*) '*** 13A. DATA FOR SNAPSHOT TALLIES'
          READ (IUNIN,6666) NSNVI
          IREAD=0  !  fix by XB
        END IF
      END IF
      NSNV = MAX(NSNV,NSNVI)

      IF (NLERG.AND.NPRNLI.LE.0) THEN
C  NO TIME HORIZON DEFINED, DESPITE NLERG=.TRUE.
C  THEREFORE: SET A DEFAULT TIME HORIZON HERE
        IF (NTIME.EQ.0) NTIME=1
        NPRNLI=100
      ENDIF
      NPRNL = MAX(NPRNL,NPRNLI)

cdr        NPRNL is only valid for writing census arrays onto fort.15
cdr  tbd:  when reading fort 15 (census), the size is determined by the
cdr        size of that file, (IPRNL) not by NPRNL

C  SKIP READING REST OF THIS BLOCK
      IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO

C
C  INPUT BLOCK 14 BEGIN
C
C  READ DATA IN INTERFACING SUBROUTINE INFCOP  1400 -- 1499
C
      WRITE (iunout,*) '*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"'
      IF (NMODE.EQ.0) THEN
        READ (IUNIN,6666) NAINI,NCOPII,NCOPIE
        NCPVI=NCOPIE
      ELSE
        NAINI=0
        NCPVI=0
        CALL EIRENE_IF0PRM(IUNIN,IUNOUT)
      ENDIF

cdr  some parameters may have gotten changed in IF0PRM, case-specific
      NAIN = MAX(NAIN,NAINI)
      NCPV = MAX(NCPV,NCPVI)

cdr  due to these changes there, also some derived storage parameters may have changed....
      NRAD=MAX(N1ST*N2ND*N3RD,NTRI*N3RD,NTETRA)+NADD+1 ! as in parmmod


      REWIND IUNIN
      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) 'AUTOMATED STORAGE SETTING (FIND_PARAM.F)'
      CALL EIRENE_LEER(1)
C
cdr  grid size
      WRITE (iunout,'(a14,i8)') 'N1ST        = ',N1ST
      WRITE (iunout,'(a14,i8)') 'N2ND        = ',N2ND
      WRITE (iunout,'(a14,i8)') 'N3RD        = ',N3RD
      WRITE (iunout,'(a14,i8)') 'NADD        = ',NADD
      WRITE (iunout,'(a14,i8)') 'NTOR        = ',NTOR
      WRITE (iunout,'(a14,i8)') 'NRTAL       = ',NRTAL
      WRITE (iunout,'(a14,i8)') 'NLIM        = ',NLIM
      WRITE (iunout,'(a14,i8)') 'NSTS        = ',NSTS
      WRITE (iunout,'(a14,i8)') 'NPLG        = ',NPLG
      WRITE (iunout,'(a14,i8)') 'NPPART      = ',NPPART
      WRITE (iunout,'(a14,i8)') 'NKNOT       = ',NKNOT
      WRITE (iunout,'(a14,i8)') 'NTRI        = ',NTRI
      WRITE (iunout,'(a14,i8)') 'NTETRA      = ',NTETRA
      WRITE (iunout,'(a14,i8)') 'NCOORD      = ',NCOORD

      WRITE (iunout,*) ' '
      WRITE (iunout,'(a14,i8)') 'NRAD        = ',NRAD
cdr  primary source
      CALL EIRENE_LEER(1)
      WRITE (iunout,'(a14,i8)') 'NSTRA       = ',NSTRA
      WRITE (iunout,'(a14,i8)') 'NSRFS       = ',NSRFS
      WRITE (iunout,'(a14,i8)') 'NSTEP       = ',NSTEP
      WRITE (iunout,'(a14,i8)') 'NPTRGT      = ',NPTRGT
cdr species
      CALL EIRENE_LEER(1)
      WRITE (iunout,'(a14,i8)') 'NATM        = ',NATM
      WRITE (iunout,'(a14,i8)') 'NMOL        = ',NMOL
      WRITE (iunout,'(a14,i8)') 'NION        = ',NION
      WRITE (iunout,'(a14,i8)') 'NPHOT       = ',NPHOT
      WRITE (iunout,'(a14,i8)') 'NPLS        = ',NPLS

      CALL EIRENE_LEER(1)
      WRITE (iunout,'(a14,i8)') 'NADV        = ',NADV
      WRITE (iunout,'(a14,i8)') 'NCLV        = ',NCLV
      WRITE (iunout,'(a14,i8)') 'NSNV        = ',NSNV
      WRITE (iunout,'(a14,i8)') 'NALV        = ',NALV

      WRITE (iunout,'(a14,i8)') 'NADS        = ',NADS
      WRITE (iunout,'(a14,i8)') 'NALS        = ',NALS

      WRITE (iunout,'(a14,i8)') 'NAIN        = ',NAIN
      WRITE (iunout,'(a14,i8)') 'NCPV        = ',NCPV
      WRITE (iunout,'(a14,i8)') 'NBGK        = ',NBGK
      WRITE (iunout,'(a14,i8)') 'NSD         = ',NSD
      WRITE (iunout,'(a14,i8)') 'NSDW        = ',NSDW
      WRITE (iunout,'(a14,i8)') 'NCV         = ',NCV

      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) 'MAX. NO. OF ATOMIC/MOLECULAR "REACTIONS"'
      WRITE (iunout,'(a14,i8)') 'NREAC       = ',NREAC
      WRITE (IUNOUT,*) 'NREC,NREI,NRCX,NREL,NRPI: DETERMINED LATER'
C     WRITE (iunout,'(a14,i8)') 'NREC        = ',NREC
C     WRITE (iunout,'(a14,i8)') 'NREI        = ',NREI
C     WRITE (iunout,'(a14,i8)') 'NRCX        = ',NRCX
C     WRITE (iunout,'(a14,i8)') 'NREL        = ',NREL
C     WRITE (iunout,'(a14,i8)') 'NRPI        = ',NRPI

      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) 'SETTING OF STORAGE OPTIMIZATION OPTIONS'
      IF (.NOT.LDEFSTOR) THEN
C  ADJUST SOME DEFAULT STORAGE OPTIMIZATION SETTING
        NOPTIM=N1ST*N2ND*N3RD+NADD     ! = NRAD ??
C       NOPTM1=   ??
      ENDIF
C  TRY TO BE INTELLIGENT:  AUTOMATIC STORAGE REDUCTION....
C  SWITCH OFF SUM OVER STRATA IF THERE IS ONLY ONE STRATUM TO BE CALCULATED
c  TO BE TESTED, E.G. CHECK VARIANCES, AND VARIANCES SUM OVER STRATA
      IF (NSTRAI == 1.AND.NSMSTRA.NE.0) THEN
         WRITE (iunout,*) 'NSMSTRA SET = 0, BECAUSE NSTRAI=1 '
         NSMSTRA = 0
      ENDIF
      CALL EIRENE_LEER(1)

C  OPTIONAL STORAGE/PERFORMANCE HANDLING FLAGS
      WRITE (iunout,'(a14,i8)') 'NOPTIM      = ',NOPTIM
      WRITE (iunout,'(a14,i8)') 'NOPTM1      = ',NOPTM1
      WRITE (iunout,'(a14,i8)') 'NGEOM_USR   = ',NGEOM_USR
      WRITE (iunout,'(a14,i8)') 'NCOUP_INPUT = ',NCOUP_INPUT
      WRITE (iunout,'(a14,i8)') 'NSMSTRA     = ',NSMSTRA
      WRITE (iunout,'(a14,i8)') 'NSTORAM     = ',NSTORAM
      WRITE (iunout,'(a14,i8)') 'NGSTAL      = ',NGSTAL
      WRITE (iunout,'(a14,i8)') 'NREAC_ADD   = ',NREAC_ADD
C
      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) 'SETTING OF CENSUS STORAGE FOR T-DEP. MODE'
cdr  time-dependent options: census array size
      WRITE (iunout,'(a14,i8)') 'NPRNL       = ',NPRNL
C
      CALL EIRENE_LEER(2)
C
      RETURN
C
 6664 FORMAT (6E12.4)
 6665 FORMAT (12(5L1,1X))
 6666 FORMAT (12I6)
C
 6999 WRITE (IUNOUT,*) 'Empty input file found!'
      WRITE (IUNOUT,*)
     . 'Either remove it or replace it with a correct file.'
      CALL EIRENE_EXIT_OWN(1)
 7999 WRITE (IUNOUT,*) 'Could not open input file!'
      CALL EIRENE_EXIT_OWN(1)
      RETURN

      END
