cdr  oct 18 :   unify reading of A&M data (reaction decks) from external file:
cdr             Formerly from block 4, block 5 ("density models") 
cdr             and block 12 (line emissivities).
cdr             Now: unified interpreter of reaction card, subr. READ_REACLINES.f.
cdr             Reading only from the list in  block 4, IR=1,NREACI.
cdr             ALL reaction data potentially needed in blocks 5 and 12
cdr             must have been read in block 4, as well as their transfer to internal reaction
cdr             reaction data structures (calls to SLREAC.f) are already carried out.
cdr             Block 5 TDMPAR (density model) option: revised
cdr sept. 18:   iopt:  ?? further optional input lines at the end of block 5?
cdr             for turning on/off input tally storage, and for gradients of
cdr             input tallies (optional)
cdr  sept.18:   XDR format options for fort.13 stream: removed.
cdr  apr. 18:   fully connected and tested: trchktm option, in block 11.
cdr  july 17 :  GR cleanup: wrmesh option split into writing and plotting
Cdr  april 17:  some cleanup (spelling, trim(character)) adopted from sols_iter version
cdr             added: logical NEXVS   (default: F. Unclear meaning, so far...)
c               added: logical TRCRNF  (traceback for random seeds for correlated sampling)
cdr  sept. 16:  extend options for extrapolations for A&M data beyond range
cdr             of tables or validity range fit expressions.

!               rename RMN and RMX to R1MN, R1MX, add R2MN, R2MX for range
!               of second variable in fit or data table
!               same with jfexmn,jfexmx (parameters to select extrapolation scheme)

!pb  June  16:  default for NPLSTI changed from 1 to NPLS
cdr             indpro(2) and indpro(4): try to synchronize the meaning, to be done
cdr  june  16:  comments, disable accidental use of HYDKIN interface,
cdr             option lhyddef. error exit. Tests of that interface options started.
!cd  jan   16:  reset census start time to time0, even for time0=0.
!cd  jan   16:  remove unused: ILE,tpb1,...,ian,ien,iab,reac
!cd  dec.  15:  jj-nlim, rather than jj-nlimi, for non.dev.std. surfaces
!    april 15:  esptcr, esptsr: sputtered particle energy flags introduced
!cd  29.10.14:  reading external file for block 4&5: allow comment lines at the beginning of file
!               (same in find_param)
!cd  22.09.14:  1D case, levgeo=2:  do not call grid(2)
!cd  22.03.14:  option 'include filname ' instead of block 4 and 5 tested and verified
!               some minor changes at transition from end of block ***3 and re-entry to block ***6
!pb  01.01.14:  options AMPTS, multiplier for ntcpu.... added (input block 7)
!dr  03.03.10:  READING A&M DATA in block 4, commented, sorted,.....
!pb  20.03.08:  allocate and nullify estiml(1) if no input block 10F
!pb             is available
!pb  01.08.07:  save NLSRON in case of coupled run with short cycle
!    20.06.07:  check NVOLPR against NSPEZV_DIM
!pb  22.05.07:  input of NPTSDEL added in block 7
!pb  20.04.07:  check of geometry flags added
!pb  20.04.07:  check for casename of NLFEM or NLTET
!pb  20.04.07:  option for reading triangle grid from Eirene input
!pb             file removed
!pb  22.03.07:  LEVGEO=6 --> LEVGEO=10
!pb  21.03.07:  logicals in block 1 cleaned up
!pb  02.03.07:  fourth secondary group in species specification cards
!pb             introduced
!pb  02.03.07:  remove ESCD2* arrays, use local variables instead
!pb             sum up ESCD1 and ESCD2 contributions
!pb  12.01.07:  one additional line in input block 4 defining the
!pb             HYDKIN model
!pb  09.01.07:  input of reaction cards (block 4) rewritten using
!pb             read_token
!pb  01.12.06:  bug fix: advance line in input for tetrahedra
!pb  09.10.06:  save NZADD for higher timesteps
!dr  20.04.06:  fort.10 added as density model in block 5.
!dr             Also other density model may now refere to test
!dr             particle tallies, for postprocessing and iteration
!pb  02.03.06:  NLRAY: switch on raytracing method for stratum
!pb  20.01.06:  line of sight for cell-based spectrum introduced
!pb  12.01.06:  flag for cell-based spectrum added in block 10F
C    24.09 05:  CALL TO IF2COP NOW LATER, AFTER ALL PLASMA DATA ARE SET
C               SO THAT IF2COP CAN BE USED WITHOUT HAVING TO USE IF1COP
C
cdr june-05:  spectrum input (10F) extended: SPC_SHIFT,.....
cdr                                SPCPLT_X,SPCPLT_Y,SPCPLT_SAME
cdr           see corresponding changes in CESTIM (ESTIML...)
cdr  28.4.04: nhsts(ispz) introduced, to select species
cdr           for trajectory plot
cdr           default: = 0: "plot trajectory for this species"
cdr           new    : =-1: "do not plot trajectory for this species"
cdr           see modifications in plt2d.f from 28.4.04
cpb sept-05:  specification of filenames and paths for external
cpb           databases in input block 1 added. (CFILE-cards)
cpb
cpb           example:
cpb           CFILE AMJUEL /home/boerner/Database/AMdata/amjuel.tex

      SUBROUTINE EIRENE_INPUT
C
C   READ INPUT DATA AND SET DEFAULT VALUES
c   IN CASE IITER.GT.1 OR ITIMV.GT.1 : SKIP READING NEW INPUT FROM IUNIN.
C                                      ONLY INPUT DATA PROCESSING (STATEMENT 4000 FF)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_COMSIG
      USE EIRMOD_CREF
      USE EIRMOD_CGRID
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_CTEXT
      USE EIRMOD_CLGIN
      USE EIRMOD_CSPEI
      USE EIRMOD_CESTIM
      USE EIRMOD_CUPD
      USE EIRMOD_PHOTON
      USE EIRMOD_TIMEA, ONLY: EIRENE_TIMEA0
      USE EIRMOD_PROFILES
      USE EIRMOD_JSON

      IMPLICIT NONE

C
      TYPE(VOLUMEP),POINTER :: VOLCUR
CC
      REAL(DP) :: VOLTOT_TAL

C  RUN TIME STATISTICS IN INITIALIZATION PHASE, WITHIN INPUT.F
cdr   REAL(DP) :: tpb1, tpb2, EIRENE_SECOND_OWN, timea
c
      REAL(DP), ALLOCATABLE :: SAREA_SAVE(:)
      INTEGER :: IADTYP(0:4)
      INTEGER :: IERROR, IUNIN_SAVE, JSTREAM, NSOPT, ILIMPS, 
     .           I, J, I1, I2, I3, JL, IO, IUSR, IFLG,
     .           ITALI, IRAD, IS, ISS, IN, INC, IRET,
     .           JPLS, IRE, JSPZ, JTRJ

      INTEGER, SAVE :: NITER0, IUSROUT=0
      LOGICAL :: NLSRON_SAVE(NSTRA)
      LOGICAL :: LRDJSON
      CHARACTER(10) :: CDATE, CTIME
      CHARACTER(420) :: ZEILE
C
C  DO NOT READ ANY INPUT, IF THIS IS NOT THE VERY FIRST ITERATION
C  STEP IN THIS RUN. IITER IS THE ACTUAL ITERATION NUMBER
C  DO NOT READ ANY INPUT, IF THIS IS NOT THE VERY FIRST TIMESTEP
C  IN THIS RUN. ITIMV IS THE ACTUAL TIMESTEP NUMBER
C
C  INITIALISE SOME DATA AND SET DEFAULTS
C
      IERROR=0

C  UNIT NUMBER FOR INPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15
!pb IUNIN set in COMPRT
!pb   IUNIN=1+ifoff
      IUNIN_SAVE = IUNIN
      IUSROUT = 0
C
C  UNIT NUMBER FOR OUTPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15 AND IUNIN
!pb IUNOUT has already been set in subroutine EIRENE
!pb   IUNOUT=6+ifoff
C
      IF (IITER.GT.1) THEN
        CALL EIRENE_MASBOX
     .   ('NEXT ITERATION STARTS, SKIP READING INPUT FILE')
        CALL EIRENE_MASJ1('IITER   ',IITER)
        GOTO 4000
      ENDIF

      IF (ITIMV.GT.1) THEN
        CALL EIRENE_MASBOX
     .   ('NEXT TIME-CYCLE STARTS, SKIP READING INPUT FILE ')
        CALL EIRENE_MASJ1('ITIMV   ',ITIMV)
        GOTO 4000
      ENDIF
C
      NAINI=0
      NCPVI=0
      NBGVI=0

      MTSURF=0 ! cdr  ????
C
C  SET DEFAULT REACTION MODELS
C
      CALL EIRENE_SETUP_DEFAULT_REACTIONS

C
C  SET DEFAULT SOURCE MODEL BLOCK 7
C
      NSTRAI=0
C
C  SET DEFAULT 'ADDITIONAL SURFACE' AND 'STANDARD SURFACE' DATA
C
      NBITS=BIT_SIZE(I)
      CALL EIRENE_SET_DEF_SURF_DATA

      ALLOCATE (SAREA_SAVE(NLIMPS))
      SAREA_SAVE = 666.

      NULLIFY(SURFLIST)
      NULLIFY(REFLIST)
C
C  SET DEFAULT DATA FOR BLOCK 13
C
      DTIMV=1.D30
      TIME0=0.
      NSNVI=0
      NTMSTP=1

C  AS DEFAULT: SWITCH OFF MOMENTUM DENSITY TALLIES.
C  FOR ACTIVATING THOSE TALLIES THEY NEED TO BE EXPLICITLY
C  SWITCHED ON IN BLOCK 11

C  LV?DEN.. IS AN ALIAS FOR AN ENTRY IN ARRAY LMISTALV
C  THEREFORE .TRUE. MEANS: TALLY IS SWITCHED OFF
      LVXDENA  = .TRUE.
      LVXDENM  = .TRUE.
      LVXDENI  = .TRUE.
      LVXDENPH = .TRUE.
      LVYDENA  = .TRUE.
      LVYDENM  = .TRUE.
      LVYDENI  = .TRUE.
      LVYDENPH = .TRUE.
      LVZDENA  = .TRUE.
      LVZDENM  = .TRUE.
      LVZDENI  = .TRUE.
      LVZDENPH = .TRUE.

C
      CALL EIRENE_LEER(2)

      CALL EIRENE_ALLOC_JSON_ARRAYS
      CALL EIRENE_ALLOC_CTEXT(3)
C
C  READ TEXT DESCRIBING THE RUN, 100--199
C

      CALL DATE_AND_TIME(CDATE,CTIME)
      READ(CDATE(1:4),*) I1
      READ(CDATE(5:6),*) I2
      READ(CDATE(7:8),*) I3
      WRITE (iunout,'(1X,A6,1X,2(I2,1X),I4)') 'DATE: ',I3,I2,I1
      READ(CTIME(1:2),*) I1
      READ(CTIME(3:4),*) I2
      READ(CTIME(5:6),*) I3
      WRITE (iunout,'(1X,A6,1X,3(I2,1X))') 'TIME: ',I1,I2,I3
      CALL EIRENE_LEER(2)

      REWIND IUNIN

      READ (IUNIN,'(A80)') ZEILE
      
      REWIND IUNIN
      LRDJSON = .FALSE.
      IF (ZEILE(1:1) == '{') THEN
	LRDJSON = .TRUE.
      ELSEIF (ZEILE(1:1) /= '*') THEN
        WRITE (IUNOUT,*) ' EIRENE INPUT FORT.1 HAS WRONG FORMAT'
        CALL EIRENE_EXIT_OWN(1)
      END IF

      IF (LRDJSON) THEN

        close (unit=1)
        CALL EIRENE_READ_JSON(NLIMPS, SAREA_SAVE, ierror)

      ELSE

     	 CALL EIRENE_READ_FIXFORM (NLIMPS, SAREA_SAVE, IERROR)

      END IF
C
      CALL EIRENE_SETUP_TIME_SURFACE (IERROR)
C     
 1500 IF (IERROR.GT.0) THEN
        WRITE (iunout,*) IERROR,' INPUT OR PARAMETER ERRORS DETECTED'
        WRITE (iunout,*)
     .   ' SEE THE ERROR MESSAGES LISTED ABOVE AND CORRECT'
        WRITE (iunout,*) ' THE ERRORS BEFORE RE-EXECUTION'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      CALL EIRENE_PAGE
C
C   MODIFICATION OF INPUT DUE TO EITHER INCONSISTENCIES OR DUE
C   TO COUPLED NEUTRAL-PLASMA (OR NEUTRAL-NEUTRAL) CALCULATIONS
C   SOME FURTHER CONSTANTS ARE SET.      STATEM. NO. 2000 --> 4000

      CALL EIRENE_CHECK_GEOM_CONSIST   
C
C  SOURCE PARAMETERS AND (REFLECTING) BOUNDARY CONDITIONS,
C  ON ADDITIONAL AND NON-DEFAULT STANDARD SURFACES
C
      CALL EIRENE_SETUP_LOC_REF_MODELS

      CALL EIRENE_SETUP_SURFACE_SWITCHES
C
C  SET NON-DEFAULT STANDARD SURFACE IDENTIFIERS INMP...
C
      CALL EIRENE_SETUP_INMP
C     
      CALL EIRENE_PREP_STRATA
C     
C
C
C  SPECIES INDEX DISTRIBUTION OF PRIMARY SOURCE PARTICLES
C  OR FOR THERMAL PARTICLE REFLECTION MODEL

      CALL EIRENE_SETUP_SPEC_IND_DISTRIB

      CALL EIRENE_SETUP_ATOMIC_WEIGHTS
C
C
C  SET SOME ARRAYS TO SPEED UP COMPUTATIONS
C
      CALL EIRENE_SETUP_ISPEZ

      IF (NPHOTI > 0) CALL EIRENE_PH_INIT(1)
      CALL EIRENE_SETAMD(0)
      CALL EIRENE_ALLOC_CTEXT(2)

      CALL EIRENE_SETTXT_INTAL
      CALL EIRENE_SETPRM_INTAL
      CALL EIRENE_SETTXT
C
C
C  ADDITIONAL INPUT FOR THIS RUN COMES FROM EITHER
C  ANOTHER CODE (DATA FILE) OR FROM AN EARLIER RUN OF EIRENE
C
C  INPUT BLOCK 14 BEGIN
C
C  READ DATA IN INTERFACING SUBROUTINE INFCOP  1400 -- 1499
C


! CALL TO ALLOC_BCKGRND MOVED HERE TO ALLOW SPECIFICATION OF VOL
! IN IF0COP

      IF (ANY(INDPRO(1:12) == 6)) CALL EIRENE_ALLOC_BCKGRND

      IF (LRDJSON) THEN

        CALL EIRENE_READ_BLK14_JSON(IERROR)

C  CHECK FOR USER SPECIFIC INPUT
        OPEN(NEWUNIT=IUSR,FILE='user_data.input',STATUS='OLD',
     .       IOSTAT=IO)
        IUSROUT = IUSR
        IUNIN_SAVE = IUNIN
        IF (IO == 0) THEN
          IUNIN = IUSROUT
          CALL EIRENE_LEER(1)
          WRITE (IUNOUT,*) 'USER SPECIFIC INPUT READ FROM ',
     .         'user_data.input'
        ELSE
          CALL EIRENE_LEER(1)
          WRITE (IUNOUT,*) 'NO FILE FOR USR SPECIFIC INPUT FOUND'
        END IF

      ELSE

        CALL EIRENE_READ_BLK14_FIXFORM(IERROR)

C  COPY USER SPECIFIC DATA TO FILE user_data.input
        JL = 0
        IO = 0
        IUNIN_SAVE = IUNIN
        DO WHILE (IO == 0)
          READ (IUNIN,'(A72)',IOSTAT=IO) ZEILE
          IF (IO == 0) THEN
            JL = JL + 1
            IF (JL == 1) THEN
              OPEN(NEWUNIT=IUSR,FILE='user_data.input')
              IUSROUT = IUSR
            END IF
            WRITE (IUSROUT,'(A)') TRIM(ZEILE)
          END IF
        END DO
        IF (JL > 0) THEN
          REWIND IUSROUT
         IUNIN = IUSROUT
        ELSE
          IUSROUT = 0
        END IF

      END IF
     
C
C  INPUT BLOCK 14 DONE
C
      CALL EIRENE_PREP_PLOTTING
      
      CALL EIRENE_CORRECT_STATS_INPUT
      
C
C  NO MODIFICATION OF INPUT VARIABLES BEYOND THIS POINT
C  WITHOUT WARNING
C  EXCEPT IN SUBROUTINE MODUSR FOR THE NEXT ITERATION STEP
C         IN SUBROUTINE TMSUSR FOR THE NEXT TIME STEP
C
C
      CALL EIRENE_INIUSR
C
C  SET DERIVED INPUT PARAMETERS, GRIDS AND PROFILES    
C
      CALL EIRENE_SET_DERIVED_INPUT_PARAMETERS(IERROR)
      

      CALL EIRENE_SET_PARMMOD(3)
      CALL EIRENE_ALLOC_CGEOM(2)
      CALL EIRENE_ALLOC_COMUSR(3)

      SAREA(1:NLIMPS) = SAREA_SAVE(1:NLIMPS)
      DEALLOCATE (SAREA_SAVE)

C
C  WRITE JSON FILE
      IF (LRDJSON) THEN
        CALL EIRENE_WRITE_JSON_AMData ('eirene_AMData.json')
        CALL EIRENE_WRITE_JSON_FILE('eirene_input_json.out')
      ELSE
        CALL EIRENE_WRITE_JSON_FILE('eirene.input.json')
      END IF

!  REMOVE LIST OF REFLECTION MODELS
!  MOVED HERE AS TO BE AVAILABLE FOR WRITING ON JSON-FILE
      
      CALL EIRENE_DEALLOC_REFLIST

      CALL EIRENE_PAGE


C
      IF (NFILEM.LE.1) THEN
C
C  SET GRIDS AND VOLUMES OF THE CELLS FOR THE STANDARD MESHES
C
C  SET RADIAL OR X GRID
        IF (NLRAD) CALL EIRENE_GRID (1)
C  SET POLOIDAL OR Y GRID
        IF (NLPOL.OR.(LEVGEO == 3)) CALL EIRENE_GRID (2)  ! NO NEED TO SET UP POLYGON GRID IN CASE OF 1D RUN, LEVGEO=2
C  SET TOROIDAL OR Z GRID
        IF (NLTOR) CALL EIRENE_GRID (3)
C
        IF (INDPRO(12).LT.4) THEN
C  INITIALISE SUBROUTINE VOLUME FOR LATER CALLS
C  TO BE WRITTEN:    CALL VOLUME(0)
C  SET VOLUMES, 1ST DIMENSION (R/X-GRID)
          IF (NLRAD) CALL EIRENE_VOLUME(1)
C  SET VOLUMES, 2ND DIMENSION (THETA/Y-GRID)
          IF (NLPOL) CALL EIRENE_VOLUME(2)
C  SET VOLUMES, 3RD DIMENSION (PHI/Z-GRID)
         IF (NLTOR) CALL EIRENE_VOLUME(3)
C  SET VOLUMES, IN ADDITIONAL CELL REGION
          IF (NLADD) CALL EIRENE_VOLUME(4)
        ELSEIF (INDPRO(12).EQ.4) THEN
C  READ VOLUMES FROM STREAM JSTREAM=VL0
          JSTREAM=INT(VL0)
          ITALI=NTALO
          CALL EIRENE_READTL(TXTPLS(1,ITALI),TXTPSP(1,ITALI),
     .                TXTPUN(1,ITALI),
     .                VOL,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .                3,JSTREAM)
C  TAKE VOLUMES FROM USER-SUPPLIED ROUTINE
        ELSEIF (INDPRO(12).EQ.5) THEN
          CALL EIRENE_PROUSR (VOL,5+5*NPLS,0._DP,0._DP,0._DP,0._DP,
     .                              0._DP,0._DP,0._DP,NSBOX)
C  TAKE VOLUMES FROM EXTERNAL FILE
        ELSEIF (INDPRO(12).EQ.6) THEN
          CALL EIRENE_PROFR (VOL,5+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
          IF (NLADD) CALL EIRENE_VOLUME(4)
        ELSEIF (INDPRO(12).EQ.7) THEN
          CALL EIRENE_PROFR (VOL,5+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        ENDIF
C
C  MULTIPLY GEOMETRICAL DATA, IF NBMLT.GT.1
C
        IF (NBMLT.GT.1) CALL EIRENE_MULTIG
C
C  SET CELL DIAMETER
C
        CALL EIRENE_SET_CELL_DIAMETER
C     
C   INCLUDE INFORMATION PROVIDED BY INPUT BLOCK 8: ADDITIONAL
C   DATA FOR SPECIFIC ZONES
C
        DO WHILE (ASSOCIATED(VOLLIST))
          VOL(VOLLIST%IN) = VOLLIST%VOL
          VOLCUR => VOLLIST
          VOLLIST => VOLLIST%NEXT
          DEALLOCATE(VOLCUR)
        ENDDO
        IF (NPHOTI > 0) CALL EIRENE_PH_INIT(2)
C
C   MODIFY SOME GEOMETRICAL DATA, USER-SUPPLIED ROUTINE
C
        CALL EIRENE_GEOUSR

        IF (LEVGEO == 4) CALL EIRENE_CUT_ADS_CELL
C
C  SET SOME DATA FOR ADDITIONAL SURFACES: INITIALISE SUBR. TIMEA
C
cdr     timea = EIRENE_SECOND_OWN()
        CALL EIRENE_TIMEA0
cdr     WRITE(iunout,*)'cpu time for timea0 ',EIRENE_SECOND_OWN()-timea
cdr     WRITE(iunout,*)
C
C   MODIFY THE BOUNDARIES OF SOME SURFACES TO AVOID ROUND-OFF
C   ERRORS
C
        CALL EIRENE_SETFIT (TRCSUR)
C
C   SET THE COEFFICIENTS OF SOME SURFACES IDENTICAL TO THOSE
C   OF SOME OTHER, TO AVOID ROUND-OFF ERRORS
C
        CALL EIRENE_SETEQ
C
C   WRITE A LIST OF CLOSED POLYGONIAL LINES, DETERMINED
C   FROM EIRENE ADDITIONAL AND NON-DEFAULT STANDARD SURFACES
C   (WITH THEIR ORIENTATION) ONTO STREAM 78+IFOFF
C   FOR FURTHER USE IN TRIANGULARISATION CODES, WHICH MAY THEN
C   PRODUCE GRIDS OF UNSTRUCTURED TRIANGLES INSIDE THESE CLOSED
C   POLYGONS (EXCLUDING THOSE AREAS WHICH ARE DESCRIBED BY
C   CLOSED POLYGONS WITH NEGATIVE ORIENTATION)
C
        IF (NLWRMSH) THEN
          CALL EIRENE_WRMESH
          CALL EIRENE_PLMESH
        ENDIF
C
C
        CALL EIRENE_INTVOL (VOL,1,1,NSBOX,VOLTOT,
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        WRITE (iunout,*) 'TOTAL VOLUME, SUM VOL(:)  ',VOLTOT
C
C  SET 'VISIBLE ADDITIONAL SURFACES' RANGES nlimii(j),nlimie(j), for each grid cell j
C  FROM INFORMATION ON IGJUM3
        CALL EIRENE_SETUP_VIS_ADD_SURF_RANGES
        
C
C  ALL GEOMETRICAL DATA (GRIDS, VOLUMES, SWITCHES) ARE DEFINED NOW
C
C   SAVE GEOMETRICAL DATA ON FILE FT12
C
        DO 8006 IRAD=1,NSBOX
          VOLG(IRAD)=VOL(IRAD)
 8006   CONTINUE
C
        DO ILIMPS=1,NLMPGS
          AREAG(ILIMPS)=SAREA(ILIMPS)
        ENDDO
C
        IF (NFILEM.EQ.1) CALL EIRENE_WRGEOM(TRCFLE)
C
      ELSEIF (NFILEM.EQ.2) THEN
C
C   RESTORE GEOMETRICAL DATA FROM FILE FT12
C
        CALL EIRENE_RGEOM(TRCFLE)
C
        DO 8010 IRAD=1,NSBOX
          VOL(IRAD)=VOLG(IRAD)
 8010   CONTINUE
C
        DO ILIMPS=1,NLMPGS
          SAREA(ILIMPS)=AREAG(ILIMPS)
        ENDDO
C
C
      ENDIF
C
      DO 8011 IS=1,NLIMPS
        IF (ILACLL(IS).NE.0.AND..NOT.NLADD) THEN
          WRITE (iunout,*) 'ADDITIONAL CELL SWITCHES DEFINED, BUT NO'
          WRITE (iunout,*) 'ADDITIONAL CELLS DEFINED'
          ISS=IS
          IF (ISS.GT.NLIM) ISS=-(ISS-NLIM)
          WRITE (iunout,*) 'SURFACE NO. IS= ',ISS
          WRITE (iunout,*) 'CHECK INPUT BLOCK 2E'
          CALL EIRENE_EXIT_OWN(1)
        ELSEIF (ILBLCK(IS).NE.0.AND..NOT.(NLMLT.OR.NLADD)) THEN
          WRITE (iunout,*) 'BLOCK SWITCHES DEFINED, BUT NEITHER BLOCKS'
          WRITE (iunout,*) 'NOR ADDITIONAL CELLS DEFINED'
          ISS=IS
          IF (ISS.GT.NLIM) ISS=-(ISS-NLIM)
          WRITE (iunout,*) 'SURFACE NO. IS= ',ISS
          WRITE (iunout,*) 'CHECK INPUT BLOCKS 2D AND 2E'
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
 8011 CONTINUE

C
 4000 CONTINUE
C

!  NOTHING IS DONE IF ARRAYS FOR BACKGROUND ARE ALREADY ALLOCATED
      IF (ANY(INDPRO(1:12) == 6)) CALL EIRENE_ALLOC_BCKGRND


      IF ((NMODE.NE.0.AND.IITER.LE.MAX(1,NITER0)) .OR.
     .    (ABS(NMODE).EQ.2)) THEN
C  READ PLASMA BACKGROUND
c  EITHER:  FROM EXTERNAL DATABASE (FT31) (NOT NLPLAS)
C  OR    :  FROM COMMON BRAEIR (NLPLAS)
        IF (ANY(INDPRO(1:12) == 6)) CALL EIRENE_ALLOC_BCKGRND
        CALL EIRENE_IF1COP
      ENDIF

C
      IF (NSTEP > 0) CALL EIRENE_ALLOC_CSTEP
C
cdr  VOL are cell volumes on the fine grid
cdr  coarser grid FOR SCORING may have been set, find cell volumes
cdr  VOLTAL on coarser grid
      VOLTAL = EPS60
      DO IN=1,NSBOX
        INC = NCLTAL(IN)
        IF (INC > 0) VOLTAL(INC) = VOLTAL(INC) + VOL(IN)
      END DO
      CALL EIRENE_INTVOL (VOLTAL,1,1,NSBOX_TAL,VOLTOT_TAL,
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      WRITE (iunout,*) 'TOTAL VOLUME, SUM VOLTAL(:) ',VOLTOT_TAL


C
cpb add nlshrt13
      IF ((NFILEL.LE.1) .OR. NLSHRT13) THEN
C
C  SET PLASMA PARAMETERS AND SOURCE PARAMETERS
C
        CALL EIRENE_PLASMA


C
C  MULTIPLY PLASMA PARAMETERS, IF NBLCKS.GT.1
C
        IF (NBLCKS.GT.1) CALL EIRENE_MULTIP
C
C   INCLUDE INFORMATION PROVIDED BY INPUT BLOCK 8: ADDITIONAL
C   DATA FOR SPECIFIC ZONES
C
        IF ((NZADD.GT.0) .AND. (IITER == 1) .AND. (ITIMV == 1))
     .    CALL EIRENE_SET_SPECIFIC_ZONES
C     
C  MODIFY SOME PLASMA DATA, USER-SUPPLIED ROUTINE
C
        CALL EIRENE_PLAUSR


C  THIS POINT IS ONLY REACHED WITH NFILEL=3 IF NLSHRT13=.T.
c  (SHORT VERSION OF FORT.13 ONLY).
C  READ ONLY PLASMA DATA OF THOSE BACKGROUND SPECIES
C  WHICH ARE NOT CONTAINED IN EXTERNAL PLASMA CODE,
C  I.E. ONLY THOSE WHICH ARE NEEDED FOR INTERNAL EIRENE CYCLING (NON-LINEARITIES)
        IF (NFILEL.EQ.3) CALL EIRENE_RPLAM(TRCFLE,0,IRET)

C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE INPUT PROFILES
C
        CALL EIRENE_PLASMA_DERIV(0)


C
C  SET ATOMIC DATA TABLES
C
        CALL EIRENE_SETAMD(1)


C
        IF (NFILEL.EQ.1) CALL EIRENE_WRPLAM(TRCFLE,0)


C
      ELSEIF (NFILEL.GE.2.AND.NFILEL.LE.4) THEN
C
C  READ PLASMA DATA, ATOMIC DATA, SOURCE DATA FROM FT13
C
        NLSRON_SAVE = NLSRON

        CALL EIRENE_ALLOC_CORNERS

        IFLG = 0
        IF ((ITIMV > 1) .AND. NLMOVIE) IFLG = 1
cdr  iflg rather than nfilel  ??
        IF ((NFILEL == 2) .OR. (NFILEL == 3)) THEN
          CALL EIRENE_RPLAM(TRCFLE,IFLG,IRET)

cdr  from now on: iflg=nfilel
        ELSEIF (NFILEL == 4) THEN
          CALL EIRENE_RPLAM(TRCFLE,NFILEL,IRET)
        END IF
        CALL EIRENE_XSECTPH

        IF (IITER > 1) NLSRON = NLSRON_SAVE
C
      ENDIF



C
C  SET UP TABLE OF CONTRIBUTIONS OF MONTE CARLO PARTICLES TO BACKGROUND SPECIES
C
      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)

      DO JPLS = 1, NPLSI
        IF ((LEN_TRIM(CDENMODEL(JPLS)) > 0) .AND.
     .      (INDEX(CDENMODEL(JPLS),'CONSTANT') == 0)) THEN
          DO IRE = 1, TDMPAR(JPLS)%TDM%NRE
            ITYP = TDMPAR(JPLS)%TDM%ITP(IRE)
            ISPZ = IADTYP(ITYP) + TDMPAR(JPLS)%TDM%ISP(IRE)
            ISPZ_BACK(ISPZ,JPLS) = 1
          END DO
        END IF
      END DO

      CALL EIRENE_LEER(2)
      WRITE (IUNOUT,*) ' LIST OF CONTRIBUTIONS TO BACKGROUND SPECIES'
      DO JSPZ = 1, NSPZ
        DO JPLS = 1, NPLSI
          IF (ISPZ_BACK(JSPZ,JPLS) > 0)
     .      WRITE (IUNOUT,*) TEXTS(JSPZ), ' CONTRIBUTES TO ',
     .                       TEXTS(NSPAMI+JPLS)
        END DO
      END DO
C
C  AT THIS POINT THE BACKGROUND MEDIUM DATA ARE ALL SET.
C
C  COMPUTE SOURCE DATA (OVERRULE SOME OF INPUT BLOCK 7)
!pb 22.10.2009
!pb   IF ((NMODE.NE.0.AND.IITER.LE.1) .OR. (IITER > NITER))  THEN
      IF (NMODE.NE.0.AND.((IITER.LE.1) .OR. (IITER > NITER)))  THEN
        DO ISTRA=1,NSTRAI
          IF (INDSRC(ISTRA).GE.0) CALL EIRENE_IF2COP(ISTRA)
        ENDDO
      ENDIF
C
C
      IF (TRCSUR) THEN
C
        CALL EIRENE_LEER(2)
        WRITE (iunout,*) 'COEFFICIENTS FOR ADDITIONAL SURFACES'
        WRITE (iunout,*)
     .    'THIS IS AFTER IF0COP, GEOUSR, SETEQ AND SETFIT ARE CALLED'
        DO 7701 J=1,NLIMI
          CALL EIRENE_LEER(2)
          WRITE (iunout,*) TXTSFL(J)
          CALL EIRENE_LEER(1)
          IF (IGJUM0(J) == 1) THEN
            WRITE (iunout,*) 'THIS SURFACE IS NOT DEFINED'
          ELSE
            WRITE (iunout,*) 'A0       ',A0LM(J)
            WRITE (iunout,*) 'A1,A2,A3 ',A1LM(J),A2LM(J),A3LM(J)
            WRITE (iunout,*) 'A4,A5,A6 ',A4LM(J),A5LM(J),A6LM(J)
            WRITE (iunout,*) 'A7,A8,A9 ',A7LM(J),A8LM(J),A9LM(J)
            WRITE (iunout,*) 'JUMLIM ',JUMLIM(J)
            WRITE (iunout,*) 'ISWICH(1),ISWICH(2),ISWICH(3),ISWICH(4),',
     .                  'ISWICH(5),ISWICH(6)'
            WRITE (iunout,*)  ISWICH(1,J),ISWICH(2,J),ISWICH(3,J),
     .                   ISWICH(4,J),ISWICH(5,J),ISWICH(6,J)
          ENDIF
 7701   CONTINUE
        CALL EIRENE_LEER(2)
        CALL EIRENE_MASIR2('IGJUM0 ',IGJUM0,1,1,1,1,NLIMPS)
        CALL EIRENE_LEER(1)
        IF (NLIMPB >= NLIMPS) THEN
          CALL EIRENE_MASIR2('IGJUM1 ',IGJUM1,0,NLIMPS,1,NLIMPS,NLIMPS)
        ELSE
          CALL
     .  EIRENE_MASBR2('IGJUM1 ',IGJUM1,0,NLIMPS,1,NLIMPS,NLIMPS,NBITS)
        END IF
        CALL EIRENE_LEER(1)
        IF (NLIMPB >= NLIMPS) THEN
          CALL EIRENE_MASIR2('IGJUM2 ',IGJUM2,0,NLIMPS,1,NLIMPS,NLIMPS)
        ELSE
          CALL
     .  EIRENE_MASBR2('IGJUM2 ',IGJUM2,0,NLIMPS,1,NLIMPS,NLIMPS,NBITS)
        END IF
        CALL EIRENE_LEER(1)
        NSOPT=MIN(NSBOX,NOPTIM)
        IF (NLIMPB >= NLIMPS) THEN
          CALL EIRENE_MASIR2('IGJUM3 ',IGJUM3,0,NOPTIM,1,NSOPT,NLIMPS)
        ELSE
          CALL
     .  EIRENE_MASBR2('IGJUM3 ',IGJUM3,0,NOPTIM,1,NSOPT,NLIMPS,NBITS)
        END IF
        CALL EIRENE_LEER(1)

        DO 7702 J=1,NSOPT
          CALL EIRENE_MASJ3
     .  ('J,NLIMII,NLIMIE          ',J,NLIMII(J),NLIMIE(J))
 7702   CONTINUE
C
      ENDIF

      CALL EIRENE_DEALLOC_BCKGRND
C
      IF (NPHOTI > 0) CALL EIRENE_PH_INIT(3)

!  allocate and initialise storage for trajectories
      IF (.NOT.ALLOCATED(TRAJ)) THEN
        ALLOCATE (TRAJ(NCHORI+NTRJ))

        DO JTRJ = 1, NCHORI+NTRJ
          ALLOCATE(TRAJ(JTRJ)%TRJ)
          TRAJ(JTRJ)%TRJ%NCOU_CELL = 0
          NULLIFY(TRAJ(JTRJ)%TRJ%CELLS)
        END DO
      END IF

      IF (NCHORI > 0) THEN
        IF (ANY(NLSTCHR)) CALL EIRENE_SETUP_CHORD_SPECTRA
      END IF

!  determine number of background spectra
      CALL EIRENE_NUM_BACKGROUND_SPECTRA(IADTYP)

      IF (.NOT.ALLOCATED(BACK_SPEC) .AND. (NBACK_SPEC > 0))
     .   ALLOCATE(BACK_SPEC(NBACK_SPEC))

      IUNIN = IUNIN_SAVE
      IF (IUSROUT /= 0) CLOSE(IUSROUT)
C
      RETURN
      END SUBROUTINE EIRENE_INPUT

c=======================================================================

      subroutine fix_logical_input(a,l)
cxb: * For gfortran: it does not accept empty field for logicals
c      Fill up empty spaces with 'f'
      implicit none
      character*(*) a
      integer l
      integer i,j
      j=1
      do i=1,l !{
        if(a(j:j).eq.' ') a(j:j)='f'
        j=j+1
        if(mod(i,5).eq.0) j=j+1
      end do !}
c      write (iunout,'(a,a)') 'zeile :',trim(a)   !###
      end

!========================================
