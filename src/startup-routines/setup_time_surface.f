      SUBROUTINE EIRENE_SETUP_TIME_SURFACE(IERROR)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CTEXT
      USE EIRMOD_CLGIN
      USE EIRMOD_CINIT
      USE EIRMOD_COMPRT
      USE EIRMOD_COMSOU
      USE EIRMOD_CCONA
      USE EIRMOD_COMNNL
      USE EIRMOD_CCONA
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: IERROR
      INTEGER :: I, ISTR_A(1), ISTR, IPNT
      REAL(DP) :: DTIMVO
      
C
      IF (NTIME.GE.1) THEN
C  DEFINE ONE MORE SURFACE
        NSTSI=NSTSI+1
C  CHECK STORAGE
        CALL EIRENE_LEER(1)
        IF (NSTSI.GT.NSTS) THEN
          CALL EIRENE_MASPRM('NSTS',4,NSTS,'NSTSI',5,NSTSI,IERROR)
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
C  SET DEFAULTS FOR "TIME HORIZON" SURFACE (ABSORBING)
C
        TXTSFL(NLIM+NSTSI)='"TIME HORIZON"                           '
        ILIIN(NLIM+NSTSI)=2
C
C
cdr     IF (NFILEJ.EQ.2.OR.NFILEJ.EQ.3) THEN
cdr functioniert noch nicht, falls mehrere timesteps, davon nur
cdr der erste: initialisierung, die anderen: fortsetzung.
cdr denn dann wird bei der fortsetzung das stratum nicht gemacht.
cdr wg. goto 4000. angefangen: "mkcens" (make stratum for census array)
C
C  DEFINE ONE MORE STRATUM
        NSTRAI=NSTRAI+1
C  CHECK STORAGE
        IF (NSTRAI.GT.NSTRA) THEN
          CALL EIRENE_MASPRM('NSTRA',5,NSTRA,'NSTRAI',6,NSTRAI,IERROR)
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
C  SET DEFAULTS FOR SOURCE DUE TO INITIAL CONDITION, VALID ONLY FOR
C  FIRST TIMESTEP. MODIFIED FOR LATER TIMESTEPS IN SUBR. TMSTEP
C
        TXTSOU(NSTRAI)='SOURCE DUE TO INITIAL CONDITION          '
C  SOURCE DISTRIBUTION SAMPLED FROM CENSUS ARRAYS RPARTC,IPARTC
        NLCNS(NSTRAI)=.TRUE.
        NLPNT(NSTRAI)=.FALSE.
        NLLNE(NSTRAI)=.FALSE.
        NLSRF(NSTRAI)=.FALSE.
        NLVOL(NSTRAI)=.FALSE.
C  DO NOT CALL IF2COP(NSTRAI)
        INDSRC(NSTRAI)=-1
C
        NLAVRP(NSTRAI)=.FALSE.
        NLAVRT(NSTRAI)=.FALSE.
        NLSYMP(NSTRAI)=.FALSE.
        NLSYMT(NSTRAI)=.FALSE.
        NPTS(NSTRAI)=0
        NINITL(NSTRAI)=2000*NINITL(NSTRAI-1)+1
        IF (NINITL_READ /= 0) NINITL(NSTRAI)=NINITL_READ
        NEMODS(NSTRAI)=1
        NAMODS(NSTRAI)=1
        FLUX(NSTRAI)=0.
        NLATM(NSTRAI)=.FALSE.
        NLMOL(NSTRAI)=.FALSE.
        NLION(NSTRAI)=.FALSE.
        NLPLS(NSTRAI)=.FALSE.
        NSPEZ(NSTRAI)=0
        NSRFSI(NSTRAI)=0
C     
        SORENI(NSTRAI)=0.
        SORENE(NSTRAI)=0.
        SORVDX(NSTRAI)=0.
        SORVDY(NSTRAI)=0.
        SORVDZ(NSTRAI)=0.
        SORCOS(NSTRAI)=0.
        SORMAX(NSTRAI)=0.
        SORCTX(NSTRAI)=0.
        SORCTY(NSTRAI)=0.
        SORCTZ(NSTRAI)=0.
C
C  NEW TIMESTEP
C
        IF (DTIMVN.LE.0.D0) THEN
          DTIMVN=DTIMV
C       ELSE
C         DTIMVN=DTIMVN
        ENDIF
C
C  OLD TIMESTEP
C
C  READ INITIAL POPULATION FROM FILE, FORT.15, OVERWRITE DEFAULTS
C
        IPRNL=0
        IF (NFILEJ.EQ.2.OR.NFILEJ.EQ.3) THEN
          CALL EIRENE_RSNAP(NSTRAI)
          DTIMVO=DTIMV
C
          WRITE (iunout,*) 'INITIAL POPULATION FOR FIRST TIMESTEP'
          WRITE (iunout,*) 'READ FROM FILE ', FORT, '15'
          WRITE (iunout,*) 'PARTICLES AND FLUX RETRIEVED FOR'
          WRITE (iunout,*) 'INITIAL DISTRIBUTION AT T0= ',TIME0
          CALL EIRENE_MASJ1('IPRNL   ',IPRNL)
          CALL EIRENE_MASR1('FLUX    ',FLUX(NSTRAI))
C
          IF (DTIMVN.NE.DTIMVO) THEN
            FLUX(NSTRAI)=FLUX(NSTRAI)*DTIMVO/DTIMVN
C
            WRITE (iunout,*) 'FLUX IS RESCALED BY DTIMV_OLD/DTIMV_NEW'
            CALL EIRENE_MASR1('FLUX    ',FLUX(NSTRAI))
            CALL EIRENE_LEER(1)
          ENDIF
C
          CALL EIRENE_LEER(2)
          IF (TIME0.GE.0.) THEN
c    reset clock of source particles from old census to time0
cdr  must be done also for time0=0.0, for otherwise flight time =0 is possible
cdr  for census source particles and resulting error exits
C Would gain performance by turning RPSTT into a pointer
cpb  changed due to optimizer problem            
cpb  try to determine the index of the element of RPSTT which TIME points to
            RPSTT(1:NPARTT)=RPARTC(1:NPARTT,1)
            TIME=HUGE(1._DP)
            IPNT = MAXLOC(RPSTT,1)
cpb  reset time stamp
            IF (RPSTT(IPNT) > EPS30) THEN
              RPARTC(IPNT,1:IPRNL) = TIME0
              TIME=TIME0
            ELSE
              DO I=1,IPRNL
                RPSTT(1:NPARTT)=RPARTC(1:NPARTT,I)
                TIME=TIME0
                RPARTC(1:NPARTT,I)=RPSTT(1:NPARTT)
              ENDDO
            END IF
            WRITE (iunout,*) 'PARTICLE CLOCK RESET'
            WRITE (iunout,*) 'FIRST TIMESTEP RUNS FROM TIM1 TO TIM2:'
            CALL EIRENE_MASR2('TIM1, TIM2      ',TIME0,TIME0+DTIMV)
            CALL EIRENE_LEER(2)
          ENDIF
C
        ENDIF
C
        DTIMV=DTIMVN
C
        IF (NPTST.EQ.0) THEN
C  SET NUMBER OF PARTICLES FOR RELAUNCH FROM CENSUS EQUAL TO THE NUMBER OF PREVIOUS SCORES ON CENSUS
C  (BUT STILL: SAMPLING WITH REPLACEMENT, BOOTSTRAPPING)
C  N.B.: WE HAVE ALREADY MADE SURE ABOVE, THAT IN CASE NLMOVIE: NPTST = -1
          NPTS(NSTRAI)=IPRNL
          NMINPTS(NSTRAI)=IPRNL   ! NMINPTS is currently not used anywhere
        ELSEIF (NPTST.GT.0) THEN
          NPTS(NSTRAI)=NPTST
        ELSEIF (NPTST.LT.0) THEN
C  ONE BY ONE RELAUNCH FROM OLD CENSUS
C  OLD CENSUS CONTAINS IPRNL ENTRIES.
          NPTS(NSTRAI)=IPRNL
          NMINPTS(NSTRAI)=IPRNL   ! NMINPTS is currently not used anywhere
        ENDIF

        CALL EIRENE_MASJ1('NPTS=    ',NPTS(NSTRAI))
C
        IF (NPTS(NSTRAI).GT.0.AND.FLUX(NSTRAI).GT.0) THEN
          NSRFSI(NSTRAI)=1
          SORWGT(1,NSTRAI)=1.D0
        ENDIF
C
      ELSEIF (NFILEJ.EQ.2.OR.NFILEJ.EQ.3) THEN
!pb        IF ( SIZE( PACK((/ (i, i = 1, NSTRA) /),NLCNS) ) == 1 ) THEN
        IF ( COUNT(NLCNS(1:NSTRA)) == 1 ) THEN
C Only read census from file if exactly one stratum is a census stratum
          ISTR_A = PACK((/ (i, I = 1, NSTRA) /),NLCNS)
          ISTR = ISTR_A(1)
          CALL EIRENE_RSNAP( ISTR )
C
          WRITE (iunout,*) 'INITIAL POPULATION READ FROM FILE ', FORT, 
     .                     '15'
          CALL EIRENE_MASJ1('IPRNL   ',IPRNL)
          CALL EIRENE_MASR1('FLUX    ',FLUX(ISTR))
C
C  ONE BY ONE RELAUNCH FROM OLD CENSUS
C  OLD CENSUS CONTAINS IPRNL ENTRIES.
          NPTS(ISTR)=IPRNL
          NMINPTS(ISTR)=IPRNL   ! NMINPTS is currently not used anywhere
          CALL EIRENE_MASJ1('NPTS=    ',NPTS(ISTR))

          IF (NPTS(ISTR).GT.0.AND.FLUX(ISTR).GT.0) THEN
            NSRFSI(ISTR)=1
            SORWGT(1,ISTR)=1.D0
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE EIRENE_SETUP_TIME_SURFACE
