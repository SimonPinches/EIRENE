CDR  OCT.14 ADDED: READ ESBPARM (PROJECTILE SURFACE BINDING ENERGY) FROM TRIM FILES.
c                  CURRENTLY NOT IN USE.
C                  CURRENTLY ALSO NOT YET READ (TO BE DONE):
c                            FIND END OF LINE AND READ DATA ONLY IF AVAILABLE.
cdr  Jan 17      : started: read INR (resolution on data file), rather than fixed INR=5
c                  tbd:     read DUMMY=INR from 1st of the 84 files.
C
C                         DEFAULT              : DUMMY=5
C                         SOME FILES READY WITH: DUMMY=10
c                         (DUMMY was reserved for sputter data in TRIM format?)
cdr  june 17:  remove variable NFLR, was same as NHD6
c
      SUBROUTINE EIRENE_RDTRIM
C
C  THIS SUBROUTINE READS, AS EXPLICITLY SELECTED, SINGLE
C  REFLECTION DATA FILES "A_ON_B" PRODUCED E.G. BY MONTE CARLO BCA CODES
C  THERE ARE NHD6 SUCH FILES IN THIS RUN
C  INPUT
C     STREAM: IUN=20+IFOFF
C     NHD6  : FROM PARMMOD
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CREF
      USE EIRMOD_CSPEI
      USE EIRMOD_CINIT, ONLY: MASTER_PATH

      IMPLICIT NONE

      REAL(DP) :: PID180
C     REAL(DP) :: DUMMY, ESBPARM  ! ESBPARM SHOULD BE ARRAY(IFILE)
      INTEGER :: I1, I2, I3, I4, I5, IUN, IFILE, I,
     .           IWWW, ITTT, INR2, IFLR, IERR
C
C
      IFLR=NHD6
      INE=12
      INW=7
      INR=5  !  this value should now come from data file itself,
             !  variable :dummy.
C
      IF (INE*INW*IFLR.GT.NH0 .OR.
     .    INE*INW*INR*IFLR.GT.NH1  .OR.
     .    INE*INW*INR*INR*IFLR.GT.NH2  .OR.
     .    INE*INW*INR*INR*INR*IFLR.GT.NH3) THEN
        WRITE (iunout,*)
     .    'ERROR IN PARAMETER STATEMENT FOR REFLECTION DATA'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      DO 7 IFILE=1,IFLR

        IUN=20+ifoff
        OPEN (UNIT=IUN,FILE=REFFIL(IFILE),ACCESS='SEQUENTIAL',
     .        FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
        if ( ierr /= 0 ) then
!         If open fails, try again interpreting the path as relative
!         to MASTER_PATH
          OPEN (UNIT=IUN,FILE=trim(MASTER_PATH)//'/'//REFFIL(IFILE),
     .          ACCESS='SEQUENTIAL',
     .          FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
          if ( ierr /= 0 ) then
            WRITE (iunout,'(a,a)')
     .       'RDTRIM.F: UNABLE TO FIND TRIM FILE '//TRIM(REFFIL(IFILE))
            CALL EIRENE_EXIT_OWN(1)
          endif
        endif
C
        READ (IUN,*)  !  READ 1 LINE: name (HEADER) of trim file: here
        READ (IUN,*)  !  comment on file
        READ (IUN,*)  !  comment on file
        READ (IUN,*)  !  comment on file
        DO 2 I1=1,INE
          DO 3 I2=1,INW
            READ (IUN,*)  !comment on file per incident angle and energy
            READ (IUN,*)  !comment on file per incident angle and energy
C  READ: PROJECTILE CHARGE AND MASS TC,TM
C  READ: WALL (TARGET) CHARGE AND MASS WC,WM
C  READ: INCIDENT ENERGY, ANGLE ENAR,WIAR
C  READ: REFLECTION PROBABILITY HFTR0, ???, PROJECTILE SURFACE BINDING ENERGY PARAMETER
            IF (I1.EQ.1.AND.I2.EQ.1) THEN
cdr
!  additional data DUMMY and ESBPARM only for the first of the 12*7=84 datasets
!  SOME COMPILERS DO NOT LIKE READING MORE DATA THAN THERE ARE IN A SINGLE LINE
!  to be done: check length of input line, and decide then whether to read
!  DUMMY and ESBPARM, or not.
              READ (IUN,*) TC(IFILE),TM(IFILE),WC(IFILE),WM(IFILE),
     .                     enar(i1),wiar(i2),HFTR0(I1,I2,IFILE)
!2 NEW PARAMETERS, MAYBE ONLY IN FIRST OF THE 84 BLOCKS i1=i2=1??  IF AT ALL?
C    .                    ,DUMMY, ESBPARM(IFILE))

            ELSE
              READ (IUN,*) TC(IFILE),TM(IFILE),WC(IFILE),WM(IFILE),
     .                     enar(i1),wiar(i2),HFTR0(I1,I2,IFILE)
            ENDIF
C  FIND NEAREST INTEGER FOR NUCLEAR MASS NUMBER OF WALL MATERIAL
            IWWW=NINT(WM(IFILE))
            WM(IFILE)=IWWW
C  FIND NEAREST INTEGER FOR NUCLEAR MASS NUMBER OF PROJECTILE
            ITTT=NINT(TM(IFILE))
            TM(IFILE)=ITTT
            READ (IUN,*)
            READ (IUN,*) (HFTR1(I1,I2,I3,IFILE),I3=1,INR)
            READ (IUN,*)
            DO 5 I3=1,INR
              READ (IUN,*) (HFTR2(I1,I2,I3,I4,IFILE),I4=1,INR)
    5       CONTINUE
            READ (IUN,*)
            DO I3=1,INR
             DO I4=1,INR
              READ (IUN,*) (HFTR3(I1,I2,I3,I4,I5,IFILE),I5=1,INR)
             END DO
            END DO
    3     CONTINUE
    2   CONTINUE
        CLOSE (UNIT=IUN)
    7 CONTINUE
C
      INEM=INE-1
      DO I=1,INEM
        DENAR(I)=1./(ENAR(I+1)-ENAR(I))
      END DO
      PID180=ATAN(1.)/45.
      DO I=1,INW
        WIAR(I)=COS(WIAR(I)*PID180)
      END DO
      INWM=INW-1
      DO I=1,INWM
        DWIAR(I)=1./(WIAR(I+1)-WIAR(I))
      END DO
      INRM=INR-1
cdr  old version: hard-wired INR=5
      RAAR(1)=0.1
      RAAR(2)=0.3
      RAAR(3)=0.5
      RAAR(4)=0.7
      RAAR(5)=0.9
cdr  new version (not ready, allow higher resolution "INR" in quantile data tables)
cdr  for INR=5: should produce the same RAAR as above.
      INR2=2*INR
      do I=1,INR
        raar(I)=float(1+2*(I-1))/float(inr2)
      enddo

      DO I=1,INRM
        DRAAR(I)=1./(RAAR(I+1)-RAAR(I))
      END DO
C
      RETURN
      END SUBROUTINE EIRENE_RDTRIM
