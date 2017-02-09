CDR  OCT.14 ADDED:  READ ESBPARM  (PROJECTILE SURFACE BINDING ENERGY) FROM TRIM FILES.  
c                   CURRENTLY NOT IN USE.
C                   CURRENTLY ALSO NOT YET READ (TO BE DONE): 
c                             FIND END OF LINE AND READ DATA ONLY IF AVAILABLE.
cdr  Jan 17      :  started: read inr (resolution on data file), rather than fixed inr=5
c                   tbd:  read   dummy=inr from 1st file. 
c                         (DUMMY was reserved for sputter data in TRIM format?)
c
      SUBROUTINE EIRENE_RDTRIM
C
C  THIS SUBROUTINE READS SELECTIVELY INDIVIDUAL
C  REFLECTION DATA FILES "A_ON_B" PRODUCED E.G. BY MONTE CARLO BCA CODES
C  THERE ARE NFLR SUCH FILES IN THIS RUN
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CREF
      USE EIRMOD_CSPEI
 
      IMPLICIT NONE
 
      REAL(DP) :: PID180, DUMMY, ESBPARM  ! ESBPARM SHOULD BE ARRAY(IFILE)
      INTEGER :: I1, I2, I3, I4, I5, IUN, IFILE, I, IWWW, ITTT, INR2
C
C
      INE=12
      INW=7
      INR=5  !  this value should now come from data file itself.
C
      IF (INE*INW*NFLR.GT.NH0 .OR.
     .    INE*INW*INR*NFLR.GT.NH1  .OR.
     .    INE*INW*INR*INR*NFLR.GT.NH2  .OR.
     .    INE*INW*INR*INR*INR*NFLR.GT.NH3) THEN
        WRITE (iunout,*)
     .    'ERROR IN PARAMETER STATEMENT FOR REFLECTION DATA'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      DO 7 IFILE=1,NFLR
!pb        IUN=20
        IUN=20+ifoff
        OPEN (UNIT=IUN,FILE=REFFIL(IFILE),ACCESS='SEQUENTIAL',
     .        FORM='FORMATTED')
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
C  READ: INCIDENT ENERGY, ANGLE
C  READ: REFLECTION PROBABILITY HFTR0, ???, PROJECTILE SURFACE BINDING ENERGY PARAMETER
            IF (I1.EQ.1.AND.I2.EQ.1) THEN    
cdr
!  additional data DUMMY and ESBPARM only for the first of the 12*7=84 datasets
!  SOME COMPILERS DON'T LIKE READING MORE DATA THAN THERE ARE IN A SINGLE LINE
!  to be done: check length of input line, and decide then whether to read
!  DUMMY and ESBPARM, or not.
              READ (IUN,*) TC(IFILE),TM(IFILE),WC(IFILE),WM(IFILE),
     .                     enar(i1),wiar(i2),HFTR0(I1,I2,IFILE) 
C    .                    ,DUMMY, ESBPARM,  !2 NEW PARAMETERS, MAYBE ONLY IN FIRST OF THE 84 BLOCKS ??  IF AT ALL?

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
5           CONTINUE
            READ (IUN,*)
            DO 6 I3=1,INR
            DO 6 I4=1,INR
              READ (IUN,*) (HFTR3(I1,I2,I3,I4,I5,IFILE),I5=1,INR)
6           CONTINUE
3         CONTINUE
2       CONTINUE
        CLOSE (UNIT=IUN)
7     CONTINUE
C
      INEM=INE-1
      DO 11 I=1,INEM
11      DENAR(I)=1./(ENAR(I+1)-ENAR(I))
      PID180=ATAN(1.)/45.
      DO 12 I=1,INW
12      WIAR(I)=COS(WIAR(I)*PID180)
      INWM=INW-1
      DO 13 I=1,INWM
13      DWIAR(I)=1./(WIAR(I+1)-WIAR(I))
      INRM=INR-1
cdr  old version: hard wired INR=5
      RAAR(1)=0.1
      RAAR(2)=0.3
      RAAR(3)=0.5
      RAAR(4)=0.7
      RAAR(5)=0.9
cdr  new version (not ready, allow higher resolution "INR" in quantile data tables)
      INR2=2*INR
      do I=1,INR
        raar(I)=float(1+2*(I-1))/float(inr2)
      enddo

      DO 15 I=1,INRM
15      DRAAR(I)=1./(RAAR(I+1)-RAAR(I))
C
      RETURN
      END
