cdr jan. 2020: remove argument IRET. Unused.
c  jan. 2019:  lgdft removed.  read (....,IOSTAT=IO)
c  oct. 2018:  iflg=0 read primary source data (incl. stepfunctions)
c              else   no reading of primary source data 
c             (also not of  stepfunctions)
c              remove redundant logical tally LGDFT
c  feb. 2018:  restructured because of switchable input tallies
c              Tests: are the same input tallies active in read and write runs?

c  sept. 05:  five more tallies added to step function, see also CSTEP.f
c  nov.  05:  add eltot and ve to step function data

C  write plasma (background) data, source distribution and atomic data
C  on unit 13.
C
c  at subroutine RPLAM:
C  read plasma (background) data, source distribution and atomic data
C  from unit 13.
C
C  trcfle:  confirm writing on printout on unit IUNOUT
C  IFLGIN  :  only for  RPLAM:
C           /= 0   do not read primary source data COMSOU
C            = 0   do also read data from COMSOU
C            = 10  read only plasma background

      SUBROUTINE EIRENE_WRPLAM_LONG(TRCFLE,IFLGIN)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CZT1
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMXS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFLGIN
      LOGICAL TRCFLE
cdr  for testing: are identical input tallies active in write and read ?
      INTEGER :: IFLG
C
      IFLG = IFLGIN 
      OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13+ifoff
#ifdef CHECKBIN
      OPEN (UNIT=113,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 113
#endif
      
C  write those input tallies which are active in the present run
      WRITE (13+ifoff) LIVTALI
#ifdef CHECKBIN
      WRITE (113,*) 'LIVTALI '
      WRITE (113,*) LIVTALI
#endif
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: LIVTALI '
      WRITE (13+ifoff) NFRSTP, NADDP
#ifdef CHECKBIN
      WRITE (113,*) 'NFRSTP, NADDP '
      WRITE (113,*) NFRSTP, NADDP
#endif
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: NFRSTP, NADDP '
      WRITE (13+ifoff) PLSTLS
#ifdef CHECKBIN
      WRITE (113,*) 'PLSTLS '
      WRITE (113,*) PLSTLS
#endif
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: input tallies PLSTLS '
C
      WRITE (13+ifoff)
C  REAL
     R           FLXOUT,SAREA,
     R           TEINL,TIINL,DEINL,DIINL,
     R           RMASSI,RMASSA,RMASSM,RMASSP,
     R           DIOD,DATD,DMLD,DPLD,DPHD,
     R           DION,DATM,DMOL,DPLS,DPHOT,
     R           TVAC,DVAC,VVAC,ALLOC,
     R           CORNER_PROFILES,
     T           TEXTS,
C  MUSR, INTEGER
     I           NSPH  ,NPHOTI,NPHOTIM,NFOLPH,NGENPH,
     I           NSPA  ,NATMI,NATMIM,NMASSA,NCHARA,NFOLA,NGENA,
     I           NSPAM ,NMOLI,NMOLIM,NMASSM,NCHARM,NFOLM,NGENM,
     I           NSPAMI,NIONI,NIONIM,NMASSI,NCHARI,NCHRGI,NFOLI,NGENI,
     I           NSPTOT,NPLSI,NPLSIM,NMASSP,NCHARP,NCHRGP,NBITS,
     I           NSNVI,NCPVI,NADVI,NBGVI,NALVI,NCLVI,NADSI,NALSI,NAINI,
     I           NPRT,ISPEZ,ISPEZI,MPLSTI,MPLSV,
C  LUSR, LOGICAL
     L           LGVAC,LSMOPRO
#ifdef CHECKBIN
      write (113,*) 'FLXOUT, SAREA ...', nlimps, nlmpgs
      WRITE (113,*)
C  REAL
     R           FLXOUT,SAREA,
     R           TEINL,TIINL,DEINL,DIINL,
     R           RMASSI,RMASSA,RMASSM,RMASSP,
     R           DIOD,DATD,DMLD,DPLD,DPHD,
     R           DION,DATM,DMOL,DPLS,DPHOT,
     R           TVAC,DVAC,VVAC,ALLOC,
     R           CORNER_PROFILES,
     T           TEXTS,
C  MUSR, INTEGER
     I           NSPH  ,NPHOTI,NPHOTIM,NFOLPH,NGENPH,
     I           NSPA  ,NATMI,NATMIM,NMASSA,NCHARA,NFOLA,NGENA,
     I           NSPAM ,NMOLI,NMOLIM,NMASSM,NCHARM,NFOLM,NGENM,
     I           NSPAMI,NIONI,NIONIM,NMASSI,NCHARI,NCHRGI,NFOLI,NGENI,
     I           NSPTOT,NPLSI,NPLSIM,NMASSP,NCHARP,NCHRGP,NBITS,
     I           NSNVI,NCPVI,NADVI,NBGVI,NALVI,NCLVI,NADSI,NALSI,NAINI,
     I           NPRT,ISPEZ,ISPEZI,MPLSTI,MPLSV,
C  LUSR, LOGICAL
     L           LGVAC,LSMOPRO
#endif
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: module EIRMOD_COMUSR.f '
      CALL EIRENE_WRITE_CMDTA
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMDTA,ICMDTA'
      CALL EIRENE_WRITE_CMAMF
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMAMF,ICMAMF'
      WRITE (13+ifoff) RCZT1,RCZT2,ZT1,ZRG
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCZT1,RCZT2,ZT1,ZRG'

c  write primary source parameters
      WRITE (13+ifoff) RCMSOU,SREC,EIO,EEL,
     .           ICMSOU,INGRDA,INGRDE,NSTRAI,
     .           LCMSOU,NLSYMP,NLSYMT
#ifdef CHECKBIN
      WRITE (113,*) 'COMSOU'
      WRITE (113,*) RCMSOU,SREC,EIO,EEL,
     .           ICMSOU,INGRDA,INGRDE,NSTRAI,
     .           LCMSOU,NLSYMP,NLSYMT
#endif
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMSOU,ICMSOU,LCMSOU'
      IF (ALLOCATED(FLSTEP))
     .  WRITE (13+ifoff) FLSTEP,ELSTEP,FLTOT,ELTOT,VF,VE,
     .             QUOT,ADD,QUOTI,ADDIV,
     .             TESTEP,TISTEP,RRSTEP,VXSTEP,VYSTEP,VZSTEP,DISTEP,
     .             FESTEP,FISTEP,SHSTEP,VPSTEP,MCSTEP,
     .             IRSTEP,IPSTEP,ITSTEP,IASTEP,IBSTEP,IGSTEP,
     .             ISTUF,NSMAX,NSPSTI,NSPSTE
#ifdef CHECKBIN
      IF (ALLOCATED(FLSTEP)) THEN
        write (113,*) 'CSTEP NGITT ', NGITT
        WRITE (113,*) FLSTEP,ELSTEP,FLTOT,ELTOT,VF,VE,
     .             QUOT,ADD,QUOTI,ADDIV,
     .             TESTEP,TISTEP,RRSTEP,VXSTEP,VYSTEP,VZSTEP,DISTEP,
     .             FESTEP,FISTEP,SHSTEP,VPSTEP,MCSTEP,
     .             IRSTEP,IPSTEP,ITSTEP,IASTEP,IBSTEP,IGSTEP,
     .             ISTUF,NSMAX,NSPSTI,NSPSTE
      END IF
#endif
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: module EIRMOD_CSTEP.f'
      CLOSE (UNIT=13+ifoff)
#ifdef CHECKBIN
      CLOSE (UNIT=113)
#endif
      END
C
c...............................................................
C
      SUBROUTINE EIRENE_RPLAM_LONG(TRCFLE,IFLGIN,IRET)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CINIT, ONLY: FORT
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CZT1
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFLGIN
      INTEGER, INTENT(INOUT) :: IRET
      LOGICAL TRCFLE
cdr, jan 2019
      INTEGER :: NFRS(NTALI), NAD(NTALI), IO, I, IFLG
      LOGICAL :: LIVT(NTALI)
      REAL(DP), ALLOCATABLE :: PTL(:,:)

      IFLG = IFLGIN
      IRET = 0
      IF (.NOT.ALLOCATED(PTL))
     .  ALLOCATE(PTL, MOLD=PLSTLS)
      IF (IFLG == 10) CALL EIRENE_ALLOC_BCKGRND    

      OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13+ifoff

      READ (13+ifoff,IOSTAT=IO) LIVT
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: LIVTALI '
      IF (IO /= 0) THEN
        DEALLOCATE (PTL)
        IF (IFLG == 10) THEN  ! plasma_bckground = 0 in EIRENE_ALLOC_BCKGRND
          IRET = IO
          RETURN
        ELSE 
          GOTO 990
        END IF
      END IF
c   verify: same active tallies as in previous write?
      DO I=1, NTALI
        IF ((LIVTALI(I).AND.LIVT(I)).OR.
     .      (.NOT.LIVTALI(I).AND..NOT.LIVT(I))) THEN
          CYCLE
        ELSE
          GOTO 991
        END IF
      END DO
      READ (13+ifoff,IOSTAT=IO) NFRS, NAD
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: NFRSTP, NADDP '
      IF (IO /= 0) GOTO 990
      IF (ANY(NFRSTP(1:NTALI) /= NFRS(1:NTALI))) GOTO 992
      IF (ANY(NADDP(1:NTALI) /= NAD(1:NTALI))) GOTO 993

!PB      READ (13+ifoff,IOSTAT=IO) PLSTLS
      READ (13+ifoff,IOSTAT=IO) PTL
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: input tallies PLSTLS '
      IF (IO /= 0) GOTO 990
      IF (IFLG == 10) THEN
        TEINTF(1:NRAD) = PTL(NADDP(1)+1,1:NRAD)
        IF (INTLOPTS(2) >= 0) 
     .    TIINTF(1:NPLSTI,1:NRAD) = PTL(NADDP(2)+1:NADDP(3),1:NRAD)       
        DIINTF(1:NPLS,1:NRAD) = PTL(NADDP(4)+1:NADDP(5),1:NRAD)
        VXINTF(1:NPLSV,1:NRAD) = PTL(NADDP(5)+1:NADDP(6),1:NRAD)
        VYINTF(1:NPLSV,1:NRAD) = PTL(NADDP(6)+1:NADDP(7),1:NRAD)
        VZINTF(1:NPLSV,1:NRAD) = PTL(NADDP(7)+1:NADDP(8),1:NRAD)
        DEALLOCATE (PTL)
        RETURN
      ELSE
        PLSTLS = PTL
        DEALLOCATE (PTL)
      END IF

      READ (13+ifoff,IOSTAT=IO)
C  REAL
     R           FLXOUT,SAREA,
     R           TEINL,TIINL,DEINL,DIINL,
     R           RMASSI,RMASSA,RMASSM,RMASSP,
     R           DIOD,DATD,DMLD,DPLD,DPHD,
     R           DION,DATM,DMOL,DPLS,DPHOT,
     R           TVAC,DVAC,VVAC,ALLOC,
     R           CORNER_PROFILES,
     T           TEXTS,
C  MUSR, INTEGER
     I           NSPH  ,NPHOTI,NPHOTIM,NFOLPH,NGENPH,
     I           NSPA  ,NATMI,NATMIM,NMASSA,NCHARA,NFOLA,NGENA,
     I           NSPAM ,NMOLI,NMOLIM,NMASSM,NCHARM,NFOLM,NGENM,
     I           NSPAMI,NIONI,NIONIM,NMASSI,NCHARI,NCHRGI,NFOLI,NGENI,
     I           NSPTOT,NPLSI,NPLSIM,NMASSP,NCHARP,NCHRGP,NBITS,
     I           NSNVI,NCPVI,NADVI,NBGVI,NALVI,NCLVI,NADSI,NALSI,NAINI,
     I           NPRT,ISPEZ,ISPEZI,MPLSTI,MPLSV,
C  LUSR, LOGICAL
     L           LGVAC,LSMOPRO
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: module EIRMOD_COMUSR.f '
      IF (IO /= 0) GOTO 990
      CALL EIRENE_READ_CMDTA
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMDTA,ICMDTA'
      CALL EIRENE_READ_CMAMF
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMAMF,ICMAMF'

      READ (13+ifoff,IOSTAT=IO) RCZT1,RCZT2,ZT1,ZRG
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCZT1,RCZT2,ZT1,ZRG'
      IF (IO /= 0) GOTO 990

      IF (IFLG == 0) THEN
c  read primary source parameters
        READ (13+ifoff,IOSTAT=IO) RCMSOU,SREC,EIO,EEL,
     .            ICMSOU,INGRDA,INGRDE,NSTRAI,
     .            LCMSOU,NLSYMP,NLSYMT
        IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMSOU,ICMSOU,LCMSOU,...'
        IF (IO /= 0) GOTO 990
        IF (ALLOCATED(FLSTEP))
     .    READ (13+ifoff,IOSTAT=IO) FLSTEP,ELSTEP,FLTOT,ELTOT,VF,VE,
     .             QUOT,ADD,QUOTI,ADDIV,
     .             TESTEP,TISTEP,RRSTEP,VXSTEP,VYSTEP,VZSTEP,DISTEP,
     .             FESTEP,FISTEP,SHSTEP,VPSTEP,MCSTEP,
     .             IRSTEP,IPSTEP,ITSTEP,IASTEP,IBSTEP,IGSTEP,
     .             ISTUF,NSMAX,NSPSTI,NSPSTE
        IF (TRCFLE) WRITE (iunout,*) 'READ 13: module EIRMOD_CSTEP.f '
        IF (IO /= 0) GOTO 990
      ELSE
        IF (TRCFLE)
     .   WRITE (iunout,*) 'SOURCE DATA NOT READ FROM ', FORT, '13'
      END IF

      IRET = IO
      CLOSE (UNIT=13+ifoff)
      RETURN

 990  CONTINUE
      WRITE (IUNOUT,*) ' ERROR READING FILE '//FORT//'13 '
      CALL EIRENE_EXIT_OWN(1)
 991  CONTINUE
      WRITE (IUNOUT,*) ' AVAILABLE INPUT TALLIES ARE DIFFERENT FROM',
     .                 ' PRIOR JOB WHICH WROTE '//FORT//'13 '
      CALL EIRENE_EXIT_OWN(1)
 992  CONTINUE
      WRITE (IUNOUT,*) ' LEADING DIMENSIONS OF INPUT TALLIES ARE',
     .                 ' DIFFERENT FROM',
     .                 ' PRIOR JOB WHICH WROTE '//FORT//'13 '
      CALL EIRENE_EXIT_OWN(1)
 993  CONTINUE
      WRITE (IUNOUT,*) ' STARTING POSITIONS OF INPUT TALLIES ',
     .                 ' IN ARRAY PLSTLS ARE DIFFERENT FROM',
     .                 ' PRIOR JOB WHICH WROTE '//FORT//'13 '
      CALL EIRENE_EXIT_OWN(1)

      END
