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
C  IFLG  :  only for  RPLAM:
C           = 0   do not read primary source data COMSOU
C          else   do also read data from COMSOU

      SUBROUTINE EIRENE_WRPLAM_LONG(TRCFLE,IFLG)
      USE EIRMOD_PARMMOD
      USE EIRMOD_CINIT, ONLY: FORT
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CZT1
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMXS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFLG
      LOGICAL TRCFLE
cdr, jan 2019
      INTEGER :: IO

C
      OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13+ifoff
      WRITE (13+ifoff)
C  REAL
     R           TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN,
     R           BXIN,BYIN,BZIN,BFIN,ADIN,EDRIFT,
     R           VOL,WGHT,BXPERP,BYPERP,
     R           EXIN,EYIN,EZIN,EFIN,
     R           FLXOUT,SAREA,
     R           TEINL,TIINL,DEINL,DIINL,BVIN,PARMOM,
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
     I           NPRT,ISPEZ,ISPEZI,
C  LUSR, LOGICAL
     L           LGVAC,LSMOPRO
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
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMSOU,ICMSOU,LCMSOU'
      IF (ALLOCATED(FLSTEP))
     .  WRITE (13+ifoff) FLSTEP,ELSTEP,FLTOT,ELTOT,VF,VE,
     .             QUOT,ADD,QUOTI,ADDIV,
     .             TESTEP,TISTEP,RRSTEP,VXSTEP,VYSTEP,VZSTEP,DISTEP,
     .             FESTEP,FISTEP,SHSTEP,VPSTEP,MCSTEP,
     .             IRSTEP,IPSTEP,ITSTEP,IASTEP,IBSTEP,IGSTEP,
     .             ISTUF,NSMAX,NSPSTI,NSPSTE
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: module EIRMOD_CSTEP.f'
      CLOSE (UNIT=13+ifoff)
      END
C
c...............................................................
C
      SUBROUTINE EIRENE_RPLAM_LONG(TRCFLE,IFLG)
      USE EIRMOD_PARMMOD
      USE EIRMOD_CINIT, ONLY: FORT
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CZT1
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMXS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFLG
      LOGICAL TRCFLE
cdr, jan 2019
      INTEGER :: IO

      OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13+ifoff
      READ (13+ifoff,IOSTAT=IO)
C  REAL
     R           TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN,
     R           BXIN,BYIN,BZIN,BFIN,ADIN,EDRIFT,
     R           VOL,WGHT,BXPERP,BYPERP,
     R           EXIN,EYIN,EZIN,EFIN,
     R           FLXOUT,SAREA,
     R           TEINL,TIINL,DEINL,DIINL,BVIN,PARMOM,
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
     I           NPRT,ISPEZ,ISPEZI,
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
        IF (TRCFLE) WRITE (iunout,*) 'SOURCE DATA NOT READ FROM FORT.13' 
      END IF


      CLOSE (UNIT=13+ifoff)
      RETURN

 990  CONTINUE
      WRITE (IUNOUT,*) ' ERROR READING FILE FORT.13 '
      CALL EIRENE_EXIT_OWN(1)
 991  CONTINUE
      WRITE (IUNOUT,*) ' AVAILABLE INPUT TALLIES ARE DIFFERENT FROM',
     .                 ' PRIOR JOB WHICH WROTE FORT.13 '
      CALL EIRENE_EXIT_OWN(1)
 992  CONTINUE
      WRITE (IUNOUT,*) ' LEADING DIMENSIONS OF INPUT TALLIES ARE',
     .                 ' DIFFERENT FROM',
     .                 ' PRIOR JOB WHICH WROTE FORT.13 '
      CALL EIRENE_EXIT_OWN(1)
 993  CONTINUE
      WRITE (IUNOUT,*) ' STARTING POSITIONS OF INPUT TALLIES ',
     .                 ' IN ARRAY PLSTLS ARE DIFFERENT FROM',
     .                 ' PRIOR JOB WHICH WROTE FORT.13 '
      CALL EIRENE_EXIT_OWN(1)

      END
