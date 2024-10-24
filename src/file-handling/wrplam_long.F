cdr jan. 2020: remove argument IRET. Unused.
c  jan. 2019:  lgdft removed. read (....,IOSTAT=IO)
c  oct. 2018:  iflg=0 read primary source data (incl. stepfunctions)
c              else   no reading of primary source data
c             (also not of stepfunctions)
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
C  from unit 13. And: allocate target: plasma_bckgrnd
C
C  trcfle: confirm writing on printout on unit IUNOUT

      SUBROUTINE EIRENE_WRPLAM_LONG(TRCFLE,CALLEDFROM)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CZT1
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMXS
      IMPLICIT NONE

      LOGICAL TRCFLE
      CHARACTER(*) CALLEDFROM
      EXTERNAL EIRENE_LEER
C
      OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13+ifoff
      write (iunout,*) 'WRPLAM_LONG called from ',calledfrom
      CALL EIRENE_LEER(1)
#ifdef CHECKBIN
      OPEN (UNIT=113,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 113
#endif

C  write those input tallies which are active in the present run
      WRITE (13+ifoff) LIVTALI
#ifdef CHECKBIN
      WRITE (113,*) 'LIVTALI'
      WRITE (113,*) LIVTALI
#endif
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: LIVTALI'
      WRITE (13+ifoff) NFRSTP, NADDP
#ifdef CHECKBIN
      WRITE (113,*) 'NFRSTP, NADDP'
      WRITE (113,*) NFRSTP, NADDP
#endif
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: NFRSTP, NADDP'
      WRITE (13+ifoff) PLSTLS
#ifdef CHECKBIN
      WRITE (113,*) 'PLSTLS'
      WRITE (113,*) PLSTLS
#endif
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: input tallies PLSTLS'
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
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: module EIRMOD_COMUSR.f'

cdr  Write A&M data. The routines WRITE_CMDTA, WRITE_CMAMF are contained
cdr  in module COMXS
      CALL EIRENE_WRITE_CMDTA
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMDTA,ICMDTA'
      CALL EIRENE_WRITE_CMAMF
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMAMF,ICMAMF'

      WRITE (13+ifoff) RCZT1,RCZT2,ZT1,ZRG
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCZT1,RCZT2,ZT1,ZRG'

c  write primary source parameters, step functions, etc.
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
      RETURN
      END SUBROUTINE EIRENE_WRPLAM_LONG
C
c...............................................................
C
      SUBROUTINE EIRENE_RPLAM_LONG(TRCFLE,IFLG,CALLEDFROM)
cdr  feb 20: commented, and: removed argument IRET ==  IO
cdr
C  IFLG  :  only for  RPLAM:
C            = 0   read background data and profiles (COMUSR)
C                  read A&M data (COMXS)
C                  read primary source data (COMSOU)
C            > 0   same, but do not read primary source data COMSOU
C            < 0   ...
cdr iflag=10:  something special, called from
cdr
cdr
cdr
C            = 10  read only plasma background (PLSTLS part of COMUSR)
C                  do not read atomic rates COMXS,
C                  do not read primary source data COMSOU
C                  do not set PLSTLS target,
C                  but
C                  set pointers DIINTF,.... for PLASMA_BCKGRND data structure.

cdr  unclear coding:
cdr  Reading COMXS, CSTEP, COMSOU and setting pointers DIINTF excludes each other ! WHY?
cdr  In former versions: not reading of COMSOU, COMXS,
cdr                                    although we are in ..._long.f
cdr  Also strange: INTLOPTS(2):  what is special about that? use livtali instead?

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CINIT, ONLY: FORT, FORT_LC
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CZT1
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: IFLG
      LOGICAL TRCFLE
      CHARACTER(*) CALLEDFROM

cdr, jan 2019
      INTEGER :: NFRS(NTALI), NAD(NTALI), IO, I
      LOGICAL :: LIVT(NTALI)
      REAL(DP), ALLOCATABLE :: PTL(:,:)
      CHARACTER*8 FILENAME
      CHARACTER*3 FILENUMBER
      LOGICAL EX
      EXTERNAL EIRENE_LEER, EIRENE_EXIT_OWN

      IF (.NOT.ALLOCATED(PTL))
     . ALLOCATE(PTL(size(PLSTLS,1),size(PLSTLS,2)))
      if (13+ifoff.ge.100) then
        WRITE(FILENUMBER,'(I3)') 13+ifoff
      else
        WRITE(FILENUMBER,'(I2)') 13+ifoff
      end if
      FILENAME = FORT_LC//trim(FILENUMBER)
      INQUIRE(file=trim(FILENAME),exist=EX)
      IF (IFLG == 10 .or. .not.EX) CALL EIRENE_ALLOC_BCKGRND
      IF (.not.EX) THEN
        DEALLOCATE(PTL)
        RETURN
      END IF

      OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13+ifoff

      write (iunout,*) 'RPLAM_LONG called from ',calledfrom
      CALL EIRENE_LEER(1)

      READ (13+ifoff,IOSTAT=IO) LIVT
      IF (IO /= 0) THEN
        GOTO 990
      END IF
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: LIVTALI'
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
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: NFRSTP, NADDP'
      IF (IO /= 0) GOTO 990
      IF (ANY(NFRSTP(1:NTALI) /= NFRS(1:NTALI))) GOTO 992
      IF (ANY(NADDP(1:NTALI) /= NAD(1:NTALI))) GOTO 993

cdr next: read background input "plasma" tallies PLSTLS
      IF (.NOT.ALLOCATED(PTL))
     . ALLOCATE(PTL(size(PLSTLS,1),size(PLSTLS,2)))
      READ (13+ifoff,IOSTAT=IO) PTL
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: input tallies PLSTLS'
      IF (IO /= 0) THEN
        DEALLOCATE (PTL)
        GOTO 990
      ENDIF
      IF (IFLG == 10) THEN
cdr set pointers (only for tallies of nlshrt13 but in rplm_long?)
        CALL EIRENE_ALLOC_BCKGRND
        TEINTF(1:NRAD) = PTL(NADDP(1)+1,1:NRAD)
cdr to be checked. intlopts should not be used after input.
cdr only livtali and lsmopro.
cdr further: intlopts(2)< 0 is prohibited anyway?
        IF (INTLOPTS(2) >= 0)
     .  TIINTF(1:NPLSTI,1:NRAD) = PTL(NADDP(2)+1:NADDP(3),1:NRAD)
cdr     DEINTF...
        DIINTF(1:NPLS,1:NRAD)   = PTL(NADDP(4)+1:NADDP(5),1:NRAD)
        VXINTF(1:NPLSV,1:NRAD)  = PTL(NADDP(5)+1:NADDP(6),1:NRAD)
        VYINTF(1:NPLSV,1:NRAD)  = PTL(NADDP(6)+1:NADDP(7),1:NRAD)
        VZINTF(1:NPLSV,1:NRAD)  = PTL(NADDP(7)+1:NADDP(8),1:NRAD)
        DEALLOCATE (PTL)
        RETURN
        if (trcfle) then
          write (iunout,*) 'set pointers to plasma_backgrnd fields'
          write (iunout,*) 'called rplam_long. But iflg=10'
          write (iunout,*) 'Rest of comusr, as well as comxs, comsou:'
          write (iunout,*) 'is not read.'
        endif
        RETURN  ! this seems to be acting a bit like rplam_shrt.
cdr               No A&M data, and no primary source data are read.
cdr               But we are in rplam_long
cdr               So this is for field particles 1:npls, not only for npls_fix+1:npls
      ELSE
cdr set full plasma background data structure PLSTLS;
cdr     but no pointers DIINTF... to target PLASMA_BCKGRND
        PLSTLS = PTL
      END IF
      DEALLOCATE (PTL)

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
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: module EIRMOD_COMUSR.f'
      IF (IO /= 0) GOTO 990

cdr  Read A&M data. The routines READ_CMDTA, READ_CMAMF are contained
cdr  in module COMXS
      CALL EIRENE_READ_CMDTA
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMDTA,ICMDTA'
      CALL EIRENE_READ_CMAMF
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMAMF,ICMAMF'

      READ (13+ifoff,IOSTAT=IO) RCZT1,RCZT2,ZT1,ZRG
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCZT1,RCZT2,ZT1,ZRG'
      IF (IO /= 0) GOTO 990

      IF (IFLG == 0) THEN
c  read primary source parameters as well, step functions, etc.
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
        IF (TRCFLE) WRITE (iunout,*) 'READ 13: module EIRMOD_CSTEP.f'
        IF (IO /= 0) GOTO 990
      ELSE
        IF (TRCFLE)
     .   WRITE (iunout,*) 'SOURCE DATA NOT READ FROM ', FORT, '13'
      END IF

      CLOSE (UNIT=13+ifoff)
      RETURN

  990 CONTINUE
      WRITE (IUNOUT,*) ' RPLAM_LONG: cannot read FILE '//FORT//'13'
      CALL EIRENE_LEER(1)
c     CALL EIRENE_EXIT_OWN(1)
      IFLG=-1  ! Try to indicate to calling routine: fort.13 not found.
c                Then there we may choose to write a "default" fort.13
c                for an initial iteration step.
      CLOSE (UNIT=13+ifoff)
      return

  991 CONTINUE
      WRITE (IUNOUT,*) ' AVAILABLE INPUT TALLIES ARE DIFFERENT FROM',
     .                 ' PRIOR JOB WHICH WROTE '//FORT//'13'
      CALL EIRENE_EXIT_OWN(1)
  992 CONTINUE
      WRITE (IUNOUT,*) ' LEADING DIMENSIONS OF INPUT TALLIES ARE',
     .                 ' DIFFERENT FROM',
     .                 ' PRIOR JOB WHICH WROTE '//FORT//'13'
      CALL EIRENE_EXIT_OWN(1)
  993 CONTINUE
      WRITE (IUNOUT,*) ' STARTING POSITIONS OF INPUT TALLIES',
     .                 ' IN ARRAY PLSTLS ARE DIFFERENT FROM',
     .                 ' PRIOR JOB WHICH WROTE '//FORT//'13'
      CALL EIRENE_EXIT_OWN(1)

      END SUBROUTINE EIRENE_RPLAM_LONG
