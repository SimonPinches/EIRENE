cdr  jan 16: started to cleanup, comment
cdr          remove redundant parameter iflg
cdr  jan 18: comments
c
c
C
C MODIFIED BY V. KOTOV  (when ?)
C
      SUBROUTINE EIRENE_WRPLAM_SHRT(TRCFLE)

cdr Only the input tallies of the last (virtual) plasma species: nfla+1,...nplsi
cdr are written/read using I/O stream fort.13

c  if NLSRT13=true : wrplam_short and (entry) rplam_short are called from WRPLAM,
c  if NLSRT13=false: wrplam_long and (entry)  rplam_long are called from WRPLAM,
c

cdr this is the SHORT version of WRPLAM.F
cdr It writes and reads (entry RPLAM_SHRT) background data onto/from fort.13
cdr Distinct from WRPLAM_long here only the background tallies are written/read,
cdr (tallies T, n, V for ipls=1,nplsi), but not the atomic data,
cdr nor the primary source sampling information.
cdr for BGK type non-linear iterations, with velocity-independent rates,
cdr this may be sufficient.
cdr Better: add here also the rates and other atomic data needed to streamline
cdr         non-linear iterations.



      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: IFOFF, NRAD
      USE EIRMOD_CINIT, ONLY: FORT
      USE EIRMOD_COMUSR, ONLY: NPLSI, TIIN, DIIN, VXIN, VYIN, VZIN
      USE EIRMOD_CCOUPL, ONLY: NFLA
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      LOGICAL,INTENT(IN) :: TRCFLE
      INTEGER IO


      OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13+ifoff
      IF(.NOT.ASSOCIATED(NFLA)) THEN
        WRITE(iunout,*)
     w         "ERROR IN WRPLAM_SHRT: NFLA WAS NOT ASSOCIATED. ",
     w         "NO DATA WILL BE STORED IN ", FORT, "13"
        RETURN
      END IF
cdr  only write plasma background data for species, which are not already
cdr  transfered via Common BRAEIR, i.e. only:  nfla+1,....nplsi
cdr  I.e. the virtual background species for non-linear iterations
cdr  have to come last in the list of all background species.
      IF(NFLA.LT.NPLSI) THEN
        WRITE (13+ifoff,IOSTAT=IO)
     w           TIIN(NFLA+1:NPLSI,1:NRAD),DIIN(NFLA+1:NPLSI,1:NRAD),
     w           VXIN(NFLA+1:NPLSI,1:NRAD),VYIN(NFLA+1:NPLSI,1:NRAD),
     w           VZIN(NFLA+1:NPLSI,1:NRAD)
      END IF

      CLOSE (UNIT=13+ifoff)

      RETURN
C .......................................................................

      ENTRY EIRENE_RPLAM_SHRT(TRCFLE)

C ........................................................................


      OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     o      STATUS='OLD',IOSTAT=IO)
      IF(IO.NE.0) THEN
        WRITE(iunout,*)
     w   'ERROR IN RPLAM_SHRT: CANNOT READ ', FORT, '13'
        RETURN
      END IF

      REWIND 13+ifoff
      IF(.NOT.ASSOCIATED(NFLA)) THEN
        WRITE(IUNOUT,*)
     w       "ERROR IN RPLAM_SHRT: NFLA IS NOT ASSOCIATED. ",
     w       "NO DATA WILL BE READ FROM ", FORT, "13"
        RETURN
      END IF

      IF(NFLA.LT.NPLSI) THEN
cdr  only read plasma background data for species, which are not already
cdr  transfered via Common BRAEIR, i.e. only: nfla+1,....nplsi
cdr  I.e. the virtual background species for non-linear iterations
cdr  have to come last in the list of all background species.  
        REWIND 13+ifoff
        READ (13+ifoff,IOSTAT=IO)
     R        TIIN(NFLA+1:NPLSI,1:NRAD),DIIN(NFLA+1:NPLSI,1:NRAD),
     R        VXIN(NFLA+1:NPLSI,1:NRAD),VYIN(NFLA+1:NPLSI,1:NRAD),
     R        VZIN(NFLA+1:NPLSI,1:NRAD)
        IF(IO.NE.0) GOTO 200
        IF (TRCFLE) WRITE (iunout,*)
     w   'RPLAM: BGK BACKGROUND IS READ FROM ', FORT, '13'
      END IF
csw   CALL READ_TABEF(TRCFLE) !VK, READS TABEF, SEE CCRM
      CLOSE (UNIT=13+ifoff)
      RETURN

  200 CONTINUE

      WRITE(iunout,*) 
     w 'ERROR IN RPLAM_SHRT: CANNOT READ ', FORT, '13: ',
     w 'ZERO BACKGROUND WILL BE ASSIGNED'
      TIIN(NFLA+1:NPLSI,1:NRAD)=0._DP
      DIIN(NFLA+1:NPLSI,1:NRAD)=0._DP
      VXIN(NFLA+1:NPLSI,1:NRAD)=0._DP
      VYIN(NFLA+1:NPLSI,1:NRAD)=0._DP
      VZIN(NFLA+1:NPLSI,1:NRAD)=0._DP

      RETURN

      END SUBROUTINE EIRENE_WRPLAM_SHRT
