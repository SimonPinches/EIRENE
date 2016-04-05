cdr  jan 16: started to cleanup, comment
cdr          remove redundant parameter iflg
c
c
c
C
C MODIFIED BY V. KOTOV  (when ??)
C
      SUBROUTINE EIRENE_WRPLAM_SHRT(TRCFLE)

cdr this is the SHORT version of WRPLAM.F 
cdr It writes and reads (entry RPLAM_SHRT) background data onto/from fort.13
cdr Distinct from WRPLAM_long here only the background tallies are written-read, 
cdr (tallies T, n, V for ipls=1,nplsi), but not the atomic data, 
cdr nor the primary source sampling information 



      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CZT1
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMXS
      USE EIRMOD_CESTIM
csw      USE EIRMOD_CCRM
      USE EIRMOD_CCOUPL
      USE EIRMOD_COMPRT,ONLY:IUNOUT !VKMPI
csw      USE IFWRITE !VK

      IMPLICIT NONE

      LOGICAL,INTENT(IN) :: TRCFLE
      INTEGER IO


        OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        REWIND 13+ifoff
        IF(.NOT.ASSOCIATED(NFLA)) THEN
         WRITE(iunout,*) 
     w         "ERROR IN WRPLAM_SHRT: NFLA WAS NOT ASSOCIATED. ",
     w         "NO DATA WILL BE STORED IN FORT.13"
         RETURN
        END IF
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
        WRITE(iunout,*) 'ERROR IN WRPLAM_SHRT: CAN NOT READ FORT.13'
        RETURN
      END IF

      REWIND 13+ifoff
      IF(.NOT.ASSOCIATED(NFLA)) THEN 
        WRITE(IUNOUT,*) 
     w       "ERROR IN WRPLAM_SHRT: NFLA IS NOT ASSOCIATED ",
     w       "NO DATA WILL BE STORED IN FORT.13"
        RETURN
      END IF

       IF(NFLA.LT.NPLSI) THEN
C FIRST TRY TO READ IN THE OLD "LONG" FORMAT
csw 02jan2012 NO! will kill DIIN coming from B2.5 by memory transfer..
csw         READ (13,IOSTAT=IO) TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN
csw         IF(IO.EQ.0) THEN
csw          WRITE(IUNOUT,*) "WARNING FROM WRPLAM: ", 
csw     w                 "THE DATA IS READ IN THE OLD (LONG) FORMAT"
csw         ELSE
C IF READING IN OLD FORMAT DOESN'T WORK, THEN TRY THE NEW ONE        
          REWIND 13+ifoff  
          READ (13+ifoff,IOSTAT=IO)
     R          TIIN(NFLA+1:NPLSI,1:NRAD),DIIN(NFLA+1:NPLSI,1:NRAD),
     R          VXIN(NFLA+1:NPLSI,1:NRAD),VYIN(NFLA+1:NPLSI,1:NRAD),
     R          VZIN(NFLA+1:NPLSI,1:NRAD)
          IF(IO.NE.0) GOTO 200
          IF (TRCFLE) WRITE (iunout,*)
     w                'WRPLAM: BGK BACKGROUND IS READ FROM FORT.13'
csw        END IF !IF(IO.EQ.0) THEN
       END IF
csw      CALL READ_TABEF(TRCFLE) !VK, READS TABEF, SEE CCRM
      CLOSE (UNIT=13+ifoff)
      RETURN

 200  CONTINUE

       WRITE(iunout,*) 'ERROR IN WRPLAM_SHRT: CAN NOT READ FORT.13',
     w                 'ZERO BACKGROUND WILL BE ASSIGNED'
       TIIN(NFLA+1:NPLSI,1:NRAD)=0._DP
       DIIN(NFLA+1:NPLSI,1:NRAD)=0._DP
       VXIN(NFLA+1:NPLSI,1:NRAD)=0._DP
       VYIN(NFLA+1:NPLSI,1:NRAD)=0._DP
       VZIN(NFLA+1:NPLSI,1:NRAD)=0._DP

       RETURN

      END SUBROUTINE EIRENE_WRPLAM_SHRT
