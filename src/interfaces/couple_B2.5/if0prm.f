cdr called from find_param.f in initialization phase,
cdr when eirene is in "coupled mode": i.e. IF(NMODE.NE.0)
C
cdr Read block 14 from interfacing routines (not from eirene_input.f)
c   this version: couple_dummy, i.e. only dummy interfacing routines.
c
c   and set storage for allocatable arrays:
c   NSTEP :
c   NPTRGT:
c   NAIN  :
c   NAOT  :
c   NCPV  :
c   NKNOT :
C   NTRII :
C   NCPVI :  no. of special couple tallies

c  also set: NDX,NDY,NFL, NDXP, NDYP

      SUBROUTINE EIRENE_IF0PRM(IUNIN)

      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CINIT
      USE EIRMOD_BRAEIR

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IUNIN
      INTEGER :: NFLA, NCUTB, NCUTL, NDXA, NDYA, IPL, NTARGI, IT, IPRT,
     .           NAINB, IAIN, NAOTB, IAOT, NRKNOT,
     .           NTRII
      INTEGER, ALLOCATABLE :: NTGPRT(:)
      CHARACTER(72) :: ZEILE

C  READ INPUT BLOCK 14
      READ (IUNIN,*)
      READ (IUNIN,'(3I6)') NFLA,NCUTB,NCUTL
      DO IPL=1,NPLS
        READ (IUNIN,*)
      END DO
C  GRID SIZE IN 2D PLASMA FLUID CODE
      READ (IUNIN,'(2I6)') NDXA,NDYA
C  NUMBER OF TARGET SOURCES ON B2 SURFACES: NTARGI
      READ (IUNIN,'(I6)') NTARGI
C  NUMBER OF PARTS PER TARGET RECYCLING SOURCE
      IF (NTARGI.GT.0) THEN
        NSTEP=MAX(NSTEP,NTARGI)
        ALLOCATE (NTGPRT(NTARGI))
        READ (IUNIN,'(12I6)') (NTGPRT(IT),IT=1,NTARGI)
        NPTRGT=SUM(NTGPRT)
        DO IT=1,NTARGI
          DO IPRT=1,NTGPRT(IT)
  331       READ (IUNIN,'(A72)') ZEILE
            IF (ZEILE(1:1).EQ.'*') THEN
              GOTO 331
            ENDIF
          END DO
        END DO
        DEALLOCATE (NTGPRT)
      END IF
      READ (IUNIN,*)
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2 INTO EIRENE
C  HERE: B2 VOLUME TALLIES
      READ (IUNIN,'(I6)') NAINB
      NAIN = MAX(NAIN,NAINB)
      DO IAIN=1,NAINB
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
      END DO
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO B2
C  HERE: EIRENE SURFACE TALLIES
      READ (IUNIN,'(I6)') NAOTB
      DO IAOT=1,NAOTB
        READ (IUNIN,*)
      END DO
C
C READING BLOCK 14 FROM FORMATTED INPUT FILE (IUNIN) FINISHED
C
C
C  DEFINE ADDITIONAL TALLIES FOR COUPLING (UPDATED IN SUBR. UPTCOP
C                                          AND IN SUBR. COLLIDE)
      NCPVI=NPLS
      NCPV = MAX(NCPV,NCPVI)
C
C SAVE SOME MORE INPUT DATA FOR SHORT CYCLE ON COMMON CCOUPL
      NDX = NDXA
      NDY = NDYA
      NFL = NFLA
      IF (ALLOCATED(DNIB)) THEN
        NDX=UBOUND(DNIB,1)-1
        NDY=UBOUND(DNIB,2)-1
        NFL=UBOUND(DNIB,3)
      END IF
      NDXP = NDX+1
      NDYP = NDY+1

      RETURN
      END
