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

      SUBROUTINE EIRENE_IF0PRM_JSON(json,p)

      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CINIT
      USE EIRMOD_BRAEIR
      USE EIRMOD_COMPRT, only: iunout      
      use json_module
     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck

      IMPLICIT NONE

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      INTEGER :: NFLA, NCUTB, NCUTL, NDXA, NDYA, IPL, NTARGI, IT, IPRT,
     .           NAINB, IAIN, NAOTB, IAOT, NRKNOT,
     .           NTRII, NR1ST
      INTEGER, ALLOCATABLE :: NTGPRT(:)
      LOGICAL :: FOUND
      CHARACTER(72) :: ZEILE

      WRITE (iunout,*) '*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"'

C  READ INPUT BLOCK 14

      call json%get(p,'NFLA',nfla,found)
      call json%get(p,'NCUTB',ncutb,found)
      call json%get(p,'NCUTL',ncutl,found)
      
C  GRID SIZE IN 2D PLASMA FLUID CODE
      call json%get(p,'NDXA',ndxa,found)
      call json%get(p,'NDYA',ndya,found)

C     NUMBER OF TARGET SOURCES ON B2 SURFACES: NTARGI
      call json%get(p,'NTARGI',ntargi,found)

C     NUMBER OF PARTS PER TARGET RECYCLING SOURCE
      IF (NTARGI.GT.0) THEN
        NSTEP=MAX(NSTEP,NTARGI)
        call json%get(p,'NTGPRT',ntgprt,found)
        NPTRGT=SUM(NTGPRT)
        DEALLOCATE (NTGPRT)
      END IF

C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2 INTO EIRENE
C  HERE: B2 VOLUME TALLIES
      call json%get(p,'NAINB',nainb,found)
      NAIN = MAX(NAIN,NAINB)

C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO B2
C  HERE: EIRENE SURFACE TALLIES
      call json%get(p,'NAOTB',naotb,found)
C
C READING BLOCK 14 FROM JSON INPUT FILE FINISHED
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
C
C
C  READ DATA FOR TRIANGULAR MESH
C
      nr1st = n1st
      call eirene_find_triang_dim(nr1st,ntri,ntrii,nknot,ngitt,' ',0)
      
      RETURN
      END
