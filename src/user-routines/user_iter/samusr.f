C
C
      SUBROUTINE EIRENE_SAMUSR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)
C
C  SAMPLE INITIAL COORDIATES X0,Y0,Z0 ON SURFACE NLSF
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CADGEO
      USE EIRMOD_COMUSR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6
      REAL(DP), INTENT(OUT) :: X0,Y0,Z0,TEWL,TIWL(*),DIWL(*),
     .                         VXWL(*),VYWL(*),VZWL(*),
     .                         EFWL(*), SHWL, WEISPZ(*)
      INTEGER, INTENT(IN) :: NLSF,is1, is2
      INTEGER, INTENT(OUT) :: IRUSR, IPUSR, ITUSR, IAUSR, IBUSR
      REAL(DP) :: X, Y, T, B0, B1, B2, Z1, Z2
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER :: IER

C  CALLED IN INITIALIZATION PHASE
C  E.G. TO INITIALIZE SAMPLING ON SURFACE 
      ENTRY EIRENE_SM0USR
     .  (is1,is2,sorad1,sorad2,sorad3,sorad4,sorad5,sorad6)
      return

C.............................................................................
C  CALLED FROM SUBR. SAMSRF, FOR SURFACE SAMPLING 
      ENTRY EIRENE_SM1USR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)


C  RETURN BIRTH POINT OF TEST PARTICLE 
      x0 = 0._dp
      y0 = 0._dp
      z0 = 0._dp
C  RETURN CELL NO. INFORMATION AT BIRTH POINT 
      irusr = 0
      ipusr = 0
      itusr = 0
      iausr = 0
      ibusr = 0
C  RETURN BACKGROUND MEDIUM PARAMETERS AT BIRTHPOINT 
      tiwl(1:nplsti) = 0._dp
      tewl = 0._dp
      diwl(1:nplsi) = 0._dp
      vxwl(1:nplsv) = 0._dp
      vywl(1:nplsv) = 0._dp
      vzwl(1:nplsv) = 1._dp
      efwl(1:nplsi) = 0._dp
      shwl = 0._dp
      weispz(1:nspz) = 0._dp
 
      RETURN
      END
