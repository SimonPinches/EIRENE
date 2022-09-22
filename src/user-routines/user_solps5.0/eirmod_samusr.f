      MODULE EIRMOD_SAMUSR
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CADGEO
      USE EIRMOD_COMUSR
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_SAMUSR, EIRENE_SAMUSR_INIT

      REAL(DP), SAVE :: SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6

      CONTAINS
      
C
C  SAMPLE INITIAL COORDINATES X0,Y0,Z0 ON SURFACE NLSF
C

C  CALLED IN INITIALIZATION PHASE
C  E.G. TO INITIALIZE SAMPLING ON SURFACE
      SUBROUTINE EIRENE_SAMUSR_INIT (is1,is2,
     .                               sorad1_in,sorad2_in,sorad3_in,
     .                               sorad4_in,sorad5_in,sorad6_in)
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: SORAD1_IN,SORAD2_IN,SORAD3_IN,
     .                        SORAD4_IN,SORAD5_IN,SORAD6_IN
      INTEGER, INTENT(IN) :: is1, is2

      SORAD1 = SORAD1_IN
      SORAD2 = SORAD2_IN
      SORAD3 = SORAD3_IN
      SORAD4 = SORAD4_IN
      SORAD5 = SORAD5_IN
      SORAD6 = SORAD6_IN

      RETURN
      END SUBROUTINE EIRENE_SAMUSR_INIT

C.............................................................................
C  CALLED FROM SUBR. SAMSRF, FOR SURFACE SAMPLING
      SUBROUTINE EIRENE_SAMUSR (NLSF,X0,Y0,Z0,
     .              SORAD1_IN,SORAD2_IN,SORAD3_IN,
     .              SORAD4_IN,SORAD5_IN,SORAD6_IN,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,ZIWL,
     .              WEISPZ)
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: SORAD1_IN,SORAD2_IN,SORAD3_IN,
     .                        SORAD4_IN,SORAD5_IN,SORAD6_IN
      REAL(DP), INTENT(OUT) :: X0,Y0,Z0,TEWL,TIWL(*),DIWL(*),
     .                         VXWL(*),VYWL(*),VZWL(*),
     .                         EFWL(*), SHWL, WEISPZ(*), ZIWL(*)
      INTEGER, INTENT(IN) :: NLSF
      INTEGER, INTENT(OUT) :: IRUSR, IPUSR, ITUSR, IAUSR, IBUSR

      SORAD1 = SORAD1_IN
      SORAD2 = SORAD2_IN
      SORAD3 = SORAD3_IN
      SORAD4 = SORAD4_IN
      SORAD5 = SORAD5_IN
      SORAD6 = SORAD6_IN

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
      ziwl(1:nplsi) = dble(nchrgp(1:nplsi))
      shwl = 0._dp
      weispz(1:nspz) = 0._dp

      RETURN
      END SUBROUTINE EIRENE_SAMUSR

      END MODULE EIRMOD_SAMUSR
