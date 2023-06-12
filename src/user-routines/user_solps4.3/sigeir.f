C 
C
      SUBROUTINE EIRENE_SIGEIR(INIT,JJJ,ZDS,DUM1,PSIG,DUM2,ARGST)
CMG 3/9/21 (MG) User routine to integrate plasma parameters,
C               e.g. atomic density    
C  INPUT:
C          IFIRST: FLAG FOR INITIALISATION
C          NCELL:   INDEX IN TALLY ARRAYS FOR CURRENT ZONE
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C               
C  OUTPUT: CONTRIB. FROM CELL NCELL AND CHORD SEGMENT JJJ TO
C          THE PARAMETER PSIG (TOTAL) and ARGST (CONTR.)
C          ISP = NSPSPZ in block 13, defined for each ICHORI      
C          PSIG(ISP=1), ARGST(ISP=1,) - nH0 ---> to be checked
C          for order of species!
C          PSIG(2), ARGST (2,) - nH2      
C      
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CGEOM
      USE EIRMOD_CTRCEI

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: INIT,JJJ
      REAL(DP), INTENT(INOUT) :: PSIG(0:),ARGST(0:,:)
      REAL(DP), INTENT(IN) :: ZDS,DUM1,DUM2
      INTEGER :: NCELC,UDIM
      INTEGER, SAVE :: ITROLD
      LOGICAL :: LARGST
!pb   SAVE

      LARGST = SIZE(ARGST,2) >= NRAD
      udim=ubound(psig,1)

      IF (INIT.EQ.0) THEN
        PSIG=0.
        IF (LARGST) ARGST=0.

        IF (IITER .NE. ITROLD) THEN
          WRITE (IUNOUT,*) ' sigeir '
           
        ENDIF
        ITROLD=IITER
        IF (TRCSIG) WRITE(IUNOUT,*) 'SIGEIR, UDIM= ',UDIM
        RETURN
      ENDIF
 

C  LINE INTEGRAL: NEUTRAL PARAMETER * CM
C  
      udim=ubound(psig,1)
      ncelc=ncltal(ncell)
      if (udim >= 1)
     .  PSIG(1)=PSIG(1)+ZDS*PDENA(1,NCELC)
      if (udim >= 2)
     .  PSIG(2)=PSIG(2)+ZDS*PDENM(1,NCELC)

      IF (LARGST) THEN
        if (udim >= 1)
     .    ARGST(1,JJJ)=PDENA(1,NCELC)
        if (udim >= 2)
     .    ARGST(2,JJJ)=PDENM(1,NCELC)
      ENDIF
      
      RETURN
      END
