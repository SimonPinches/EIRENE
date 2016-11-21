!pb APR  16: piods -> pioei
cdr Nov 16: nmdsi --> nmeii,  and comments

C
C
      SUBROUTINE EIRENE_CONDENSE
C
C  CONDENSE COLLISION KERNEL, IF SOME SECONDARIES ARE NOT FOLLOWED
C  BY EIRENE, I.E., IF NFOLA(IATM), NFOLM(IMOL), NFOLI(IION) LT 0
C  FOR SOME TEST PARTICLE SPECIES
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CTEXT
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
 
      IMPLICIT NONE
 
<<<<<<< HEAD
      INTEGER :: ISP, ISP0, IION, ICOL, IATM, IMOL, IREI, IMEI
=======
      INTEGER :: ISP, IION, ICOL, IATM, IMOL, IREI, IMEI
>>>>>>> generation-limit
 
      DO 10 IATM=1,NATMI
C  NRCA=0 ?
        DO 100 ICOL=1,NRCA(IATM)
100     CONTINUE
10    CONTINUE
C
C
      DO 20 IMOL=1,NMOLI
C  currently: only electron impact collisions on molecules
        DO 200 IMEI=1,NMEII(IMOL)
          IREI=LGMEI(IMOL,IMEI)
<<<<<<< HEAD
          ISP0=NSPA+IMOL
=======
C  electron impact process no. irei, on molecules imol
c  search for secondaries, that are not followed:
c  atom secondaries
c         DO   NATMI 
c  molecule secondaries
C         DO   NMOLI
c  photonic secondaries
C         DO   NPHOTI
C  test ion secondaries:  
>>>>>>> generation-limit
          DO 220 IION=1,NIONI
            ISP=NSPAM+IION
            IF (PIOEI(IREI,IION).GT.0) THEN
              IF (NFOLI(IION).LT.0) THEN
                WRITE (iunout,*) 'TEST ION ',TEXTS(ISP),
     .                           'BORN FROM MOLECULE ', TEXTS(ISP0),
     .                           'CAN BE CONDENSED'
              ENDIF
            ENDIF
220       CONTINUE
200     CONTINUE
20    CONTINUE
C
      RETURN
      END
