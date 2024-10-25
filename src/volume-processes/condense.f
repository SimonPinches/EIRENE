!pb APR 16: piods -> pioei
cdr Nov 16: nmdsi -> nmeii,  and comments

C
C
      SUBROUTINE EIRENE_CONDENSE
C
C  CONDENSE COLLISION KERNEL, IF SOME SECONDARIES ARE NOT FOLLOWED
C  BY EIRENE, I.E., IF NFOLA(IATM), NFOLM(IMOL), NFOLI(IION) LT 0
C  FOR SOME TEST PARTICLE SPECIES
C
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

      INTEGER :: ISP, ISP0, IION, ICOL, IATM, IMOL, IREI, IMEI
      EXTERNAL :: EIRENE_LEER

      DO 10 IATM=1,NATMI
C  NRCA=0 ?
        DO 100 ICOL=1,NRCA(IATM)
  100   CONTINUE
   10 CONTINUE
C
C
      DO 20 IMOL=1,NMOLI
C  currently: only electron impact collisions on molecules
        DO 200 IMEI=1,NMEII(IMOL)
          IREI=LGMEI(IMOL,IMEI)
          ISP0=NSPA+IMOL
C  electron impact process no. irei, on molecules imol
c  search for secondaries, that are not followed:
c  atom secondaries
c         DO   NATMI
c         ENDDO
c  molecule secondaries
C         DO   NMOLI
c         ENDDO
c  photonic secondaries
C         DO   NPHOTI
c         ENDDO
C  test ion secondaries:
          DO IION=1,NIONI
            ISP=NSPAM+IION
            IF (PIOEI(IREI,IION).GT.0) THEN
              IF (NFOLI(IION).LT.0) THEN
                WRITE (iunout,*) 'TEST ION ',TEXTS(ISP),
     .                           ' BORN FROM MOLECULE ',TEXTS(ISP0),
     .                           ' CAN BE CONDENSED'
              ENDIF
            ENDIF
          ENDDO
  200   CONTINUE
C  here the same for heavy particle impact collisions on molecules
c  tbd.
   20 CONTINUE
      call eirene_leer(1)
C
      RETURN
      END SUBROUTINE EIRENE_CONDENSE
