      MODULE EIRMOD_UPTBGK
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_UPTBGK, EIRENE_uptbgk_reinit

      INTEGER, SAVE :: IFIRST=0

!$OMP THREADPRIVATE(IFIRST)

      CONTAINS

cdr Aug. 2015: revisited:  comments,...
c
c  code segment: bgk
c
c  only needed, if some test particle species are labeled as "bgk species"
c               with one or more elastic nonlinear self-interactions
c               to be treated by iteration.
c               This segment contains a routine UPTBGK which updates the tallies
c               required for iteration (carried out in MODBGK).
c
C  CURRENTLY:  3 TALLIES ARE SCORED PER BGK COLLISION SPECIES,IBGK_SP, IBGK_SP=1,NRBGI/3
c              On input: npbgk= npbgka(iatm), or npbgkm(imol), or npbgki(iion)
c              ibgk_sp=npbgk, and update three tallies for bgk species no. ibgk_sp.
c
c  do not confuse: ibgk is the bgk reaction number, the bgk reactions
c                  form a subset of the elastic reactions, IREL=1,NREL.
c
c                  ibgk_sp is the counter for the number of those test particle species
c                  which have at least one BGK collision.
c                  For each test particle species ibgk_sp there are currently
c                  three so-called additional "bgk tallies" scored
c                  (by default: the transport flux vector components).

c  Note:  for velocity-dependent BGK collision rates probably 5 tallies per bgk collision (ibgk)
c        need to be scored, rather than the three tallies per bgk species (ibgk_sp),
c         to enforce the 5 collision invariants by iteration.
c  Note:  for ES-BGK models (correct Prandtl number models) more than 3 bgk tallies
c         are needed per BGK species ibgk_sp (non-diagonal pressure tensor elements)
c
c  A routine (MODBGK) carries out the iterations at the end of an iteration step.
c
c
      SUBROUTINE EIRENE_UPTBGK(XSTOR2,XSTORV2,WV,NPBGK)
C
C  UPDATE BGK-SPECIFIC TALLIES, TRACKLENGTH ESTIMATORS
C
C  INPUT: NPBGK IDENTIFIER FOR THE BGK SPECIES
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CUPD
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT
      USE EIRMOD_CTEXT
      USE EIRMOD_CSDVI
      USE EIRMOD_COMXS
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                      XSTORV2(NSTORV,N2ND+N3RD)
      REAL(DP), INTENT(IN) :: WV
      INTEGER, INTENT(IN) :: NPBGK
      REAL(DP) :: DIST, WTRV, WTRVX, WTRVY, WTRVZ
      INTEGER :: I, NMTSP, IUPD2, IUPD3, 
     .           IRD, NSBGK, IBGK_SP,
     .           IML, IIO, IUPD1, ITP, ISP, IAT, IRDO
      CHARACTER(8) :: TXT

!$OMP THREADPRIVATE(DIST,WTRV,WTRVX,WTRVY,WTRVZ, 
!$OMP&              I, NMTSP, IUPD2, IUPD3, IRD, NSBGK, IBGK_SP,
!$OMP&              IML,IIO,IUPD1,ITP,ISP,IAT,IRDO,TXT)      

      SAVE
C
      IF (IFIRST.EQ.0) THEN

C  FIND TEST PARTICLE SPECIES FLAG (TYPE ITP, TEXT 'TXT') FOR BGK SPECIES NO. IBGK_SP
        IFIRST=1
C  NUMBER OF (ADDITIONAL) BGK TALLIES: NRBGI
C  NUMBER OF BGK SPECIES: NSBGK

        NSBGK=NRBGI/3
        DO IBGK_SP=1,NSBGK
          ITP=0
          DO ISP=1,NATMI
            IF (NPBGKA(ISP).EQ.IBGK_SP) THEN
C  ISP IS ONE OF THE ATOMIC TEST SPECIES WHICH HAVE AT LEAST ONE BGK COLLISION
              ITP=1
              IAT=ISP
              TXT=TEXTS(NSPH+IAT)
              GOTO 1
            ENDIF
          ENDDO
          DO ISP=1,NMOLI
            IF (NPBGKM(ISP).EQ.IBGK_SP) THEN
C  ISP IS ONE OF THE MOLECULAR TEST SPECIES WHICH HAVE AT LEAST ONE BGK COLLISION
              ITP=2
              IML=ISP
              TXT=TEXTS(NSPA+IML)
              GOTO 1
            ENDIF
          ENDDO
          DO ISP=1,NIONI
            IF (NPBGKI(ISP).EQ.IBGK_SP) THEN
C  ISP IS ONE OF THE TEST ION SPECIES WHICH HAVE AT LEAST ONE BGK COLLISION
              ITP=3
              IIO=ISP
              TXT=TEXTS(NSPAM+IIO)
              GOTO 1
            ENDIF
          ENDDO
C  PHOTONIC BGK COLLISIONS:  TO BE DONE ??

          WRITE (iunout,*) 'SPECIES ERROR IN UPTBGK'
          CALL EIRENE_EXIT_OWN(1)
    1     CONTINUE
C
C  BGK SPECIES NO. IBGK_SP
cdr  tbd. find corresponding bgk reaction irbg, which
cdr       has a vel. dep. reaction rate.
cdr  not ready:  this next code is partially
cdr              from obsolete old (vel. indep.) procedure still.
          IUPD1=(IBGK_SP-1)*3+1
          IUPD2=(IBGK_SP-1)*3+2
          IUPD3=(IBGK_SP-1)*3+3
cym these are all shared variables updated with private TXT & ITP
cym atomic  not usable with strings ?
!$OMP CRITICAL
          TXTTAL(IUPD1,NTALB)='BGK TALLY: FLUX DENSITY IN X DIRECTION '
          TXTTAL(IUPD2,NTALB)='BGK TALLY: FLUX DENSITY IN Y DIRECTION '
          TXTTAL(IUPD3,NTALB)='BGK TALLY: FLUX DENSITY IN Z DIRECTION '

          TXTUNT(IUPD1,NTALB)='#/CM**3*CM/S            '
          TXTUNT(IUPD2,NTALB)='#/CM**3*CM/S            '
          TXTUNT(IUPD3,NTALB)='#/CM**3*CM/S            '

          TXTSPC(IUPD1,NTALB)=TXT
          TXTSPC(IUPD2,NTALB)=TXT
          TXTSPC(IUPD3,NTALB)=TXT
!$OMP END CRITICAL
!$OMP ATOMIC WRITE
          IBGVE(IUPD1)=1
!$OMP ATOMIC WRITE
          IBGVE(IUPD2)=1
!$OMP ATOMIC WRITE
          IBGVE(IUPD3)=1
!$OMP ATOMIC WRITE
          IBGRC(IUPD1)=ITP
!$OMP ATOMIC WRITE
          IBGRC(IUPD2)=ITP
!$OMP ATOMIC WRITE          
          IBGRC(IUPD3)=ITP
        ENDDO

cdr  Species index increment for bgk tallies, used for LMETSP arrays.
cdr: This species index increment should be set in input.f,
cdr  like all the others.
cdr  Sequence: test species, bulk species, add tallies, alg. tallies, collest tallies,
cdr             cop tallies, bgk tallies.
        NMTSP=NPHOTI+NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI+NCLVI+NCPVI
C
C  END OF IFIRST BLOCK
      ENDIF
C
c  BGK tallies scoring starts here

C  UPDATE BGK TALLIES FOR THE NPBGK "BGK SPECIES"
C  PRESENTLY: UPDATE TRANSPORT FLUX VECTOR ON BGKV TALLY,
C  THREE TALLIES PER BGK SPECIES CONTRIBUTING IN BGK PROCESSES.
C
      IBGK_SP=NPBGK
C  FROM CALLING PROGRAM: IBGK_SP.NE.0, I.E. FOR THIS TEST PARTICLE (IATM, IMOL OR IION)
C  THE BGK TALLIES NO. IUPD1,IUPD2,IUPD3 NEED TO BE SCORED.
      IUPD1=(IBGK_SP-1)*3+1
      IUPD2=(IBGK_SP-1)*3+2
      IUPD3=(IBGK_SP-1)*3+3
      LMETSP(NMTSP+IUPD1)=.TRUE.
      LMETSP(NMTSP+IUPD2)=.TRUE.
      LMETSP(NMTSP+IUPD3)=.TRUE.
      DO 51 I=1,NCOU
        DIST=CLPD(I)
        WTRV=WV*DIST*VEL
        WTRVX=WTRV*VELX
        WTRVY=WTRV*VELY
        WTRVZ=WTRV*VELZ
! cdr: May 2017.
! BGKV output tallies are now identical with the default
!      VXDEN..., VYDEN..., VZDEN... tallies. Compare UPDATE.f, identical code!
!
        IRDO=NRCELL+NUPC(I)*NR1P2+NBLCKA
        IRD=NCLTAL(IRDO)
!$OMP ATOMIC
        BGKV(IUPD1,IRD)=BGKV(IUPD1,IRD)+WTRVX
!$OMP ATOMIC
        BGKV(IUPD2,IRD)=BGKV(IUPD2,IRD)+WTRVY
!$OMP ATOMIC
        BGKV(IUPD3,IRD)=BGKV(IUPD3,IRD)+WTRVZ
   51 CONTINUE

      END SUBROUTINE EIRENE_UPTBGK

csw 19apr07
      SUBROUTINE EIRENE_uptbgk_reinit
      IMPLICIT NONE
      ifirst=0
      return
csw
      END SUBROUTINE EIRENE_uptbgk_reinit

      END MODULE EIRMOD_UPTBGK
