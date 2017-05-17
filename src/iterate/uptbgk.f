cdr Aug. 2015: revisited:  comments,...
c
c  code segment: bgk
c
c  only needed, if some test species are labeled as bgk-species
c               with one or more non-linear self interactions
c               this segment contains a routine which updates the tallies
c               required for iteration (UPTBGK).
c
C  CURRENTLY:  3 TALLIES ARE SCORED PER BGK COLLISION IBGK_SP, IBGK_SP=1,NRBGI/3
c              on input: npbgk= npbgka(iatm), or npbgkm(imol), npbgki(iion) 
c              ibgk=npbgk, and update three tallies for bgk collision no. ibgk_sp.
c
c  no not confuse: ibgk is the bgk reaction number, the bgk reactions form a subset of the elastic reactions
c                  ibgk_sp is the counter for the number of test particle species which have at least one bgk collision.
c                  For each test particle species ibgk_sp there are currently three so called additional "bgk tallies" scored
c                  (by default: the transport flux vector components). 
c  Note:  for velocity dependent bgk collision rates probably 5 tallies per bgk collision (ibgk)
c         need to be scored, rather than the three per bgk species (ibgk_sp),
c         to enforce the 5 collision invariants by iteration.
c
c               A routine (MODBGK) carries out the iterations at the end of an
c               iteration.

c               The standard deviations for the "bgk-tallies" are
c               computed in subroutine STATIS_BGK   ??? why  ???
C
c
c
      SUBROUTINE EIRENE_UPTBGK(XSTOR2,XSTORV2,WV,NPBGK,IFLAG)
C
C  UPDATE BGK-SPECIFIC TALLIES, TRACKLENGTH ESTIMATORS
C
C  INPUT:  NPBGK IDENTIFIER FOR THE BGK-SPECIES 
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
      INTEGER, INTENT(IN) :: NPBGK, IFLAG
      REAL(DP) :: DIST, WTRV, WTRVX, WTRVY, WTRVZ
      INTEGER :: I, NMTSP, IUPD2, IUPD3, IRD, IFIRST, NSBGK, IBGK_SP,
     .           IML, IIO, IUPD1, ITP, ISP, IAT, IRDD
      CHARACTER(8) :: TXT
      DATA IFIRST/0/
      SAVE
C
      IF (IFIRST.EQ.0) THEN

C  FIND TEST PARTICLE SPECIES FLAG (TYPE ITP, TEXT 'TXT') FOR BGK SPECIES NO. IBGK_SP
        IFIRST=1
C  NUMBER OF (ADDITIONAL) BGK TALLIES: NRBGI
C  NUMBER OF BGK-SPECIES:  NSBGK
        NSBGK=NRBGI/3
        DO IBGK_SP=1,NSBGK
          ITP=0
          DO ISP=1,NATMI
            IF (NPBGKA(ISP).EQ.IBGK_SP) THEN
C  ISP IS ONE OF THE ATOMIC TEST SPECIES WHICH HAVE AT LEAST ONE BKG COLLISION
              ITP=1
              IAT=ISP
              TXT=TEXTS(NSPH+IAT)
              GOTO 1
            ENDIF
          ENDDO
          DO ISP=1,NMOLI
            IF (NPBGKM(ISP).EQ.IBGK_SP) THEN
C  ISP IS ONE OF THE MOLECULAR TEST SPECIES WHICH HAVE AT LEAST ONE BKG COLLISION
              ITP=2
              IML=ISP
              TXT=TEXTS(NSPA+IML)
              GOTO 1
            ENDIF
          ENDDO
          DO ISP=1,NIONI
            IF (NPBGKI(ISP).EQ.IBGK_SP) THEN
              ITP=3
              IIO=ISP
              TXT=TEXTS(NSPAM+IIO)
              GOTO 1
            ENDIF
          ENDDO
C  PHOTONIC BGK COLLISIONS:  TO BE DONE ??

          WRITE (iunout,*) 'SPECIES ERROR IN UPTBGK'
          CALL EIRENE_EXIT_OWN(1)
1         CONTINUE
C
C  BGK-SPECIES NO. IBGK_SP
          IUPD1=(IBGK_SP-1)*3+1
          IUPD2=(IBGK_SP-1)*3+2
          IUPD3=(IBGK_SP-1)*3+3
          TXTTAL(IUPD1,NTALB)='BGK TALLY: FLUX DENSITY IN X DIRECTION '
          TXTTAL(IUPD2,NTALB)='BGK TALLY: FLUX DENSITY IN Y DIRECTION '
          TXTTAL(IUPD3,NTALB)='BGK TALLY: FLUX DENSITY IN Z DIRECTION '
          TXTUNT(IUPD1,NTALB)='#/CM**3*CM/S            '
          TXTUNT(IUPD2,NTALB)='#/CM**3*CM/S            '
          TXTUNT(IUPD3,NTALB)='#/CM**3*CM/S            '
          TXTSPC(IUPD1,NTALB)=TXT
          TXTSPC(IUPD2,NTALB)=TXT
          TXTSPC(IUPD3,NTALB)=TXT
          IBGVE(IUPD1)=1
          IBGVE(IUPD2)=1
          IBGVE(IUPD3)=1
          IBGRC(IUPD1)=ITP
          IBGRC(IUPD2)=ITP
          IBGRC(IUPD3)=ITP
        ENDDO
 
        NMTSP=NPHOTI+NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI+NCLVI+NCPVI
C
C  END OF IFIRST BLOCK
      ENDIF
C
C  UPDATE BGK TALLIES FOR THE NPBGK "BGK-SPECIES"
C  PRESENTLY: UPDATE TRANSPORT FLUX VECTOR ON BGKV-TALLY, THREE TALLIES PER BGK-SPECIES
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
!pb 05.02.2013  BGKV is output tally thus only defined for (COARSE) NRTAL cells
        IRDD=NRCELL+NUPC(I)*NR1P2+NBLCKA
        IRD=NCLTAL(IRDD)
        BGKV(IUPD1,IRD)=BGKV(IUPD1,IRD)+WTRVX
        BGKV(IUPD2,IRD)=BGKV(IUPD2,IRD)+WTRVY
        BGKV(IUPD3,IRD)=BGKV(IUPD3,IRD)+WTRVZ
51    CONTINUE
 
      RETURN
csw 19apr07
      entry EIRENE_uptbgk_reinit
      ifirst=0
      return
csw
      END
