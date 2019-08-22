c   modbgk  FZJ-MASTER, eirene_git, sept. 2014

C       ??     INDIRECT SPECIES INDEXING IPLSV, IPLSTI INCLUDED, FOR TIIN and VXIN,VYIN,VZIN
CDR            This is risky, because bgk collisions
cdr            may then overwrite temperatures/velocities for non-bgk background
CDR            depending on setting of MPLSTI(ipls) and MPLSV(ipls) arrays.
cdr            in case of bgk: always multiple temperatures and multiple flow velocities.

cdr  tbd:  implement corresponding check.

c
c  to be done: cross-temperature correction (Kotov)  (already done in solps-iter part)
C              HERE: ALLOCATION, DEALLOCATION: DONE, but not used
c
!pb  24.11.06: flag for shifting of first parameter of rate coeff introduced
!pb  24.11.06: BZIN initialized with 1
C
!pb  05.04.11: BFIN initialized with 1
!pb  18.02.13: take index transformation via NCLTAL into account


cdr 29.08.15 :  CROSSTEMP WAS INTRODUCED BY VK SO THAT different
c               EFFECTIVE COLLISION RATE IS EVALUATED
c               WITH THIS EFFECTIVE TEMPERATURE
c               this appears to be another option for the 5th (free)
c               bgk model parameter, to match certain transport coefficients.
c               (the other 4 are fixed by conservation laws.
c
cdr  june 2017: revisited:
c  multigrid option (NCLTAL):
c  if not NCLTAL (I) = I everywhere, then we have two grids, grid structures
C  NCLTAL(I-FINE):  CELL I-FINE IS ONLY A PART OF COARSER (SCORING) GRID CELL NCELL,
C                   NCELL=NCLTAL(I-FINE)
C                   SCORING OF VOLUME-AVERAGED TALLIES IS ON COARSE GRID CELLS NCELL ONLY.
cdr  nov. 17:  PLS added to call xstel. (strictly not needed here until now,
cdr                but for other EL rates enhanced by CR effects)
c
C
      SUBROUTINE EIRENE_MODBGK

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_BGK
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI

      IMPLICIT NONE
C
      REAL(DP), ALLOCATABLE :: PDEN(:),  EDEN(:),
     .                         PDEN2(:), EDEN2(:), ENERGY(:,:),
     .                         CROSSTEMP(:,:)
!pb 05.02.2013
      REAL(DP), ALLOCATABLE :: GBGKV(:,:)  ! TALLIES SCORED FOR BGK RELAXATION
cdr:  Nov. 17:for sync between xstel, xstpi, etc...
      REAL(DP), ALLOCATABLE :: PLS(:)

CDR THESE TALLIES bgkv ARE SCORED IN (COARSER) SCORING GRID.
C   ITERATION IS ON TALLIES DEFINED ON FINER GRID.
C   IF ONLY ONE GRID IS USED: NCLTAL(ICELL)=ICELL ALWAYS, AND
C   GBGKV == BGKV EVERYWHERE

      REAL(DP) :: RATM(3)
      REAL(DP) :: VXIN1, VXIN2, VYIN1, VYIN2, VZIN1, VZIN2, VYMIX, T1,
     .            T2, ED1, ED2, VXMIX, VZMIX, DELX, DELY, DELZ, VX, VY,
     .            VZ, EOLD, ED, RM, FACTKK, TMIX, EBULK, TS1, DS1,
     .            FACT2, RMAS2, A_ROBIN, FACT1, RMAS1,
     .            RESE, RESM, TBEL, DOLD, DEIMIN,
     .            DEL, TII, CNDYN, RRN, RATE, RESN, RATN, RRE, RRM,
     .            RM1, RM2, TCSUM, !VK
     .            EIRENE_RATE_COEFF
      INTEGER ::  ITYP1(NPLS), ITYP2(NPLS), ISPZ1(NPLS), ISPZ2(NPLS),
     .            IREL1(NPLS), INRC1(NPLS),CROSSINDEX(NPLS),NCROSS,
     .            ICROSS1,ICROSS2,ICROSS
      INTEGER ::  J, IESTM, IBGV,
     .            ISCDE, IAIN, IPLS1, NRWK1, I, ISP, IPLS2, IATM2, IIEL,
     .            IAEL, IMEL, IUP12, IUP22, IION2, IBGK2, IMOL2, IUP2,
     .            IUP3, IUP1, IBGK1, IP, NRC, IR, IT, IREL,
     .            KK, NXM, NYM, NZM, IUP32, IPLSTI, IPLSTI1, IPLSTI2,
     .            IPLSV, IRD, I_FINE, IRAD, IFLG
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      LOGICAL :: LMARK(NPLS)
      LOGICAL :: TRCSAV
C
C
      CALL EIRENE_LEER(3)
      WRITE (iunout,*) 'MODBGK CALLED AFTER ITERATION IITER= ',IITER
      CALL EIRENE_LEER(3)
C
cdr  FOR STOCH. APPROX. UNDER-RELAXTION.  NOT IN USE
      A_ROBIN=1.D0/REAL(IITER,KIND(1.D0))

c  iterate on tallies from "sum over strata".
      ISTRA=0
      IF (NSTRAI.EQ.1.AND.IESTR.EQ.1) ISTRA=1
c
c read output tally data, sum over strata
      IF (ISTRA.EQ.IESTR) THEN
C  NOTHING TO BE DONE, DATA ARE ALREADY FOR "SUM OVER STRATA"
      ELSEIF (NFILEN.NE.0) THEN
        IESTR=0
        CALL EIRENE_RSTRT(0,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSIGI_SPC,TRCFLE)
      ELSE
        WRITE (iunout,*) 'ERROR IN MODBGK: DATA FOR STRATUM ISTRA= ',
     .                    ISTRA
        WRITE (iunout,*) 'ARE NOT AVAILABLE. EXIT CALLED'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      if (NBMLT.gt.1) then
        write (iunout,*) 'NBMLT-option not ready in MODBGK '
        call EIRENE_exit_own(1)
      endif
      if (any(MPLSTI(1:npls) /= (/(i,i=1,npls)/))) then
        write (iunout,*) 'MPLSTI-option not ready in MODBGK '
        call EIRENE_exit_own(1)
      endif
      if (any(MPLSV(1:npls) /=  (/(i,i=1,npls)/))) then
        write (iunout,*) 'MPLSV-option not ready in MODBGK '
        call EIRENE_exit_own(1)
      endif
      if (any(NCLTAL(1:nrad) /=  (/(i,i=1,nrad)/))) then
        write (iunout,*) 'NCLTAL-option perhaps ready in MODBGK '
c       call EIRENE_exit_own(1)
      endif
C
      NBLCKA=0
      ALLOCATE (PDEN(NRAD))
      ALLOCATE (EDEN(NRAD))
      ALLOCATE (PDEN2(NRAD))
      ALLOCATE (EDEN2(NRAD))
      ALLOCATE (ENERGY(NPLS,NRAD))

cdr  PLS:  ELECTRON DENSITY PARAMETER in CR MODELS
cdr       (NOT TO BE CONFUSED WITH THE DENSITY FACTOR BETWEEN RATES AND RATE COEFF.)
cdr: set hard-wired lower density for H.4, H.10 type fits from AMJUEL: 1e8 cm**-3
cdr: at this lower limit density the fits are produced such
cdr: that they collapse to the corona limit values.
      ALLOCATE (PLS(NSTORDR))
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO J=1,NSBOX
          PLS(J)=MAX(DEIMIN,DEINL(J))
        ENDDO
      END IF

CVK  CALCULATES NUMBER OF CROSS-COLLISION PROCESSES 'ncross' AND ALLOCATES MEMORY FOR
C    "CROSS-COLLISION TEMPERATURE" CORRECTION
      NCROSS=0
      DO I=1,NPLS
        IF(NPBGKP(I,2).NE.0) NCROSS=NCROSS+1
      END DO
      IF(NCROSS.GT.0) THEN
        ALLOCATE (CROSSTEMP(NCROSS,NRAD))
        CROSSTEMP=TVAC
      END IF
      ICROSS=0 !VK COUNTER
CVK END

cdr  set bgk volume-averaged tallies BGKV, which MAY HAVE been scored on goarser grid only,
CDR  now also on fine grid, because input tallies may be needed on finer grid.
cdr  Assume: constant "extensive" fine grid values IN ALL I_fine cells
cdr  within one coarse grid cell IRD

      ALLOCATE (GBGKV(NBGV,NRAD))
      GBGKV = 0._DP
      DO IBGV = 1, NBGVI
        DO I_fine=1,NRAD  ! fine grid cell loop
          IRD = NCLTAL(I_fine) ! coarse scoring cell IRD contains cell I_fine
          IF (IRD == 0) CYCLE
C  ird=0: additional cells for averaging over all fine grid cells.
c         no value assigned.
          GBGKV(IBGV,I_fine)=BGKV(IBGV,IRD)
c  COARSER INPUT CELL VALUES ARE TURNED INTO FINER INPUT CELL VALUES
C  IN USER-SPECIFIC PARTS.
C  FINE GRID INPUT CELLS ARE AVERAGED BACK TO COARSER CELLS VALUES,
c  DEPENDING ON... (EXTENSIVE, INTENSIVE QUANTITIES) SEE OUTPLA.
        END DO
      END DO

c
C  LOOP OVER THOSE BACKGROUND ION SPECIES, WHICH ARE ARTIFICIAL
C  SPECIES FOR (NON-LINEAR) ITERATIONS
C .........................................................................
      DO 1000 IPLS=1,NPLSI
C .........................................................................
C
C  IS IPLS AN ARTIFICIAL "BGK BACKGROUND SPECIES"?
C
        ITYP1(IPLS)=-1
        ISPZ1(IPLS)=0

        IPLSTI = MPLSTI(IPLS)
        IPLSV = MPLSV(IPLS)
C
        IF (NPBGKP(IPLS,1).EQ.0) GOTO 1000
C
C  YES. FIND CORRESPONDING INCIDENT TEST PARTICLE COLLISION PARTNER: ITYP, ISPZ, IREL

        IBGK1=NPBGKP(IPLS,1)
        IUP1=(IBGK1-1)*3+1
        IUP2=(IBGK1-1)*3+2
        IUP3=(IBGK1-1)*3+3
C
C  TRY ATOMS
        DO IATM=1,NATMI
          IF (NPBGKA(IATM).EQ.IBGK1) THEN
            ITYP1(IPLS)=1
            ISPZ1(IPLS)=IATM
            FACT1=CVRSSA(IATM)
            RMAS1=RMASSA(IATM)
            DO I_fine=1,NRAD
!pb 05.02.2013  take care of coarse scoring cells
c    fine cell i_fine belongs to coarser cell ird. Scoring was done on coarse cell ird only
cdr  This is not obvious for additional cells for averaging, e.g. IR=NR1ST, etc..
cdr  Those IRD should not appear here  (hopefully)
              IRD = NCLTAL(I_fine)
              IF (IRD == 0) CYCLE
              PDEN(I_fine)=PDENA(IATM,IRD)
              EDEN(I_fine)=EDENA(IATM,IRD)
            ENDDO
C  FIND INDEX NRC
            DO NRC=1,NRCA(IATM)
              IP=EIRENE_IDEZ(IBULKA(IATM,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 10
              ENDIF
            ENDDO
            GOTO 995
C  FIND INDEX IREL
   10       DO IAEL=1,NAELI(IATM)
              IF (LGAEL(IATM,IAEL,1).EQ.IPLS) THEN
                IREL1(IPLS)=LGAEL(IATM,IAEL,0)
                GOTO 50
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
C  TRY MOLECULES
        DO IMOL=1,NMOLI
          IF (NPBGKM(IMOL).EQ.IBGK1) THEN
            ITYP1(IPLS)=2
            ISPZ1(IPLS)=IMOL
            FACT1=CVRSSM(IMOL)
            RMAS1=RMASSM(IMOL)
            DO I_fine=1,NRAD
!pb 05.02.2013  take care of coarse scoring cells
c     fine cell i_fine belongs to coarse cell ird. scoring was done on coarse cell ird
cdr  This is not obvious for additional cells for averaging, e.g. IR=NR1ST, etc..
cdr  Those IRD should not appear here  (hopefully)
              IRD = NCLTAL(I_fine)
              IF (IRD == 0) CYCLE
              PDEN(I_fine)=PDENM(IMOL,IRD)
              EDEN(I_fine)=EDENM(IMOL,IRD)
            ENDDO
C  FIND INDEX NRC
            DO NRC=1,NRCM(IMOL)
              IP=EIRENE_IDEZ(IBULKM(IMOL,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 20
              ENDIF
            ENDDO
            GOTO 995
C  FIND INDEX IREL
   20       DO IMEL=1,NMELI(IMOL)
              IF  (LGMEL(IMOL,IMEL,1).EQ.IPLS) THEN
                IREL1(IPLS)=LGMEL(IMOL,IMEL,0)
                GOTO 50
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
C  TRY TEST IONS
        DO IION=1,NIONI
          IF (NPBGKI(IION).EQ.IBGK1) THEN
            ITYP1(IPLS)=3
            ISPZ1(IPLS)=IION
            FACT1=CVRSSI(IION)
            RMAS1=RMASSI(IION)
            DO I_fine=1,NRAD
!pb 05.02.2013  take care of coarse scoring cells
c    fine cell i_fine belongs to coarse cell ird. scoring was done on coarse cell ird
cdr  This is not obvious for additional cells for averaging, e.g. IR=NR1ST, etc..
cdr  Those IRD should not appear here  (hopefully)
              IRD = NCLTAL(I_fine)
              IF (IRD == 0) CYCLE
              PDEN(I_fine)=PDENI(IION,IRD)
              EDEN(I_fine)=EDENI(IION,IRD)
            ENDDO
C  FIND INDEX NRC
            DO NRC=1,NRCI(IION)
              IP=EIRENE_IDEZ(IBULKI(IION,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 30
              ENDIF
            ENDDO
C  FIND INDEX IREL
   30       DO IIEL=1,NIELI(IION)
              IF  (LGIEL(IION,IIEL,1).EQ.IPLS) THEN
                IREL1(IPLS)=LGIEL(IION,IIEL,0)
                GOTO 50
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
        GOTO 995

C  AT THIS POINT: COLLISION PARTNER AMONGST TEST PARTICLES HAS BEEN IDENTIFIED
C  (ITYP1, IATM1, IMOL1, IION1)
C  AS WELL AS THE NUMBER OF COLLISION PROCESS IREL1
C  DENSITY AND ENERGY DENSITY TALLIES OF COLLISION PARTNER "PDEN,EDEN"
C  HAVE NOW BEEN SET ON FINE GRID
C
C  SELF-COLLISION OR CROSS-COLLISION
C
   50   CONTINUE
C
        IF (NPBGKP(IPLS,2).EQ.0) THEN
          ITYP2(IPLS)=-1
          ISPZ2(IPLS)=0
C
        ELSEIF (NPBGKP(IPLS,2).NE.0) THEN
C
C  CROSS-COLLISION, FIND SECOND COLLISION PARTNER
C  THIS IS NOT THE INGOING COLLIDING TEST PARTICLE, WHICH WE HAVE ALREADY IDENTIFIED,
C  (AND WHICH, E.G., DETERMINES MASS AND DENSITY OF ARTIFICIAL BACKGROUND PARTICLE IPLS)
C  BUT, INSTEAD, IT IS THE TEST PARTICLE WHICH PLAYS THE ROLE
c  OF THE "SECOND" PARTICLE, AMONGST THE TEST PARTICLES
C
          ITYP2(IPLS)=EIRENE_IDEZ(NPBGKP(IPLS,2),1,3)
          ISPZ2(IPLS)=EIRENE_IDEZ(NPBGKP(IPLS,2),3,3)
C
          IF (ITYP2(IPLS).EQ.1) THEN
            IATM2=ISPZ2(IPLS)
            FACT2=CVRSSA(IATM2)
            RMAS2=RMASSA(IATM2)
            DO I_fine=1,NRAD
!pb 05.02.2013  take care of scoring cells
cdr  This is not obvious for additional cells for averaging, e.q. IR=NR1ST, etc..
cdr  Those IRD should not appear here  (hopefully)
              IRD = NCLTAL(I_fine)
              IF (IRD == 0) CYCLE
              PDEN2(I_fine)=PDENA(IATM2,IRD)
              EDEN2(I_fine)=EDENA(IATM2,IRD)
            ENDDO
            IBGK2=NPBGKA(IATM2)
          ELSEIF (ITYP2(IPLS).EQ.2) THEN
            IMOL2=ISPZ2(IPLS)
            FACT2=CVRSSM(IMOL2)
            RMAS2=RMASSM(IMOL2)
            DO I_fine=1,NRAD
!pb 05.02.2013  take care of scoring cells
cdr  This is not obvious for additional cells for averaging, e.g. IRD=NR1ST, etc..
cdr  Those IRD should not appear here  (hopefully)
              IRD = NCLTAL(I_fine)
              IF (IRD == 0) CYCLE
              PDEN2(I_fine)=PDENM(IMOL2,IRD)
              EDEN2(I_fine)=EDENM(IMOL2,IRD)
            ENDDO
            IBGK2=NPBGKM(IMOL2)
          ELSEIF (ITYP2(IPLS).EQ.3) THEN
            IION2=ISPZ2(IPLS)
            FACT2=CVRSSI(IION2)
            RMAS2=RMASSI(IION2)
            DO I_fine=1,NRAD
!pb 05.02.2013  take care of scoring cells
cdr  This is not obvious for additional cells for averaging, e.g. IR=NR1ST, etc..
cdr  Those IRD should not appear here  (hopefully)
              IRD = NCLTAL(I_fine)
              IF (IRD == 0) CYCLE
              PDEN2(I_fine)=PDENI(IION2,IRD)
              EDEN2(I_fine)=EDENI(IION2,IRD)
            ENDDO
            IBGK2=NPBGKI(IION2)
          ENDIF
          IUP12=(IBGK2-1)*3+1
          IUP22=(IBGK2-1)*3+2
          IUP32=(IBGK2-1)*3+3
C
        ENDIF
c
cdr:  Parameters of virtual background species are set,
cdr   their densities, energy densities.
cdr   are on 1D arrays: pden,eden, pden2, eden2.
cdr   momentum densities are still on 2d array GBGKV(ibgk, icell).
C
        IF (RMASSP(IPLS).NE.RMAS1) THEN
          RM=RMAS1
          WRITE (iunout,*) 'MODBGK: INCONSISTENT MASS FOR IPLS= ',IPLS
          WRITE (iunout,*) '        RMASSP(IPLS),RM= ',RMASSP(IPLS),RM
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
        CNDYN=AMUA*RMAS1
C
        RRN=0.
        RRE=0.
        RRM=0.
        RATN=0.
        RATE=0.
        RATM=0.
        RESN=0.
        RESE=0.
        RESM=0.
C
        IREL=IREL1(IPLS)
C
        IF (ITYP2(IPLS).EQ.-1) THEN
C
        IF (TRCMOD) THEN
          WRITE (iunout,*) 'MODBGK: SELF-COLLISION WITH, IPLS',IPLS
          WRITE (iunout,*) 'ITYP,ISPZ,IBGK_SP,IREL ',ITYP1(IPLS),
     .                                ISPZ1(IPLS),IBGK1,IREL1(IPLS)
        ENDIF
C


C  this cell loop is referring to the underlying 'fine' grid, not the coarse grid
C                              on which the eirene tallies had been updated.

        NXM=MAX(1,NR1STM)
        NYM=MAX(1,NP2NDM)
        NZM=MAX(1,NT3RDM)

        DO IR=1,NXM
          DO IP=1,NYM
            DO IT=1,NZM
              IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA  !irad=i_fine
C
              TBEL=0.
              IF (LGVAC(IRAD,IPLS)) GOTO 81
              IF (NSTORDR >= NRAD) THEN
                TBEL = TABEL3(IREL,IRAD,1)
              ELSE
cdr             TBEL=EIRENE_FTABEL3(IREL,IRAD)  ! this should replace the next three cards
                KK=NREAEL(IREL)
                TII=TIINL(IPLSTI,IRAD)+ADDEL(IREL,IPLS)
                TBEL = EIRENE_RATE_COEFF(KK,IRAD,TII,0._DP,.TRUE.,0)
     .                 *DIIN(IPLS,IRAD)*FACREL(IREL,1)
              END IF
   81         CONTINUE
C DELTA_N
              DOLD=DIIN(IPLS,IRAD)
              DEL=DOLD-PDEN(IRAD)

c  rate of particle exchange: 1/s, per cell I_fine
              RATN=RATN+TBEL*DEL*VOL(IRAD)
c  L1 norm of particle residual
              RESN=RESN+TBEL*ABS(DEL)*VOL(IRAD)

c             IF (TRCMOD) THEN
c               WRITE (iunout,*) 'IR,T,TBEL,RATN ',
c    .                            IRAD,TIIN(IPLSTI,IRAD),
c    .                               TBEL,TBEL*DEL*VOL(IRAD)*ELCHA
c               CALL EIRENE_MASR3('DOLD,DNEW,DEL           ',
c    .                      DOLD,PDEN(IRAD),DEL)
c             ENDIF

C DELTA_E
c             EOLD=(1.5*TIIN(IPLSTI,IRAD)+EDRIFT(IPLS,IRAD))*
c    .             DIIN(IPLS,IRAD)
!pb              EOLD=(1.5*TIIN(IPLSTI,IRAD)+EDRIFT(IPLS,IRAD))*
!pb     .             PDEN(IRAD)
              EOLD=1.5*TIIN(IPLSTI,IRAD)*PDEN(IRAD)
              IF (LEDRIFT) EOLD=EOLD+EDRIFT(IPLS,IRAD)*PDEN(IRAD)
              DEL=EOLD-EDEN(IRAD)

c  rate of energy exchange: (eV)/s, per cell I_fine
              RATE=RATE+TBEL*DEL*VOL(IRAD)
c  L1 norm of energy residual
              RESE=RESE+TBEL*ABS(DEL)*VOL(IRAD)
C DELTA_V
CDR  next lines: gbgkv (=bgkv) is redundant, because momentum density vector
cdr              components have become default volume-averaged output tallies
              DELX=GBGKV(IUP1,IRAD)-VXIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
              DELY=GBGKV(IUP2,IRAD)-VYIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
              DELZ=GBGKV(IUP3,IRAD)-VZIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)

c  cdr  rate of momentum exchange: (g*cm/s)/s, per cell I_fine
              RATM(1)=RATM(1)+TBEL*DELX*VOL(IRAD)
              RATM(2)=RATM(2)+TBEL*DELY*VOL(IRAD)
              RATM(3)=RATM(3)+TBEL*DELZ*VOL(IRAD)
C NEW V
              VX=GBGKV(IUP1,IRAD)/(PDEN(IRAD)+EPS60)
              VY=GBGKV(IUP2,IRAD)/(PDEN(IRAD)+EPS60)
              VZ=GBGKV(IUP3,IRAD)/(PDEN(IRAD)+EPS60)
              VXIN(IPLSV,IRAD)=VX
              VYIN(IPLSV,IRAD)=VY
              VZIN(IPLSV,IRAD)=VZ
C NEW T
              ED=(VX**2+VY**2+VZ**2)*FACT1
              TIIN(IPLSTI,IRAD)=(EDEN(IRAD)/(PDEN(IRAD)+EPS60)-ED)/1.5
C NEW N
              DIIN(IPLS,IRAD)=PDEN(IRAD)

C  For  normalization of global residuals:
c  RRN: total particle [1]
c  RRM: total momentum content
c  RRE: total energy content (eV)
              RRN=RRN+PDEN(IRAD)*VOL(IRAD)
C             RRM=?
              RRE=RRE+EDEN(IRAD)*VOL(IRAD)
            END DO
          END DO
        END DO
c
c  same as do loop above, for additional cell region
c
        DO 90 IRAD=NSURF+1,NSURF+NRADD
C
          TBEL=0.
          IF (LGVAC(IRAD,IPLS)) GOTO 91
          IF (NSTORDR >= NRAD) THEN
            TBEL = TABEL3(IREL,IRAD,1)
          ELSE
cdr         TBEL=EIRENE_FTABEL3(IREL,IRAD)  ! this should replace the next three cards
            KK=NREAEL(IREL)
            TII=TIINL(IPLSTI,IRAD)+ADDEL(IREL,IPLS)
            TBEL = EIRENE_RATE_COEFF(KK,IRAD,TII,0._DP,.TRUE.,0)
     .             *DIIN(IPLS,IRAD)*FACREL(IREL,1)
          END IF
   91     CONTINUE
C DELTA_N
          DOLD=DIIN(IPLS,IRAD)
          DEL=DOLD-PDEN(IRAD)
c  rate of particle exchange: 1/s, per cell I_fine
          RATN=RATN+TBEL*DEL*VOL(IRAD)
c  L1 norm of particle residual
          RESN=RESN+TBEL*ABS(DEL)*VOL(IRAD)

c         IF (TRCMOD) THEN
c           WRITE (iunout,*) 'IR,T,TBEL,RATN ',IRAD,TIIN(IPLSTI,IRAD),
c    .                           TBEL,TBEL*DEL*VOL(IRAD)*ELCHA
c           CALL EIRENE_MASR3('DOLD,DNEW,DEL           ',
c    .                  DOLD,PDEN(IRAD),DEL)
c         ENDIF

C DELTA_E
c         EOLD=(1.5*TIIN(IPLSTI,IRAD)+EDRIFT(IPLS,IRAD))*
c    .          DIIN(IPLS,IRAD)
!pb          EOLD=(1.5*TIIN(IPLSTI,IRAD)+EDRIFT(IPLS,IRAD))*
!pb     .          PDEN(IRAD)
          EOLD=1.5*TIIN(IPLSTI,IRAD)*PDEN(IRAD)
          IF (LEDRIFT) EOLD=EOLD+EDRIFT(IPLS,IRAD)*PDEN(IRAD)
          DEL=EOLD-EDEN(IRAD)
c  cdr  rate of energy exchange: (eV)/s, per cell I_fine
          RATE=RATE+TBEL*DEL*VOL(IRAD)
c  L1 norm of energy residual
          RESE=RESE+TBEL*ABS(DEL)*VOL(IRAD)
C DELTA_V
          DELX=GBGKV(IUP1,IRAD)-VXIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
          DELY=GBGKV(IUP2,IRAD)-VYIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
          DELZ=GBGKV(IUP3,IRAD)-VZIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)

c  cdr  rate of momentum exchange: (g*cm/s)/s, per cell I_fine
          RATM(1)=RATM(1)+TBEL*DELX*VOL(IRAD)
          RATM(2)=RATM(2)+TBEL*DELY*VOL(IRAD)
          RATM(3)=RATM(3)+TBEL*DELZ*VOL(IRAD)
C NEW V
          VX=GBGKV(IUP1,IRAD)/(PDEN(IRAD)+EPS60)
          VY=GBGKV(IUP2,IRAD)/(PDEN(IRAD)+EPS60)
          VZ=GBGKV(IUP3,IRAD)/(PDEN(IRAD)+EPS60)
          VXIN(IPLSV,IRAD)=VX
          VYIN(IPLSV,IRAD)=VY
          VZIN(IPLSV,IRAD)=VZ
C NEW T
          ED=(VX**2+VY**2+VZ**2)*FACT1
          TIIN(IPLSTI,IRAD)=(EDEN(IRAD)/(PDEN(IRAD)+EPS60)-ED)/1.5
C NEW N
          DIIN(IPLS,IRAD)=PDEN(IRAD)

C FOR GLOBAL BALANCES
          RRN=RRN+PDEN(IRAD)*VOL(IRAD)
C         RRM=?
          RRE=RRE+EDEN(IRAD)*VOL(IRAD)
   90   CONTINUE
C
C   IF IPLS IS AN ARTIFICIAL BACKGROUND SPECIES FOR A SELF-COLLISION  :  DONE !
C
C   NEXT, IF IPLS IS AN ARTIFICIAL BACKGROUND SPECIES FOR A CROSS-COLLISION BETWEEN TWO DISTINCT
C        SPECIES OF TEST PARTICLES
        ELSE
C
          IF (TRCMOD) THEN
            WRITE (iunout,*) 'MODBGK: CROSS-COLLISION, IPLS ',IPLS
            WRITE (iunout,*) 'ITYP1,ISPZ1,IBGK1,IREL1 ',
     .                        ITYP1(IPLS),ISPZ1(IPLS),
     .                        IBGK1,IREL1(IPLS)
            WRITE (iunout,*) 'ITYP2,ISPZ2,IBGK2       ',
     .                        ITYP2(IPLS),ISPZ2(IPLS),IBGK2
          ENDIF

CVK FOR "CROSS-COLLISION TEMPERATURE"
          ICROSS=ICROSS+1
          CROSSINDEX(IPLS)=ICROSS
          RM1=RMAS1/(RMAS1+RMAS2)
          RM2=RMAS2/(RMAS1+RMAS2)
          RM=2.D0*RMAS1*RMAS2/(RMAS1+RMAS2)**2
CVK END

        DO IR=1,NXM
          DO IP=1,NYM
            DO IT=1,NZM
              IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
C
              TBEL=0
              IF (LGVAC(IRAD,IPLS)) GOTO 181
              IF (NSTORDR >= NRAD) THEN
                TBEL = TABEL3(IREL,IRAD,1)
              ELSE
cdr             TBEL=EIRENE_FTABEL3(IREL,IRAD)  ! this should replace the next three cards
                KK=NREAEL(IREL)
                TII=TIINL(IPLSTI,IRAD)+ADDEL(IREL,IPLS)
                TBEL = EIRENE_RATE_COEFF(KK,IRAD,TII,0._DP,.TRUE.,0)
     .                 *DIIN(IPLS,IRAD)*FACREL(IREL,1)
              END IF
  181         CONTINUE
c             EOLD=(1.5*TIIN(IPLSTI,IRAD)+EDRIFT(IPLS,IRAD))*
c    .              DIIN(IPLS,IRAD)
!pb              EOLD=(1.5*TIIN(IPLSTI,IRAD)+EDRIFT(IPLS,IRAD))*
!pb     .              PDEN(IRAD)
              EOLD=1.5*TIIN(IPLSTI,IRAD)*PDEN(IRAD)
              IF (LEDRIFT) EOLD=EOLD+EDRIFT(IPLS,IRAD)*PDEN(IRAD)
              DOLD=DIIN(IPLS,IRAD)
              DEL=EOLD-EDEN(IRAD)
              RATE=RATE+TBEL*DEL*VOL(IRAD)
              DELX=GBGKV(IUP1,IRAD)-VXIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
              DELY=GBGKV(IUP2,IRAD)-VYIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
              DELZ=GBGKV(IUP3,IRAD)-VZIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
              RATM(1)=RATM(1)+TBEL*DELX*VOL(IRAD)
              RATM(2)=RATM(2)+TBEL*DELY*VOL(IRAD)
              RATM(3)=RATM(3)+TBEL*DELZ*VOL(IRAD)
C
              VXIN1=GBGKV(IUP1 ,IRAD)/(PDEN (IRAD)+EPS60)
              VXIN2=GBGKV(IUP12,IRAD)/(PDEN2(IRAD)+EPS60)
              VXMIX=(RMAS1*VXIN1+RMAS2*VXIN2)/(RMAS1+RMAS2)
              VYIN1=GBGKV(IUP2 ,IRAD)/(PDEN (IRAD)+EPS60)
              VYIN2=GBGKV(IUP22,IRAD)/(PDEN2(IRAD)+EPS60)
              VYMIX=(RMAS1*VYIN1+RMAS2*VYIN2)/(RMAS1+RMAS2)
              VZIN1=GBGKV(IUP3 ,IRAD)/(PDEN (IRAD)+EPS60)
              VZIN2=GBGKV(IUP32,IRAD)/(PDEN2(IRAD)+EPS60)
              VZMIX=(RMAS1*VZIN1+RMAS2*VZIN2)/(RMAS1+RMAS2)
C  SET NEW VELOCITY OF ARTIFICIAL BACKGROUND SPECIES
              VXIN(IPLSV,IRAD)=VXMIX
              VYIN(IPLSV,IRAD)=VYMIX
              VZIN(IPLSV,IRAD)=VZMIX

              ED1=(VXIN1**2+VYIN1**2+VZIN1**2)*FACT1
              ED2=(VXIN2**2+VYIN2**2+VZIN2**2)*FACT2
              T1=(EDEN(IRAD)/(PDEN(IRAD)+EPS60)-ED1)/1.5
              T2=(EDEN2(IRAD)/(PDEN2(IRAD)+EPS60)-ED2)/1.5
              CROSSTEMP(ICROSS,IRAD)=T1*RM2+T2*RM1 !VK THE "COLLISION RATE" TEMPERATURE

              TMIX=T1+RM*(T2-T1+
     .             FACT2/3.D0*((VXIN1-VXIN2)**2+(VYIN1-VYIN2)**2+
     .                         (VZIN1-VZIN2)**2))
              TIIN(IPLSTI,IRAD)=TMIX

              DEL=DIIN(IPLS,IRAD)-PDEN(IRAD)
              RATN=RATN+TBEL*DEL*VOL(IRAD)
              RESN=RESN+TBEL*ABS(DEL)*VOL(IRAD)

C  SET NEW DENSITY OF ARTIFICIAL BACKGROUND SPECIES
              DIIN(IPLS,IRAD)=PDEN(IRAD)

              RRN=RRN+PDEN(IRAD)*VOL(IRAD)
C             RRM=?
              RRE=RRE+EDEN(IRAD)*VOL(IRAD)
            END DO
          END DO
        END DO
c
c  same as do loop above, for additional cell region
c
        DO 190 IRAD=NSURF+1,NSURF+NRADD
C
          TBEL=0
          IF (LGVAC(IRAD,IPLS)) GOTO 191
          IF (NSTORDR >= NRAD) THEN
            TBEL = TABEL3(IREL,IRAD,1)
          ELSE
cdr         TBEL=EIRENE_FTABEL3(IREL,IRAD)  ! this should replace the next three cards
            KK=NREAEL(IREL)
            TII=TIINL(IPLSTI,IRAD)+ADDEL(IREL,IPLS)
            TBEL = EIRENE_RATE_COEFF(KK,IRAD,TII,0._DP,.TRUE.,0)
     .             *DIIN(IPLS,IRAD)*FACREL(IREL,1)
          END IF
  191     CONTINUE
C         EOLD=(1.5*TIIN(IPLSTI,IRAD)+EDRIFT(IPLS,IRAD))*
C    .          DIIN(IPLS,IRAD)
!pb          EOLD=(1.5*TIIN(IPLSTI,IRAD)+EDRIFT(IPLS,IRAD))*
!pb     .          PDEN(IRAD)
          EOLD=1.5*TIIN(IPLSTI,IRAD)*PDEN(IRAD)
          IF (LEDRIFT) EOLD=EOLD+EDRIFT(IPLS,IRAD)*PDEN(IRAD)
          DOLD=DIIN(IPLS,IRAD)
          DEL=EOLD-EDEN(IRAD)
          RATE=RATE+TBEL*DEL*VOL(IRAD)
          DELX=GBGKV(IUP1,IRAD)-VXIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
          DELY=GBGKV(IUP2,IRAD)-VYIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
          DELZ=GBGKV(IUP3,IRAD)-VZIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
          RATM(1)=RATM(1)+TBEL*DELX*VOL(IRAD)
          RATM(2)=RATM(2)+TBEL*DELY*VOL(IRAD)
          RATM(3)=RATM(3)+TBEL*DELZ*VOL(IRAD)
C
          VXIN1=GBGKV(IUP1 ,IRAD)/(PDEN (IRAD)+EPS60)
          VXIN2=GBGKV(IUP12,IRAD)/(PDEN2(IRAD)+EPS60)
          VXMIX=(RMAS1*VXIN1+RMAS2*VXIN2)/(RMAS1+RMAS2)
          VYIN1=GBGKV(IUP2 ,IRAD)/(PDEN (IRAD)+EPS60)
          VYIN2=GBGKV(IUP22,IRAD)/(PDEN2(IRAD)+EPS60)
          VYMIX=(RMAS1*VYIN1+RMAS2*VYIN2)/(RMAS1+RMAS2)
          VZIN1=GBGKV(IUP3 ,IRAD)/(PDEN (IRAD)+EPS60)
          VZIN2=GBGKV(IUP32,IRAD)/(PDEN2(IRAD)+EPS60)
          VZMIX=(RMAS1*VZIN1+RMAS2*VZIN2)/(RMAS1+RMAS2)

C  SET NEW VELOCITY OF ARTIFICIAL BACKGROUND SPECIES
          VXIN(IPLSV,IRAD)=VXMIX
          VYIN(IPLSV,IRAD)=VYMIX
          VZIN(IPLSV,IRAD)=VZMIX

          ED1=(VXIN1**2+VYIN1**2+VZIN1**2)*FACT1
          ED2=(VXIN2**2+VYIN2**2+VZIN2**2)*FACT2
          T1=(EDEN(IRAD)/(PDEN(IRAD)+EPS60)-ED1)/1.5
          T2=(EDEN2(IRAD)/(PDEN2(IRAD)+EPS60)-ED2)/1.5
          CROSSTEMP(ICROSS,IRAD)=T1*RM2+T2*RM1 !VK THE "COLLISION RATE" TEMPERATURE

          TMIX=T1+RM*(T2-T1+
     .         FACT2/3.D0*((VXIN1-VXIN2)**2+(VYIN1-VYIN2)**2+
     .                     (VZIN1-VZIN2)**2))
          TIIN(IPLSTI,IRAD)=TMIX

          DEL=DIIN(IPLS,IRAD)-PDEN(IRAD)
          RATN=RATN+TBEL*DEL*VOL(IRAD)
          RESN=RESN+TBEL*ABS(DEL)*VOL(IRAD)
C  SET NEW DENSITY OF ARTIFICIAL BACKGROUND SPECIES
          DIIN(IPLS,IRAD)=PDEN(IRAD)

          RRN=RRN+PDEN(IRAD)*VOL(IRAD)
C         RRM=?
          RRE=RRE+EDEN(IRAD)*VOL(IRAD)
C
  190   CONTINUE
c
        ENDIF
c
        CALL EIRENE_LEER(2)
        WRITE (iunout,*) 'PARTICLE, MOMENTUM AND ENERGY EXCHANGE RATES'
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASR1('RATN=   ',RATN*ELCHA)
        WRITE (iunout,'(1X,A8,3X,3(1PE12.4))') 'RATM=   ',
     .                    RATM(1)*ELCHA*CNDYN,
     .                    RATM(2)*ELCHA*CNDYN,RATM(3)*ELCHA*CNDYN
        CALL EIRENE_MASR1('RATE=   ',RATE*ELCHA)
        CALL EIRENE_LEER(2)
        WRITE (iunout,*) 'RESIDUA (1/SEC)'
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASR1('RESN=   ',RESN/(RRN+EPS60))
C       CALL EIRENE_MASR3('RESM=                   ',RATM(1)/(RRM+EPS60),
C    .                   ,RATM(2)/(RRM+EPS60),RATM(3)/(RRM+EPS60))
        CALL EIRENE_MASR1('RESE=   ',RESE/(RRE+EPS60))
        CALL EIRENE_LEER(2)
C
C.................................................................
 1000 CONTINUE
C.................................................................

      CALL EIRENE_LEER(2)
C
C  SAVE OVERHEAD, IF GEOMETRY DATA ALREADY AVAILABLE ON FILE
C
      IF (NFILEM.EQ.1) NFILEM=2
C
C  SET INDPRO=7, AND
C  WRITE PLASMA DATA ONTO PLASMA_BCKGRND FOR CALL TO SUBR. PLASMA BELOW
C  TI,NI AND (VX,VY,VZ) FOR IPLS=1,NPLSI
C  PLAY SAVE: WRITE WHOLE PLASMA_BCKGRND ARRAY.
C

      DO 500 I=1,6
        INDPRO(I)=7
  500 CONTINUE
! STORAGE FOR INPUT TALLIES 1 (TEIN) TO 13 (ADIN), WITHOUT NO.3 (DEIN)
      NRWK1=6+NPLS+NPLSTI+3*NPLSV+NAIN  ! STORAGE FOR INPUT TALLIES 1 TO 13, WITHOUT NO.3
!pb      IF (NIDV < NRWK1) THEN
      IF (NIDC < NRWK1) THEN
        WRITE (iunout,*) ' PLASMA_BCKGRND ARRAY IS TOO SMALL TO HOLD'
        WRITE (iunout,*) ' PLASMA DATA'
        WRITE (iunout,*) ' CHECK PARAMETER NSMSTRA'
        CALL EIRENE_EXIT_OWN(1)
      END IF
      CALL EIRENE_ALLOC_BCKGRND
      PLASMA_BCKGRND(1:NRWK1,:) = 0.D0

cdr
cdr initialize BXIN=0
      PLASMA_BCKGRND(1+1*NPLS+NPLSTI+3*NPLSV+1,:)= 0._DP  
cdr initialize BYIN=0
      PLASMA_BCKGRND(2+1*NPLS+NPLSTI+3*NPLSV+1,:)= 0._DP  
!pb initialize BZIN=1
      PLASMA_BCKGRND(3+1*NPLS+NPLSTI+3*NPLSV+1,:)= 1._DP
!pb initialize BFIN=1
      PLASMA_BCKGRND(4+1*NPLS+NPLSTI+3*NPLSV+1,:)= 1._DP

      DO IR=1,NXM
        DO IP=1,NYM
          DO IT=1,NZM
            IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
            PLASMA_BCKGRND  (0+0*NPLS+1,IRAD)= TEIN(IRAD)
            DO IPLSTI=1,NPLSTI
              PLASMA_BCKGRND(1+0*NPLS+IPLSTI,IRAD)= TIIN(IPLSTI,IRAD)
            END DO
            DO IPLS=1,NPLS
              PLASMA_BCKGRND(1+0*NPLS+NPLSTI+IPLS,IRAD)= DIIN(IPLS,IRAD)
            ENDDO
            DO IPLSV=1,NPLSV
              PLASMA_BCKGRND(1+1*NPLS+NPLSTI+0*NPLSV+IPLSV,IRAD)=
     .               VXIN(IPLSV,IRAD)
              PLASMA_BCKGRND(1+1*NPLS+NPLSTI+1*NPLSV+IPLSV,IRAD)=
     .               VYIN(IPLSV,IRAD)
              PLASMA_BCKGRND(1+1*NPLS+NPLSTI+2*NPLSV+IPLSV,IRAD)=
     .               VZIN(IPLSV,IRAD)
            END DO
            IF (LBXIN)
     .        PLASMA_BCKGRND(1+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BXIN(IRAD)
            IF (LBYIN)
     .        PLASMA_BCKGRND(2+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BYIN(IRAD)
            IF (LBZIN)
     .        PLASMA_BCKGRND(3+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BZIN(IRAD)
            IF (LBFIN)
     .        PLASMA_BCKGRND(4+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BFIN(IRAD)

            IF (LVOL)
     .        PLASMA_BCKGRND(5+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= VOL(IRAD)

            IF (LADIN) THEN
              DO IAIN=1,NAINI
                PLASMA_BCKGRND(6+1*NPLS+NPLSTI+3*NPLSV+IAIN,IRAD)=
     .               ADIN(IAIN,IRAD)
              ENDDO
            END IF
          END DO
        END DO
      END DO
C
c  same as do loop above, for additional cell region
c
      DO IRAD=NSURF+1,NSURF+NRADD
        PLASMA_BCKGRND  (0+0*NPLS+1   ,IRAD)= TEIN(IRAD)
        DO IPLSTI=1,NPLSTI
          PLASMA_BCKGRND(1+0*NPLS+IPLSTI,IRAD)= TIIN(IPLSTI,IRAD)
        END DO
        DO IPLS=1,NPLSI
          PLASMA_BCKGRND(1+0*NPLS+NPLSTI+IPLS,IRAD)= DIIN(IPLS,IRAD)
        ENDDO
        DO IPLSV=1,NPLSV
          PLASMA_BCKGRND(1+1*NPLS+NPLSTI+0*NPLSV+IPLSV,IRAD)=
     .               VXIN(IPLSV,IRAD)
          PLASMA_BCKGRND(1+1*NPLS+NPLSTI+1*NPLSV+IPLSV,IRAD)=
     .               VYIN(IPLSV,IRAD)
          PLASMA_BCKGRND(1+1*NPLS+NPLSTI+2*NPLSV+IPLSV,IRAD)=
     .               VZIN(IPLSV,IRAD)
        END DO

        IF (LBXIN)
     .    PLASMA_BCKGRND(1+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BXIN(IRAD)
        IF (LBYIN)
     .    PLASMA_BCKGRND(2+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BYIN(IRAD)
        IF (LBZIN)
     .    PLASMA_BCKGRND(3+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BZIN(IRAD)
        IF (LBFIN)
     .    PLASMA_BCKGRND(4+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BFIN(IRAD)

        IF (LVOL)
     .    PLASMA_BCKGRND(5+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= VOL(IRAD)

        IF (LADIN) THEN
          DO IAIN=1,NAINI
            PLASMA_BCKGRND(6+1*NPLS+NPLSTI+3*NPLSV+IAIN,IRAD)=
     .               ADIN(IAIN,IRAD)
          END DO
        END IF
      ENDDO
C
      CALL EIRENE_PLASMA_DERIV(0)
C
C .........................................................................
C  NOW: NEW COLLISION RATES MUST BE SET FOR THE NEXT ITERATION
C .........................................................................
C
C  IN CASE OF CROSS-COLLISION, SOME MODIFICATIONS OF THE
C  BACKGROUND PARAMETERS ARE REQUIRED TEMPORARILY TO ENFORCE
C  A SPECIFIC RELATION BETWEEN TAU_1,2 AND TAU_2,1
C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE PROFILES
C  E.G.: EDRIFT, NEEDED FOR EPLEL3
C
C
      TRCSAV=TRCAMD
      TRCAMD=.FALSE.
C
      CALL EIRENE_LEER(1)
      DO IPLS=1,NPLSI
        LMARK(IPLS)=.FALSE.
      ENDDO
      DO IPLS1=1,NPLSI
        IF (NPBGKP(IPLS1,2).NE.0) THEN
C  IPLS1 IS A CROSS-COLLISION TALLY
C  FIND CORRESPONDING 2ND CROSS-COLLISION TALLY
          IPLS2=0
          DO IPLS=1,NPLSI
            IF (ITYP1(IPLS).EQ.ITYP2(IPLS1).AND.
     .          ISPZ1(IPLS).EQ.ISPZ2(IPLS1).AND.
     .          ITYP2(IPLS).EQ.ITYP1(IPLS1).AND.
     .          ISPZ2(IPLS).EQ.ISPZ1(IPLS1)) IPLS2=IPLS
          ENDDO
          IF (IPLS2.EQ.0) GOTO 800
          CALL EIRENE_LEER(1)
          IF (TRCMOD) THEN
            WRITE (iunout,*)
     .        'MODBGK: CORRESPONDING CROSS-COLLISION SPECIES'
            WRITE (iunout,*) 'IPLS1,IPLS2 ',IPLS1,IPLS2
          ENDIF
          IF (LMARK(IPLS1).OR.LMARK(IPLS2)) GOTO 800
C  IPLS2 IS THE SECOND CROSS-COLLISION TALLY
          IF (TRCMOD) THEN
            WRITE (iunout,*)
     .        'MODBGK: MODIFY PARAMETERS FOR CROSS-COLLISIONALITIES'
            WRITE (iunout,*) 'IPLS1,IPLS2 ',IPLS1,IPLS2
            CALL EIRENE_LEER(1)
          ENDIF

          LMARK(IPLS1)=.TRUE.
          LMARK(IPLS2)=.TRUE.
CVK
          ICROSS1=CROSSINDEX(IPLS1)
          ICROSS2=CROSSINDEX(IPLS2)
          IF(ABS(CROSSTEMP(ICROSS1,1)-CROSSTEMP(ICROSS2,1)).GT.EPS12)
     .      WRITE(iunout,*) "WARNING FROM MODBGK: ",
     .                 "CROSS-COLLISION TEMPERATURES ARE WRONG ",
     .                 "IPLS1,IPLS2,CROSSIND1,CROSSIND2",
     .                  IPLS1,IPLS2,CROSSINDEX(IPLS1),CROSSINDEX(IPLS2)
          ICROSS=ICROSS1
          TCSUM=0
CVK END
          IPLSTI1 = MPLSTI(IPLS1)
          IPLSTI2 = MPLSTI(IPLS2)
          DO IRAD=1,NSBOX
            DS1=DIIN(IPLS1,IRAD)
            DIIN(IPLS1,IRAD)=DIIN(IPLS2,IRAD)
            DIIN(IPLS2,IRAD)=DS1
CVK PLASMA_BACKGROUND CORRECTION
            DIINTF(IPLS1,IRAD)=DIIN(IPLS1,IRAD)
            DIINTF(IPLS2,IRAD)=DIIN(IPLS2,IRAD)
CSW 02jan2012 CHECK FOLLOWING TWO LINES !!!! ipls1/ipls2 mixed up?
            TIINTF(IPLSTI1,IRAD)=TIIN(IPLSTI2,IRAD) ! HAS BEEN CHANGED IPLS1->IPLS2
            TIINTF(IPLSTI2,IRAD)=TIIN(IPLSTI1,IRAD) !
CVK END
C
C
            ENERGY(IPLS1,IRAD)=1.5*TIIN(IPLSTI1,IRAD)
            IF (LEDRIFT) ENERGY(IPLS1,IRAD)=ENERGY(IPLS1,IRAD)
     .                                      +EDRIFT(IPLS1,IRAD)

            ENERGY(IPLS2,IRAD)=1.5*TIIN(IPLSTI2,IRAD)
            IF (LEDRIFT) ENERGY(IPLS2,IRAD)=ENERGY(IPLS2,IRAD)
     .                                      +EDRIFT(IPLS2,IRAD)
C
C EFFECTIVE TEMPERATURE FOR THE CALCULATION OF THE MUTUAL REACTION RATES
C [L.H. Holway, Phys. Fluids, vol. 9,pp.1658-1673,1966]
C           TS1=0.5*(TIIN(IPLS1,IRAD)+TIIN(IPLS2,IRAD))
C
            TS1=CROSSTEMP(ICROSS,IRAD) !VK
            TIIN(IPLSTI1,IRAD)=TS1
            TIIN(IPLSTI2,IRAD)=TS1

            TCSUM=TCSUM+TS1*VOL(IRAD)   !VK FOR DIAGNOSTIC
          ENDDO
          IF (TRCMOD)
     f     WRITE(iunout,*) "MODBGK: AVERAGE CROSS-COLLISION TEMPERATURE"
     .                     ," IPLSTI1, IPLSTI2",
     .                        IPLSTI1, IPLSTI2,TCSUM/SUM(VOL)

  800     CONTINUE
        ENDIF
      ENDDO
      CALL EIRENE_LEER(2)
C
C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE MODIFIED PROFILES
C
      CALL EIRENE_PLASMA_DERIV(0)
C
C  RESET BGK ATOMIC AND MOLECULAR DATA ARRAYS
C
      DO IPLS=1,NPLSI
        IF (NPBGKP(IPLS,1).NE.0) THEN
          ITYP=ITYP1(IPLS)
          ISPZ=ISPZ1(IPLS)
          IREL=IREL1(IPLS)
          NRC=INRC1(IPLS)
          IF (ITYP.EQ.1) THEN
            ISP=NSPH+ISPZ
            KK=IREACA(ISPZ,NRC)
            EBULK=EBULKA(ISPZ,NRC)
            ISCDE=ISCDEA(ISPZ,NRC)
            IESTM=IESTMA(ISPZ,NRC)
            FACTKK=FREACA(ISPZ,NRC)
          ELSEIF (ITYP.EQ.2) THEN
            ISP=NSPA+ISPZ
            KK=IREACM(ISPZ,NRC)
            EBULK=EBULKM(ISPZ,NRC)
            ISCDE=ISCDEM(ISPZ,NRC)
            IESTM=IESTMM(ISPZ,NRC)
            FACTKK=FREACM(ISPZ,NRC)
          ELSEIF (ITYP.EQ.3) THEN
            ISP=NSPAM+ISPZ
            KK=IREACI(ISPZ,NRC)
            EBULK=EBULKI(ISPZ,NRC)
            ISCDE=ISCDEI(ISPZ,NRC)
            IESTM=IESTMI(ISPZ,NRC)
            FACTKK=FREACI(ISPZ,NRC)
          ENDIF
          IF (FACTKK.EQ.0.D0) FACTKK=1.D0
C  BGK COLLISION, RESET TABEL3, EPLEL3
          CALL EIRENE_XSTEL(IREL,ISP,IPLS,EBULK,ISCDE,IESTM,
     .                      KK,FACTKK,PLS)
          IF (NPBGKP(IPLS,2).NE.0) THEN
C  CROSS-COLLISION, RESET EPLEL3 FOR TRACKLENGTH ESTIMATOR
            IF (NSTORDR >= NRAD) THEN
              DO J=1,NSBOX
                EPLEL3(IREL,J,1)=ENERGY(IPLS,J)
              ENDDO
            ELSE
              NELREL(IREL)=-3
            END IF
          ENDIF
        ENDIF
      ENDDO
C
      TRCAMD=TRCSAV
C
C
C  RESTORE PLASMA DATA FROM PLASMA_BCKGRND ARRAY
C
      CALL EIRENE_PLASMA
      CALL EIRENE_PLASMA_DERIV(0)
C
C  SAVE PLASMA DATA AND ATOMIC DATA ON FORT.13
C
      NFILEL=3
      IFLG=0
      CALL EIRENE_WRPLAM(TRCFLE,IFLG)
C
      DEALLOCATE (PDEN)
      DEALLOCATE (EDEN)
      DEALLOCATE (PDEN2)
      DEALLOCATE (EDEN2)
      DEALLOCATE (ENERGY)
      DEALLOCATE (PLS)
      IF (ALLOCATED(CROSSTEMP)) DEALLOCATE(CROSSTEMP)
      DEALLOCATE (GBGKV)

C
      RETURN
C
  995 CONTINUE
      WRITE (iunout,*) 'SPECIES ERROR IN MODBGK'
      CALL EIRENE_EXIT_OWN(1)
      END
