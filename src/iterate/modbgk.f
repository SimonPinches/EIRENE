cdr Nov  19:   From testing modbgk in internal eirene iterations
cdr            i.e. not via B2, or other plasma codes:
cdr
cdr            Bug fixes: alloc of plasma_bckground was too late.
cdr            In ca. 2008, when merging SOLPS4.3 and the Master EIRENE
cdr            version at FZ Juelich further 
cdr            bugs regarding cross-collisions got implemented:
cdr            Cross-collisions led to erroneous Ti and Ni profiles,
cdr            due to wrong setting of TIINTF, DIINTF pointers.
cdr            Runs were perhaps still correct, (by luck?) 
cdr            but convergence monitoring was entirely wrong since. 
cdr            New further convergence diagnostics added, and code cleanup,
cdr            documentation.
cdr Sept 19:   problems with NIDC, reuse of plasma_bckgrnd structure
cdr            only Diin, Tiin, V_xyzIN of bgk virt. species should be
cdr            modified.
cdr            Storage on plasma_bckgrnd ?
cdr            In EMC3 case MODBGK is not used at all, instead the
cdr            affected virt. plasma data
cdr            are directly modified in infcop_emc3.


c   modbgk  FZJ-MASTER, eirene_git, sept. 2014

C       ??     INDIRECT SPECIES INDEXING IPLSV, IPLSTI INCLUDED, FOR TIIN and VXIN,VYIN,VZIN
CDR            This is risky, because bgk collisions
cdr            may then overwrite temperatures/velocities for non-bgk background
CDR            depending on setting of MPLSTI(ipls) and MPLSV(ipls) arrays.
cdr            In case of bgk: always use fully multiple temperatures and multiple flow velocities.

cdr  implemented corresponding check.

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
c               WITH THIS EFFECTIVE TEMPERATURE.
c               This appears to be another option for the 5th (free)
c               bgk model parameter, to match certain transport coefficients.
c               (the other 4 are fixed by other criteria)
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
      USE EIRMOD_COMUSR  !TEXTS, NSPAMI
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_CSPEZ ! logatm, logmol, logion

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
     .            T2, ED1, ED2, VXMIX, VZMIX, DELX, DELY, DELZ,
     .            VX, VY, VZ, EOLD, TMIX,
     .            ED, RM, FACTKK, EBULK, TS1, DS1,
     .            FACT1, FACT2, FACTNEW, RMAS1, RMAS2, A_ROBIN,
     .            RESE, RESM, TBEL, DOLD, DEIMIN,
     .            DEL, TII, CNDYN, RRN, RATE, RESN, RATN, RRE, RRM,
     .            RM1, RM2, TCSUM,
     .            RESDEN, RESENE, RESMOM,
     .            RATDEN, RATENE, RATMOM
      INTEGER ::  ITYP1(NPLS), ITYP2(NPLS), ISPZ1(NPLS), ISPZ2(NPLS),
     .            IREL1(NPLS), INRC1(NPLS),CROSSINDEX(NPLS),NCROSS,
     .            ICROSS1,ICROSS2,ICROSS
      INTEGER ::  J, JATM, JMOL, JION, JPLS, IESTM, IBGV,
     .            ISCDE, IPLS1, NRWK1, I, ISP, IPLS2, IATM2, IIEL,
     .            IAEL, IMEL, IUP12, IUP22, IION2, IBGK2, IMOL2, IUP2,
     .            IUP3, IUP1, IBGK1, IP, NRC, IR, IT, IREL,
     .            KK, NXM, NYM, NZM, IUP32, IPLSTI, IPLSTI1, IPLSTI2,
     .            IPLSV, IRD, I_FINE, IRAD, IFLG
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      LOGICAL :: LMARK(NPLS)
      LOGICAL :: TRCSAV, LSKIP
C
C
      CALL EIRENE_LEER(3)
      WRITE (iunout,*) 'MODBGK CALLED AFTER ITERATION IITER= ',IITER
      CALL EIRENE_LEER(3)

      IF (NFILEL.EQ.0) THEN
        WRITE (0,*) 'BGK COLLISIONS ARE REQUESTED, BUT'
        WRITE (0,*) 'INFORMATION FROM '//FORT//'13 FILE IS NOT PROVIDED'
        WRITE (0,*) 'THE BGK REACTIONS WILL HAVE NO EFFECT!'
        WRITE (0,*) 'CHECK SETTINGS OF NFILE IN INPUT FILE'
      END IF
CC
cdr  FOR STOCH. APPROX. UNDER-RELAXATION.  NOT IN USE
      A_ROBIN=1.D0/REAL(IITER,DP)

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
c       call EIRENE_exit_own(1)
      endif
      if (any(MPLSTI(1:npls) /= (/(i,i=1,npls)/))) then
        write (iunout,*) 'MPLSTI-option not ready in MODBGK '
        write (iunout,*) 'NPLS, NPLSTI ',NPLS, NPLSTI
        call EIRENE_exit_own(1)
      endif
      if (any(MPLSV(1:npls) /=  (/(i,i=1,npls)/))) then
        write (iunout,*) 'MPLSV-option not ready in MODBGK '
        write (iunout,*) 'NPLS, NPLSV ',NPLS, NPLSV
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
      ENERGY = 0._DP

cdr  PLS:  ELECTRON DENSITY PARAMETER in CR MODELS
cdr       (NOT TO BE CONFUSED WITH THE DENSITY FACTOR BETWEEN RATES AND RATE COEFF.)
cdr: set hard-wired lower density for H.4, H.10 type fits from AMJUEL: 1e8 cm**-3.
cdr: At this lower limit density the fits are produced such
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

cdr  Set bgk volume-averaged tallies BGKV, which MAY HAVE been scored on coarser grid only,
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

      NXM=MAX(1,NR1STM)
      NYM=MAX(1,NP2NDM)
      NZM=MAX(1,NT3RDM)
c
C  LOOP OVER THOSE BACKGROUND ION SPECIES, WHICH ARE ARTIFICIAL
C  SPECIES FOR (NONLINEAR) ITERATIONS
C .........................................................................
      DO 1000 JPLS=1,NPLSI
        IPLS=JPLS
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
        DO JATM=1,NATMI
          IF (NPBGKA(JATM).EQ.IBGK1) THEN
            ITYP1(IPLS)=1
            ISPZ1(IPLS)=JATM
            FACT1=CVRSSA(JATM)
            RMAS1=RMASSA(JATM)
            DO I_fine=1,NRAD
!pb 05.02.2013  take care of coarse scoring cells
c    fine cell i_fine belongs to coarser cell ird. Scoring was done on coarse cell ird only
cdr  This is not obvious for additional cells for averaging, e.g. IR=NR1ST, etc..
cdr  Those IRD should not appear here  (hopefully)
              IRD = NCLTAL(I_fine)
              IF (IRD == 0) CYCLE
              PDEN(I_fine)=PDENA(JATM,IRD)
              EDEN(I_fine)=EDENA(JATM,IRD)
            ENDDO
C  FIND INDEX NRC
            DO NRC=1,NRCA(JATM)
              IP=EIRENE_IDEZ(IBULKA(JATM,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 10
              ENDIF
            ENDDO
            GOTO 995
C  FIND INDEX IREL
   10       DO IAEL=1,NAELI(JATM)
              IF (LGAEL(JATM,IAEL,1).EQ.IPLS) THEN
cdr elastic collision between jatm and ipls
                IREL1(IPLS)=LGAEL(JATM,IAEL,0)
                GOTO 50
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
C  TRY MOLECULES
        DO JMOL=1,NMOLI
          IF (NPBGKM(JMOL).EQ.IBGK1) THEN
            ITYP1(IPLS)=2
            ISPZ1(IPLS)=JMOL
            FACT1=CVRSSM(JMOL)
            RMAS1=RMASSM(JMOL)
            DO I_fine=1,NRAD
!pb 05.02.2013  take care of coarse scoring cells
c     fine cell i_fine belongs to coarse cell ird. scoring was done on coarse cell ird
cdr  This is not obvious for additional cells for averaging, e.g. IR=NR1ST, etc..
cdr  Those IRD should not appear here  (hopefully)
              IRD = NCLTAL(I_fine)
              IF (IRD == 0) CYCLE
              PDEN(I_fine)=PDENM(JMOL,IRD)
              EDEN(I_fine)=EDENM(JMOL,IRD)
            ENDDO
C  FIND INDEX NRC
            DO NRC=1,NRCM(JMOL)
              IP=EIRENE_IDEZ(IBULKM(JMOL,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 20
              ENDIF
            ENDDO
            GOTO 995
C  FIND INDEX IREL
   20       DO IMEL=1,NMELI(JMOL)
              IF  (LGMEL(JMOL,IMEL,1).EQ.IPLS) THEN
                IREL1(IPLS)=LGMEL(JMOL,IMEL,0)
                GOTO 50
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
C  TRY TEST IONS
        DO JION=1,NIONI
          IF (NPBGKI(JION).EQ.IBGK1) THEN
            ITYP1(IPLS)=3
            ISPZ1(IPLS)=JION
            FACT1=CVRSSI(JION)
            RMAS1=RMASSI(JION)
            DO I_fine=1,NRAD
!pb 05.02.2013  take care of coarse scoring cells
c    fine cell i_fine belongs to coarse cell ird. scoring was done on coarse cell ird
cdr  This is not obvious for additional cells for averaging, e.g. IR=NR1ST, etc..
cdr  Those IRD should not appear here  (hopefully)
              IRD = NCLTAL(I_fine)
              IF (IRD == 0) CYCLE
              PDEN(I_fine)=PDENI(JION,IRD)
              EDEN(I_fine)=EDENI(JION,IRD)
            ENDDO
C  FIND INDEX NRC
            DO NRC=1,NRCI(IION)
              IP=EIRENE_IDEZ(IBULKI(JION,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 30
              ENDIF
            ENDDO
C  FIND INDEX IREL
   30       DO IIEL=1,NIELI(JION)
              IF  (LGIEL(JION,IIEL,1).EQ.IPLS) THEN
                IREL1(IPLS)=LGIEL(JION,IIEL,0)
                GOTO 50
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
        GOTO 995

C  AT THIS POINT: COLLISION PARTNER AMONGST TEST PARTICLES
C  interacting with virtual background IPLS
C  HAS BEEN IDENTIFIED: (ITYP1, JATM, JMOL, JION)
C  AS WELL AS THE LABEL NUMBER OF COLLISION PROCESS IREL1.
C  DENSITY AND ENERGY DENSITY TALLIES OF COLLISION PARTNER "PDEN,EDEN"
C  HAVE NOW BEEN SET ON FINE GRID
C  THE RELEVANT FURTHER BGK TALLIES ARE: IUP1, IUP2, IUP3.
C
   50   CONTINUE
C
        IREL=IREL1(IPLS)
C
C  Skip iteration processes for IPLS,
C  when corresponding test particle did not exist in this run
        LSKIP=.FALSE.
        SELECT CASE (ITYP1(IPLS))
          CASE (1)
            IF (.NOT.LOGATM(JATM,ISTRA)) LSKIP=.TRUE.
          CASE (2)
            IF (.NOT.LOGMOL(JMOL,ISTRA)) LSKIP=.TRUE.
          CASE (3)
            IF (.NOT.LOGION(JION,ISTRA)) LSKIP=.TRUE.
        END SELECT
        IF (LSKIP) goto 2000
C
C
C  SELF-COLLISION OR CROSS-COLLISION
C
        IF (NPBGKP(IPLS,2).EQ.0) THEN

C  SELF-COLLISION: ALL DONE
          ITYP2(IPLS)=-1
          ISPZ2(IPLS)=0
C
        ELSEIF (NPBGKP(IPLS,2).NE.0) THEN
C
C  CROSS-COLLISION, FIND SECOND COLLISION PARTNER

C  THIS IS NOT THE INGOING COLLIDING TEST PARTICLE, WHICH WE HAVE ALREADY IDENTIFIED,
C  (AND WHICH, E.G., DETERMINES MASS AND DENSITY OF ARTIFICIAL BACKGROUND PARTICLE IPLS)
C  BUT, INSTEAD, IT IS THE TEST PARTICLE WHICH PLAYS THE ROLE
C  OF THE "SECOND" PARTICLE, AMONGST THE TEST PARTICLES
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

cdr  relevant further BGK tallies for this cross-collision
cdr  parameters of second involved species.
          IUP12=(IBGK2-1)*3+1
          IUP22=(IBGK2-1)*3+2
          IUP32=(IBGK2-1)*3+3
C
        ENDIF  ! Cross-collision, two different test species
c
cdr   Parameters of virtual background species are set,
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
        RESDEN=0.
        RESENE=0.
        RESMOM=0.
        RATDEN=0.
        RATENE=0.
        RATMOM=0.
C
        IF (ITYP2(IPLS).EQ.-1) THEN
C
        IF (TRCMOD) THEN
          WRITE (iunout,*) 'MODBGK: SELF-COLLISION WITH, IPLS',IPLS
          WRITE (iunout,*) 'ITYP,ISPZ,IBGK_SP,IREL ',ITYP1(IPLS),
     .                                ISPZ1(IPLS),IBGK1,IREL1(IPLS)
        ENDIF
C
C
C  this cell loop is referring to the underlying 'fine' grid, not the coarse grid
C                              on which the eirene tallies had been updated.

        DO IR=1,NXM
          DO IP=1,NYM
            DO IT=1,NZM
              IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA  !irad=i_fine
              CALL EIRENE_SELF_COLLISION
            END DO
          END DO
        END DO
c
c  same as do loop above, for additional cell region
c
        DO IRAD=NSURF+1,NSURF+NRADD
          CALL EIRENE_SELF_COLLISION
        END DO  ! IRAD
cdr
c  check proper mass factor:
        FACTNEW=CVRSSP(IPLS)
        IF (FACT1.NE.FACTNEW) THEN
          write (IUNOUT,*) 'Incorrect mass in MODBGK, IPLS: ',IPLS
          call eirene_exit_own(1)
        END IF
C
C   IF IPLS IS AN ARTIFICIAL BACKGROUND SPECIES FOR A SELF-COLLISION  :  DONE !
C
C   NEXT, IF IPLS IS AN ARTIFICIAL BACKGROUND SPECIES FOR A CROSS-COLLISION BETWEEN TWO DISTINCT
C        SPECIES OF TEST PARTICLES "1" AND "2".
C
      ELSE
C
          IF (TRCMOD) THEN
            WRITE (iunout,*) 'MODBGK: CROSS-COLLISION, IPLS ',IPLS
            WRITE (iunout,*) 'ITYP1,ISPZ1,IBGK1,IREL1 ',
     .                        ITYP1(IPLS),ISPZ1(IPLS),IBGK1,IREL1(IPLS)
            WRITE (iunout,*) 'ITYP2,ISPZ2,IBGK2       ',
     .                        ITYP2(IPLS),ISPZ2(IPLS),IBGK2
          ENDIF

Cdr MASS FACTOR FOR "CROSS-COLLISION TEMPERATURE",
cdr i.e. for the effective temperature entering IJ and JI cross-collision terms.
cdr  1-2 symmetric.
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
              CALL EIRENE_CROSS_COLLISION
            END DO
          END DO
        END DO
c
c  same as do loop above, for additional cell region
c
        DO IRAD=NSURF+1,NSURF+NRADD
          CALL EIRENE_CROSS_COLLISION
        END DO  ! IRAD
c
        ENDIF   !  SELF- OR CROSS-COLLISION
c

 2000   CONTINUE
        CALL EIRENE_LEER(2)
        WRITE (iunout,*) 'Virt. Species IPLS ',TEXTS(NSPAMI+IPLS)
        WRITE (iunout,*) 'Reaction IREL: ',IREL
        IF (NPBGKP(IPLS,2).EQ.0) then
          WRITE (IUNOUT,*) 'SELF-COLLISION'
        ELSE
          WRITE (IUNOUT,*) 'CROSS-COLLISION'
        ENDIF

        IF (LSKIP) THEN
          WRITE (iunout,*) 'BGK iteration skipped for this IPLS'
          WRITE (iunout,*) 'No related test particle in this run'
          CALL EIRENE_LEER(2)
          CYCLE
        ENDIF

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
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'TOTAL PARAMETER CHANGES:'
        WRITE (IUNOUT,*) 'TOTAL NO. OF IPLS PARTICLES [1]'
        CALL EIRENE_MASR2('RESDEN: TOT, REL',RESDEN, RESDEN/(RRN+EPS30))
        WRITE (IUNOUT,*) 'TOTAL ENERGY CONTENT OF IPLS PARTICLES [eV]'
        CALL EIRENE_MASR2('RESENE: TOT, REL',RESENE, RESENE/(RRE+EPS30))
        WRITE (IUNOUT,*) 'TOTAL MOMENTUM OF IPLS PARTICLES [g cm/s]'
        CALL EIRENE_MASR1('RESMOM= ',RESMOM)

        CALL EIRENE_LEER(2)
C
C.................................................................
 1000 CONTINUE  ! IPLS LOOP
C.................................................................

      CALL EIRENE_LEER(2)
C
C  SAVE OVERHEAD, IF GEOMETRY DATA ALREADY AVAILABLE ON FILE
C
      IF (NFILEM.EQ.1) NFILEM=2
C
cdr
C    SET INDPRO(2:6)=7, AND
C    WRITE PLASMA DATA ONTO PLASMA_BCKGRND FOR CALL TO SUBR. PLASMA BELOW
C    TI (IPLS=1,NPLSTI) ,
C    NI(IPLS=1,NPLS) AND
C   (VX,VY,VZ) (IPLS=1,NPLSV)
cdr try to leave everything else untouched
C   PLAY SAFE: WRITE ENTIRE PLASMA_BCKGRND ARRAY.
C

      DO 500 I=1,6   ! should be: 2:6, leave TEIN untouched?
        INDPRO(I)=7
  500 CONTINUE
! STORAGE FOR INPUT TALLIES 1 (TEIN) TO 13 (ADIN), WITHOUT NO.3 (DEIN)
      NRWK1=6+NPLS+NPLSTI+3*NPLSV+NAIN
cdr Oct.2019
cdr   NIINTF is set in eirene_alloc_bckgrnd. Careful: "Hidden link"?
cdr   If alloc_backgrnd has not been called, then also NIINTF is not known.
      IF (.NOT. ALLOCATED(PLASMA_BCKGRND)) CALL EIRENE_ALLOC_BCKGRND



      IF (NIINTF < NRWK1) THEN
        WRITE (iunout,*) ' PLASMA_BCKGRND ARRAY IS TOO SMALL TO HOLD'
        WRITE (iunout,*) ' PLASMA DATA in MODBGK'
        WRITE (iunout,*) ' CHECK PARAMETER NSMSTRA'
        CALL EIRENE_EXIT_OWN(1)
      END IF
      PLASMA_BCKGRND(1:NRWK1,:) = 0.D0

cdr
cdr There should be no need to fiddle around with B field here at all.
cdr initialize BXIN=0
      PLASMA_BCKGRND(1+1*NPLS+NPLSTI+3*NPLSV+1,:)= 0._DP  
cdr initialize BYIN=0
      PLASMA_BCKGRND(2+1*NPLS+NPLSTI+3*NPLSV+1,:)= 0._DP  
!pb initialize BZIN=1
      PLASMA_BCKGRND(3+1*NPLS+NPLSTI+3*NPLSV+1,:)= 1._DP
!pb initialize BFIN=1
      PLASMA_BCKGRND(4+1*NPLS+NPLSTI+3*NPLSV+1,:)= 1._DP


cdr  set new (virtual) background data species data onto plasma_bckgrnd,
cdr  for next iteration
      DO IR=1,NXM
        DO IP=1,NYM
          DO IT=1,NZM
            IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
            CALL EIRENE_SET_VIRT_SPECIES
          END DO
        END DO
      END DO
C
c  same as do loop above, for additional cell region
c
      DO IRAD=NSURF+1,NSURF+NRADD
        CALL EIRENE_SET_VIRT_SPECIES
      ENDDO
C  the target plasma_bckgrnd, as well as the fields DIIN, TIIN, VXIN, VYIN, VZIN
C  are now set consistently for the virtual species in the next iteration.
C
      CALL EIRENE_PLASMA_DERIV(0)
C
C .........................................................................
C  NOW: NEW COLLISION RATES MUST BE SET FOR THE NEXT ITERATION
C .........................................................................
C
C  IN CASE OF CROSS-COLLISION, SOME MODIFICATIONS OF THE
C  BACKGROUND PARAMETERS ARE REQUIRED TEMPORARILY, TO ENFORCE
C  A SPECIFIC RELATION BETWEEN TAU_1,2 AND TAU_2,1, LGVAC, etc...
C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE PROFILES
C
C
      TRCSAV=TRCAMD
      TRCAMD=.FALSE.
C
c     write (iunout,*) 'before diin '
c     do ipls=1,npls
c       write (iunout,*) ipls,diin(ipls,1)
c     enddo

      CALL EIRENE_LEER(1)
      DO JPLS=1,NPLSI
        LMARK(JPLS)=.FALSE.
      ENDDO
      DO IPLS1=1,NPLSI
        IF (NPBGKP(IPLS1,2).NE.0) THEN
C  IPLS1 IS A CROSS-COLLISION SPECIES
C  FIND CORRESPONDING 2ND CROSS-COLLISION SPECIES
          IPLS2=0
          DO JPLS=1,NPLSI
            IF (ITYP1(JPLS).EQ.ITYP2(IPLS1).AND.
     .          ISPZ1(JPLS).EQ.ISPZ2(IPLS1).AND.
     .          ITYP2(JPLS).EQ.ITYP1(IPLS1).AND.
     .          ISPZ2(JPLS).EQ.ISPZ1(IPLS1)) IPLS2=JPLS
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
cdr crosstemp are the temperature for evaluating the new cross-collision rates.
cdr These must be symmetric to make also the collision rates symmetic.
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
cdr  exchange the two densities, to get that TABEL rates right.
cdr  For the f_12 term the density is n1, but the rate must be scaled with n2.
          DO IRAD=1,NSBOX
            DS1=DIIN(IPLS1,IRAD)
            DIIN(IPLS1,IRAD)=DIIN(IPLS2,IRAD)
            DIIN(IPLS2,IRAD)=DS1
CVK PLASMA BACKGROUND CORRECTION
cdr This must be wrong: DIINTF is a pointer to the plasma_bckgrnd target.
cdr
cdr         DIINTF(IPLS1,IRAD)=DIIN(IPLS1,IRAD)
cdr         DIINTF(IPLS2,IRAD)=DIIN(IPLS2,IRAD)

CSW 02jan2012 CHECK FOLLOWING TWO LINES !!!! ipls1/ipls2 mixed up?
cdr Nov. 2019: No, not mixed up, but entirely wrong.
cdr tbd: check the proper handling of TI in cross-collisions,
cdr depending on the model chosen.
cdr Compare with original paper: Reiter et al., J Nucl. Mat. 241-243 (1997) 342
cdr and the modifications (VK, Juel-report,2007) from the
cdr originally implemented Morse/Hamel form now to the Holway form.
cdr         TIINTF(IPLSTI1,IRAD)=TIIN(IPLSTI2,IRAD) ! HAS BEEN CHANGED IPLS1->IPLS2
cdr         TIINTF(IPLSTI2,IRAD)=TIIN(IPLSTI1,IRAD) !
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
C  DIIN,TIIN  (such as lgvac, tiinl, diinl,...)
C
      CALL EIRENE_PLASMA_DERIV(0)
C
C  RESET BGK ATOMIC AND MOLECULAR DATA ARRAYS TO NEW BACKGROUND
cdr  careful: the pointers diintf, tiintf,...etc. are used. 
C
      DO JPLS=1,NPLSI
        IF (NPBGKP(JPLS,1).NE.0) THEN
          ITYP=ITYP1(JPLS)
          ISPZ=ISPZ1(JPLS)
          IREL=IREL1(JPLS)
          NRC=INRC1(JPLS)
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
          CALL EIRENE_XSTEL(IREL,ISP,JPLS,EBULK,ISCDE,IESTM,
     .                      KK,FACTKK,PLS)
          IF (NPBGKP(JPLS,2).NE.0) THEN
C  CROSS-COLLISION, RESET EPLEL3 FOR TRACKLENGTH ESTIMATOR
            IF (NSTORDR >= NRAD) THEN
              DO J=1,NSBOX
                EPLEL3(IREL,J,1)=ENERGY(JPLS,J)
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
C  TABEL3, EPLEL3 are set now with intermediate ni, Ti, Vi.

C  RESTORE PLASMA DATA FROM PLASMA_BCKGRND ARRAY
cdr
cdr error ?  diint tiint are pointers to target plasma_bckground.
cdr          so they are still changed !
cdr error !
C
      CALL EIRENE_PLASMA
      CALL EIRENE_PLASMA_DERIV(0)
C
C  SAVE PLASMA DATA AND ATOMIC DATA ON FORT.13
C
      NFILEL=3
      IFLG=0
      CALL EIRENE_WRPLAM(TRCFLE,IFLG)

c     write (iunout,*) 'after diin '
c     do ipls=1,npls
c       write (iunout,*) ipls,diin(ipls,1)
c     enddo

C
cdr  remove temporary storage
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

      CONTAINS

      SUBROUTINE EIRENE_CROSS_COLLISION
      IMPLICIT NONE
cdr  set parameters for cross-collision Maxwellian f_ij in the f_i BGK collision term.
cdr  i=ipls corresponds to the incident test particle, colliding with particle "j".

CDR Also return: CROSSTEMP: the temperature to be used in post collision rates TABEL, ...
cdr Note: this rules out "on the fly" A&M data evaluation, because then crosstemp is not known.

cdr Also returns: cumulated residuals

      REAL(DP) :: EIRENE_RATE_COEFF
      REAL(DP) :: ENEW, EDMIX

      TBEL=0
      IF (LGVAC(IRAD,IPLS)) GOTO 191
cdr elastic collision rate tabel3, at old ni, Ti, Vi
        IF (NSTORDR >= NRAD) THEN
          TBEL = TABEL3(IREL,IRAD,1)
        ELSE
cdr       TBEL=EIRENE_FTABEL3(IREL,IRAD)  ! this should replace the next three cards
          KK=NREAEL(IREL)
          TII=TIINL(IPLSTI,IRAD)+ADDEL(IREL,IPLS)
          TBEL = EIRENE_RATE_COEFF(KK,IRAD,TII,0._DP,.TRUE.,0)
     .           *DIIN(IPLS,IRAD)*FACREL(IREL,1)
        END IF
  191 CONTINUE

      DEL=DIIN(IPLS,IRAD)-PDEN(IRAD)

cdr  all residuals are evaluated
      RATN=RATN+TBEL*DEL*VOL(IRAD)
      RESN=RESN+TBEL*ABS(DEL)*VOL(IRAD)
c  density residuals; UNITS: [1] (even if TBEL=0)
      RATDEN=RATDEN+DEL*VOL(IRAD)
      RESDEN=RESDEN+ABS(DEL)*VOL(IRAD)

c  prepare for momentum residuals, and mixed V
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

c  prepare for energy residuals, and mixed Ti
      EOLD=1.5*TIIN(IPLSTI,IRAD)*DIIN(IPLS,IRAD)
      IF (LEDRIFT) EOLD=EOLD+EDRIFT(IPLS,IRAD)*DIIN(IPLS,IRAD)

cdr  in case of cross-collisions: use new mixed energy density to compare with

      ED1=(VXIN1**2+VYIN1**2+VZIN1**2)*FACT1
      ED2=(VXIN2**2+VYIN2**2+VZIN2**2)*FACT2
      T1=(EDEN(IRAD)/(PDEN(IRAD)+EPS60)-ED1)/1.5
      T2=(EDEN2(IRAD)/(PDEN2(IRAD)+EPS60)-ED2)/1.5

! THE TEMPERATURE FOR EVALUATING THE NEW CROSS-COLLISION RATE for i+j cross-collisions.
cdr  clearly: must be the same as that for j+i cross-collisions.
      CROSSTEMP(ICROSS,IRAD)=T1*RM2+T2*RM1

cdr temperature in f_ij, not necessarily the same as temperature in f_ji
      TMIX=T1+RM*(T2-T1+
     .         FACT2/3.D0*((VXIN1-VXIN2)**2+(VYIN1-VYIN2)**2+
     .                     (VZIN1-VZIN2)**2))

      FACTNEW=CVRSSP(IPLS)
      EDMIX=(VXMIX**2+VYMIX**2+VZMIX**2)*FACTNEW
      ENEW=(1.5*TMIX + EDMIX)*PDEN(IRAD)

      DEL=EOLD-ENEW
      RATE=RATE+TBEL*DEL*VOL(IRAD)
      RESE=RESE+TBEL*ABS(DEL)*VOL(IRAD)
c  energy residuals; UNITS: [eV] (even if TBEL=0)
      RATENE=RATENE+DEL*VOL(IRAD)
      RESENE=RESENE+ABS(DEL)*VOL(IRAD)

c  Next: set new virtual background parameters: V.IN, TIIN, DIIN

C NEW T
      TIIN(IPLSTI,IRAD)=TMIX

C NEW V
      VXIN(IPLSV,IRAD)=VXMIX
      VYIN(IPLSV,IRAD)=VYMIX
      VZIN(IPLSV,IRAD)=VZMIX

C NEW N
      DIIN(IPLS,IRAD)=PDEN(IRAD)

C  TOTALS: PARTICLES, ENERGY; MOMENTUM
      RRN=RRN+PDEN(IRAD)*VOL(IRAD)  ! New particle content
C     RRM=?
      RRE=RRE+ENEW*VOL(IRAD)  ! New energy content

      RETURN
      END SUBROUTINE EIRENE_CROSS_COLLISION


      SUBROUTINE EIRENE_SELF_COLLISION
      IMPLICIT NONE
      REAL(DP) :: EIRENE_RATE_COEFF
      REAL(DP) :: VNEW, VOLD

      TBEL=0.
      IF (LGVAC(IRAD,IPLS)) GOTO 81
      IF (NSTORDR >= NRAD) THEN
        TBEL = TABEL3(IREL,IRAD,1)
      ELSE
cdr     TBEL=EIRENE_FTABEL3(IREL,IRAD)  ! this should replace the next three cards
        KK=NREAEL(IREL)
        TII=TIINL(IPLSTI,IRAD)+ADDEL(IREL,IPLS)
        TBEL = EIRENE_RATE_COEFF(KK,IRAD,TII,0._DP,.TRUE.,0)
     .         *DIIN(IPLS,IRAD)*FACREL(IREL,1)
      END IF
   81 CONTINUE
C DELTA_N
      DOLD=DIIN(IPLS,IRAD)
      DEL=DOLD-PDEN(IRAD)

c  rate of particle exchange: [1/s], per cell I_fine, for reaction IREL
      RATN=RATN+TBEL*DEL*VOL(IRAD)
c  L1 norm of particle residual; INTEGRATED.
      RESN=RESN+TBEL*ABS(DEL)*VOL(IRAD)
c  density residuals; UNITS: [1] (even if TBEL=0)
      RATDEN=RATDEN+DEL*VOL(IRAD)
      RESDEN=RESDEN+ABS(DEL)*VOL(IRAD)

c     IF (TRCMOD) THEN
c       CALL EIRENE_MASJR3('IR,T,TBEL,RATN                  ',
c    .   IRAD,TIIN(IPLSTI,IRAD),TBEL,TBEL*DEL*VOL(IRAD)*ELCHA)
c       CALL EIRENE_MASR3('DOLD,DNEW,DEL           ',
c    .                     DOLD,PDEN(IRAD),DEL)
c     ENDIF

C DELTA_E
      EOLD=1.5*TIIN(IPLSTI,IRAD)*DIIN(IPLS,IRAD)
      IF (LEDRIFT) EOLD=EOLD+EDRIFT(IPLS,IRAD)*DIIN(IPLS,IRAD)
      DEL=EOLD-EDEN(IRAD)

c  rate of energy exchange: (eV)/s, per cell I_fine, for reaction IREL
      RATE=RATE+TBEL*DEL*VOL(IRAD)
c  L1 norm of energy residual
      RESE=RESE+TBEL*ABS(DEL)*VOL(IRAD)
c  energy residuals; UNITS: [eV] (even if TBEL=0)
      RATENE=RATENE+DEL*VOL(IRAD)
      RESENE=RESENE+ABS(DEL)*VOL(IRAD)

C DELTA_V
CDR  next lines: gbgkv (=bgkv) is redundant, because momentum density vector
cdr              components have become default volume-averaged output tallies
cdr
      DELX=GBGKV(IUP1,IRAD)-VXIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
      DELY=GBGKV(IUP2,IRAD)-VYIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)
      DELZ=GBGKV(IUP3,IRAD)-VZIN(IPLSV,IRAD)*DIIN(IPLS,IRAD)

      VOLD=0.
      VNEW=0.
      DEL=VOLD-VNEW  ! to be written

c  cdr  rate of momentum exchange: (cm/s)/s, per cell I_fine, for reaction IREL
      RATM(1)=RATM(1)+TBEL*DELX*VOL(IRAD)
      RATM(2)=RATM(2)+TBEL*DELY*VOL(IRAD)
      RATM(3)=RATM(3)+TBEL*DELZ*VOL(IRAD)

c  momentum residuals; UNITS: [cm/s] (even if TBEL=0)
      RATMOM=RATMOM+DEL*VOL(IRAD)
      RESMOM=RESMOM+ABS(DEL)*VOL(IRAD)


c  Next: set new virtual background parameters: V.IN, TIIN, DIIN
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

C  For normalization of global residuals:
c  RRN: total particle [1], new
c  RRM: total momentum content, new
c  RRE: total energy content (eV), new
      RRN=RRN+PDEN(IRAD)*VOL(IRAD)
C     RRM=?
      RRE=RRE+EDEN(IRAD)*VOL(IRAD)

      RETURN
      END SUBROUTINE EIRENE_SELF_COLLISION

      SUBROUTINE EIRENE_SET_VIRT_SPECIES
      IMPLICIT NONE
      INTEGER IPLSTI, IPLSV, IAIN, JPLS
cdr  set virt. species IPLS parameters NI(IPLS), TI(IPLS), VX(IPLS),VY(IPLS),VZ(IPLS)
cdr  in calling program: 2 identical loops, IRAD

      PLASMA_BCKGRND  (0+0*NPLS+1   ,IRAD)= TEIN(IRAD)
      DO IPLSTI=1,NPLSTI
        PLASMA_BCKGRND(1+0*NPLS+IPLSTI,IRAD)= TIIN(IPLSTI,IRAD)
      END DO
      DO JPLS=1,NPLS
        PLASMA_BCKGRND(1+0*NPLS+NPLSTI+JPLS,IRAD)= DIIN(JPLS,IRAD)
      ENDDO
      DO IPLSV=1,NPLSV
        PLASMA_BCKGRND(1+1*NPLS+NPLSTI+0*NPLSV+IPLSV,IRAD)=
     .             VXIN(IPLSV,IRAD)
        PLASMA_BCKGRND(1+1*NPLS+NPLSTI+1*NPLSV+IPLSV,IRAD)=
     .             VYIN(IPLSV,IRAD)
        PLASMA_BCKGRND(1+1*NPLS+NPLSTI+2*NPLSV+IPLSV,IRAD)=
     .             VZIN(IPLSV,IRAD)
      END DO

      IF (LBXIN)
     .  PLASMA_BCKGRND(1+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BXIN(IRAD)
      IF (LBYIN)
     .  PLASMA_BCKGRND(2+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BYIN(IRAD)
      IF (LBZIN)
     .  PLASMA_BCKGRND(3+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BZIN(IRAD)
      IF (LBFIN)
     .  PLASMA_BCKGRND(4+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= BFIN(IRAD)

      IF (LVOL)
     .  PLASMA_BCKGRND(5+1*NPLS+NPLSTI+3*NPLSV+1,IRAD)= VOL(IRAD)

      IF (LADIN) THEN
        DO IAIN=1,NAINI
          PLASMA_BCKGRND(6+1*NPLS+NPLSTI+3*NPLSV+IAIN,IRAD)=
     .             ADIN(IAIN,IRAD)
        END DO
      END IF

      RETURN
      END SUBROUTINE EIRENE_SET_VIRT_SPECIES



      END SUBROUTINE EIRENE_MODBGK
