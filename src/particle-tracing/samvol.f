cdr Jan   18 : only notational change, to distuingish surface substrata from volume substrata
cdr  5.14.15 : vecusr called with ncell, and 0,0,0 (center of gravity)
cdr  2.11.14 : new function eirene_brems: bremsstrahlung in W per ion
cdr            replaces explicit expression.
cdr 21.10.14 : bug fix: spectral cut off flag ICCT set to zero for default vol.rec (KK=0)
cdr           -->now runs again on eirene default vol.rec model.
cdr 30.10.14 :  lplssr true even if npts=0, to allow setting up volume source tallies,
cdr             even if npts=0 for the vol-rec stratum

cdr  1111.07: "istep out of range" error message removed once again.
!pb  2203.07: LEVGEO=6 --> LEVGEO=10
!pb  2710.06: use flux set by user defined sampling routine
!pb  1001.06: ENTRY SAMVOL_REINIT added for reinitialsation of Eirene
!pb  1812.06: calculate bremsstrahlung
!pb  2408.06: set output values for DIWL and SHWL
cdr  2008.06: tiwl(*), ... instead of tiwl(npls),... to unify code.
cdr  0604.06: check "istep out of range" moved to correct place
c    0311.05: iplsti moved after check of validity of ipls, to produce legal exit
c             rather than code crash
C  JET 2005, PATCH 1: NEW ARGUMENTS EFWL AND SHWL IN PARAMETER LIST
c                     AT ENTRY SMVOL1 AND SMUSR1
C
      SUBROUTINE EIRENE_SAMVOL
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT
      USE EIRMOD_COMSOU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_CTRIG
      USE EIRMOD_CTETRA
      USE EIRMOD_PHOTON
      IMPLICIT NONE
C
 
      REAL(DP), INTENT(OUT) :: TEWL, SHWL, TIWL(*), DIWL(*),
     .                         VXWL(*), VYWL(*), VZWL(*),
     .                         EFWL(*), WEISPZ(*)
      INTEGER, INTENT(IN) :: NVLM
      REAL(DP), ALLOCATABLE, SAVE :: FREC(:,:,:), VSOURC(:,:), VSMXI(:)
      REAL(DP), ALLOCATABLE, SAVE :: RQ21(:), PS21(:)
      REAL(DP), ALLOCATABLE, SAVE :: ASIMP(:,:)
      INTEGER, ALLOCATABLE, SAVE  :: ISOURC(:,:), ICMX(:),
     .                               IFREC(:)
      REAL(DP) :: ZEP1, X1, Y1, X2, Y2, X3, Y3, RR, RRI, WINK,
     .            ZRM1, CNORM, EPR, ELR, RRD, RRN, ZZ, X01, Z1, Z2, Z3,
     .            REC, BX, BY, BZ, ADD, EIRENE_FTABRC1, CDYN, 
     .            VX, VY, VZ, VPARA, EELRC, 
     .            EIRENE_FEELRC1, SUMM, EISUMM, EISUM, SUM,
     .            X4, Y4, Z4, MOMPARA, BREMS, TOT_BREMS(NPLS), Z, BF,
     .            EIRENE_BREMS,XC,YC,ZC
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER :: IC1, IC2, ICELL, IAUSR, IBUSR, IRUSR, IPUSR,
     .           ITUSR, IN, IIRC, IRC, IRRC, J, IT1, IT2, ISTEP, IFRC,
     .           IR2, IP1, IP2, IND, IR, IP, IT, IVOLSI, I,
     .           ICC, IR1, IVL, ISTR, IL, IU, IM, MXREC, MXPLS, IFPLS,
     .           IPLSTI, IPLSV, KK, ICCT
      INTEGER, SAVE :: ISTROLD=-1
      LOGICAL, ALLOCATABLE, SAVE :: LPLSSR(:)
C
C  AT ENTRY SAMVL0:
C    DEFINE THE CUMULATIVE DISTRIBUTION FUNCTION
C    FREC(IPLS,IRRC,ICELL) FOR EACH VOLUME SOURCE DISTRIBUTION, FOR SAMPLING
C    THE CELL INDEX ICELL OF THE VOLUME SOURCE PARTICLE.
C
C    A FEW GEOMETRICAL CONSTANTS FOR RANDOM SAMPLING
C    OF THE STARTING POINTS IN EACH CELL ARE PRE-COMPUTED
C
C    THE SOURCE STRENGTH FLUX(ISTRA) IS MODIFIED FOR THE
C    STRATA WITH NLVOL(ISTRA)=.TRUE.
C
C  AT ENTRY SAMVL1:
C    THE INITIAL CO-ORDINATES OF A TEST FLIGHT ARE SAMPLED,
C    AND THE CELL NUMBERS ARE COMPUTED
C
      ENTRY EIRENE_SAMVL0
C
 
      IF (.NOT.ALLOCATED(FREC)) THEN
 
C  LPLSSR(IPLS):
C  IDENTIFY THOSE IPLS WHICH NEED A VOLUME SOURCE DISTRIBUTION
 
        ALLOCATE (LPLSSR(NPLSI))
        LPLSSR = .FALSE.
        DO ISTR=1,NSTRAI
          IF (NLVOL(ISTR) .AND. NLPLS(ISTR)
     .        .AND. (FLUX(ISTR) > 0._DP)) THEN
            IPLS = NSPEZ(ISTR)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) THEN
c  nspez out of range: Set volumetric sources for ALL species 
              LPLSSR = .TRUE.
            ELSE
              LPLSSR(IPLS) = .TRUE.
            END IF
          END IF
        END DO

        MXREC=MAXVAL(NPRCI(1:NPLSI))
        MXPLS=COUNT(LPLSSR(1:NPLSI))
        ALLOCATE (FREC(0:MXPLS,0:MXREC,0:NRAD))
        ALLOCATE (VSOURC(NSRFS,0:NRAD))
        ALLOCATE (VSMXI(NSRFS))
        ALLOCATE (RQ21(N1ST))
        ALLOCATE (PS21(N2ND))
        IF (LEVGEO == 3) ALLOCATE (ASIMP(2,NRAD))
        ALLOCATE (ISOURC(NSRFS,0:NRAD))
        ALLOCATE (ICMX(NSRFS))
        ALLOCATE (IFREC(NPLS))
      END IF
C
      FREC=0.
      SREC=0.
      IFREC=0
C
      IFPLS=0
      DO 2 IPLS=1,NPLSI
        IF (LGPRC(IPLS,0).EQ.0) GOTO 2
        IF (.NOT.LPLSSR(IPLS)) GOTO 2
        IFPLS=IFPLS+1
        IFREC(IPLS)=IFPLS
        DO 3 IIRC=1,NPRCI(IPLS)
          IRRC=LGPRC(IPLS,IIRC)
          KK=NREARC(IRRC)
          ICCT=0
C  SPECTRAL CUT OFF FOR SOURCE RATE: ONLY FOR PHOTONS SO FAR.
          IF (KK.GT.0) THEN 
            ICCT=NREACT(KK)
          ENDIF
          DO 3 J=1,NSBOX
            ADD=0.
C  EXCLUDE DEAD CELLS (GRID CUTS, ISOLATED CELLS FROM COUPLE_.., ETC)
C  EXCLUDE IPLS-VACUUM CELLS

C  tabrc1 is a volumetric rate per ion  (1/s)
c  turn this into rate in (amp per ion), factor 'elch'
c  and then into source rate amp per cell, factor di * vol
            IF (NSTGRD(J).EQ.0.AND..NOT.LGVAC(J,IPLS)) THEN
              IF (NSTORDR >= NRAD) THEN
                ADD=TABRC1(IRRC,J)*DIIN(IPLS,J)*VOL(J)*ELCHA
              ELSE
                ADD=EIRENE_FTABRC1(IRRC,J)*DIIN(IPLS,J)*VOL(J)*ELCHA
              END IF
            END IF
C  SPECTRAL CUT OFF FOR SOURCE RATE (ONLY USED FOR PHOTONS SO FAR)
            IF (ICCT > 0)
     .        ADD = ADD*(XINTLEFT(ICCT,J) +
     .                   XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))

            FREC(IFPLS,IIRC,J)  =FREC(IFPLS,IIRC,J-1)+ADD
            SREC(IPLS,IRRC)     =SREC(IPLS,IRRC)+ADD
3       CONTINUE
2     CONTINUE
 
C  SUM OVER SPECIES AND RECOMBINATION TYPE INDICES
      DO 4 IPLS=1,NPLSI
        IF (LGPRC(IPLS,0).EQ.0) GOTO 4
        IF (.NOT.LPLSSR(IPLS)) GOTO 4
        IFPLS=IFREC(IPLS)
        DO 5 IIRC=1,NPRCI(IPLS)
          IRRC=LGPRC(IPLS,IIRC)
          SREC(IPLS,0)=SREC(IPLS,0)+SREC(IPLS,IRRC)
          SREC(0,IRRC)=SREC(0,IRRC)+SREC(IPLS,IRRC)
          SREC(0,0)   =SREC(0,0)   +SREC(IPLS,IRRC)
          DO 5 J=1,NSBOX
            FREC(IFPLS,0,J)=FREC(IFPLS,0,J)+FREC(IFPLS,IIRC,J)
5         CONTINUE
4     CONTINUE
C
C
      IF (TRCSOU.AND.IFPLS.GT.0) THEN
        EIO=0.
        EEL=0.
        MOM=0.
C
        DO 7 IPLS=1,NPLSI
          CDYN=CNDYNP(IPLS)
          IF (LGPRC(IPLS,0).EQ.0) GOTO 7
          IF (.NOT.LPLSSR(IPLS)) GOTO 7
          IFPLS=IFREC(IPLS)
          IPLSTI = MPLSTI(IPLS)
          IPLSV = MPLSV(IPLS)
          DO 6 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            KK=NREARC(IRRC)

            ICCT=0
            IF (KK.GT.0) THEN 
              ICCT=NREACT(KK)
            ENDIF
            DO 6 J=1,NSBOX
              IF (NSTGRD(J).EQ.0.AND..NOT.LGVAC(J,IPLS)) THEN
c  FREC is in Amp, so ADD is in: eV * Amp = Watt
                REC=FREC(IFPLS,IIRC,J)-FREC(IFPLS,IIRC,J-1)
                IF (REC.LE.0.D0) GOTO 6
                ADD=(1.5*TIIN(IPLSTI,J)+EDRIFT(IPLS,J))*REC
C  SPECTRAL CUT OFF, CURRENTLY ONLY FOR PHOTONS
                IF (ICCT > 0)
     .            ADD = ADD*(XINTLEFT(ICCT,J) +
     .                       XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))

                EIO(IPLS,IRRC)=EIO(IPLS,IRRC)-ADD
                EIO(IPLS,0)   =EIO(IPLS,0   )-ADD
                 
CDR  position x0,y0,z0 is not yet known here
cdr  take center of gravity in cell, if needed (last parameter (logical) in bfield.f
                xc=0.
                yc=0.
                zc=0.                
                CALL EIRENE_BFIELD (J, XC,YC,ZC, BX,BY,BZ, BF,.FALSE.)
                IF (INDPRO(4) == 8) THEN
                  CALL EIRENE_VECUSR(2,J,XC,YC,ZC,VX,VY,VZ,IPLS,.FALSE.)
                  VPARA=VX*BX+VY*BY+VZ*BZ
                  MOMPARA=VPARA*CDYN*SIGN(1._DP,VPARA)
                ELSE IF (INDPRO(5) == 8) THEN
                  VX = VXIN(IPLSV,J)
                  VY = VYIN(IPLSV,J)
                  VZ = VZIN(IPLSV,J)
                  VPARA=VX*BX+VY*BY+VZ*BZ
                  MOMPARA=VPARA*CDYN*SIGN(1._DP,VPARA)
                ELSE
                  MOMPARA=PARMOM(IPLS,J)
                ENDIF
                ADD=MOMPARA*REC

                IF (ICCT > 0)
     .            ADD = ADD*(XINTLEFT(ICCT,J) +
     .                       XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))

                MOM(IPLS,IRRC)=MOM(IPLS,IRRC)-ADD
                MOM(IPLS,0)   =MOM(IPLS,0   )-ADD
              ENDIF
6         CONTINUE
C
C  associated electron cooling/heating rate: eelrc: EV *CM**3/S
          DO 8 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            KK=NREARC(IRRC)

            ICCT=0
            IF (KK.GT.0) THEN 
              ICCT=NREACT(KK)
            ENDIF
            DO J=1,NSBOX
              ADD=0.D0
              IF (NSTGRD(J).EQ.0.AND..NOT.LGVAC(J,IPLS)) THEN
                IF (NSTORDR >= NRAD) THEN
                  EELRC = EELRC1(IRRC,J)
                ELSE
                  EELRC = EIRENE_FEELRC1(IRRC,J)
                END IF
c  Turn eV/s/particle into Watt/cell
                ADD=EELRC*DIIN(IPLS,J)*VOL(J)*ELCHA
C  SPECTRAL CUT OFF (PHOTONS ONLY)
                IF (ICCT > 0)
     .            ADD = ADD*(XINTLEFT(ICCT,J) +
     .                       XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))

              ENDIF
              
              EEL(IPLS,IRRC)=EEL(IPLS,IRRC)+ADD
              EEL(IPLS,0   )=EEL(IPLS,0   )+ADD
            ENDDO   !  nsbox loop 
8         CONTINUE  !  irrc loop 

cdr  testing internal CR model, using amjuel and h_colrad rates, nrrc=2,
cdr  with scaling factor 0.5 each. ....TEST OK, FEB 18, out again.
cdr   if (ipls.eq.1) then
c          do j=1,nsbox
c            write (iunout,*) j,eelrc1(1,j),eelrc1(2,j),tein(j),
c    .         dein(j),lgvac(j,1),nstgrd(j)
c          enddo
cdr   endif

7       CONTINUE    !  npls loop

C  BREMSSTRAHLUNG ORIGINATING FROM IONS IPLS, CHARGE Z=NCHRGP(IPLS) 
        TOT_BREMS = 0._DP
        DO IPLS=1,NPLSI
          IF (NCHRGP(IPLS) == 0) CYCLE
          Z = NCHRGP(IPLS)
          DO J = 1, NSBOX
            IF (LGVAC(J,NPLS+1).OR.LGVAC(J,IPLS)) CYCLE
            BREMS=EIRENE_BREMS(TEIN(J),DEIN(J),Z)  ! Watt per ion
            TOT_BREMS(IPLS) = TOT_BREMS(IPLS) +
     .                        BREMS*DIIN(IPLS,J)*VOL(J) ! Watt per cell
          END DO
        END DO
C
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'DIAGNOSTICS FROM SUBR. SAMVOL: '
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'VOLUME RECOMBINATION RATES INTEGRATED OVER'
        WRITE (iunout,*) 'ENTIRE COMPUTATIONAL GRID '
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'RECOMBINATION ION PARTICLE LOSS (AMP): '
        ITYP=4
        DO 10 IPLS=1,NPLSI
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 11 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL EIRENE_MASAJR('IPLS,IRRC, SREC         ',
     .                   TEXTS(ISPZ),IRRC,-SREC(IPLS,IRRC))
11        CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL EIRENE_MASAJR('IPLS,TOT., SREC(IPLS,0) ',
     .                   TEXTS(ISPZ),0   ,-SREC(IPLS,0))
          ENDIF
10      CONTINUE
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'RECOMBINATION ION ENERGY LOSS (WATT): '
        DO 12 IPLS=1,NPLSI
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 13 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL EIRENE_MASAJR('IPLS,IRRC,EIO           ',
     .                   TEXTS(ISPZ),IRRC,EIO(IPLS,IRRC))
13        CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL EIRENE_MASAJR('IPLS,TOT.,EIO(IPLS,0)   ',
     .                   TEXTS(ISPZ),0   ,EIO(IPLS,0))
          ENDIF
12      CONTINUE
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'RECOMBINATION ELECTRON ENERGY LOSS (WATT): '
        DO 14 IPLS=1,NPLSI
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 15 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL EIRENE_MASAJR('IPLS,IRRC,EEL           ',
     .                   TEXTS(ISPZ),IRRC,EEL(IPLS,IRRC))
15        CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL EIRENE_MASAJR('IPLS,TOT.,EEL(IPLS,0)   ',
     .                   TEXTS(ISPZ),0   ,EEL(IPLS,0))
          ENDIF
14      CONTINUE
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'RECOMBINATION PARALLEL MOMENTUM LOSS : '
        DO 16 IPLS=1,NPLSI
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 17 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL EIRENE_MASAJR('IPLS,IRRC,MOM           ',
     .                   TEXTS(ISPZ),IRRC,MOM(IPLS,IRRC))
17        CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL EIRENE_MASAJR('IPLS,TOT.,MOM(IPLS,0)   ',
     .                   TEXTS(ISPZ),0   ,MOM(IPLS,0))
          ENDIF
16      CONTINUE
        CALL EIRENE_LEER(1)
 
        WRITE (iunout,*) 'BREMSSTRAHLUNG (WATT): '
        DO IPLS=1,NPLSI
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          CALL EIRENE_MASAJR('IPLS,TOT.BREMSSTRAHLUNG ',
     .                 TEXTS(ISPZ),0   ,TOT_BREMS(IPLS))
        END DO
 
      ENDIF    !trcsou
C
C  SET TOTAL SOURCE STRENGTH FOR STRATA WITH NLVOL(ISTRA)=.TRUE.,
C
      DO 50 ISTRA=1,NSTRAI
        IF (NLVOL(ISTRA).AND.NLPLS(ISTRA)) THEN
          IPLS=NSPEZ(ISTRA)
          IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) THEN
            WRITE (iunout,*) 'SOURCE SPECIES INDEX NSPEZ OUT OF RANGE'
            WRITE (iunout,*) 'ISTRA, NSPEZ(ISTRA) ',ISTRA,NSPEZ(ISTRA)
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
          IPLSTI = MPLSTI(IPLS)
          SUMM=0.D0
          EISUMM=0.D0
C  VOLUMETRIC SUB-STRATA  
          DO 53 IVOLSI=1,NSRFSI(ISTRA)
            IVL=IVOLSI
            SUM=0.D0
            EISUM=0.D0
            IF (SORLIM(IVL,ISTRA).LT.0) THEN
C  INITIALIZE SAMPLING DISTRIBUTIONS FOR USER SPECIFIED VOLUME SOURCE
              CALL EIRENE_SM0USR(IVL,ISTRA,
     .                    SORAD1(IVL,ISTRA),SORAD2(IVL,ISTRA),
     .                    SORAD3(IVL,ISTRA),SORAD4(IVL,ISTRA),
     .                    SORAD5(IVL,ISTRA),SORAD6(IVL,ISTRA))
!pb assume flux is set in samusr
              SUMM=FLUX(ISTRA)
            ELSE
C  INITIALIZE SAMPLING DISTRIBUTIONS FOR DEFAULT VOLUME RECOMBINATION SOURCES
C  ACCOUNT FOR INGRDA(IVOLSI,ISTRA,...), INGRDE(IVOLSI,ISTRA,...)
              I=ISTRA
              ICC=0
              IRC=-1
              IF (NR1ST.GT.1) THEN
              IF (INGRDA(IVL,I,1).LE.0..OR.INGRDE(IVL,I,1).LE.0.D0) THEN
                CALL EIRENE_LEER(1)
                WRITE (iunout,*) 'WARNING FROM SAMVL0, ISTRA= ',ISTRA
                WRITE (iunout,*)
     .            'NEW INPUT FOR INGRDA(.,.,1),INGRDE(.,.,1)'
                WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT '
                INGRDA(IVL,I,1)=1
                INGRDE(IVL,I,1)=MAX0(1,NR1ST)
                CALL EIRENE_LEER(1)
              ENDIF
              ENDIF
              IF (NP2ND.GT.1) THEN
              IF (INGRDA(IVL,I,2).LE.0..OR.INGRDE(IVL,I,2).LE.0.D0) THEN
                CALL EIRENE_LEER(1)
                WRITE (iunout,*) 'WARNING FROM SAMVL0, ISTRA= ',ISTRA
                WRITE (iunout,*)
     .            'NEW INPUT FOR INGRDA(.,.,2),INGRDE(.,.,2)'
                WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT '
                INGRDA(IVL,I,2)=1
                INGRDE(IVL,I,2)=MAX0(1,NP2ND)
                CALL EIRENE_LEER(1)
              ENDIF
              ENDIF
              IF (NT3RD.GT.1) THEN
              IF (INGRDA(IVL,I,3).LE.0..OR.INGRDE(IVL,I,3).LE.0.D0) THEN
                CALL EIRENE_LEER(1)
                WRITE (iunout,*) 'WARNING FROM SAMVL0, ISTRA= ',ISTRA
                WRITE (iunout,*)
     .            'NEW INPUT FOR INGRDA(.,.,3),INGRDE(.,.,3)'
                WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT '
                INGRDA(IVL,I,3)=1
                INGRDE(IVL,I,3)=MAX0(1,NT3RD)
                CALL EIRENE_LEER(1)
              ENDIF
              ENDIF
              IF (NPRCI(IPLS).EQ.0) THEN
                WRITE (iunout,*) 'NO DEFAULT VOLUME SOURCE DISTRIBUTION'
                WRITE (iunout,*) 'DEFINED. SUBSTRATUM TURNED OFF'
                WRITE (iunout,*) 'IPLS,IVOLSI,ISTRA ',IPLS,IVOLSI,ISTRA
                SORWGT(IVL,ISTRA)=0.D0
                GOTO 53
              ENDIF
              IF (NLRAD) THEN
                IR1=MAX0(1,INGRDA(IVL,ISTRA,1))
                IR2=MIN0(NR1ST,INGRDE(IVL,ISTRA,1))
              ELSE
                IR1=1
                IR2=2
              ENDIF
              IF (NLPOL) THEN
                IP1=MAX0(1,INGRDA(IVL,ISTRA,2))
                IP2=MIN0(NP2ND,INGRDE(IVL,ISTRA,2))
              ELSE
                IP1=1
                IP2=2
              ENDIF
              IF (NLTOR) THEN
                IT1=MAX0(1,INGRDA(IVL,ISTRA,3))
                IT2=MIN0(NT3RD,INGRDE(IVL,ISTRA,3))
              ELSE
                IT1=1
                IT2=2
              ENDIF
 
              ISTEP=SORIND(IVL,ISTRA)
              IFPLS=IFREC(IPLS)
              DO 52 IIRC=1,NPRCI(IPLS)
                IRRC=LGPRC(IPLS,IIRC)
                IF (ISTEP.EQ.IRRC) THEN
C  RECOMBINATION PROCESS IRRC IDENTIFIED
                  IRC=IRRC
                  IFRC=IIRC
                ELSEIF (ISTEP.EQ.0) THEN
C  SUM OVER ALL RECOMBINATION PROCESSES FOR SPECIES IPLS
                  IRC=0
                  IFRC=0
                  ISTEP=-1
                ELSE
C  TRY OTHER RECOMBINATION PROCESS ASSIGNED TO IPLS
                  GOTO 52
                ENDIF
                DO 51 IR=IR1,IR2-1
                  DO 51 IP=IP1,IP2-1
                    DO 51 IT=IT1,IT2-1
                      NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                      REC=FREC(IFPLS,IFRC,NCELL)-
     .                    FREC(IFPLS,IFRC,NCELL-1)
C  INDIRECT ADDRESSING
                      IF (REC.GT.0.D0) THEN
                        ICC=ICC+1
                        SUM=SUM+REC
                        EISUM=EISUM-
     .                   (1.5*TIIN(IPLSTI,NCELL)+EDRIFT(IPLS,NCELL))*REC
                      ENDIF
51              CONTINUE
52            CONTINUE   ! summing over irrc
c
              IF (SUM.EQ.0.D0) THEN
                WRITE (IUNOUT,*) 'NO VOL. RECOMBINATION SOURCE FOR: '
                WRITE (IUNOUT,*) 'ISTRA, IVOLSI, IPLS, ISTEP ',
     .                            ISTRA, IVL   , IPLS, ISTEP
                WRITE (IUNOUT,*) 'EITHER: ISTEP OUT OF RANGE IN SAMVOL'
                WRITE (IUNOUT,*) 'OR:  DENSITY OF RECOMBINING IPLS = 0 '
                SORWGT(IVL,ISTRA)=0.D0
                GOTO 53
              ENDIF
              SORWGT(IVL,ISTRA)=SUM
              CALL EIRENE_LEER(1)
              WRITE (iunout,*) 'SUB-STRATUM WEIGHT REDEFINED '
              CALL EIRENE_MASJ2R
     .          ('IVOLSI,ISTRA,SORWGT     ',IVOLSI,ISTRA,SUM)
              IF (TRCSOU) THEN
                CALL EIRENE_MASJ3 ('IRRC,IPLS,ICMX          ',
     .                              IRC ,IPLS,ICC)
                CALL EIRENE_LEER(1)
              ENDIF
              SUMM=SUMM+SUM
              EISUMM=EISUMM+EISUM
            ENDIF
53        CONTINUE
C
          IF (SUMM.GT.0.D0) THEN
            FLUX(ISTRA)=SUMM
            WRITE (iunout,*) 'SOURCE STRENGTH REDEFINED '
            CALL EIRENE_MASJR2('ISTRA, FLUX, EIFLUX     ',
     .                   ISTRA,FLUX(ISTRA),EISUMM)
            CALL EIRENE_LEER(1)
          ELSE
            FLUX(ISTRA)=0.D0
            WRITE (iunout,*) 'SOURCE ISTRA= ',ISTRA,' TURNED OFF '
            CALL EIRENE_LEER(1)
          ENDIF
        ENDIF
50    CONTINUE
C
C  PREPARE SOME GEOMETRICAL CONSTANTS FOR RANDOM SAMPLING IN STANDARD MESH CELLS
      select case (LEVGEO)
      case (2)
        IF (NLPOL) THEN
          DO 54 IP=1,NP2NDM
            PS21(IP)=PSURF(IP+1)-PSURF(IP)
54        CONTINUE
        ENDIF
        DO 55 IR=1,NR1STM
          RQ21(IR)=RQ(IR+1)-RQ(IR)
55      CONTINUE

      case (3)
c  split quadrangle into two triangles, 
c  then 1st: sample triangle according to its relative area, 
c  then 2nd: sample uniform within this triangle
        IT=1
        DO 56 IR=1,NR1ST-1
        DO 56 IP=1,NP2ND-1
          IND=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          X1=XPOL(IR,IP)
          X2=XPOL(IR,IP+1)
          X3=XPOL(IR+1,IP+1)
          Y1=YPOL(IR,IP)
          Y2=YPOL(IR,IP+1)
          Y3=YPOL(IR+1,IP+1)
          ASIMP(1,IND)=0.5*(X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2))
          X1=XPOL(IR+1,IP)
          X2=XPOL(IR,IP)
          X3=XPOL(IR+1,IP+1)
          Y1=YPOL(IR+1,IP)
          Y2=YPOL(IR,IP)
          Y3=YPOL(IR+1,IP+1)
          ASIMP(2,IND)=0.5*(X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2))
56      CONTINUE
      end select
C
      RETURN
C
C  AT THIS POINT: CALLED FROM PARTICLE LOOP TO INITIALIZE TEST FLIGHT
C
      ENTRY EIRENE_SAMVL1
     .      (NVLM,TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)
C  USER SUPPLIED SOURCE
C
      IF (SORLIM(NVLM,ISTRA).LT.0) THEN
        CALL EIRENE_SM1USR(NVLM,X0,Y0,Z0,
     .              SORAD1(NVLM,ISTRA),SORAD2(NVLM,ISTRA),
     .              SORAD3(NVLM,ISTRA),SORAD4(NVLM,ISTRA),
     .              SORAD5(NVLM,ISTRA),SORAD6(NVLM,ISTRA),
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)
        NRCELL=IRUSR
        NPCELL=IPUSR
        NTCELL=ITUSR
        NACELL=IAUSR
        NBLOCK=IBUSR
        NBLCKA=NSTRD*(NBLOCK-1)+NACELL
        NCELL=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
C
        MTSURF=0
        NLSRFZ=.FALSE.
        MPSURF=0
        NLSRFY=.FALSE.
        MRSURF=0
        NLSRFX=.FALSE.
        EFWL(1:NPLSI)=0._DP
        RETURN
      ENDIF
C
C  VOLUME RECOMBINATION SOURCE
C
C  TENTATIVELY ASSUME: A BULK ION WILL BE GENERATED
      LGPART=.TRUE.
      ITYP=4
C
      IF (.NOT.NLPLS(ISTRA)) GOTO 999
 
      IF (ISTROLD /= ISTRA) THEN
        ISTROLD=ISTRA
        IPLS=NSPEZ(ISTRA)
        DO IVOLSI=1,NSRFSI(ISTRA)
          IVL=IVOLSI
          ICC=0
          VSOURC(IVL,0)=0.D0
          IF (NLRAD) THEN
            IR1=MAX0(1,INGRDA(IVL,ISTRA,1))
            IR2=MIN0(NR1ST,INGRDE(IVL,ISTRA,1))
          ELSE
            IR1=1
            IR2=2
          ENDIF
          IF (NLPOL) THEN
            IP1=MAX0(1,INGRDA(IVL,ISTRA,2))
            IP2=MIN0(NP2ND,INGRDE(IVL,ISTRA,2))
          ELSE
            IP1=1
            IP2=2
          ENDIF
          IF (NLTOR) THEN
            IT1=MAX0(1,INGRDA(IVL,ISTRA,3))
            IT2=MIN0(NT3RD,INGRDE(IVL,ISTRA,3))
          ELSE
            IT1=1
            IT2=2
          ENDIF
 
          ISTEP=SORIND(IVL,ISTRA)
          IFPLS=IFREC(IPLS)
          DO IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            IF (ISTEP.EQ.IRRC) THEN
              IRC=IRRC
              IFRC=IIRC
            ELSEIF (ISTEP.EQ.0) THEN
C  SUM OVER ALL RECOMBINATION PROCESSES FOR SPECIES IPLS
              IRC=0
              IFRC=0
              ISTEP=-1
            ELSE
              CYCLE
            ENDIF
            DO IR=IR1,IR2-1
              DO IP=IP1,IP2-1
                DO IT=IT1,IT2-1
                  NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  ADD=FREC(IFPLS,IFRC,NCELL)-
     .                FREC(IFPLS,IFRC,NCELL-1)
C  INDIRECT ADDRESSING
                  IF (ADD.GT.0.D0) THEN
                    ICC=ICC+1
                    ISOURC(IVL,ICC)=NCELL
                    VSOURC(IVL,ICC)=VSOURC(IVL,ICC-1)+ADD
                  ENDIF
                END DO
              END DO
            END DO
          END DO ! IIRC
          ICMX(IVL)=ICC
          VSMXI(IVL) = 1._DP / VSOURC(IVL,ICC)
        END DO ! IVOLSI
C
      END IF
C
C  FIND CELL NUMBER: NCELL
C
      IF (INDIM(NVLM,ISTRA) .GE. 0) THEN
cdr analog sampling, no weighting
!PB  choose cell according to cell contribution to total source strength
        IC1=0
        IC2=ICMX(NVLM)
        ZEP1=RANF_EIRENE()*VSOURC(NVLM,IC2)
 
        IL=0
        IU=IC2
 
c  binary search
        DO WHILE (IU-IL.gt.1)
          IM=(IU+IL)*0.5
          IF(ZEP1.GE.VSOURC(NVLM,IM)) THEN
            IL=IM
          ELSE
            IU=IM
          ENDIF
        END DO
c
        ICELL=IU
 
        NCELL=ISOURC(NVLM,ICELL)
      ELSE
 
cdr non-analog sampling.  
cdr Here use uniform distribution of cell indices and weighting
cdr tbd: correlation sampling: use previous (reference) distribution and weighting
cdr      rather than uniform sampling. 
        IC1=0
        IC2=ICMX(NVLM)
 
        ICELL = MIN(INT(1+RANF_EIRENE()*(IC2-1)),IC2)
        NCELL = ISOURC(NVLM,ICELL)
 
        WEIGHT=(VSOURC(NVLM,ICELL)-VSOURC(NVLM,ICELL-1))*VSMXI(NVLM)*IC2
 
      END IF
C
      IF (NCELL.GT.NSURF) GOTO 991
C
C  A CELL NUMBER NCELL HAS NOW BEEN SAMPLED
C
      CALL EIRENE_NCELLN(NCELL,NRCELL,NPCELL,NTCELL,NACELL,NBLOCK,
     .            NR1ST,NP2ND,NT3RD,NBMLT,NLRAD,NLPOL,NLTOR)
C
C  FIND TOROIDAL CO-ORDINATE IN NTCELL
C
      IF (.NOT.NLTOR) THEN
C       NTCELL=1
        IPERID=1
        IF (NLTRZ) THEN
          Z0=0.
        ELSEIF (NLTRA) THEN
C  TACTICALLY ASSUME: PARTICLE STARTS IN LOCAL TOR. BASIS CELL NO.1
          ZRM1=ZSURF(1)
          PHI=ZRM1+RANF_EIRENE()*ZFULL
          IPERID=1
C         Z0=??, TO BE FOUND FROM X01,PHI LATER
C         IPERID=LEARCA(PHI,ZSURF,1,NTTRA,1,'SAMVOL      ')
        ELSEIF (NLTRT) THEN
          GOTO 999
        ENDIF
      ELSEIF (NLTOR) THEN
        IPERID=NTCELL
C  SAMPLE IN CELL NTCELL
        IF (NLTRZ) THEN
          Z0=ZSURF(NTCELL)+RANF_EIRENE()*(ZSURF(NTCELL+1)-ZSURF(NTCELL))
        ELSEIF (NLTRT) THEN
          PHI=ZSURF(NTCELL)+RANF_EIRENE()*
     .        (ZSURF(NTCELL+1)-ZSURF(NTCELL))
C         Z0=??, TO BE FOUND FROM X01,PHI LATER
        ELSEIF (NLTRA) THEN
          ZRM1=ZFULL*(NTCELL-1)
          PHI=ZRM1+RANF_EIRENE()*ZFULL
C         Z0=??, TO BE FOUND FROM X01,PHI LATER
        ENDIF
      ENDIF
C
C  FIND RADIAL AND POLOIDAL CO-ORDINATE
C
      select case (LEVGEO)
      case (1)
        X0=RSURF(NRCELL)+RANF_EIRENE()*(RSURF(NRCELL+1)-RSURF(NRCELL))
        IF (NLPOL) THEN
          Y0=PSURF(NPCELL)+RANF_EIRENE()*(PSURF(NPCELL+1)-PSURF(NPCELL))
        ELSE
          Y0=YIA+RANF_EIRENE()*(YAA-YIA)
        END IF
C..........................................................................
      case (2)
        IF (NLCRC) THEN
C  POLOIDAL CO-ORDINATE
          IF (NLPOL) THEN
            WINK=PSURF(NPCELL)+RANF_EIRENE( )*PS21(NPCELL)
          ELSEIF (.NOT.NLPOL) THEN
            WINK=RANF_EIRENE( )*PI2A
          ENDIF
C  RADIAL CO-ORDINATE
          RR=SQRT(RQ(NRCELL)+RANF_EIRENE( )*RQ21(NRCELL))
C
          X0=RR*COS(WINK)
          Y0=RR*SIN(WINK)
        ELSEIF (NLELL) THEN
CDR NOT READY. STRICKLY, THETA AND R ARE CORRELATED. USE
CDR            MARGINAL AND CONDITIONAL DISTRIBUTION F1(R) AND
CDR            F2(PHI, GIVEN R)
C  POLOIDAL CO-ORDINATE
          IF (NLPOL) THEN
            WINK=PSURF(NPCELL)+RANF_EIRENE( )*PS21(NPCELL)
          ELSEIF (.NOT.NLPOL) THEN
            WINK=RANF_EIRENE( )*PI2A
          ENDIF
C  RADIAL CO-ORDINATE
          RR=SQRT(RQ(NRCELL)+RANF_EIRENE( )*RQ21(NRCELL))
C
          RRI=RSURF(NRCELL)
          RRD=RSURF(NRCELL+1)-RRI
          RRN=(RR-RRI)/RRD
C
          ELR=ELL(NRCELL)+RRN*(ELL(NRCELL+1)-ELL(NRCELL))
          EPR=EP1(NRCELL)+RRN*(EP1(NRCELL+1)-EP1(NRCELL))
          X0=RR*COS(WINK)+EPR
          Y0=RR*SIN(WINK)*ELR
        ELSEIF (NLTRI) THEN
          GOTO 999
        ENDIF
C...................................................................  
      case (3)
        IF (.NOT.NLPOL) THEN
          GOTO 999
        ENDIF
        IN = NRCELL + (NPCELL-1)*NR1ST
        ZEP1=AREA(IN)*RANF_EIRENE()
        IF (ZEP1.LE.ASIMP(1,IN)) THEN
C   POINT TO BE SAMPLED WITHIN TRIANGLE 1
          X1=XPOL(NRCELL,NPCELL)
          X2=XPOL(NRCELL,NPCELL+1)
          X3=XPOL(NRCELL+1,NPCELL+1)
          Y1=YPOL(NRCELL,NPCELL)
          Y2=YPOL(NRCELL,NPCELL+1)
          Y3=YPOL(NRCELL+1,NPCELL+1)
        ELSE
C   POINT TO BE SAMPLED WITHIN TRIANGLE 2
          X1=XPOL(NRCELL+1,NPCELL)
          X2=XPOL(NRCELL,NPCELL)
          X3=XPOL(NRCELL+1,NPCELL+1)
          Y1=YPOL(NRCELL+1,NPCELL)
          Y2=YPOL(NRCELL,NPCELL)
          Y3=YPOL(NRCELL+1,NPCELL+1)
        ENDIF
        IPOLG=NPCELL
        Z1=0.
        Z2=0.
        Z3=0.
        CALL EIRENE_FPOLYT_3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X0,Y0,ZZ)

C...................................................................  
      case (4)
        X1=XTRIAN(NECKE(1,NCELL))
        X2=XTRIAN(NECKE(2,NCELL))
        X3=XTRIAN(NECKE(3,NCELL))
        Y1=YTRIAN(NECKE(1,NCELL))
        Y2=YTRIAN(NECKE(2,NCELL))
        Y3=YTRIAN(NECKE(3,NCELL))
        Z1=0.
        Z2=0.
        Z3=0.
        CALL EIRENE_FPOLYT_3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X0,Y0,ZZ)
C.................................................................
      case (5)
        X1=XTETRA(NTECK(1,NCELL))
        Y1=YTETRA(NTECK(1,NCELL))
        Z1=ZTETRA(NTECK(1,NCELL))
        X2=XTETRA(NTECK(2,NCELL))
        Y2=YTETRA(NTECK(2,NCELL))
        Z2=ZTETRA(NTECK(2,NCELL))
        X3=XTETRA(NTECK(3,NCELL))
        Y3=YTETRA(NTECK(3,NCELL))
        Z3=ZTETRA(NTECK(3,NCELL))
        X4=XTETRA(NTECK(4,NCELL))
        Y4=YTETRA(NTECK(4,NCELL))
        Z4=ZTETRA(NTECK(4,NCELL))
        CALL
     .  EIRENE_FPOLYT_4(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,X0,Y0,Z0)
C....................................................................
      case (10)
        WRITE (iunout,*) 'ERROR EXIT FROM SAMVOL. LEVGEO ',LEVGEO
        WRITE (iunout,*) 'TO BE DONE: RETURN CENTER OF GRAVITY IN NCELL'
        CALL EIRENE_EXIT_OWN(1)
      end select
C
      IF (NLTRA) THEN
C  FIND Z0 FROM X01,PHI IN LOCAL TOROIDAL CELL NTCELL
        X01=X0+RMTOR
        CALL EIRENE_FZRTRI(X0,Z0,NTCELL,X01,PHI,NTCELL)
      ENDIF
C
      MTSURF=0
      NLSRFZ=.FALSE.
      MPSURF=0
      NLSRFY=.FALSE.
      MRSURF=0
      NLSRFX=.FALSE.
C
C  NEXT: ANALOG SPECIES INDEX DISTRIBUTION: WEISPZ(IPL)
C
      DO 630 ISPZ=1,NSPZ
        WEISPZ(ISPZ)=-1.
630   CONTINUE
C
C  NOT IN USE ANYMORE
C  CURRENTLY: ONLY SINGLE SPECIES VOLUME SOURCES POSSIBLE
C  MULTI SPECIES VOL-SOURCES HAVE TO BE TREATED BY STRATIFIED SAMPLING
C     IF (NSPEZ(ISTRA).LE.0) THEN
C       IF (NCELL.EQ.1) THEN
C         DO 640 IPL=1,NPLSI
C           IREC=0
C           IFPLS=IFREC(IPLS)
C           WEISPZ(IPL)=(FREC(IFPLS,0,1))/
C    .                  (FREC(0,  0,1))
C           IF (WEISPZ(IPL).LT.0) GOTO 991
640       CONTINUE
C       ELSE
C         DO 645 IPL=1,NPLSI
C           IFPLS=IFREC(IPLS)
C           WEISPZ(IPL)=(FREC(IFPLS,0,NCELL)-FREC(IFPLS,0,NCELL-1))/
C     .                 (FREC(0,    0,NCELL)-FREC(0,    0,NCELL-1))
C           IF (WEISPZ(IPL).LT.0) GOTO 991
645       CONTINUE
C       ENDIF
C     ENDIF
C
      CRTX=SORAD4(NVLM,ISTRA)
      CRTY=SORAD5(NVLM,ISTRA)
      CRTZ=SORAD6(NVLM,ISTRA)
      CNORM=SQRT(CRTX**2+CRTY**2+CRTZ**2)+EPS60
      CRTX=CRTX/CNORM
      CRTY=CRTY/CNORM
      CRTZ=CRTZ/CNORM
!PB
      TEWL=TEIN(NCELL)
      TIWL(1:NPLSI)=TIIN(MPLSTI(1:NPLSI),NCELL)
      DIWL(1:NPLSI)=DIIN(1:NPLSI,NCELL)
      VXWL(1:NPLSI)=VXIN(MPLSV(1:NPLSI),NCELL)
      VYWL(1:NPLSI)=VYIN(MPLSV(1:NPLSI),NCELL)
      VZWL(1:NPLSI)=VZIN(MPLSV(1:NPLSI),NCELL)
      EFWL(1:NPLSI)=0._DP
      SHWL=0._DP
C
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMVOL'
      CALL EIRENE_EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'SAMPLING ERROR IN SAMVOL'
      WRITE (iunout,*) 'NCELL,NSURF,NSBOX ',NCELL,NSURF,NSBOX
      CALL EIRENE_EXIT_OWN(1)
997   CONTINUE
      WRITE (iunout,*) 'SORIND (=IRRC) OUT OF RANGE IN SAMVOL'
      WRITE (iunout,*) 'IRRC,NREC ',IRRC,NREC
      CALL EIRENE_EXIT_OWN(1)
999   CONTINUE
      WRITE (iunout,*) 'UNWRITTEN OPTION IN SAMVOL'
      CALL EIRENE_EXIT_OWN(1)
 
C     the following ENTRY is for reinitialization of EIRENE (DMH)
 
      ENTRY EIRENE_SAMVOL_REINIT
 
      ISTROLD = -1
 
      DEALLOCATE (LPLSSR)
      DEALLOCATE (FREC)
      DEALLOCATE (VSOURC)
      DEALLOCATE (VSMXI)
      DEALLOCATE (RQ21)
      DEALLOCATE (PS21)
      IF (LEVGEO == 3) DEALLOCATE (ASIMP)
      DEALLOCATE (ISOURC)
      DEALLOCATE (ICMX)
      DEALLOCATE (IFREC)
 
      return
 
      END
