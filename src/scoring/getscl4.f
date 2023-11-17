C  oct.2014  only comments added
cdr  aug 2016  comments only:
cdr            potati, potmli, potioi and potphi (OUTGOING FLUXES)
cdr            include also the fluxes onto census.
cdr  hence: in t-dep mode and nlscl=T, total particle balances should be exact
cdr  to be done: apply that scaling also to flux and weights in census for re-sampling
cdr  in subr. tmstep (to be done)
cdr  Aug. 2020: if a type of particle (e.g., test ions), is very unlikely to be sampled,
cdr             but still has a finite (non-zero) source rate from tracklength estimator,
cdr             then the particle balance matrix may become singular (source, but no sinks).
cdr             A safety has been added now in one place.
cdr             Perhaps still needed in
cdr             more cases (icol,irow) in this routine?
C
      SUBROUTINE EIRENE_GETSCL4 (ISTRA,FA,FM,FI,FPH)
C
C  FIND SCALING FACTORS TO ENFORCE PARTICLE BALANCE,
c  BECAUSE OF NON-CONSERVATIVE PROPERTY OF TRACKLENGTH ESTIMATORS
C  SIMPLE VERSION: NOT SPLIT BY SPECIES, ONLY BY TYPE.
C  THE OUTPUT SCALING FACTORS FC=(FA,FM,FI,FPH) ARE "PER STRATUM"
c
c  build a 4 x 4 matrix P of 4 global particle balance equations,
c  one for atoms, molecules, test ions and photons each.
c
c  solve for the 4 factors FC(1)...FC(4), such that
c
c    P * FC = B
c
c  B is the external (fixed, given) source strength for each type of particle,
c  per stratum, B = B(1),....,B(4).
c
c  The factor FC(1) is then to be applied to all tallies (volumetric or surface fluxes)
c  which scale linearly with the external source B(1) for particles of type 1 (i.e. for "atoms").
c  Similarly for the other types, FC(2),... etc...

c  most of the programming below deals with possible zeroes in rows and columns,
c  i.e. with cases that some type of particle (e.g. photons) may not be present.
c
c
c
C  MODIFIED JAN/95: INCLUDE SURFACE TALLIES IN MATRIX, NOT IN
C  INHOMOGENEITY
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_COUTAU
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CTRCEI, ONLY: TRCSCL

      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: ISTRA
      REAL(DP), INTENT(OUT) :: FA, FM, FI, FPH
      REAL(DP) :: FC(4), P(4,5), B(4), PP(4,4), FFC(3), BB(3)
      REAL(DP) :: DTB1, DTB2, DTA, DTB3, FNEN, DTB4,
     .            EIRENE_DETER3X3, EIRENE_DETER4X4
      REAL(DP) :: P11, P12, P13, P21, P22, P23, P31, P32, P33,
     .            B1, B2, B3,
     .            ap3ma, ap3mi, ap3m, ap3n3, ap3n2
      INTEGER :: ICOL, IROW, I, J, I1, I2, J1, J2, IC, JOUT, IOUT,
     .           ICOL2, IROW2, I3
      LOGICAL :: LCOLM(4), LROW(4), lcol2(3), lrow2(3)
C
C
      FC(1)=1.
      FC(2)=1.
      FC(3)=1.
      FC(4)=1.
C
C P(..,1)*FC(1)
C P(..,2)*FC(2)
C P(..,3)*FC(3)
C P(..,4)*FC(4)
C
cdr  atomic sinks (1,1) and sources from other types
cdr  surface tallies: net sink: potati+prfaai
      P(1,1)=PAATI(0,ISTRA)+POTATI(0,ISTRA)+PRFAAI(0,ISTRA)+
     .       PGENAI(0,ISTRA)
      P(1,2)=PMATI(0,ISTRA)+PRFMAI(0,ISTRA)
      P(1,3)=PIATI(0,ISTRA)+PRFIAI(0,ISTRA)
      P(1,4)=0._DP
cdr  molecular sinks (2,2) and sources from other types
      P(2,1)=PAMLI(0,ISTRA)+PRFAMI(0,ISTRA)
cdr  surface tallies: net sink: potmli+prfmmi
      P(2,2)=PMMLI(0,ISTRA)+POTMLI(0,ISTRA)+PRFMMI(0,ISTRA)+
     .       PGENMI(0,ISTRA)
      P(2,3)=PIMLI(0,ISTRA)+PRFIMI(0,ISTRA)
      P(2,4)=0._DP
cdr  test-ion sinks(3,3) and sources from other types
      P(3,1)=PAIOI(0,ISTRA)+PRFAII(0,ISTRA)
      P(3,2)=PMIOI(0,ISTRA)+PRFMII(0,ISTRA)
cdr  surface tallies: net sink: potioi+prfiii
      P(3,3)=PIIOI(0,ISTRA)+POTIOI(0,ISTRA)+PRFIII(0,ISTRA)+
     .       PGENII(0,ISTRA)
      P(3,4)=0._DP
Cdr PHOTONIC TALLIES ARE CURRENTLY NOT INCLUDED IN RESCALING. TO BE DONE
CDR FC(4) SHOULD ALWAYS TURN OUT TO BE EXACTLY 1.0
cdr  photon sinks(4,4) and sources from other types
      P(4,1)=0._DP
      P(4,2)=0._DP
      P(4,3)=0._DP
cdr  surface tallies: net sink: potphti+prfphphti
      P(4,4)=1._DP

      IF (TRCSCL) THEN
        write (iunout,'(A,ES14.7)') 'PAATI   ',PAATI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PRFAAI  ',PRFAAI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'POTATI  ',POTATI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PGENAI  ',PGENAI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PMATI   ',PMATI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PRFMAI  ',PRFMAI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PIATI   ',PIATI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PRFIAI  ',PRFIAI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PAMLI   ',PAMLI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PRFAMI  ',PRFAMI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PMMLI   ',PMMLI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PRFMMI  ',PRFMMI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'POTMLI  ',POTMLI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PGENMI  ',PGENMI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PIMLI   ',PIMLI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PRFIMI  ',PRFIMI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PAIOI   ',PAIOI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PRFAII  ',PRFAII(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PMIOI   ',PMIOI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PRFMII  ',PRFMII(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PIIOI   ',PIIOI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PRFIII  ',PRFIII(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'POTIOI  ',POTIOI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PGENII  ',PGENII(0,ISTRA)
      END IF
C
cdr  (direct) primary sources, and secondaries from primary bulk particles
      B(1)=-(PPATI(0,ISTRA)+WTOTA(0,ISTRA))
      B(2)=-(PPMLI(0,ISTRA)+WTOTM(0,ISTRA))
      B(3)=-(PPIOI(0,ISTRA)+WTOTI(0,ISTRA))
      B(4)=1._DP
      IF (TRCSCL) THEN
        write (iunout,'(A,ES14.7)') 'PPATI   ',PPATI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'WTOTA   ',WTOTA(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PPMLI   ',PPMLI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'WTOTM   ',WTOTM(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'PPIOI   ',PPIOI(0,ISTRA)
        write (iunout,'(A,ES14.7)') 'WTOTI   ',WTOTI(0,ISTRA)
      END IF

      IF (TRCSCL) THEN
        CALL EIRENE_MASAGE('GETSCL4                             ')
        CALL EIRENE_MASRR1('MATRIX    :',P,4*4,4)
        CALL EIRENE_MASRR1('RHS VECTOR:',B,4,4)
        CALL EIRENE_LEER(1)
      END IF
C
      ICOL=0
      IROW=0
      DO 1 I=1,4
        LROW(I)=P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,4)**2.GT.EPS30
        LCOLM(I)=P(1,I)**2+P(2,I)**2+P(3,I)**2+P(4,I)**2.GT.EPS30
        IF (LROW(I)) IROW=IROW+1
        IF (LCOLM(I)) ICOL=ICOL+1
    1 CONTINUE
C
      IF (IROW.EQ.0) THEN
C  NO ROW IS NONZERO, I.E. NO TEST PARTICLES (OF ANY TYPE) FOLLOWED
        GOTO 1000
C
      ELSEIF (IROW.EQ.1) THEN
C  ONLY ONE ROW (NO. I) IS NONZERO, I.E., ONLY ATOMS, ONLY MOLECULES
C                          ONLY TEST IONS OR ONLY PHOTONS ARE FOLLOWED
C                          FOR THE PRESENT STRATUM.
         DO 10 I=1,4
           IF (LROW(I)) THEN
             IF (LCOLM(1)) THEN
               FC(1)=(B(I)-P(I,2)-P(I,3)-P(I,4))/P(I,1)
             ELSEIF (LCOLM(2)) THEN
               FC(2)=(B(I)-P(I,3)-P(I,4))/P(I,2)
             ELSEIF (LCOLM(3)) THEN
               FC(3)=(B(I)-P(I,4))/P(I,3)
             ELSE
               FC(4)=B(I)/P(I,4)
             ENDIF
           ENDIF
   10    CONTINUE
C
C
      ELSEIF (IROW.EQ.2) THEN
C  TWO ROWS ARE NONZERO
        J1=0
        J2=0
C  DETERMINE THE INDICES FOR THE NONZERO ROWS
        DO 20 J=1,4
          IF (LROW(J)) THEN
            IF (J1.EQ.0) THEN
              J1=J
            ELSE
              J2=J
            ENDIF
          ENDIF
   20   CONTINUE
C
        IF (ICOL.EQ.1) THEN
C  ONLY ONE COLUMN IS NONZERO
          DO 30 I=1,4
            IF (LCOLM(I)) FC(I)=B(J1)/P(J1,I)
   30     CONTINUE
C
        ELSE
C  MORE THAN ONE COLUMN IS NONZERO
          I1=0
          I2=0
C  DETERMINE THE INDICES FOR THE FIRST TWO NONZERO COLUMNS
          DO 40 I=1,4
            IF (LCOLM(I)) THEN
              IF (I1.EQ.0) THEN
                I1=I
              ELSE
                I2=I
                EXIT
              ENDIF
            ENDIF
   40     CONTINUE
C
! only columns to the right of column i2 can have non-zero entries
          DO I=I2+1,4
            B(J1)=B(J1)-P(J1,I)
            B(J2)=B(J2)-P(J2,I)
          END DO
C
          FNEN=P(J1,I1)*P(J2,I2)-P(J2,I1)*P(J1,I2)
          IF (ABS(FNEN).GT.EPS12) THEN
            FC(I1)=(B(J1)*P(J2,I2)-B(J2)*P(J1,I2))/FNEN
            IF (ABS(P(J1,I2)).GT.EPS12) THEN
              FC(I2)=(B(J1)-P(J1,I1)*FC(I1))/P(J1,I2)
            ELSEIF (ABS(P(J2,I2)).GT.EPS12) THEN
              FC(I2)=(B(J2)-P(J2,I1)*FC(I1))/P(J2,I2)
            ENDIF
          ENDIF
        ENDIF
C
C
      ELSEIF (IROW.EQ.3) THEN
C  THREE ROWS ARE NONZERO
C
C  DETERMINE THE ROW WHICH IS COMPLETELY ZERO
        JOUT = 0
        DO J=1,4
          IF (.NOT.LROW(J)) JOUT=J
        END DO

        IF (ICOL.EQ.1) THEN
C  ONLY ONE COLUMN IS NONZERO
          J1=1
          IF (JOUT.EQ.1) J1=2
          DO  I=1,4
            IF (LCOLM(I)) FC(I)=B(J1)/P(J1,I)
          END DO

        ELSEIF (ICOL.EQ.2) THEN
C  TWO COLUMNS ARE NONZERO
          J1=1
          IF (JOUT.EQ.J1) J1=J1+1
          J2=J1+1
          IF (JOUT.EQ.J2) J2=J2+1

C  DETERMINE THE INDICES FOR THE FIRST TWO NONZERO COLUMNS
          I1=0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IF (I1.EQ.0) THEN
                I1=I
              ELSE
                I2=I
                EXIT
              ENDIF
            ENDIF
          END DO

          FNEN=P(J1,I1)*P(J2,I2)-P(J2,I1)*P(J1,I2)
          IF (ABS(FNEN).GT.EPS12) THEN
            FC(I1)=(B(J1)*P(J2,I2)-B(J2)*P(J1,I2))/FNEN
            IF (ABS(P(J1,I2)).GT.EPS12) THEN
              FC(I2)=(B(J1)-P(J1,I1)*FC(I1))/P(J1,I2)
            ELSEIF (ABS(P(J2,I2)).GT.EPS12) THEN
              FC(I2)=(B(J2)-P(J2,I1)*FC(I1))/P(J2,I2)
            ENDIF
          ENDIF

        ELSE

C  AT LEAST THREE COLUMNS ARE NONZERO
          IF (ICOL.EQ.4) THEN
            DO J=1,4
              B(J)=B(J)-P(J,4)
            END DO
            IOUT=4
          ELSE
            IOUT=0
            DO I=1,4
              IF (.NOT.LCOLM(I)) IOUT=I
            END DO
          END IF

          J1=0
          DO J=1,4
            IF (J.EQ.JOUT) CYCLE
            J1=J1+1
            I1=0
            BB(J1)=B(J)
            DO I=1,4
              IF (I.EQ.IOUT) CYCLE
              I1=I1+1
              PP(J1,I1)=P(J,I)
            END DO
          END DO

cdr  We are in IROW=3 case, so we need to solve a 3x3 linear eq. system
cdr  We use the explicit Cramer Rule.
          P11=PP(1,1)
          P21=PP(2,1)
          P31=PP(3,1)
          P12=PP(1,2)
          P22=PP(2,2)
          P32=PP(3,2)
          P13=PP(1,3)
          P23=PP(2,3)
          P33=PP(3,3)
          B1=BB(1)
          B2=BB(2)
          B3=BB(3)
          dta=EIRENE_deter3x3(p11,p21,p31,
     .                        p12,p22,p32,
     .                        p13,p23,p33)
cdr  tbd:  check if determinant=0
cdr        see below, same as for IROW=4/ICOL=3 case

          dtb1=EIRENE_deter3x3(b1,b2,b3,
     .                         p12,p22,p32,
     .                         p13,p23,p33)
          dtb2=EIRENE_deter3x3(p11,p21,p31,
     .                         b1,b2,b3,
     .                         p13,p23,p33)
          dtb3=EIRENE_deter3x3(p11,p21,p31,
     .                         p12,p22,p32,
     .                         b1,b2,b3)
          Ffc(1)=dtb1/(dta+1.d-30)
          Ffc(2)=dtb2/(dta+1.d-30)
          Ffc(3)=dtb3/(dta+1.d-30)

          IC = 0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IC = IC + 1
              FC(I) = FFC(IC)
              IF (IC == 3) EXIT
            END IF
          END DO
        END IF
C
C
      ELSEIF (IROW.EQ.4) THEN
C  ALL ROWS ARE NONZERO
C
        IF (ICOL.EQ.1) THEN
          DO 50 I=1,4
            IF (LCOLM(I)) FC(I)=B(I)/P(1,I)
   50     CONTINUE
C
        ELSEIF (ICOL.EQ.2) THEN
          I1=0
          I2=0
C  DETERMINE THE INDICES FOR THE NONZERO COLUMNS
          DO 60 I=1,4
            IF (LCOLM(I)) THEN
              IF (I1.EQ.0) THEN
                I1=I
              ELSE
                I2=I
              ENDIF
            ENDIF
   60     CONTINUE
C
          FNEN=P(1,I1)*P(2,I2)-P(2,I1)*P(1,I2)
          IF (ABS(FNEN).GT.EPS12) THEN
            FC(I1)=(B(1)*P(2,I2)-B(2)*P(1,I2))/FNEN
            IF (ABS(P(1,I2)).GT.EPS12) THEN
              FC(I2)=(B(1)-P(1,I1)*FC(I1))/P(1,I2)
            ELSEIF (ABS(P(2,I2)).GT.EPS12) THEN
              FC(I2)=(B(2)-P(2,I1)*FC(I1))/P(2,I2)
            ENDIF
          ENDIF
C
        ELSEIF (ICOL.EQ.3) THEN

          IC=0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IC = IC + 1
              PP(1:3,IC) = P(1:3,I)
            END IF
          END DO

          P11=PP(1,1)
          P21=PP(2,1)
          P31=PP(3,1)
          P12=PP(1,2)
          P22=PP(2,2)
          P32=PP(3,2)
          P13=PP(1,3)
          P23=PP(2,3)
          P33=PP(3,3)
          B1=B(1)
          B2=B(2)
          B3=B(3)
          dta=EIRENE_deter3x3(p11,p21,p31,
     .                        p12,p22,p32,
     .                        p13,p23,p33)

cdr build a certain norm (L1) of the 3x3 matrix, to compare with determinant
          ap3ma=maxval(pp(1:3,1:3))
          ap3mi=minval(pp(1:3,1:3))
cdr largest absolute value of an element in matrix
          ap3m=max(abs(ap3ma),abs(ap3mi))
cdr 3x3 matrix, hence: **3, to compare with determinant
          ap3n3=ap3m**3
cdr check determinant=0, relative to L1 norm
          if (abs(dta)/ap3n3 > eps10) then
            dtb1=EIRENE_deter3x3(b1,b2,b3,
     .                           p12,p22,p32,
     .                           p13,p23,p33)
            dtb2=EIRENE_deter3x3(p11,p21,p31,
     .                           b1,b2,b3,
     .                           p13,p23,p33)
            dtb3=EIRENE_deter3x3(p11,p21,p31,
     .                           p12,p22,p32,
     .                           b1,b2,b3)
            Ffc(1)=dtb1/(dta+1.d-30)
            Ffc(2)=dtb2/(dta+1.d-30)
            Ffc(3)=dtb3/(dta+1.d-30)
          else
            write (iunout,*) 'SINGULAR MATRIX ENCOUNTERED IN GETSCL4'
            write (iunout,*) 'for ISTRA = ',ISTRA
            write (iunout,*) 'Possible reason: one type of particles'
            write (iunout,*) 'is created analytically (tracklength)'
            write (iunout,*) 'but not by random sampling'
cdr  proper norm for comparison with row or column norm squared
            ap3n2=ap3m**2
            ICOL2=0
            IROW2=0
            DO I=1,3
              LROW2(I)=
     .          ((PP(I,1)**2+PP(I,2)**2+PP(I,3)**2)/ap3n2.GT.EPS30)
              LCOL2(I)=
     .          ((PP(1,I)**2+PP(2,I)**2+PP(3,I)**2)/ap3n2.GT.EPS30)
              IF (LROW(I)) IROW=IROW+1
              IF (LCOL2(I)) ICOL=ICOL+1
            END DO
C  DETERMINE THE INDICES FOR THE FIRST TWO NONZERO COLUMNS,
C  and use a 2x2 matrix instead.
            I1=0
            I2=0
            DO I=1,3
              IF (LCOL2(I)) THEN
                IF (I1.EQ.0) THEN
                  I1=I
                ELSE
                  I2=I
                  EXIT
                ENDIF
              ENDIF
            END DO
            J1=0
            J2=0
            DO I=1,3
              IF (LROW2(I)) THEN
                IF (J1.EQ.0) THEN
                  J1=I
                  B1=B(I)
                ELSE
                  J2=I
                  B2=B(I)
                  EXIT
                ENDIF
              ENDIF
            END DO

            I3 = 6 - I1 - I2
            FFC(I3) = 1._DP

            DTA=PP(J1,I1)*PP(J2,I2)-PP(J2,I1)*PP(J1,I2)
cdr now we should check det=0 for the 2x2 determinant
cdr DTA should be compared, again, with a proper norm rather than with EPS12.
cdr build: ap2mi, ap2ma,...etc.
cdr tbd.
            IF (ABS(DTA).GT.EPS12) THEN
              FFC(I1)=(B1*PP(J2,I2)-B2*PP(J1,I2))/DTA
              IF (ABS(PP(J1,I2)).GT.EPS12) THEN
                FFC(I2)=(B1-PP(J1,I1)*FFC(I1))/PP(J1,I2)
              ELSEIF (ABS(PP(J2,I2)).GT.EPS12) THEN
                FFC(I2)=(B2-PP(J2,I1)*FFC(I1))/PP(J2,I2)
              ENDIF
cdr         else
cdr  ?
            ENDIF
          ENDIF

          IC = 0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IC = IC + 1
              FC(I) = FFC(IC)
            END IF
          END DO

        ELSE
C  THE ENTIRE 4 x 4 MATRIX IS TO BE USED
cdr unfinished, we should check for det (=dta) = 0. relative to
cdr a proper L1 norm of the matrix, and in that case reduce to the
cdr 3x3 matrix case, same as done above for the 3x3 matrix.
          pp(1:4,1:4) = p(1:4,1:4)
          dta=EIRENE_deter4x4(pp)

          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,1) = b(1:4)
          dtb1=EIRENE_deter4x4(pp)

          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,2) = b(1:4)
          dtb2=EIRENE_deter4x4(pp)

          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,3) = b(1:4)
          dtb3=EIRENE_deter4x4(pp)

          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,4) = b(1:4)
          dtb4=EIRENE_deter4x4(pp)

          fc(1)=dtb1/(dta+1.d-30)
          fc(2)=dtb2/(dta+1.d-30)
          fc(3)=dtb3/(dta+1.d-30)
          fc(4)=dtb4/(dta+1.d-30)
        ENDIF
      ENDIF
C
 1000 CONTINUE

!  FOR THE TIME BEING

      CALL EIRENE_LEER(1)
      WRITE (iunout,'(1X,A)')
     .  'EIRENE RECOMMENDED RESCALING OF VOLUME-AVERAGED TALLIES'
      WRITE (iunout,'(1X,A)')
     .  'DUE TO STATISTICAL ERRORS IN BALANCE'
      CALL EIRENE_MASR4 ('FATM,FMOL,FION,FPHOT            ',
     .             FC(1),FC(2),FC(3),FC(4))
      CALL EIRENE_LEER(2)

      IF (TRCSCL) THEN
        CALL EIRENE_MASAGE('GETSCL4                             ')
        CALL EIRENE_MASRR1('MATRIX    :',P,4*4,4)
        CALL EIRENE_MASRR1('RHS VECTOR:',B,4,4)
        CALL EIRENE_LEER(1)
      END IF
C
CNR IF ANY OF THESE FACTORS ARE NEGATIVE, KEEP AT 1.
CNR CAN HAPPEN IF BALANCES ARE BADLY BROKEN BY
CNR PARTICLES WITH ENORMOUS WEIGHTS...
      IF (ANY(FC < 0._DP)) THEN
        WRITE(*,*) '/!\ NEGATIVE RESCALING FOUND !',
     .    'BALANCES ARE BROKEN ! RESCALING RESET AT 1.'
        FA=1._DP
        FM=1._DP
        FI=1._DP
        FPH=1._DP
      ELSE
        FA=FC(1)
        FM=FC(2)
        FI=FC(3)
        FPH=FC(4)
      END IF
C
      RETURN
      END SUBROUTINE EIRENE_GETSCL4


      SUBROUTINE EIRENE_GETSCLS (ISTRA,NATMI,NMOLI,NIONI,NPHOTI,
     .                           FATM,FMOL,FION,FPHOT)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
     , ,ONLY : NATM, NMOL, NION, NPHOT, NATMP
      USE EIRMOD_COMSOU
      USE EIRMOD_COUTAU
      USE EIRMOD_CCOUPL
      USE EIRMOD_COMPRT
     , ,ONLY : IUNOUT
      USE EIRMOD_CLOGAU
     , ,ONLY : NLSPCSCL, NLSPCSCL_ATM, NLSPCSCL_MOL, NLSPCSCL_ION,
     ,         NLSPCSCL_PHOT, NLSPCSCL_ON
      USE EIRMOD_CESTIM
     , ,ONLY : LB_ATM, LB_MOL, LB_ION, LB_PHOT
      USE EIRMOD_CTRCEI
     , ,ONLY : TRCSCL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ISTRA, NATMI, NMOLI, NIONI, NPHOTI
      REAL(DP) :: FATM(0:NATMI), FMOL(0:NMOLI), FION(0:NIONI),
     .           FPHOT(0:NPHOTI)
      INTEGER :: NNP, IATM, IMOL, IION, IPHOT,
     .           JSP, JMOL, JION, JPHOT
      INTEGER :: IER
      INTEGER, ALLOCATABLE :: IW(:)
      REAL(DP), ALLOCATABLE :: PP(:,:), B(:)

      IF (NATMI.LE.1.AND.NMOLI.LE.1.AND.
     .    NIONI.LE.1.AND.NPHOTI.LE.1) THEN
        IF (NATMI .EQ.1) FATM(1) = FATM(0)
        IF (NMOLI .EQ.1) FMOL(1) = FMOL(0)
        IF (NIONI .EQ.1) FION(1) = FION(0)
        IF (NPHOTI.EQ.1) FPHOT(1)= FPHOT(0)
      ELSE IF (.NOT.NLSPCSCL) THEN
        IF (NATMI .GE.1) FATM (1:NATMI) = FATM(0)
        IF (NMOLI .GE.1) FMOL (1:NMOLI) = FMOL(0)
        IF (NIONI .GE.1) FION (1:NIONI) = FION(0)
        IF (NPHOTI.GE.1) FPHOT(1:NPHOTI)= FPHOT(0)
      ELSE
        NNP = NATMI+NMOLI+NIONI+NPHOTI
        ALLOCATE(PP(NNP,NNP),B(NNP))
        PP = 0._DP
        B = 0._DP

        PAATI2(0:NATM,LB_ATM:NATM) => PAATI(:,ISTRA)
        PAMLI2(0:NMOL,LB_ATM:NATM) => PAMLI(:,ISTRA)
        PAIOI2(0:NION,LB_ATM:NATM) => PAIOI(:,ISTRA)
        PAPHTI2(0:NPHOT,LB_ATM:NATM) => PAPHTI(:,ISTRA)
        PMATI2(0:NATM,LB_MOL:NMOL) => PMATI(:,ISTRA)
        PMMLI2(0:NMOL,LB_MOL:NMOL) => PMMLI(:,ISTRA)
        PMIOI2(0:NION,LB_MOL:NMOL) => PMIOI(:,ISTRA)
        PMPHTI2(0:NPHOT,LB_MOL:NMOL) => PMPHTI(:,ISTRA)
        PIATI2(0:NATM,LB_ION:NION) => PIATI(:,ISTRA)
        PIMLI2(0:NMOL,LB_ION:NION) => PIMLI(:,ISTRA)
        PIIOI2(0:NION,LB_ION:NION) => PIIOI(:,ISTRA)
        PIPHTI2(0:NPHOT,LB_ION:NION) => PIPHTI(:,ISTRA)
        PPHATI2(0:NATM,LB_PHOT:NPHOT) => PPHATI(:,ISTRA)
        PPHMLI2(0:NMOL,LB_PHOT:NPHOT) => PPHMLI(:,ISTRA)
        PPHIOI2(0:NION,LB_PHOT:NPHOT) => PPHIOI(:,ISTRA)
        PPHPHTI2(0:NPHOT,LB_PHOT:NPHOT) => PPHPHTI(:,ISTRA)

        IF (TRCSCL) THEN
          write (iunout,*) 'paati in getscl'
          do iatm=0,natmi
            write (iunout,'(A,i3,(1x,5es15.7))') 'iatm = ',
     .           iatm,paati2(iatm,lb_atm:natmi)
          end do
          write (iunout,*) 'pamli in getscl'
          do imol=0,nmoli
            write (iunout,'(A,i3,(1x,5es15.7))') 'imol = ',
     .           imol,pamli2(imol,lb_atm:natmi)
          end do
          write (iunout,*) 'paioi in getscl'
          do iion=0,nioni
            write (iunout,'(A,i3,(1x,5es15.7))') 'iion = ',
     .           iion,paioi2(iion,lb_atm:natmi)
          end do
          write (iunout,*) 'pmati in getscl'
          do iatm=0,natmi
            write (iunout,'(A,i3,(1x,5es15.7))') 'iatm = ',
     .           iatm,pmati2(iatm,lb_mol:nmoli)
          end do
          write (iunout,*) 'pmmli in getscl'
          do imol=0,nmoli
            write (iunout,'(A,i3,(1x,5es15.7))') 'imol = ',
     .           imol,pmmli2(imol,lb_mol:nmoli)
          end do
          write (iunout,*) 'pmioi in getscl'
          do iion=0,nioni
            write (iunout,'(A,i3,(1x,5es15.7))') 'iion = ',
     .           iion,pmioi2(iion,lb_mol:nmoli)
          end do
          write (iunout,*) 'piati in getscl'
          do iatm=0,natmi
            write (iunout,'(A,i3,(1x,5es15.7))') 'iatm = ',
     .           iatm,piati2(iatm,lb_ion:nioni)
          end do
          write (iunout,*) 'pimli in getscl'
          do imol=0,nmoli
            write (iunout,'(A,i3,(1x,5es15.7))') 'imol = ',
     .           imol,pimli2(imol,lb_ion:nioni)
          end do
          write (iunout,*) 'piioi in getscl'
          do iion=0,nioni
            write (iunout,'(A,i3,(1x,5es15.7))') 'iion = ',
     .           iion,piioi2(iion,lb_ion:nioni)
          end do
        END IF

        PRFAAI2(0:NATM,LB_ATM:NATM) => PRFAAI(:,ISTRA)
        PRFAMI2(0:NMOL,LB_ATM:NATM) => PRFAMI(:,ISTRA)
        PRFAII2(0:NION,LB_ATM:NATM) => PRFAII(:,ISTRA)
        PRFAPHTI2(0:NPHOT,LB_ATM:NATM) => PRFAPHTI(:,ISTRA)
        PRFMAI2(0:NATM,LB_MOL:NMOL) => PRFMAI(:,ISTRA)
        PRFMMI2(0:NMOL,LB_MOL:NMOL) => PRFMMI(:,ISTRA)
        PRFMII2(0:NION,LB_MOL:NMOL) => PRFMII(:,ISTRA)
        PRFMPHTI2(0:NPHOT,LB_MOL:NMOL) => PRFMPHTI(:,ISTRA)
        PRFIAI2(0:NATM,LB_ION:NION) => PRFIAI(:,ISTRA)
        PRFIMI2(0:NMOL,LB_ION:NION) => PRFIMI(:,ISTRA)
        PRFIII2(0:NION,LB_ION:NION) => PRFIII(:,ISTRA)
        PRFIPHTI2(0:NPHOT,LB_ION:NION) => PRFIPHTI(:,ISTRA)
        PRFPHAI2(0:NATM,LB_PHOT:NPHOT) => PRFPHAI(:,ISTRA)
        PRFPHMI2(0:NMOL,LB_PHOT:NPHOT) => PRFPHMI(:,ISTRA)
        PRFPHII2(0:NION,LB_PHOT:NPHOT) => PRFPHII(:,ISTRA)
        PRFPHPHTI2(0:NPHOT,LB_PHOT:NPHOT) => PRFPHPHTI(:,ISTRA)

        IF (TRCSCL) THEN
          write (iunout,*) 'prfaai in getscl'
          do iatm=0,natmi
            write (iunout,'(A,i3,(1x,5es15.7))') 'iatm = ',
     .           iatm,prfaai2(iatm,lb_atm:natmi)
          end do
          write (iunout,*) 'prfami in getscl'
          do imol=0,nmoli
            write (iunout,'(A,i3,(1x,5es15.7))') 'imol = ',
     .           imol,prfami2(imol,lb_atm:natmi)
          end do
          write (iunout,*) 'prfaii in getscl'
          do iion=0,nioni
            write (iunout,'(A,i3,(1x,5es15.7))') 'iion = ',
     .           iion,prfaii2(iion,lb_atm:natmi)
          end do
          write (iunout,*) 'prfmai in getscl'
          do iatm=0,natmi
            write (iunout,'(A,i3,(1x,5es15.7))') 'iatm = ',
     .           iatm,prfmai2(iatm,lb_mol:nmoli)
          end do
          write (iunout,*) 'prfmmi in getscl'
          do imol=0,nmoli
            write (iunout,'(A,i3,(1x,5es15.7))') 'imol = ',
     .           imol,prfmmi2(imol,lb_mol:nmoli)
          end do
          write (iunout,*) 'prfmii in getscl'
          do iion=0,nioni
            write (iunout,'(A,i3,(1x,5es15.7))') 'iion = ',
     .           iion,prfmii2(iion,lb_mol:nmoli)
          end do
          write (iunout,*) 'prfiai in getscl'
          do iatm=1,natmi
            write (iunout,'(A,i3,(1x,5es15.7))') 'iatm = ',
     .           iatm,prfiai2(iatm,1:nioni)
          end do
          write (iunout,*) 'prfimi in getscl'
          do imol=1,nmoli
            write (iunout,'(A,i3,(1x,5es15.7))') 'imol = ',
     .           imol,prfimi2(imol,1:nioni)
          end do
          write (iunout,*) 'prfiii in getscl'
          do iion=0,nioni
            write (iunout,'(A,i3,(1x,5es15.7))') 'iion = ',
     .           iion,prfiii2(iion,lb_ion:nioni)
          end do
        END IF

        DO JSP = 1, NATMI+NMOLI+NIONI+NPHOTI
          IF(JSP.LE.NATMI) THEN
            DO IATM = 1, NATMI
              PP(IATM,JSP) = PAATI2(IATM,JSP)+
     .                       PRFAAI2(IATM,JSP)
            END DO
            PP(JSP,JSP) = PP(JSP,JSP)+
     .                    POTATI(JSP,ISTRA)+PGENAI(JSP,ISTRA)
            DO IMOL = 1, NMOLI
              JMOL = NATMI+IMOL
              PP(JMOL,JSP) = PAMLI2(IMOL,JSP)+
     .                      PRFAMI2(IMOL,JSP)
            END DO
            DO IION = 1, NIONI
              JION = NATMI+NMOLI+IION
              PP(JION,JSP) = PAIOI2(IION,JSP)+
     .                      PRFAII2(IION,JSP)
            END DO
            DO IPHOT = 1, NPHOTI
              JPHOT = NATMI+NMOLI+NIONI+IPHOT
              PP(JPHOT,JSP) = PAPHTI2(IPHOT,JSP)+
     .                      PRFAPHTI2(IPHOT,JSP)
            END DO
          ELSE IF (JSP.GT.NATMI.AND.JSP.LE.NATMI+NMOLI) THEN
            JMOL = JSP-NATMI
            DO IATM = 1, NATMI
              PP(IATM,JSP) = PMATI2(IATM,JMOL)+
     .                      PRFMAI2(IATM,JMOL)
            END DO
            DO IMOL = 1, NMOLI
              JMOL = NATMI+IMOL
              PP(JMOL,JSP) = PMMLI2(IMOL,JSP-NATMI)+
     .                      PRFMMI2(IMOL,JSP-NATMI)
            END DO
            JMOL = JSP-NATMI
            PP(JSP,JSP) = PP(JSP,JSP)+
     .                    POTMLI(JMOL,ISTRA)+PGENMI(JMOL,ISTRA)
            DO IION = 1, NIONI
              JION = NATMI+NMOLI+IION
              PP(JION,JSP) = PMIOI2(IION,JMOL)+
     .                      PRFMII2(IION,JMOL)
            END DO
            DO IPHOT = 1, NPHOTI
              JPHOT = NATMI+NMOLI+NIONI+IPHOT
              PP(JPHOT,JSP) = PMPHTI2(IPHOT,JMOL)+
     .                      PRFMPHTI2(IPHOT,JMOL)
            END DO
          ELSE IF (JSP.GT.NATMI+NMOLI.AND.JSP.LE.NATMI+NMOLI+NIONI) THEN
            JION = JSP-NATMI-NMOLI
            DO IATM = 1, NATMI
              PP(IATM,JSP) = PIATI2(IATM,JION)+
     .                      PRFIAI2(IATM,JION)
            END DO
            DO IMOL = 1, NMOLI
              JMOL = NATMI+IMOL
              PP(JMOL,JSP) = PIMLI2(IMOL,JION)+
     .                      PRFIMI2(IMOL,JION)
            END DO
            DO IION = 1, NIONI
              JION = NATMI+NMOLI+IION
              PP(JION,JSP) = PIIOI2(IION,JSP-NATMI-NMOLI)+
     .                      PRFIII2(IION,JSP-NATMI-NMOLI)
            END DO
            JION = JSP-NATMI-NMOLI
            PP(JSP,JSP) = PP(JSP,JSP)+
     .                    POTIOI(JION,ISTRA)+PGENII(JION,ISTRA)
            DO IPHOT = 1, NPHOTI
              JPHOT = NATMI+NMOLI+NIONI+IPHOT
              PP(JPHOT,JSP) = PIPHTI2(IPHOT,JION)+
     .                      PRFIPHTI2(IPHOT,JION)
            END DO
          ELSE IF (JSP.GT.NATMI+NMOLI+NIONI.AND.
     .             JSP.LE.NATMI+NMOLI+NIONI+NPHOTI) THEN
            JPHOT = JSP-NATMI-NMOLI-NIONI
            DO IATM = 1, NATMI
              PP(IATM,JSP) = PPHATI2(IATM,JPHOT)+
     .                      PRFPHAI2(IATM,JPHOT)
            END DO
            DO IMOL = 1, NMOLI
              JMOL = NATMI+IMOL
              PP(JMOL,JSP) = PPHMLI2(IMOL,JPHOT)+
     .                      PRFPHMI2(IMOL,JPHOT)
            END DO
            DO IION = 1, NIONI
              JION = NATMI+NMOLI+IION
              PP(JION,JSP) = PPHIOI2(IION,JPHOT)+
     .                      PRFPHII2(IION,JPHOT)
            END DO
            DO IPHOT = 1, NPHOTI
              JPHOT = NATMI+NMOLI+NIONI+IPHOT
              PP(JPHOT,JSP) =
     .            PPHPHTI2(IPHOT,JSP-NATMI-NMOLI-NIONI)+
     .          PRFPHPHTI2(IPHOT,JSP-NATMI-NMOLI-NIONI)
            END DO
            JPHOT = JSP-NATMI-NMOLI-NIONI
            PP(JSP,JSP) = PP(JSP,JSP)+
     .                    POTPHTI(JPHOT,ISTRA)+PGENPHI(JPHOT,ISTRA)
          END IF
        END DO
        DO IATM = 1, NATMI
          B(IATM) =-(PPATI(IATM,ISTRA)+WTOTA(IATM,ISTRA))
        END DO
        DO IMOL = 1, NMOLI
          JMOL = NATMI+IMOL
          B(JMOL) =-(PPMLI(IMOL,ISTRA)+WTOTM(IMOL,ISTRA))
        END DO
        DO IION = 1, NIONI
          JION = NATMI+NMOLI+IION
          B(JION) =-(PPIOI(IION,ISTRA)+WTOTI(IION,ISTRA))
        END DO
        DO IPHOT = 1, NPHOTI
          JPHOT = NATMI+NMOLI+NIONI+IPHOT
          B(JPHOT) =-(PPPHTI(IPHOT,ISTRA)+WTOTPH(IPHOT,ISTRA))
        END DO
        DO JSP = 1, NATMI+NMOLI+NIONI+NPHOTI
          IF (SUM(PP(1:NNP,JSP)).EQ.0._DP.AND.
     .        SUM(PP(JSP,1:NNP)).EQ.0._DP.AND.B(JSP).EQ.0._DP) THEN
            PP(JSP,JSP) = 1._DP
            B(JSP) = 1._DP
          END IF
        END DO
C
        IF (TRCSCL) THEN
          CALL EIRENE_MASRR1('MATRIX    :',PP,NNP*NNP,NNP)
          CALL EIRENE_MASRR1('RHS VECTOR:',B,NNP,NNP)
          CALL EIRENE_MASRR1('PAATI2    :',PAATI(:,ISTRA),
     .                                     NATMP*NATMI,NATMI)
        END IF
!! Solve the matrix
        IER=0
        ALLOCATE(IW(NNP))
        CALL EIRENE_GALPD_M(PP,NNP,NNP,B,1,1,IW,IER)
        IF(IER.EQ.0) THEN
          IF(MINVAL(B).GT.0._DP) THEN
            IF(NATMI.GE.1) FATM(1:NATMI) = B(1:NATMI)
            IF(NMOLI.GE.1) FMOL(1:NMOLI) = B(NATMI+1:NATMI+NMOLI)
            IF(NIONI.GE.1) FION(1:NIONI) = B(NATMI+NMOLI+1:
     .                                       NATMI+NMOLI+NIONI)
            IF(NPHOTI.GE.1)FPHOT(1:NPHOTI)=B(NATMI+NMOLI+NIONI+1:
     .                                       NATMI+NMOLI+NIONI+NPHOTI)
          ELSE
            CALL EIRENE_MASAGE('MATRIX SOLVE YIELDED NEGATIVE VALUES')
            CALL EIRENE_MASAGE('REVERTING BACK TO DEFAULT BEHAVIOUR')
            IF (NATMI .GE.1) FATM (1:NATMI) = FATM(0)
            IF (NMOLI .GE.1) FMOL (1:NMOLI) = FMOL(0)
            IF (NIONI .GE.1) FION (1:NIONI) = FION(0)
            IF (NPHOTI.GE.1) FPHOT(1:NPHOTI)= FPHOT(0)
          END IF
        ELSE IF (IER.GE.1) THEN
          CALL EIRENE_MASAGE('RESCALING MATRIX IS SINGULAR')
          CALL EIRENE_MASAGE('REVERTING BACK TO DEFAULT BEHAVIOUR')
          IF (NATMI .GE.1) FATM (1:NATMI) = FATM(0)
          IF (NMOLI .GE.1) FMOL (1:NMOLI) = FMOL(0)
          IF (NIONI .GE.1) FION (1:NIONI) = FION(0)
          IF (NPHOTI.GE.1) FPHOT(1:NPHOTI)= FPHOT(0)
        ENDIF
        DEALLOCATE(IW)
        DEALLOCATE(PP,B)
      END IF

      IF (.NOT.NLSPCSCL_ON) THEN
        FATM (1:NATMI) = FATM(0)
        FMOL (1:NMOLI) = FMOL(0)
        FION (1:NIONI) = FION(0)
        FPHOT(1:NPHOTI)= FPHOT(0)
      END IF

      RETURN
      END SUBROUTINE EIRENE_GETSCLS
