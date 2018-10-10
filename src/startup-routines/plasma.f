C
C
!  6.12.05  bugfix: avoid calculation of B-field in dead cells
!                   because geometrical parameters may not be known there.
!  6.8. 06  bugfix of bugfix: avoid calculation of B-field in dead cells, but
!                             still make sure to set B-field in 1D cases.
!  15.12.06 bug fix: index error corrected in call to prousr when called for ADIN
!  10.06.08 new:  default BFIN=1 T, rather than 0 T
!  10.06.08 new option: profile type 3 (profs): set BFIN using B2 and B3 parameters
!  22.09.14 bug fix re. this ind=3 option in case of type (=ind) = 1,2 .
!                       help2 was undefined --> zero b-field
!
cdr try to re-unify treatment of 1st dimension (species index) in parameters
cdr n,T,V for background (bulk) velocity distribution: not finished.
!  sept. 16 change variable names ipls --> iplsti, (for TI)
!                                 ipls --> iplsv,  (for VX,VY,VZ)
!  oct. 16  comments, one minor bug fix (VZIN(IPLSV) in one (unused) option)
!  nov. 16  nlpitch option added, for orientation of B-field in 1D and 2D runs

cpb: add parameter ndim: special treatment of Ti fields species index.
cpb: reading tiin from profr:  set 1st dimension of tiin array.

cdr: check under which conditions can nplsti be different from npls, and is that still needed?
cdr: why is that not needed for V and n profiles?
!
!
      SUBROUTINE EIRENE_PLASMA
C  SET DENSITY, TEMPERATURE AND MACH NUMBER PROFILES, B AND E FIELDS, 
C  ON: 
C  INDPRO=1,2,3    1D MESH "RHOZNE(J)", 1,NR1STM, CELL CENTERED
C                  B-FIELD (INDPRO(5)) SET ON 1:NSURF 
C  INDPRO=4        READ FROM EXTERNAL FILE ISTREAM, EVERYWHERE, 1,NSBOX, 
C  INDPRO=5        PROUSR: ONLY IN STANDARD GRID, 1:NSURF
C  INDPRO=6        PROFR : ONLY IN STANDARD GRID, 1:NSURF
C  INDPRO=7        PROFR : EVERYWHERE, 1,NSBOX
C  INDPRO=8        ??
C  INDPRO=9        INPUT TALLIES ARE ALREADY SET ELSEWHERE, 1:NSURF ??
C
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CTEXT
      USE EIRMOD_COMPRT
 
      IMPLICIT NONE
 
      REAL(DP), ALLOCATABLE :: HELP(:), HELP2(:)
      REAL(DP) :: PUX, PUY, EL, EP, PN, BD, B, BVAC, FACT
      INTEGER :: IB, IAIN, K, JJ, ITALI, ICELL, IND, ISTREAM,
     .           IP, IT, IA, J, IR, IPLSTI, IPLSV, NDIM
C
C  INDPRO=9 MEANS: THESE ARRAYS ARE ALREADY SET IN COUPLE_... (SUBR. INFCOP)
      IF (INDPRO(1) /= 9) TEIN = 0.D0
      IF (INDPRO(2) /= 9) TIIN = 0.D0
                          DEIN = 0.D0 ! DERIVED TALLY DEIN IS SET IN PLASMA_DERIV
      IF (INDPRO(3) /= 9) DIIN = 0.D0
      IF (INDPRO(4) /= 9) VXIN = 0.D0
      IF (INDPRO(4) /= 9) VYIN = 0.D0
      IF (INDPRO(4) /= 9) VZIN = 0.D0
c  magnetic field
      IF (LBXIN .AND. (INDPRO(5) /= 9)) BXIN = 0.D0
      IF (LBYIN .AND. (INDPRO(5) /= 9)) BYIN = 0.D0
      IF (LBZIN .AND. (INDPRO(5) /= 9)) BZIN = 0.D0
      IF (LBFIN .AND. (INDPRO(5) /= 9)) BFIN = 0.D0

      IF (LADIN .AND. (INDPRO(6) /= 9)) ADIN = 0.D0
c  electric field
      IF (LEXIN .AND. (INDPRO(7) /= 9)) EXIN = 0.D0
      IF (LEYIN .AND. (INDPRO(7) /= 9)) EYIN = 0.D0
      IF (LEZIN .AND. (INDPRO(7) /= 9)) EZIN = 0.D0
      IF (LEFIN .AND. (INDPRO(7) /= 9)) EFIN = 0.D0
      IF (LPOT .AND.  (INDPRO(7) /= 9)) POT  = 0.D0

      ALLOCATE (HELP(NRAD))
      ALLOCATE (HELP2(NRAD))
      HELP=0.
      HELP2=0.
C
C  SET EIRENE VACUUM BACKGROUND MODEL DATA. I.E. IF TEMPERATURES ARE
C  LESS THAN TVAC OR THE BACKGROUND DENSITY IS LESS THAN DVAC,
C  THEN THIS ZONE IS CONSIDERED TO BE AN "EIRENE VACUUM ZONE",
C  FOR A PARTICULAR BACKGROUND SPECIES:
C  PARTICLE MEAN FREE PATHES IN SUCH ZONES           ARE SET EQUAL TO 1.D10 (CM)
C  AND ALL REACTION RATES WRT: TO THIS BULK PARTICLE ARE SET EQUAL TO ZERO (1/S)
      TVAC=0.02_dp
      DVAC=1.E2_dp
      VVAC=0._dp
      BVAC=1._dp  ! dr:  B field must not be "vacuum". check use of BVAC
C

C
C  ELECTRON TEMPERATURE
      IND=INDPRO(1)
      GOTO (101,102,103,104,105,106,107,110,110),IND
101     CALL EIRENE_PROFN (TEIN,TE0,TE1,TE2,TE3,TE4,TE5,TVAC)
        GOTO 110
102     CALL EIRENE_PROFE (TEIN,TE0,TE1,TE2,TE4,TE5,TVAC)
        GOTO 110
103     CALL EIRENE_PROFS (TEIN,TE0,TE1,TE5,TVAC)
        GOTO 110
104     CONTINUE
c  INDPRO=4:  read tally from stream TEO
        ISTREAM=TE0
        ITALI=1
        CALL EIRENE_READTL(TXTPLS(1,ITALI),TXTPSP(1,ITALI),
     .              TXTPUN(1,ITALI),
     .              TEIN,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .              3,ISTREAM)
        GOTO 110
c  INDPRO=5:  tally from PROUSR, indx=0
105     CALL EIRENE_PROUSR (TEIN,0,TE0,TE1,TE2,TE3,TE4,TE5,TVAC,NSURF)
        GOTO 110

c  INDPRO=6:  tally from PROFR,  1:NSURF
106     CALL EIRENE_PROFR (TEIN,0,1,1,NSURF)
        GOTO 110
c  INDPRO=7:  tally from PROFR,  1:NSBOX=NSURF+NRADD
107     CALL EIRENE_PROFR (TEIN,0,1,1,NSBOX)
        GOTO 110
110   CONTINUE
 
 
 
C  ION TEMPERATURE
      IND=INDPRO(2)
      DO 120 IPLSTI=1,NPLSTI
cdr one profile iplsti set at a time
        GOTO (111,112,113,114,115,116,117,120,120),IND
111     CALL EIRENE_PROFN (HELP,TI0(IPLSTI),TI1(IPLSTI),TI2(IPLSTI),
     .                     TI3(IPLSTI),TI4(IPLSTI),TI5(IPLSTI),TVAC)
        TIIN(IPLSTI,1:NR1ST)=HELP(1:NR1ST)
        GOTO 120
112     CALL EIRENE_PROFE (HELP,TI0(IPLSTI),TI1(IPLSTI),TI2(IPLSTI),
     .                                 TI4(IPLSTI),TI5(IPLSTI),TVAC)
        TIIN(IPLSTI,1:NR1ST)=HELP(1:NR1ST)
        GOTO 120
113     CALL EIRENE_PROFS (HELP,TI0(IPLSTI),TI1(IPLSTI),
     .                     TI5(IPLSTI),TVAC)
        TIIN(IPLSTI,1:NR1ST)=HELP(1:NR1ST)
        GOTO 120
114     CONTINUE
c  INDPRO=4:  read tally from stream TIO(IPLSTI)
        ISTREAM=TI0(IPLSTI)
        ITALI=2
        CALL EIRENE_READTL(TXTPLS(IPLSTI,ITALI),TXTPSP(IPLSTI,ITALI),
     .              TXTPUN(IPLSTI,ITALI),
     .              HELP,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .              3,ISTREAM)
        TIIN(IPLSTI,1:NSBOX)=HELP(1:NSBOX)
        GOTO 120
c  INDPRO=5:  tally from PROUSR, indx=1, but NPLSTI calls, one for each IPLSTI
115     CALL EIRENE_PROUSR (HELP,1+0*NPLS,TI0(IPLSTI),TI1(IPLSTI),
     .                      TI2(IPLSTI),TI3(IPLSTI),TI4(IPLSTI),
     .                      TI5(IPLSTI),TVAC,NSURF)
        TIIN(IPLSTI,1:NSURF)=HELP(1:NSURF)
        GOTO 120
120   CONTINUE
      GOTO 1120

cdr all nplsti profiles set in a single call

c  INDPRO=6:  tally from PROFR, indx=1, all TIIN fields in one single call
cdr first dimension of arrays:  NDIM .ne NPLSTI possible ?
!pb Jan 17: 116     CALL EIRENE_PROFR (TIIN,1+0*NPLS,NPLSTI,NPLSTI,NSURF)
116   NDIM = SIZE(TIIN,DIM=1)
      CALL EIRENE_PROFR (TIIN,1+0*NPLS,NPLSTI,NDIM,NSURF)
      GOTO 1120
c  INDPRO=7:  tally from PROFR, indx=1, all TIIN fields in one single call
cdr first dimension of arrays:  NDIM .ne NPLSTI possible ?
!pb Jan 17 117     CALL EIRENE_PROFR (TIIN,1+0*NPLS,NPLSTI,NPLSTI,NSBOX)
117   NDIM = SIZE(TIIN,DIM=1)
      CALL EIRENE_PROFR (TIIN,1+0*NPLS,NPLSTI,NDIM,NSBOX)
      GOTO 1120
1120  CONTINUE
 
 
 
C  ION DENSITY
      IND=INDPRO(3)
      DO 130 IPLS=1,NPLSI
        IF (LEN_TRIM(CDENMODEL(IPLS)) > 0) CYCLE
        GOTO (121,122,123,124,125,126,127,130,130),IND
cdr one profile ipls set at a time
121     CALL EIRENE_PROFN (HELP,DI0(IPLS),DI1(IPLS),DI2(IPLS),
     .                   DI3(IPLS),DI4(IPLS),DI5(IPLS),DVAC)
        DIIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 130
122     CALL EIRENE_PROFE (HELP,DI0(IPLS),DI1(IPLS),DI2(IPLS),
     .                             DI4(IPLS),DI5(IPLS),DVAC)
        DIIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 130
123     CALL EIRENE_PROFS (HELP,DI0(IPLS),DI1(IPLS),DI5(IPLS),DVAC)
        DIIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 130
124     CONTINUE
c  INDPRO=4:  read tally from stream DIO(IPLS)
        ISTREAM=DI0(IPLS)
        ITALI=4
        CALL EIRENE_READTL(TXTPLS(IPLS,ITALI),TXTPSP(IPLS,ITALI),
     .              TXTPUN(IPLS,ITALI),
     .              HELP,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .              3,ISTREAM)
        DIIN(IPLS,1:NSBOX)=HELP(1:NSBOX)
        GOTO 130
c  INDPRO=5:  tally from PROUSR, indx=1+1*NPLS, but NPLSI calls, one for each IPLS
125     CALL EIRENE_PROUSR (HELP,1+1*NPLS,DI0(IPLS),DI1(IPLS),DI2(IPLS),
     .                    DI3(IPLS),DI4(IPLS),DI5(IPLS),DVAC,NSURF)
        DIIN(IPLS,1:NSURF)=HELP(1:NSURF)
        GOTO 130
130   CONTINUE
      GOTO 1130

cdr all nplsi profiles set in a single call

c  INDPRO=6:
cdr first dimension of arrays:  always NPLS
126   CALL EIRENE_PROFR (DIIN,1+0*NPLS+NPLSTI,NPLSI,NPLS,NSURF)
      GOTO 1130
c  INDPRO=7:
cdr first dimension of arrays:  always NPLS
127   CALL EIRENE_PROFR (DIIN,1+0*NPLS+NPLSTI,NPLSI,NPLS,NSBOX)
      GOTO 1130
1130  CONTINUE
 
 
 
C  DRIFT VELOCITY
      IND=INDPRO(4)
      DO 140 IPLSV=1,NPLSV
        GOTO (131,132,133,134,135,136,137,140,140),IND
cdr one vector component profile (vx,vy,vz) iplsv set at a time
131     CALL EIRENE_PROFN (HELP,VX0(IPLSV),VX1(IPLSV),VX2(IPLSV),
     .                   VX3(IPLSV),VX4(IPLSV),VX5(IPLSV),VVAC)
        VXIN(IPLSV,1:NR1ST)=HELP(1:NR1ST)
        CALL EIRENE_PROFN (HELP,VY0(IPLSV),VY1(IPLSV),VY2(IPLSV),
     .                   VY3(IPLSV),VY4(IPLSV),VY5(IPLSV),VVAC)
        VYIN(IPLSV,1:NR1ST)=HELP(1:NR1ST)
        CALL EIRENE_PROFN (HELP,VZ0(IPLSV),VZ1(IPLSV),VZ2(IPLSV),
     .                   VZ3(IPLSV),VZ4(IPLSV),VZ5(IPLSV),VVAC)
        VZIN(IPLSV,1:NR1ST)=HELP(1:NR1ST)
        GOTO 140
132     CALL EIRENE_PROFE (HELP,VX0(IPLSV),VX1(IPLSV),VX2(IPLSV),
     .                   VX4(IPLSV),VX5(IPLSV),VVAC)
        VXIN(IPLSV,1:NR1ST)=HELP(1:NR1ST)
        CALL EIRENE_PROFE (HELP,VY0(IPLSV),VY1(IPLSV),VY2(IPLSV),
     .                   VY4(IPLSV),VY5(IPLSV),VVAC)
        VYIN(IPLSV,1:NR1ST)=HELP(1:NR1ST)
        CALL EIRENE_PROFE (HELP,VZ0(IPLSV),VZ1(IPLSV),VZ2(IPLSV),
     .                   VZ4(IPLSV),VZ5(IPLSV),VVAC)
        VZIN(IPLSV,1:NR1ST)=HELP(1:NR1ST)
        GOTO 140
133     CALL EIRENE_PROFS (HELP,VX0(IPLSV),VX1(IPLSV),VX5(IPLSV),VVAC)
        VXIN(IPLSV,1:NR1ST)=HELP(1:NR1ST)
        CALL EIRENE_PROFS (HELP,VY0(IPLSV),VY1(IPLSV),VY5(IPLSV),VVAC)
        VYIN(IPLSV,1:NR1ST)=HELP(1:NR1ST)
        CALL EIRENE_PROFS (HELP,VZ0(IPLSV),VZ1(IPLSV),VZ5(IPLSV),VVAC)
        VZIN(IPLSV,1:NR1ST)=HELP(1:NR1ST)
        GOTO 140
134     CONTINUE
C  INDPRO=4: read tallies from streams VXO(IPLSV), VY0(IPLSV),VZ0(IPLSV)
c            NOT READY FOR FLOW FIELDS
        GOTO 140
c  INDPRO=5:  tally from PROUSR,
c             VX: indx=1+2*NPLS, but NPLSV calls, one for each IPLSV
c             VY: indx=1+3*NPLS, but NPLSV calls, one for each IPLSV
c             VZ: indx=1+4*NPLS, but NPLSV calls, one for each IPLSV
135     CALL EIRENE_PROUSR (HELP,1+2*NPLS,VX0(IPLSV),VX1(IPLSV),
     .                      VX2(IPLSV),VX3(IPLSV),
     .                      VX4(IPLSV),VX5(IPLSV),VVAC,NSURF)
        VXIN(IPLSV,1:NSURF)=HELP(1:NSURF)
        CALL EIRENE_PROUSR (HELP,1+3*NPLS,VY0(IPLSV),VY1(IPLSV),
     .                      VY2(IPLSV),VY3(IPLSV),
     .                      VY4(IPLSV),VY5(IPLSV),VVAC,NSURF)
        VYIN(IPLSV,1:NSURF)=HELP(1:NSURF)
        CALL EIRENE_PROUSR (HELP,1+4*NPLS,VZ0(IPLSV),VZ1(IPLSV),
     .                      VZ2(IPLSV),VZ3(IPLSV),
     .                      VZ4(IPLSV),VZ5(IPLSV),VVAC,NSURF)
        VZIN(IPLSV,1:NSURF)=HELP(1:NSURF)
        GOTO 140
140   CONTINUE

C  SCALE FROM MACH NUMBER PROFILE TO CM/SEC PROFILE?
C  USE ISOTHERMAL ACCOUSTIC SPEED OF ION IPLS.
      IF (NLMACH) THEN
        DO 1141 IPLS=1,NPLSI
          IPLSTI=MPLSTI(IPLS)
          IPLSV=MPLSV(IPLS)
          DO 1142 ICELL=1,NSURF
            FACT=CVEL2A*SQRT((TIIN(IPLSTI,ICELL)+
     .                        TEIN(ICELL))/RMASSP(IPLS))
            VXIN(IPLSV,ICELL)=VXIN(IPLSV,ICELL)*FACT
            VYIN(IPLSV,ICELL)=VYIN(IPLSV,ICELL)*FACT
            VZIN(IPLSV,ICELL)=VZIN(IPLSV,ICELL)*FACT
1142      CONTINUE
1141    CONTINUE
      ENDIF
      GOTO 1140


cdr all nplsv vector component profiles set in a single call

c  read tally from external data structure, all V.IN fields in one single call
cdr first dimension of arrays:  always NPLSV
136   CALL EIRENE_PROFR (VXIN,1+1*NPLS+NPLSTI+0*NPLSV,NPLSV,NPLSV,NSURF)
      CALL EIRENE_PROFR (VYIN,1+1*NPLS+NPLSTI+1*NPLSV,NPLSV,NPLSV,NSURF)
      CALL EIRENE_PROFR (VZIN,1+1*NPLS+NPLSTI+2*NPLSV,NPLSV,NPLSV,NSURF)
      GOTO 1140

c  read tally from external data structure, all V.IN fields in one single call
cdr first dimension of arrays:  always NPLSV
137   CALL EIRENE_PROFR (VXIN,1+1*NPLS+NPLSTI+0*NPLSV,NPLSV,NPLSV,NSBOX)
      CALL EIRENE_PROFR (VYIN,1+1*NPLS+NPLSTI+1*NPLSV,NPLSV,NPLSV,NSBOX)
      CALL EIRENE_PROFR (VZIN,1+1*NPLS+NPLSTI+2*NPLSV,NPLSV,NPLSV,NSBOX)
      GOTO 1140
1140  CONTINUE
C
C
C  MAGNETIC FIELD UNIT VECTOR
C
      IF (.NOT.(LBXIN.AND.LBYIN.AND.LBZIN.AND.LBFIN)) GOTO 154

C  FOR IND=4,5,6,7 OR 9: ALSO THE ABSOLUTE B-FIELD STRENGTH BF CAN BE SET
      IND=INDPRO(5)
C  DEFAULT: 1 TESLA BFIELD IN Z-DIRECTION, IE., PITCH=0
      IF (IND /= 9) THEN
        BXIN=0.
        BYIN=0.
        BZIN=1.
        BFIN=1.
      END IF

      GOTO (141,142,143,144,145,146,147,150,150),IND
C  HELP IS FIELD LINE PITCH ANGLE: B_POL/B_TOT

C  INDPRO(5)=1:
141     CALL EIRENE_PROFN (HELP,B0,B1,B2,B3,B4,B5,BVAC)
        GOTO 1400
C  INDPRO(5)=2:
142     CALL EIRENE_PROFE (HELP,B0,B1,B2,B4,B5,BVAC)
        GOTO 1400
C  INDPRO(5)=3:
143     CALL EIRENE_PROFS (HELP,B0,B1,B5,BVAC)
        CALL EIRENE_PROFS (HELP2,B2,B3,B5,BVAC) ! new (2008) set constant B profile, ONLY INDPRO(5)=3
        GOTO 1400
C  INDPRO(5)=4: read tally from stream B0:  NOT IN USE
144     CONTINUE
        GOTO 150

C  CONVERT PITCH ANGLE INTO B-FIELD UNIT VECTOR
1400    CONTINUE
C  AT THIS POINT: INDPRO= 1,2, OR =3. HELP2 IS KNOWN ONLY IN CASE INDPRO=3
        IF (LEVGEO.EQ.1) THEN
          DO 1401 J=1,NSURF
            CALL EIRENE_NCELLN(J,IR,IP,IT,IA,IB,
     .                  NR1ST,NP2ND,NT3RD,NBMLT,NLRAD,NLPOL,NLTOR)
            IF (IR.GE.NR1ST) GOTO 1401
            IF ((NP2ND.GT.1).AND.(IP.GE.NP2ND)) GOTO 1401
C
            IF (.NOT.NLPITCH) THEN ! OLD DEFAULT: B-FIELD IS parallel TO Y,Z
              BXIN(J)=0.0
              BYIN(J)=HELP(IR)
            ELSEIF (NLPITCH) THEN  ! NEW OPTION : B-FIELD IS parallel TO X,Z
              BXIN(J)=HELP(IR)
              BYIN(J)=0.0
            ENDIF
C
            BZIN(J)=SQRT(1.-HELP(IR)*HELP(IR))
C
            IF (IND.EQ.3) THEN
              BFIN(J)=HELP2(IR)
            ELSE
              BFIN(J)=1.
            ENDIF
1401      CONTINUE
        ELSEIF (LEVGEO.EQ.2.AND.NLPOL) THEN
          DO 1402 J=1,NSURF
            CALL EIRENE_NCELLN(J,IR,IP,IT,IA,IB,
     .                  NR1ST,NP2ND,NT3RD,NBMLT,NLRAD,NLPOL,NLTOR)
            IF (IR.GE.NR1ST) GOTO 1402
            IF ((NP2ND.GT.1).AND.(IP.GE.NP2ND)) GOTO 1402
            EP=0.5*(EP1(IR)+EP1(IR+1))
            EL=0.5*(ELL(IR)+ELL(IR+1))
            JJ=IR+(IP-1)*NR1ST
            PUY= XCOM(JJ)-EP
            PUX=-YCOM(JJ)/EL/EL
            PN=SQRT(PUX*PUX+PUY*PUY+EPS60)
            BXIN(J)=HELP(IR)*PUX/PN
            BYIN(J)=HELP(IR)*PUY/PN
            BZIN(J)=SQRT(1.-HELP(IR)*HELP(IR))
            IF (IND.EQ.3) THEN
              BFIN(J)=HELP2(IR)
            ELSE
              BFIN(J)=1.
            ENDIF
1402      CONTINUE
        ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
          DO 1403 J=1,NSURF
            IF (NSTGRD(J) /= 0) CYCLE
            CALL EIRENE_NCELLN(J,IR,IP,IT,IA,IB,
     .                  NR1ST,NP2ND,NT3RD,NBMLT,NLRAD,NLPOL,NLTOR)
            IF (IR.GE.NR1ST) GOTO 1403
            IF ((NP2ND.GT.1).AND.(IP.GE.NP2ND)) GOTO 1403
            PUX=VPLX(IR,IP)+VPLX(IR+1,IP)
            PUY=VPLY(IR,IP)+VPLY(IR+1,IP)
            PN=SQRT(PUX*PUX+PUY*PUY+EPS60)
            BXIN(J)=HELP(IR)*PUX/PN
            BYIN(J)=HELP(IR)*PUY/PN
            BZIN(J)=SQRT(1.-HELP(IR)*HELP(IR))
            IF (IND.EQ.3) THEN
              BFIN(J)=HELP2(IR)
            ELSE
              BFIN(J)=1.
            ENDIF
1403      CONTINUE
        ELSE
          CALL EIRENE_LEER(1)
          WRITE (iunout,*)
     .      'NO MAGNETIC PITCH PROFILE COULD BE DEFINED FOR THIS CASE'
          WRITE (iunout,*)
     .      'DEFAULT MAGNETIC FIELD (IN Z-DIRECTION) IS USED'
          CALL EIRENE_LEER(1)
        ENDIF
        GOTO 150
c  INDPRO(5)=5:  call prousr
145     CONTINUE
        CALL EIRENE_PROUSR (BXIN,1+1*NPLS+NPLSTI+3*NPLSV,
     .                      B0,B1,B2,B3,B4,B5,0._DP,NSURF)
        CALL EIRENE_PROUSR (BYIN,2+1*NPLS+NPLSTI+3*NPLSV,
     .                      B0,B1,B2,B3,B4,B5,0._DP,NSURF)
        CALL EIRENE_PROUSR (BZIN,3+1*NPLS+NPLSTI+3*NPLSV,
     .                      B0,B1,B2,B3,B4,B5,1._DP,NSURF)
        CALL EIRENE_PROUSR (BFIN,4+1*NPLS+NPLSTI+3*NPLSV,
     .                      B0,B1,B2,B3,B4,B5,1._DP,NSURF)
        GOTO 150
c  INDPRO(5) =6:  call profr (information comes from interfacing code)
c                 default (vacuum) parameters in additional cells
146     CALL EIRENE_PROFR (BXIN,1+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        CALL EIRENE_PROFR (BYIN,2+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        CALL EIRENE_PROFR (BZIN,3+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        CALL EIRENE_PROFR (BFIN,4+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        GOTO 150
c  INDPRO(5) =7:  call profr (information comes from interfacing code)
c                 include also additional cells
147     CALL EIRENE_PROFR (BXIN,1+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        CALL EIRENE_PROFR (BYIN,2+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        CALL EIRENE_PROFR (BZIN,3+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        CALL EIRENE_PROFR (BFIN,4+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        GOTO 150
150   CONTINUE
C
C  CHECK FOR ZERO MAGNETIC FIELD IN ANY CELL (INCL. ADD. CELL REGION)
      DO 153 JJ=1,NSBOX
        IF (BXIN(JJ)**2+BYIN(JJ)**2+BZIN(JJ)**2.LE.EPS30) THEN
          WRITE (iunout,*)
     .       'ZERO B-FIELD UNIT VECTOR IN STANDARD CELL JJ= ',JJ
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
        B=SQRT(BXIN(JJ)**2+BYIN(JJ)**2+BZIN(JJ)**2)
        IF (ABS(B-1.D0).GT.EPS12) THEN
          WRITE (iunout,*)
     .       'B-FIELD UNIT VECTOR IN STANDARD CELL JJ= ',JJ,B
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
        IF (ABS(BFIN(JJ)).LT.EPS12) THEN
          WRITE (iunout,*) 'MAGNETIC FIELD STRENGTH = ZERO, JJ= ',JJ
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
153   CONTINUE

154   CONTINUE  !  BFIELD SPECIFIED AT ALL ??
C
C  ADDITIONAL INPUT TALLIES

      IF (.NOT.LADIN) GOTO 1160
      IND=INDPRO(6)
      DO 160 K=1,NAINI
        GOTO (151,151,151,151,155,156,157,160,160),IND
C  DEFAULT: ZERO, only options ind=5,6,7 are available
c          (transfer from problem specific codes or external data structures)
151     CONTINUE
        DO 1151 J=1,NR1ST
          ADIN(K,J)=0.
1151    CONTINUE
        GOTO 160
155     CALL EIRENE_PROUSR (HELP,6+1*NPLS+NPLSTI+3*NPLSV,
     .                      BD,BD,BD,BD,BD,BD,0._DP,NSURF)
        ADIN(K,1:NSURF)=HELP(1:NSURF)
        GOTO 160
160   CONTINUE
      GOTO 1160
156   CALL EIRENE_PROFR (ADIN,6+1*NPLS+NPLSTI+3*NPLSV,NAINI,NAIN,NSURF)
      GOTO 1160
157   CALL EIRENE_PROFR (ADIN,6+1*NPLS+NPLSTI+3*NPLSV,NAINI,NAIN,NSBOX)
      GOTO 1160
C
1160  CONTINUE
C
C  ELECTRIC FIELD
      IF (.NOT.(LEXIN.AND.LEYIN.AND.LEZIN.AND.LEFIN)) GOTO 170
      IND=INDPRO(7)
      GOTO (170,170,170,170,175,176,177,170,170),IND
C  DEFAULT: ZERO, only options ind=5,6,7
c          (transfer from problem specific codes or external data structures)
      goto 170
175     CALL EIRENE_PROUSR (EXIN,7+1*NPLS+NPLSTI+3*NPLSV,
     .                      EF0,EF1,EF2,EF3,EF4,EF5,0._DP,NSURF)
        CALL EIRENE_PROUSR (EYIN,8+1*NPLS+NPLSTI+3*NPLSV,
     .                      EF0,EF1,EF2,EF3,EF4,EF5,0._DP,NSURF)
        CALL EIRENE_PROUSR (EZIN,9+1*NPLS+NPLSTI+3*NPLSV,
     .                      EF0,EF1,EF2,EF3,EF4,EF5,0._DP,NSURF)
        CALL EIRENE_PROUSR (EFIN,10+1*NPLS+NPLSTI+3*NPLSV,
     .                      EF0,EF1,EF2,EF3,EF4,EF5,0._DP,NSURF)
        GOTO 170
176     CALL EIRENE_PROFR (EXIN,7+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        CALL EIRENE_PROFR (EYIN,8+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        CALL EIRENE_PROFR (EZIN,9+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        CALL EIRENE_PROFR (EFIN,10+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        GOTO 170
177     CALL EIRENE_PROFR (EXIN,7+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        CALL EIRENE_PROFR (EYIN,8+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        CALL EIRENE_PROFR (EZIN,9+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        CALL EIRENE_PROFR (EFIN,10+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        GOTO 170
c
170   CONTINUE
C
CDR
C   SET VACUUM DATA IN ADDITIONAL REGIONS OUTSIDE THE
C   THE STANDARD MESH. 
C   EXCLUDE: INDPRO=4: ADDITIONAL CELL REGION FROM FILE ISTREAM
C   EXCLUDE: INDPRO=7: ADDITIONAL CELL REGION DATA FROM EXTERNAL CODE (PROFR) 
C   EXCLUDE: INDPRO=8: ??

cdr tbd: indpro=4:  are data set on 1:nsurf, or on 1:nsbox=nsurf+nradd ?  
C
      IF (INDPRO(1).LE.6 .OR. INDPRO(1).EQ.9) THEN
        DO J=NSURF+1,NSURF+NRADD
          TEIN(J)=TVAC
        ENDDO
      ENDIF
      IF (INDPRO(2).LE.6 .OR. INDPRO(2).EQ.9) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 17 IPLS=1,NPLSI
            IPLSTI=MPLSTI(IPLS)
            TIIN(IPLSTI,J)=TVAC
17        CONTINUE
        ENDDO
      ENDIF
      IF (INDPRO(3).LE.6 .OR. INDPRO(3).EQ.9) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 18 IPLS=1,NPLSI
            DIIN(IPLS,J)=DVAC
18        CONTINUE
        ENDDO
      ENDIF
      IF (INDPRO(4).LE.6 .OR. INDPRO(4).EQ.9) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 19 IPLS=1,NPLSI
            IPLSV=MPLSV(IPLS)
            VXIN(IPLSV,J)=VVAC
            VYIN(IPLSV,J)=VVAC
            VZIN(IPLSV,J)=VVAC
19        CONTINUE
        ENDDO
      ENDIF
      IF (LBXIN.AND.LBYIN.AND.LBZIN.AND.LBFIN) THEN
        IF (INDPRO(5).LE.6 .OR. INDPRO(5).EQ.9) THEN
          DO J=NSURF+1,NSURF+NRADD
            BXIN(J)=0.
            BYIN(J)=0.
            BZIN(J)=1.
            BFIN(J)=1.
          ENDDO
        ENDIF
      ENDIF
      IF (LADIN .AND. (INDPRO(6).LE.6 .OR. INDPRO(6).EQ.9)) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 20 IAIN=1,NAINI
            ADIN(IAIN,J)=0.
20        CONTINUE
        ENDDO
      ENDIF
      IF (LEXIN.AND.LEYIN.AND.LEZIN.AND.LEFIN) THEN
        IF (INDPRO(7) == 5 .OR. INDPRO(5).EQ.6
     .                     .OR. INDPRO(5).EQ.9) THEN
          DO J=NSURF+1,NSURF+NRADD
            EXIN(J)=0.
            EYIN(J)=0.
            EZIN(J)=0.
            EFIN(J)=1.
          ENDDO
        ENDIF
      ENDIF
 
      DEALLOCATE(HELP)
      DEALLOCATE(HELP2)
C
      RETURN
      END
