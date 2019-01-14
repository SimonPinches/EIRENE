

C
C      SUBROUTINE SETPRM_INTAL
C
      SUBROUTINE EIRENE_SETPRM_INTAL
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COUTAU
      USE EIRMOD_CCONA, ONLY: EPS10
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CTEXT
      USE EIRMOD_CTRCEI, ONLY: TRCTAL
 
      IMPLICIT NONE
 
      INTEGER :: NTESTP, J, ITAL, NLSTTL, INDGRAD, INDTL
      INTEGER :: NPLPRM_TEST
      REAL(DP):: TSAVE
C
 
C  LIVTALI: SWITCH OFF SOME INPUT TALLIES AUTOMATICALLY;
C           Finally set in COMUSR.f
c  Structure similar to LIVTALV(..) FOR OUTPUT (SCORED) VOLUME TALLIES
c            (which is finally set in CESTIM.f

c default setting:
      LIVTALI = .FALSE.

c  plasma background
      LIVTALI(1)   = .TRUE.       ! TE
      LIVTALI(2)   = .TRUE.       ! TI
      LIVTALI(3)   = .TRUE.       ! DEIN
      LIVTALI(4)   = .TRUE.       ! DIIN

C  background flow velocities
      LIVTALI(5)   = NPLSV>0      ! VXIN
      LIVTALI(6)   = NPLSV>0      ! VYIN 
      LIVTALI(7)   = NPLSV>0      ! VZIN 

c  magnetic field
      LIVTALI(8)   = .TRUE.       ! BXIN, else: =0.0
      LIVTALI(9)   = .TRUE.       ! BYIN, else: =0.0
      LIVTALI(10)  = .TRUE.       ! BZIN, else: =1.0
      LIVTALI(11)  = .TRUE.       ! BFIN, else: =1.0

      LIVTALI(12)  = NAIN>0       ! ADIN

      LIVTALI(13)  = .TRUE.       ! EDRIFT, else: =0.0
      LIVTALI(14)  = .TRUE.       ! VOL

      LIVTALI(23)  = NPLSV>0      ! BVIN, else: sign(1.,bvin)=1.0   
      LIVTALI(24)  = NPLS>0       ! PARMOM, else: = 0.0

C  CURRENTLY THE LAST INPUT TALLY IS TALLY NO. 24 (=NTALG)

C  INTLOPT < 0  : SWITCH OFF TALLY
C          = 0  : KEEP DEFAULT
C          = 1  : EXPLICITELY SWITCH ON TALLY
C          = 2  : PREPARE FOR SMOOTHING, INTERPOLATE TO CORNERPOINTS
C          = 3  : SWITCH ON GRADIENTS

      DO ITAL = 1, NTALI
        IF (INTLOPTS(ITAL) < 0) THEN
C  SWITCH OFF TALLY
          LIVTALI(ITAL) = .FALSE.
        ELSE IF (INTLOPTS(ITAL) > 0) THEN
C  EXPLICITELY SWITCH ON TALLY
          LIVTALI(ITAL) = .TRUE.
          IF (ITAL <= NTALG) THEN
C  SWITCH ON INTERPOLATION TO CELL VERTICES  
            IF (INTLOPTS(ITAL) >= 2) LSMOPRO(ITAL) = .TRUE.
C  SWITCH ON GRADIENT TALLIES d(TL)/dX, d(TL)/dY,  d(TL)/dZ 
            IF (INTLOPTS(ITAL) == 3) THEN
              INDGRAD = NTALG + (ITAL-1)*3
              LIVTALI(INDGRAD+1 : INDGRAD+3) = .TRUE.
            END IF
          ELSE   ! ITAL > NTALG
C  TALLY ITAL IS A GRADIENT TALLY; 
C  ENSURE THAT INTERPOLATION TO CELL VERTICES IS SWITCHED ON FOR
C  CORRESPONDING INPUT TALLY
            INDTL = (ITAL - NTALG) / 3 
            IF (MOD(ITAL-NTALG,3) > 0) INDTL = INDTL + 1
C  SWITCH ON INTERPOLATION TO CELL VERTICES  
            LSMOPRO(INDTL) = .TRUE.
          END IF
        END IF
      END DO 

C  ENSURE THAT THE PRIMARY INPUT TALLIES FOR THE BACKGROUND PLASMA
C  ARE SWITCHED ON

      IF (.NOT.LTEIN) THEN
        WRITE (IUNOUT,*) ' SWITCHING OFF OF ELECTRON TEMPERATURE' //
     .                   ' IS PROHIBITED'
        WRITE (IUNOUT,*) ' TALLY IS SWITCHED ON AGAIN '
        LTEIN = .TRUE.
        INTLOPTS(1) = 0
      END IF

      IF (.NOT.LTIIN) THEN
        WRITE (IUNOUT,*) ' SWITCHING OFF OF ION TEMPERATURE' //
     .                   ' IS PROHIBITED'
        WRITE (IUNOUT,*) ' SET NPLSTI = 1 AND TI = TE '
        LTIIN = .TRUE.
        INTLOPTS(2) = 0
        NPLSTI = 1
        MPLSTI = 1
      END IF

      IF (.NOT.(LDEIN .AND. LDIIN)) THEN
        WRITE (IUNOUT,*) ' SWITCHING OFF OF BACKGROUND DENSITIES' //
     .                   ' IS PROHIBITED'
        WRITE (IUNOUT,*) ' TALLIES ARE SWITCHED ON AGAIN '
        LDEIN = .TRUE.
        INTLOPTS(3) = MAX(0,INTLOPTS(3))
        LDIIN = .TRUE.
        INTLOPTS(4) = MAX(0,INTLOPTS(4))
      END IF

      IF (.NOT.(LVXIN .AND. LVYIN .AND. LVZIN)) THEN
        WRITE (IUNOUT,*) ' SWITCHING OFF OF PLASMA DRIFT VELOCITY' //
     .                   ' IS PROHIBITED'
        WRITE (IUNOUT,*) ' TALLIES ARE SWITCHED ON AGAIN '
        LVXIN = .TRUE.
        INTLOPTS(5) = 0
        LVYIN = .TRUE.
        INTLOPTS(6) = 0
        LVZIN = .TRUE.
        INTLOPTS(7) = 0
      END IF

      IF (.NOT.LVOL) THEN
        WRITE (IUNOUT,*) ' SWITCHING OFF OF CELL VOLUME' //
     .                   ' IS PROHIBITED'
        WRITE (IUNOUT,*) ' TALLY IS SWITCHED ON AGAIN '
        LVOL = .TRUE.
        INTLOPTS(14) = 0
      END IF

C  ENSURE THAT CONNECTED TALLIES HAVE THE SAME SETTING

      IF (.NOT.(LBXPERP .AND. LBYPERP)) THEN
        LBXPERP = .FALSE.
        LBYPERP = .FALSE.
      END IF

      IF (.NOT.(LBXIN .AND. LBYIN .AND. LBZIN .AND. LBFIN)) THEN
        LBXIN = .FALSE.
        LBYIN = .FALSE.
        LBZIN = .FALSE.
        LBFIN = .FALSE.
        LBXPERP = .FALSE.
        LBYPERP = .FALSE.
      END IF

      IF (.NOT.(LEXIN .AND. LEYIN .AND. LEZIN .AND. LEFIN)) THEN
        LEXIN = .FALSE.
        LEYIN = .FALSE.
        LEZIN = .FALSE.
        LEFIN = .FALSE.
        LPOT  = .FALSE.
      END IF

      
C  18 primary input tallies plus 6 derived background tallies, (# ...)
c     unfortunately mixed 
C  --> 24 rather than 18 background tallies
C  --> 30 include 6 free slots
      NFRSTP(1)=0
      NFRSTP(2)=NPLSTI
      NFRSTP(3)=0       ! # DEIN,  DERIVED QUANTITY
      NFRSTP(4)=NPLS    ! DIIN
      NFRSTP(5)=NPLSV   
      NFRSTP(6)=NPLSV
      NFRSTP(7)=NPLSV
      NFRSTP(8)=0       ! BX
      NFRSTP(9)=0       ! BY
      NFRSTP(10)=0      ! BZ
      NFRSTP(11)=0      ! BF
      NFRSTP(12)=NAIN   ! ADIN
      NFRSTP(13)=NPLS   ! # EDRIFT,  DERIVED QUANTITY
      NFRSTP(14)=0      ! VOL
      NFRSTP(15)=NSPZMC ! WEIGHT WINDOW, UNUSED  
      NFRSTP(16)=0      ! # BX_PERP,  DERIVED QUANTITY
      NFRSTP(17)=0      ! # BY_PERP,  DERIVED QUANTITY 
      NFRSTP(18)=0      ! EX
      NFRSTP(19)=0      ! EY 
      NFRSTP(20)=0      ! EZ
      NFRSTP(21)=0      ! EF
      NFRSTP(22)=0      ! POT
      NFRSTP(23)=NPLSV  ! # BVIN
      NFRSTP(24)=NPLS   ! # PARMON
      NFRSTP(25)=0      ! # PSI

      NFRSTP(26)=0      ! # FREE26
      NFRSTP(27)=0      ! # FREE27
      NFRSTP(28)=0      ! # FREE28
      NFRSTP(29)=0      ! # FREE29
      NFRSTP(30)=0      ! # FREE30

c  from here on: derivatives (gradients) of input tallies
      NFRSTP(31)=0
      NFRSTP(32)=0
      NFRSTP(33)=0
      NFRSTP(34)=NPLSTI
      NFRSTP(35)=NPLSTI
      NFRSTP(36)=NPLSTI
      NFRSTP(37)=0       ! # DEIN,  DERIVED QUANTITY
      NFRSTP(38)=0       ! # DEIN,  DERIVED QUANTITY
      NFRSTP(39)=0       ! # DEIN,  DERIVED QUANTITY
      NFRSTP(40)=NPLS    ! DIIN
      NFRSTP(41)=NPLS    ! DIIN
      NFRSTP(42)=NPLS    ! DIIN
      NFRSTP(43)=NPLSV
      NFRSTP(44)=NPLSV
      NFRSTP(45)=NPLSV
      NFRSTP(46)=NPLSV
      NFRSTP(47)=NPLSV
      NFRSTP(48)=NPLSV
      NFRSTP(49)=NPLSV
      NFRSTP(50)=NPLSV
      NFRSTP(51)=NPLSV
      NFRSTP(52)=0       ! BX
      NFRSTP(53)=0       ! BX
      NFRSTP(54)=0       ! BX
      NFRSTP(55)=0       ! BY
      NFRSTP(56)=0       ! BY
      NFRSTP(57)=0       ! BY
      NFRSTP(58)=0      ! BZ
      NFRSTP(59)=0      ! BZ
      NFRSTP(60)=0      ! BZ
      NFRSTP(61)=0      ! BF
      NFRSTP(62)=0      ! BF
      NFRSTP(63)=0      ! BF
      NFRSTP(64)=NAIN   ! ADIN
      NFRSTP(65)=NAIN   ! ADIN
      NFRSTP(66)=NAIN   ! ADIN
      NFRSTP(67)=NPLS   ! # EDRIFT,  DERIVED QUANTITY
      NFRSTP(68)=NPLS   ! # EDRIFT,  DERIVED QUANTITY
      NFRSTP(68)=NPLS   ! # EDRIFT,  DERIVED QUANTITY
      NFRSTP(70)=0      ! VOL
      NFRSTP(71)=0      ! VOL
      NFRSTP(72)=0      ! VOL
      NFRSTP(73)=NSPZMC ! WEIGHT WINDOW, UNUSED  
      NFRSTP(74)=NSPZMC ! WEIGHT WINDOW, UNUSED  
      NFRSTP(75)=NSPZMC ! WEIGHT WINDOW, UNUSED  
      NFRSTP(76)=0      ! # BX_PERP,  DERIVED QUANTITY
      NFRSTP(77)=0      ! # BX_PERP,  DERIVED QUANTITY
      NFRSTP(78)=0      ! # BX_PERP,  DERIVED QUANTITY
      NFRSTP(79)=0      ! # BY_PERP,  DERIVED QUANTITY 
      NFRSTP(80)=0      ! # BY_PERP,  DERIVED QUANTITY 
      NFRSTP(81)=0      ! # BY_PERP,  DERIVED QUANTITY 
      NFRSTP(82)=0      ! EX
      NFRSTP(83)=0      ! EX
      NFRSTP(84)=0      ! EX
      NFRSTP(85)=0      ! EY 
      NFRSTP(86)=0      ! EY 
      NFRSTP(87)=0      ! EY 
      NFRSTP(88)=0      ! EZ
      NFRSTP(89)=0      ! EZ
      NFRSTP(90)=0      ! EZ
      NFRSTP(91)=0      ! EF
      NFRSTP(92)=0      ! EF
      NFRSTP(93)=0      ! EF
      NFRSTP(94)=0      ! POT
      NFRSTP(95)=0      ! POT
      NFRSTP(96)=0      ! POT
      NFRSTP(97)=NPLSV  ! BVIN
      NFRSTP(98)=NPLSV  ! BVIN
      NFRSTP(99)=NPLSV  ! BVIN
      NFRSTP(100)=NPLS   ! PARMON
      NFRSTP(101)=NPLS   ! PARMON
      NFRSTP(102)=NPLS   ! PARMON
      NFRSTP(103)=0      ! PSI
      NFRSTP(104)=0      ! PSI
      NFRSTP(105)=0      ! PSI
      NFRSTP(106)=0      ! FREE26
      NFRSTP(107)=0      ! FREE26
      NFRSTP(108)=0      ! FREE26
      NFRSTP(109)=0      ! FREE27
      NFRSTP(110)=0      ! FREE27
      NFRSTP(111)=0      ! FREE27
      NFRSTP(112)=0      ! FREE28
      NFRSTP(113)=0      ! FREE28
      NFRSTP(114)=0      ! FREE28
      NFRSTP(115)=0      ! FREE29
      NFRSTP(116)=0      ! FREE29
      NFRSTP(117)=0      ! FREE29
      NFRSTP(118)=0      ! FREE30
      NFRSTP(119)=0      ! FREE30
      NFRSTP(120)=0      ! FREE30
C
C  NTALI=96?  number of input tallies  (19 PRIMARY + 5 DERIVED + 24 GRADIENTS)
 
cdr Since primary and derived input tallies got mixed up anyway, 
cdr add magnetic flux (vector potential). 
cdr Be careful:
cdr in some places in the code currently the numbering  of input tallies is hard coded.
cdr (ALGTAL, OUTPLA,....) 
C
      DO 5 J=1,NTALI
        IF (LIVTALI(J)) NFRSTP(J)=MAX0(1,NFRSTP(J))
5     CONTINUE
C
C  SET CUMULATED FIRST INDICES OF ACTIVE INPUT TALLIES: NADDP(J)
C  THE LAST ACTIVE INPUT TALLY IS TALLY NO. NLSTTL
      NADDP(1)=0
      DO 6 J=2,NTALI
        IF (LIVTALI(J-1)) THEN
          NADDP(J)=NADDP(J-1)+NFRSTP(J-1)
          NLSTTL=J-1
        ELSE
          NADDP(J)=NADDP(J-1)
        END IF
 6    CONTINUE

C  TOTAL NUMBER OF INPUT TALLIES, ALSO COUNTING FIRST INDICES
      NINPTL = NADDP(NTALI)

CDR  CORRECT FOR LAST TALLY:
      IF (LIVTALI(NTALI)) THEN 
        NLSTTL = NTALI
C  TOTAL NUMBER OF INPUT TALLIES, LAST TALLY NTALI WAS NOT YET IN NADDP(..)
        NINPTL = NINPTL+NFRSTP(NLSTTL)
      ENDIF
      
c  allocate and initialize plstls
      CALL EIRENE_ALLOC_COMUSR(2)
c  set pointers:  input tallies tein, tiin,....parmom,.....on plttls
      CALL EIRENE_ASSOCIATE_COMUSR

cdr  hard coded test, in case parmom is the last active input tally
cdr   tsave=parmom(npls,nrad)
cdr   parmom(npls,nrad)=1.2345678
cdr   if (abs(plstls(ninptl,nrad)-1.2345678).gt.eps10) then
cdr     write (iunout,*) 'error ninptl, assuming parmom is last tally'
cdr     call eirene_exit_own(1)
cdr   endif
cdr   parmom(npls,nrad)=tsave

!  CHECK VALUE ON LAST CELL IN LAST ACTIVE TALLY
!  THIS TEST CAN NOT BE PERFORMED DUE TO SWITCHING OFF OF INPUT TALLIES
!     NTESTP=NINPTL*NRAD  ! STORAGE POSITION OF LAST CELL IN LAST TALLY


c............................................................................. 
 
      IF (TRCTAL) THEN        
        CALL EIRENE_LEER(2)
        WRITE(IUNOUT,*) 'INPUT TALLIES USED IN THIS RUN'
        CALL EIRENE_LEER(1)
        WRITE(IUNOUT,'(A6,1X,A)') 'NO.','DESCRIPTION'
        DO ITAL=1,NTALI
          IF (LIVTALI(ITAL))
     .      WRITE (IUNOUT,'(I6,1X,A72)') ITAL,TXTPLS(1,ITAL)
        END DO
 
        IF (.NOT.ALL(LIVTALI)) THEN
          CALL EIRENE_LEER(2)
          WRITE(IUNOUT,*) 'INPUT TALLIES NOT AVAILABLE ',
     .                    'IN THIS RUN'
          CALL EIRENE_LEER(1)
          WRITE(IUNOUT,'(A6,1X,A)') 'NO.','DESCRIPTION'
          DO ITAL=1,NTALI
            IF (.NOT.LIVTALI(ITAL))
     .        WRITE (IUNOUT,'(I6,1X,A72)') ITAL,TXTPLS(1,ITAL)
          END DO
        END IF
 
        IF (ANY(INTLOPTS < 0)) THEN
          CALL EIRENE_LEER(2)
          WRITE(IUNOUT,*) 'INPUT TALLIES EXPLICITLY ',
     .                'SWITCHED OFF VIA INPUT FILE '
          CALL EIRENE_LEER(1)
          WRITE(IUNOUT,'(A6,1X,A)') 'NO.','DESCRIPTION'
          DO ITAL=1,NTALI
            IF (INTLOPTS(ITAL) < 0)
     .        WRITE (IUNOUT,'(I6,1X,A72)') ITAL,TXTPLS(1,ITAL)
          END DO
        END IF
  
        CALL EIRENE_LEER(1)
 
      END IF
C
      RETURN
      END
 
 
 
 
