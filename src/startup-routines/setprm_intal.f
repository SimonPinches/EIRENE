

C
C      SUBROUTINE SETPRM_INTAL
C
      SUBROUTINE EIRENE_SETPRM_INTAL
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COUTAU
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CTEXT
      USE EIRMOD_CTRCEI, ONLY: TRCTAL
 
      IMPLICIT NONE
 
      INTEGER :: NTESTP, J, ITAL, NLSTTL, INDGRAD, INDTL
      INTEGER :: NPLPRM_TEST
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

C  CURRENTLY THE LAST INPUT TALLY IS TALLY NO. 24 (NTALG)

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

      
C  19 primary input tallies plus 5 derived background tallies,
c     unfortunately mixed in
C  --> 24 rather than 18 background tallies
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
      NFRSTP(23)=NPLSV  ! BVIN
      NFRSTP(24)=NPLS   ! PARMON

c  from here on: derivatives (gradients) of input tallies
      NFRSTP(25)=0
      NFRSTP(26)=0
      NFRSTP(27)=0
      NFRSTP(28)=NPLSTI
      NFRSTP(29)=NPLSTI
      NFRSTP(30)=NPLSTI
      NFRSTP(31)=0       ! # DEIN,  DERIVED QUANTITY
      NFRSTP(32)=0       ! # DEIN,  DERIVED QUANTITY
      NFRSTP(33)=0       ! # DEIN,  DERIVED QUANTITY
      NFRSTP(34)=NPLS    ! DIIN
      NFRSTP(35)=NPLS    ! DIIN
      NFRSTP(36)=NPLS    ! DIIN
      NFRSTP(37)=NPLSV   
      NFRSTP(38)=NPLSV   
      NFRSTP(39)=NPLSV   
      NFRSTP(40)=NPLSV
      NFRSTP(41)=NPLSV
      NFRSTP(42)=NPLSV
      NFRSTP(43)=NPLSV
      NFRSTP(44)=NPLSV
      NFRSTP(45)=NPLSV
      NFRSTP(46)=0       ! BX
      NFRSTP(47)=0       ! BX
      NFRSTP(48)=0       ! BX
      NFRSTP(49)=0       ! BY
      NFRSTP(50)=0       ! BY
      NFRSTP(51)=0       ! BY
      NFRSTP(52)=0      ! BZ
      NFRSTP(53)=0      ! BZ
      NFRSTP(54)=0      ! BZ
      NFRSTP(55)=0      ! BF
      NFRSTP(56)=0      ! BF
      NFRSTP(57)=0      ! BF
      NFRSTP(58)=NAIN   ! ADIN
      NFRSTP(59)=NAIN   ! ADIN
      NFRSTP(60)=NAIN   ! ADIN
      NFRSTP(61)=NPLS   ! # EDRIFT,  DERIVED QUANTITY
      NFRSTP(62)=NPLS   ! # EDRIFT,  DERIVED QUANTITY
      NFRSTP(63)=NPLS   ! # EDRIFT,  DERIVED QUANTITY
      NFRSTP(64)=0      ! VOL
      NFRSTP(65)=0      ! VOL
      NFRSTP(66)=0      ! VOL
      NFRSTP(67)=NSPZMC ! WEIGHT WINDOW, UNUSED  
      NFRSTP(68)=NSPZMC ! WEIGHT WINDOW, UNUSED  
      NFRSTP(69)=NSPZMC ! WEIGHT WINDOW, UNUSED  
      NFRSTP(70)=0      ! # BX_PERP,  DERIVED QUANTITY
      NFRSTP(71)=0      ! # BX_PERP,  DERIVED QUANTITY
      NFRSTP(72)=0      ! # BX_PERP,  DERIVED QUANTITY
      NFRSTP(73)=0      ! # BY_PERP,  DERIVED QUANTITY 
      NFRSTP(74)=0      ! # BY_PERP,  DERIVED QUANTITY 
      NFRSTP(75)=0      ! # BY_PERP,  DERIVED QUANTITY 
      NFRSTP(76)=0      ! EX
      NFRSTP(77)=0      ! EX
      NFRSTP(78)=0      ! EX
      NFRSTP(79)=0      ! EY 
      NFRSTP(80)=0      ! EY 
      NFRSTP(81)=0      ! EY 
      NFRSTP(82)=0      ! EZ
      NFRSTP(83)=0      ! EZ
      NFRSTP(84)=0      ! EZ
      NFRSTP(85)=0      ! EF
      NFRSTP(86)=0      ! EF
      NFRSTP(87)=0      ! EF
      NFRSTP(88)=0      ! POT
      NFRSTP(89)=0      ! POT
      NFRSTP(90)=0      ! POT
      NFRSTP(91)=NPLSV  ! BVIN
      NFRSTP(92)=NPLSV  ! BVIN
      NFRSTP(93)=NPLSV  ! BVIN
      NFRSTP(94)=NPLS   ! PARMON
      NFRSTP(95)=NPLS   ! PARMON
      NFRSTP(96)=NPLS   ! PARMON
C
C  NTALI=24?  number of input tallies  (19 PRIMARY + 5 DERIVED)
cdr there are many more derived input tallies. 
cdr since primary and derived input tallies got mixed up anyway, 
cdr to do: change ntali, add other derived input tallies, here, and in settxt.
cdr be careful:
cdr in some places in code the numbering  of input tallies is hard coded.
cdr (ALGTAL, OUTPLA,....) 
C
      DO 5 J=1,NTALI
        IF (LIVTALI(J)) NFRSTP(J)=MAX0(1,NFRSTP(J))
5     CONTINUE
C
      NADDP(1)=0
      DO 6 J=2,NTALI
        IF (LIVTALI(J-1)) THEN
          NADDP(J)=NADDP(J-1)+NFRSTP(J-1)
          NLSTTL=J-1
        ELSE
          NADDP(J)=NADDP(J-1)
        END IF
 6    CONTINUE

      IF (LIVTALI(NTALI)) NLSTTL = NTALI
C
C  TOTAL NUMBER OF INPUT TALLIES
      NINPTL = NADDP(NTALI)+NFRSTP(NLSTTL)

      CALL EIRENE_ALLOC_COMUSR(2)
      CALL EIRENE_ASSOCIATE_COMUSR

!  THIS TEST CAN NOT BE PERFORMED DUE TO SWITCHING OFF OF INPUT TALLIES
!      NTESTP=NADDP(NTALI)+NFRSTP(NLSTTL)
!      NTESTP=NINPTL
!      NTESTP=NTESTP*NRAD

cdr  correct for the derived tallies mixed into primary input tallies.  
!      NPLPRM_TEST=NPLPRM + (2+NPLS)*NRAD
!      IF (NTESTP.NE.NPLPRM_TEST) THEN
!        WRITE (iunout,*) 'PARAMETER ERROR DETECTED IN SETPRM_INTAL: ' //
!     .                   'NPLPRM'
!        WRITE (iunout,*) 'NTESTP, NPLPRM ',NTESTP,NPLPRM_TEST
!        CALL EIRENE_EXIT_OWN(1)
!      ENDIF
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
 
 
 
 
