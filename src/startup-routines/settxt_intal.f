cdr  Aug 18: formerly this was part of settxt.f.
cdr  Generalisation of input tallies (and gradients thereof)
cdr  Aug 19: still some cleanup and documentation needed
cdr  Feb 22:  still inconsistent use of NPLS vs. NPLSI, etc...

      SUBROUTINE EIRENE_SETTXT_INTAL
c  Set default texts (volume tallies: name, species, units),
C  similar to SETTXT.f, but for INPUT TALLIES rather than output tallies.
C  Set first (leading) dimension of input tally arrays: NFSTPI.

c
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU

      IMPLICIT NONE

      INTEGER :: IPLS, ISPZ, I, J,
     .           KK, ICO   ! species index handling in vectorial tallies
      CHARACTER(24) :: TEXT24
      CHARACTER(72) :: TEXT72
C
C  TEXT FOR INPUT (BACKGROUND) TALLIES

cdr  the numbering is "a bit" illogical, due to historic reasons.
c    primary and derived tallies are mixed here.
C
      TXTPLS(1,1)='PLASMA TEMPERATURE (ELECTRONS)                   '
      TXTPLS(1,2)='PLASMA TEMPERATURE (BULK PARTICLES)              '
      TXTPLS(1,3)='PLASMA DENSITY (ELECTRONS)                       '
      TXTPLS(1,4)='PLASMA DENSITY (BULK PARTICLES)                  '
      TXTPLS(1,5)='DRIFT VELOCITY IN X-DIRECTION (BULK IONS)        '
      TXTPLS(1,6)='DRIFT VELOCITY IN Y-DIRECTION (BULK IONS)        '
      TXTPLS(1,7)='DRIFT VELOCITY IN Z-DIRECTION (BULK IONS)        '
      TXTPLS(1,8)='MAGN. FIELD UNIT VECTOR, X DIRECTION             '
      TXTPLS(1,9)='MAGN. FIELD UNIT VECTOR, Y DIRECTION             '
      TXTPLS(1,10)='MAGN. FIELD UNIT VECTOR, Z DIRECTION             '
      TXTPLS(1,11)='MAGN. FIELD STRENGTH                             '
cdr to be added here TXTPLS(1,xx)='Magn. POTENTIAL, e.g. PSI fct.'
cdr 2019: now added as tally 25 below.
      TXTPLS(1,12)='ADDITIONAL INPUT TALLIES, OPTIONAL'
      TXTPLS(1,13)='BULK ION KINETIC DRIFT ENERGY                    '
      TXTPLS(1,14)='ZONE VOLUMES                                     '

      TXTPLS(1,15)='SPACE-SPECIES WEIGHT FUNCTION                    '
      TXTPLS(1,16)='PERP. MAGN. FIELD VECTOR, X DIRECTION            '
      TXTPLS(1,17)='PERP. MAGN. FIELD VECTOR, Y DIRECTION            '

      TXTPLS(1,18)='ELEC. FIELD UNIT VECTOR, X DIRECTION             '
      TXTPLS(1,19)='ELEC. FIELD UNIT VECTOR, Y DIRECTION             '
      TXTPLS(1,20)='ELEC. FIELD UNIT VECTOR, Z DIRECTION             '
      TXTPLS(1,21)='ELEC. FIELD STRENGTH                             '
      TXTPLS(1,22)='ELECTR. POTENTIAL                                '
      TXTPLS(1,23)='FLOW VELOCITY WRT. B FIELD                       '
      TXTPLS(1,24)='MOMENTUM FLOW WRT. B FIELD                       '
      TXTPLS(1,25)='PSI                                              '
      TXTPLS(1,26)='ZI                                               '

      TXTPLS(1,27)='FREE27                                           '
      TXTPLS(1,28)='FREE28                                           '
      TXTPLS(1,29)='FREE29                                           '
      TXTPLS(1,30)='FREE30                                           '

C  gradient tallies:  derivatives wrt. cartesian coordinates
c 1
      TXTPLS(1,31)='dTE/dX                                           '
      TXTPLS(1,32)='dTE/dY                                           '
      TXTPLS(1,33)='dTE/dZ                                           '
c 2
      TXTPLS(1,34)='dTI/dX                                           '
      TXTPLS(1,35)='dTI/dY                                           '
      TXTPLS(1,36)='dTI/dZ                                           '
c 3
      TXTPLS(1,37)='dNE/dX                                           '
      TXTPLS(1,38)='dNE/dY                                           '
      TXTPLS(1,39)='dNE/dZ                                           '
c 4
      TXTPLS(1,40)='dNI/dX                                           '
      TXTPLS(1,41)='dNI/dY                                           '
      TXTPLS(1,42)='dNI/dZ                                           '

      TXTPLS(1,43)='dVX/dX                                           '
      TXTPLS(1,44)='dVX/dY                                           '
      TXTPLS(1,45)='dVX/dZ                                           '

      TXTPLS(1,46)='dVY/dX                                           '
      TXTPLS(1,47)='dVY/dY                                           '
      TXTPLS(1,48)='dVY/dZ                                           '

      TXTPLS(1,49)='dVZ/dX                                           '
      TXTPLS(1,50)='dVZ/dY                                           '
      TXTPLS(1,51)='dVZ/dZ                                           '
c 8
      TXTPLS(1,52)='dBX/dX                                           '
      TXTPLS(1,53)='dBX/dY                                           '
      TXTPLS(1,54)='dBX/dZ                                           '
      TXTPLS(1,55)='dBY/dX                                           '
      TXTPLS(1,56)='dBY/dY                                           '
      TXTPLS(1,57)='dBY/dZ                                           '
      TXTPLS(1,58)='dBZ/dX                                           '
      TXTPLS(1,59)='dBZ/dY                                           '
      TXTPLS(1,60)='dBZ/dZ                                           '

      TXTPLS(1,61)='dBF/dX                                           '
      TXTPLS(1,62)='dBF/dY                                           '
      TXTPLS(1,63)='dBF/dZ                                           '
c 12
      TXTPLS(1,64)='dADIN/dX                                         '
      TXTPLS(1,65)='dADIN/dY                                         '
      TXTPLS(1,66)='dADIN/dZ                                         '

      TXTPLS(1,67)='dEDRIFT/dX                                       '
      TXTPLS(1,68)='dEDRIFT/dY                                       '
      TXTPLS(1,69)='dEDRIFT/dZ                                       '

      TXTPLS(1,70)='dVOL/dX                                          '
      TXTPLS(1,71)='dVOL/dY                                          '
      TXTPLS(1,72)='dVOL/dZ                                          '

      TXTPLS(1,73)='dWGHT/dX                                         '
      TXTPLS(1,74)='dWGHT/dY                                         '
      TXTPLS(1,75)='dWGHT/dZ                                         '

      TXTPLS(1,76)='dBXPERP/dX                                       '
      TXTPLS(1,77)='dBXPERP/dY                                       '
      TXTPLS(1,78)='dBXPERP/dZ                                       '

      TXTPLS(1,79)='dBYPERP/dX                                       '
      TXTPLS(1,80)='dBYPERP/dY                                       '
      TXTPLS(1,81)='dBYPERP/dZ                                       '

      TXTPLS(1,82)='dEX/dX                                           '
      TXTPLS(1,83)='dEX/dY                                           '
      TXTPLS(1,84)='dEX/dZ                                           '
      TXTPLS(1,85)='dEY/dX                                           '
      TXTPLS(1,86)='dEY/dY                                           '
      TXTPLS(1,87)='dEY/dZ                                           '
      TXTPLS(1,88)='dEZ/dX                                           '
      TXTPLS(1,89)='dEZ/dY                                           '
      TXTPLS(1,90)='dEZ/dZ                                           '
      TXTPLS(1,91)='dEF/dX                                           '
      TXTPLS(1,92)='dEF/dY                                           '
      TXTPLS(1,93)='dEF/dZ                                           '
      TXTPLS(1,94)='dPOT/dX                                          '
      TXTPLS(1,95)='dPOT/dY                                          '
      TXTPLS(1,96)='dPOT/dZ                                          '
      TXTPLS(1,97)='dBV_VEC/dX                                       '
      TXTPLS(1,98)='dBV_VEC/dY                                       '
      TXTPLS(1,99)='dBV_VEC/dZ                                       '

      TXTPLS(1,100)='dPMOM_VEC/dX                                     '
      TXTPLS(1,101)='dPMOM_VEC/dY                                     '
      TXTPLS(1,102)='dPMOM_VEC/dZ                                     '
c 25 psi-function
      TXTPLS(1,103)='dPSI/dX                                          '
      TXTPLS(1,104)='dPSI/dY                                          '
      TXTPLS(1,105)='dPSI/dZ                                          '

cnh 19.08.2020 Average charge
      TXTPLS(1,106)='dZI/dX                                           '
      TXTPLS(1,107)='dZI/dY                                           '
      TXTPLS(1,108)='dZI/dZ                                           '

      TXTPLS(1,109)='dFREE27/dX                                       '
      TXTPLS(1,110)='dFREE27/dY                                       '
      TXTPLS(1,111)='dFREE27/dZ                                       '
      TXTPLS(1,112)='dFREE28/dX                                       '
      TXTPLS(1,113)='dFREE28/dY                                       '
      TXTPLS(1,114)='dFREE28/dZ                                       '
      TXTPLS(1,115)='dFREE29/dX                                       '
      TXTPLS(1,116)='dFREE29/dY                                       '
      TXTPLS(1,117)='dFREE29/dZ                                       '
      TXTPLS(1,118)='dFREE30/dX                                       '
      TXTPLS(1,119)='dFREE30/dY                                       '
      TXTPLS(1,120)='dFREE30/dZ                                       '
C
c  currently: ntali=4*30=120
      DO J=1,NTALI
        IF (J.EQ.12) CYCLE  ! retain individual tally names
                            ! for adin tally no. 12
cdr vectorial input tallies 23, 24:
        IF (J.EQ.23) THEN
          DO I=1,NPLS
            TEXT72='FLOW VELOCITY PARALLEL TO B FIELD                '
            TXTPLS(I,J)=TEXT72
          ENDDO
          DO I=1,NPLS
            TEXT72='FLOW VELOCITY PERP.    TO B FIELD                '
            TXTPLS(NPLS+I,J)=TEXT72
          ENDDO
          DO I=1,NPLS
            TEXT72='FLOW VELOCITY DIAMAGN. TO B FIELD                '
            TXTPLS(2*NPLS+I,J)=TEXT72
          ENDDO
        ELSEIF (J.EQ.24) THEN
          DO I=1,NPLS
            TEXT72='MOMENTUM FLOW PARALLEL TO B FIELD                '
            TXTPLS(I,J)=TEXT72
          ENDDO
          DO I=1,NPLS
            TEXT72='MOMENTUM FLOW PERP.    TO B FIELD                '
            TXTPLS(NPLS+I,J)=TEXT72
          ENDDO
          DO I=1,NPLS
            TEXT72='MOMENTUM FLOW DIAMAGN. TO B FIELD                '
            TXTPLS(2*NPLS+I,J)=TEXT72
          ENDDO
        ELSE
          TEXT72=TXTPLS(1,J)
          DO I=2,N1MX
            TXTPLS(I,J)=TEXT72
          ENDDO
        ENDIF
      ENDDO
C
      TXTPUN(1,1)='EV                      '
      TXTPUN(1,2)='EV                      '
      TXTPUN(1,3)='CM**-3                  '
      TXTPUN(1,4)='CM**-3                  '
      TXTPUN(1,5)='CM/SEC                  '
      TXTPUN(1,6)='CM/SEC                  '
      TXTPUN(1,7)='CM/SEC                  '
      TXTPUN(1,8)=' ---                    '
      TXTPUN(1,9)=' ---                    '
      TXTPUN(1,10)=' ---                    '
      TXTPUN(1,11)='TESLA                   '
cdr  here might come: magn. potential (at least: tor. component?): PSI fct. TESLA*CM
CDR  JUNE 2019: SEE BELOW, TALLY 25
      TXTPUN(1,12)='ADDITIONAL TALLY UNITS  '
! EDRIFT  --> DERIVED QUANTITY
      TXTPUN(1,13)='EV                      '
      TXTPUN(1,14)='CM**3                   '  ! VOL
      TXTPUN(1,15)=' ---                    '  ! WEIGHT WINDOW
! BX_PERP --> DERIVED QUANTITY
      TXTPUN(1,16)=' ---                    '
! BY_PERP --> DERIVED QUANTITY
      TXTPUN(1,17)=' ---                    '
      TXTPUN(1,18)=' ---                    '
      TXTPUN(1,19)=' ---                    '
      TXTPUN(1,20)=' ---                    '
      TXTPUN(1,21)='V/CM                    '  ! EF
      TXTPUN(1,22)='V                       '  ! POT
      TXTPUN(1,23)='CM/S                    '  ! BV_VEC
      TXTPUN(1,24)='G*CM/S                  '  ! PMOM_VEC
! PSI   ! CHECK FOR FACTOR 2Pi
      TXTPUN(1,25)='TESLA*CM                '
      TXTPUN(1,26)=' ---                    '  ! ZI

      TXTPUN(1,27)=' ---                    '  ! FREE27
      TXTPUN(1,28)=' ---                    '  ! FREE28
      TXTPUN(1,29)=' ---                    '  ! FREE29
      TXTPUN(1,30)=' ---                    '  ! FREE30

cdr derivatives in cart. coordinates (gradient vector)
      TXTPUN(1,31)='EV/CM                   '  ! grad(Te)
      TXTPUN(1,32)='EV/CM                   '
      TXTPUN(1,33)='EV/CM                   '
      TXTPUN(1,34)='EV/CM                   '  ! grad(Ti)
      TXTPUN(1,35)='EV/CM                   '
      TXTPUN(1,36)='EV/CM                   '
      TXTPUN(1,37)='CM**-3/CM               '  ! grad(ne)
      TXTPUN(1,38)='CM**-3/CM               '
      TXTPUN(1,39)='CM**-3/CM               '
      TXTPUN(1,40)='CM**-3/CM               '  ! grad(ni)
      TXTPUN(1,41)='CM**-3/CM               '
      TXTPUN(1,42)='CM**-3/CM               '
      TXTPUN(1,43)='CM/SEC/CM               '
      TXTPUN(1,44)='CM/SEC/CM               '
      TXTPUN(1,45)='CM/SEC/CM               '
      TXTPUN(1,46)='CM/SEC/CM               '
      TXTPUN(1,47)='CM/SEC/CM               '
      TXTPUN(1,48)='CM/SEC/CM               '
      TXTPUN(1,49)='CM/SEC/CM               '
      TXTPUN(1,50)='CM/SEC/CM               '
      TXTPUN(1,51)='CM/SEC/CM               '
      TXTPUN(1,52)=' ---                    '
      TXTPUN(1,53)=' ---                    '
      TXTPUN(1,54)=' ---                    '
      TXTPUN(1,55)=' ---                    '
      TXTPUN(1,56)=' ---                    '
      TXTPUN(1,57)=' ---                    '
      TXTPUN(1,58)=' ---                    '
      TXTPUN(1,59)=' ---                    '
      TXTPUN(1,60)=' ---                    '
      TXTPUN(1,61)='TESLA/CM                '
      TXTPUN(1,62)='TESLA/CM                '
      TXTPUN(1,63)='TESLA/CM                '
C     TXTPUN(1,64)='TO BE READ, ADIN        '
C     TXTPUN(1,65)='TO BE READ, ADIN        '
C     TXTPUN(1,66)='TO BE READ, ADIN        '
! EDRIFT  --> DERIVED QUANTITY
      TXTPUN(1,67)='EV/CM                   '
! EDRIFT  --> DERIVED QUANTITY
      TXTPUN(1,68)='EV/CM                   '
! EDRIFT  --> DERIVED QUANTITY
      TXTPUN(1,69)='EV/CM                   '
      TXTPUN(1,70)='CM**3/CM                '  ! VOL
      TXTPUN(1,71)='CM**3/CM                '  ! VOL
      TXTPUN(1,72)='CM**3/CM                '  ! VOL
      TXTPUN(1,73)=' ---                    '  ! WEIGHT WINDOW
      TXTPUN(1,74)=' ---                    '  ! WEIGHT WINDOW
      TXTPUN(1,75)=' ---                    '  ! WEIGHT WINDOW
! BX_PERP --> DERIVED QUANTITY
      TXTPUN(1,76)=' ---                    '
! BX_PERP --> DERIVED QUANTITY
      TXTPUN(1,77)=' ---                    '
! BX_PERP --> DERIVED QUANTITY
      TXTPUN(1,78)=' ---                    '
! BY_PERP --> DERIVED QUANTITY
      TXTPUN(1,79)=' ---                    '
! BY_PERP --> DERIVED QUANTITY
      TXTPUN(1,80)=' ---                    '
! BY_PERP --> DERIVED QUANTITY
      TXTPUN(1,81)=' ---                    '
      TXTPUN(1,82)=' ---                    '
      TXTPUN(1,83)=' ---                    '
      TXTPUN(1,84)=' ---                    '
      TXTPUN(1,85)=' ---                    '
      TXTPUN(1,86)=' ---                    '
      TXTPUN(1,87)=' ---                    '
      TXTPUN(1,88)=' ---                    '
      TXTPUN(1,89)=' ---                    '
      TXTPUN(1,90)=' ---                    '
      TXTPUN(1,91)='V/CM/CM                 '  ! EF
      TXTPUN(1,92)='V/CM/CM                 '  ! EF
      TXTPUN(1,93)='V/CM/CM                 '  ! EF
      TXTPUN(1,94)='V/CM                    '  ! POT
      TXTPUN(1,96)='V/CM                    '  ! POT
      TXTPUN(1,96)='V/CM                    '  ! POT

      TXTPUN(1,97)='CM/S/CM                 '  ! BV_VEC
      TXTPUN(1,98)='CM/S/CM                 '  ! BV_VEC
      TXTPUN(1,99)='CM/S/CM                 '  ! BV_VEC
      TXTPUN(1,100)='G*CM/S/CM               '  ! PMOM_VEC
      TXTPUN(1,101)='G*CM/S/CM               '  ! PMOM_VEC
      TXTPUN(1,102)='G*CM/S/CM               '  ! PMOM_VEC
c  grad PSI
      TXTPUN(1,103)='TESLA                   '  ! PSI
      TXTPUN(1,104)='TESLA                   '  ! PSI
      TXTPUN(1,105)='TESLA                   '  ! PSI

      TXTPUN(1,106)='/CM                     '  ! ZI
      TXTPUN(1,107)='/CM                     '  ! ZI
      TXTPUN(1,108)='/CM                     '  ! ZI

      TXTPUN(1,109)=' ---                    '  ! FREE27
      TXTPUN(1,110)=' ---                    '  ! FREE27
      TXTPUN(1,111)=' ---                    '  ! FREE27
      TXTPUN(1,112)=' ---                    '  ! FREE28
      TXTPUN(1,113)=' ---                    '  ! FREE28
      TXTPUN(1,114)=' ---                    '  ! FREE28
      TXTPUN(1,115)=' ---                    '  ! FREE29
      TXTPUN(1,116)=' ---                    '  ! FREE29
      TXTPUN(1,117)=' ---                    '  ! FREE29
      TXTPUN(1,118)=' ---                    '  ! FREE30
      TXTPUN(1,119)=' ---                    '  ! FREE30
      TXTPUN(1,120)=' ---                    '  ! FREE30
C
      DO J=1,NTALI
c  except for additional input tally: repeat units, for all species in a tally
        IF (J.EQ.12) CYCLE
        TEXT24=TXTPUN(1,J)
        DO I=2,N1MX
          TXTPUN(I,J)=TEXT24
        ENDDO
      ENDDO
C
C  True first tally dimension (physical species) for printout, plotting, etc..
C  needed for weighted averaging across several cells for intensive quantities.
C
C   Possibly reduced storage, if indirect addressing is applied: NPLSV, NPLSTI, etc.
C   is coded via the quite similar array NFRSTP set in setprm_intal.f
C
      NFSTPI(1)=1
      NFSTPI(2)=NPLSI
      NFSTPI(3)=1
      NFSTPI(4)=NPLSI
      NFSTPI(5)=NPLSI
      NFSTPI(6)=NPLSI
      NFSTPI(7)=NPLSI
      NFSTPI(8)=1
      NFSTPI(9)=1
      NFSTPI(10)=1
      NFSTPI(11)=1
      NFSTPI(12)=NAIN    ! use NAIN here, as NAINI is not yet known
      NFSTPI(13)=NPLSI
      NFSTPI(14)=1
cdr  weight window for test particles.
cdr  inconsistent with using N1MX as first dimension
      NFSTPI(15)=NATMI+NMOLI+NIONI  ! + nphoti
      NFSTPI(16)=1
      NFSTPI(17)=1
      NFSTPI(18)=1
      NFSTPI(19)=1
      NFSTPI(20)=1
      NFSTPI(21)=1
      NFSTPI(22)=1
      NFSTPI(23)=3*NPLSI
      NFSTPI(24)=3*NPLSI

      NFSTPI(25)=1

      NFSTPI(26)=NPLSI
      NFSTPI(27)=1
      NFSTPI(28)=1
      NFSTPI(29)=1
      NFSTPI(30)=1

cdr  gradients of input tallies
      NFSTPI(31)=1
      NFSTPI(32)=1
      NFSTPI(33)=1
      NFSTPI(34)=NPLSI
      NFSTPI(35)=NPLSI
      NFSTPI(36)=NPLSI
      NFSTPI(37)=1
      NFSTPI(38)=1
      NFSTPI(39)=1
      NFSTPI(40)=NPLSI
      NFSTPI(41)=NPLSI
      NFSTPI(42)=NPLSI
      NFSTPI(43)=NPLSI
      NFSTPI(44)=NPLSI
      NFSTPI(45)=NPLSI
      NFSTPI(46)=NPLSI
      NFSTPI(47)=NPLSI
      NFSTPI(48)=NPLSI
      NFSTPI(49)=NPLSI
      NFSTPI(50)=NPLSI
      NFSTPI(51)=NPLSI
      NFSTPI(52)=1
      NFSTPI(53)=1
      NFSTPI(54)=1
      NFSTPI(55)=1
      NFSTPI(56)=1
      NFSTPI(57)=1
      NFSTPI(58)=1
      NFSTPI(59)=1
      NFSTPI(60)=1
      NFSTPI(61)=1
      NFSTPI(62)=1
      NFSTPI(63)=1
      NFSTPI(64)=NAIN    ! use NAIN here as NAINI is not yet known
      NFSTPI(65)=NAIN    ! use NAIN here as NAINI is not yet known
      NFSTPI(66)=NAIN    ! use NAIN here as NAINI is not yet known
      NFSTPI(67)=NPLSI
      NFSTPI(68)=NPLSI
      NFSTPI(69)=NPLSI
      NFSTPI(70)=1
      NFSTPI(71)=1
      NFSTPI(72)=1
cdr  weight windows, mostly unused.
cdr  the next 3 lines are inconsistent with N1MX being used as 1st dimension
      NFSTPI(73)=NATMI+NMOLI+NIONI  ! NPHOTI ?
      NFSTPI(74)=NATMI+NMOLI+NIONI
      NFSTPI(75)=NATMI+NMOLI+NIONI
cdr
      NFSTPI(76)=1
      NFSTPI(77)=1
      NFSTPI(78)=1
      NFSTPI(79)=1
      NFSTPI(80)=1
      NFSTPI(81)=1
      NFSTPI(82)=1
      NFSTPI(83)=1
      NFSTPI(84)=1
      NFSTPI(85)=1
      NFSTPI(86)=1
      NFSTPI(87)=1
      NFSTPI(88)=1
      NFSTPI(89)=1
      NFSTPI(90)=1
      NFSTPI(91)=1
      NFSTPI(92)=1
      NFSTPI(93)=1
      NFSTPI(94)=1
      NFSTPI(95)=1
      NFSTPI(96)=1
c  grad BV_VEC: strictly: a 3x3 dyadic now
      NFSTPI(97)=3*NPLSI
      NFSTPI(98)=3*NPLSI
      NFSTPI(99)=3*NPLSI
c  grad PMOM_VEC: strictly: a 3x3 dyadic now
      NFSTPI(100)=3*NPLSI
      NFSTPI(101)=3*NPLSI
      NFSTPI(102)=3*NPLSI
c  grad PSI
      NFSTPI(103)=1
      NFSTPI(104)=1
      NFSTPI(105)=1

      NFSTPI(106)=NPLS
      NFSTPI(107)=NPLS
      NFSTPI(108)=NPLS

cdr  grad of free tallies 27...30
      NFSTPI(109)=1
      NFSTPI(110)=1
      NFSTPI(111)=1
      NFSTPI(112)=1
      NFSTPI(113)=1
      NFSTPI(114)=1
      NFSTPI(115)=1
      NFSTPI(116)=1
      NFSTPI(117)=1
      NFSTPI(118)=1
      NFSTPI(119)=1
      NFSTPI(120)=1
C
      TXTPSP(1,1)='ELECTRONS               '
      TXTPSP(1,3)='ELECTRONS               '
      TXTPSP(1,8)=' ---                    '
      TXTPSP(1,9)=' ---                    '
      TXTPSP(1,10)=' ---                    '
      TXTPSP(1,11)=' ---                    '
C     TXTPSP(IAIN,12)='TO BE READ            '
      TXTPSP(1,14)=' ---                    '
      TXTPSP(1,16)=' ---                    '
      TXTPSP(1,17)=' ---                    '
      TXTPSP(1,18)=' ---                    '
      TXTPSP(1,19)=' ---                    '
      TXTPSP(1,20)=' ---                    '
      TXTPSP(1,21)=' ---                    '
      TXTPSP(1,22)=' ---                    '

      TXTPSP(1,25)=' ---                    '

      TXTPSP(1,26)=' ---                    '
      TXTPSP(1,27)=' ---                    '
      TXTPSP(1,28)=' ---                    '
      TXTPSP(1,29)=' ---                    '
      TXTPSP(1,30)=' ---                    '
c
c  next: gradient vector of 1D input tallies
c   grad Te
      TXTPSP(1,31)='ELECTRONS               '
      TXTPSP(1,32)='ELECTRONS               '
      TXTPSP(1,33)='ELECTRONS               '
c  grad ne
      TXTPSP(1,37)='ELECTRONS               '
      TXTPSP(1,38)='ELECTRONS               '
      TXTPSP(1,39)='ELECTRONS               '

      TXTPSP(1,52)=' ---                    '
      TXTPSP(1,53)=' ---                    '
      TXTPSP(1,54)=' ---                    '
      TXTPSP(1,55)=' ---                    '
      TXTPSP(1,56)=' ---                    '
      TXTPSP(1,57)=' ---                    '
      TXTPSP(1,58)=' ---                    '
      TXTPSP(1,59)=' ---                    '
      TXTPSP(1,60)=' ---                    '
      TXTPSP(1,61)=' ---                    '
      TXTPSP(1,62)=' ---                    '
      TXTPSP(1,63)=' ---                    '
c  grad ADIN
C     TXTPSP(IAIN,64)='TO BE READ            '
C     TXTPSP(IAIN,65)='TO BE READ            '
C     TXTPSP(IAIN,66)='TO BE READ            '


      TXTPSP(1,70)=' ---                    '
      TXTPSP(1,71)=' ---                    '
      TXTPSP(1,72)=' ---                    '
      TXTPSP(1,76)=' ---                    '
      TXTPSP(1,77)=' ---                    '
      TXTPSP(1,78)=' ---                    '
      TXTPSP(1,79)=' ---                    '
      TXTPSP(1,80)=' ---                    '
      TXTPSP(1,81)=' ---                    '
      TXTPSP(1,82)=' ---                    '
      TXTPSP(1,83)=' ---                    '
      TXTPSP(1,84)=' ---                    '
      TXTPSP(1,85)=' ---                    '
      TXTPSP(1,86)=' ---                    '
      TXTPSP(1,87)=' ---                    '
      TXTPSP(1,88)=' ---                    '
      TXTPSP(1,89)=' ---                    '
      TXTPSP(1,90)=' ---                    '
      TXTPSP(1,91)=' ---                    '
      TXTPSP(1,92)=' ---                    '
      TXTPSP(1,93)=' ---                    '
      TXTPSP(1,94)=' ---                    '
      TXTPSP(1,95)=' ---                    '
      TXTPSP(1,96)=' ---                    '

      TXTPSP(1,103)=' ---                    '
      TXTPSP(1,104)=' ---                    '
      TXTPSP(1,105)=' ---                    '
      TXTPSP(1,106)=' ---                    '
      TXTPSP(1,107)=' ---                    '
      TXTPSP(1,108)=' ---                    '
      TXTPSP(1,109)=' ---                    '
      TXTPSP(1,110)=' ---                    '
      TXTPSP(1,111)=' ---                    '
      TXTPSP(1,112)=' ---                    '
      TXTPSP(1,113)=' ---                    '
      TXTPSP(1,114)=' ---                    '
      TXTPSP(1,115)=' ---                    '
      TXTPSP(1,116)=' ---                    '
      TXTPSP(1,117)=' ---                    '
      TXTPSP(1,118)=' ---                    '
      TXTPSP(1,119)=' ---                    '
      TXTPSP(1,120)=' ---                    '
c
c  next: 2D input tallies, = tallies with an additional  "species index"
C
cdr  txtpsp is dimensioned with N1MX
cdr  This may be 3*NPLS due to vectorial tallies, but still be smaller than NSPAMI
cdr  Most likely "weight-window tallies" 15, 73,74,75 do not exist in any run
cdr  anyway. But still corrections are needed here.
      DO 50 ISPZ=1,NSPAMI
        if (ispz.gt.n1mx) exit
        TXTPSP(ISPZ,15)=TEXTS(ISPZ)
        TXTPSP(ISPZ,73)=TEXTS(ISPZ)
        TXTPSP(ISPZ,74)=TEXTS(ISPZ)
        TXTPSP(ISPZ,75)=TEXTS(ISPZ)
   50 CONTINUE
C
      DO 80 IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
        TXTPSP(IPLS,2)=TEXTS(ISPZ)
        TXTPSP(IPLS,4)=TEXTS(ISPZ)
        TXTPSP(IPLS,5)=TEXTS(ISPZ)
        TXTPSP(IPLS,6)=TEXTS(ISPZ)
        TXTPSP(IPLS,7)=TEXTS(ISPZ)
        TXTPSP(IPLS,13)=TEXTS(ISPZ)
        DO ICO=0,2
c  grad bv_vec. ICO=0  parallel component in older code versions
          KK=ICO*NPLSI
          TXTPSP(KK+IPLS,23)=TEXTS(ISPZ)
c  grad pmom_vec. ICO=0  parallel component in older code versions
          KK=ICO*NPLSI
          TXTPSP(KK+IPLS,24)=TEXTS(ISPZ)
        ENDDO

c  dTi/dx, dTi/dy, dTi/dz
        TXTPSP(IPLS,34)=TEXTS(ISPZ)
        TXTPSP(IPLS,35)=TEXTS(ISPZ)
        TXTPSP(IPLS,36)=TEXTS(ISPZ)
c grad ni
        TXTPSP(IPLS,40)=TEXTS(ISPZ)
        TXTPSP(IPLS,41)=TEXTS(ISPZ)
        TXTPSP(IPLS,42)=TEXTS(ISPZ)
c grad vx
        TXTPSP(IPLS,43)=TEXTS(ISPZ)
        TXTPSP(IPLS,44)=TEXTS(ISPZ)
        TXTPSP(IPLS,45)=TEXTS(ISPZ)
c grad vy
        TXTPSP(IPLS,46)=TEXTS(ISPZ)
        TXTPSP(IPLS,47)=TEXTS(ISPZ)
        TXTPSP(IPLS,48)=TEXTS(ISPZ)
c grad vz
        TXTPSP(IPLS,49)=TEXTS(ISPZ)
        TXTPSP(IPLS,50)=TEXTS(ISPZ)
        TXTPSP(IPLS,51)=TEXTS(ISPZ)

        TXTPSP(IPLS,67)=TEXTS(ISPZ)
        TXTPSP(IPLS,68)=TEXTS(ISPZ)
        TXTPSP(IPLS,69)=TEXTS(ISPZ)
        DO ICO=0,2
c  grad BV_VEC
          KK=ICO*NPLSI
          TXTPSP(KK+IPLS,97)=TEXTS(ISPZ)
          TXTPSP(KK+IPLS,98)=TEXTS(ISPZ)
          TXTPSP(KK+IPLS,99)=TEXTS(ISPZ)
c  grad PMOM_VEC
          KK=ICO*NPLSI
          TXTPSP(KK+IPLS,100)=TEXTS(ISPZ)
          TXTPSP(KK+IPLS,101)=TEXTS(ISPZ)
          TXTPSP(KK+IPLS,102)=TEXTS(ISPZ)
        ENDDO
   80 CONTINUE
C
      RETURN
      END SUBROUTINE EIRENE_SETTXT_INTAL
