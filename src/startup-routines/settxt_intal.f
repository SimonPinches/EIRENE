cdr  Aug 18: formerly this was part of settxt.f.
cdr  Generalisation of input tallies (and gradients thereof)


      SUBROUTINE EIRENE_SETTXT_INTAL
c  Set default texts  (volume tallies: name, species, units), 
C  similar to SETTXT.f, but for INPUT TALLIES rather than output tallies. 
C  Set first (leading) dimension of tally arrays: nfstpi.

c  
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU
 
      IMPLICIT NONE
 
      INTEGER :: IATM, IION, IPLS, IMOL, ISPZ, IPHOT, I, J
      CHARACTER(24) :: TEXT24
      CHARACTER(72) :: TEXT72
C
C  TEXT FOR INPUT (BACKGROUND) TALLIES

cdr  the numbering is "a bit" illogical, due to historic reasons.
c    primary and derived tallies are mixed here. 
C
      TXTPLS(1,1)='PLASMA TEMPERATURE                               '
      TXTPLS(1,2)='PLASMA TEMPERATURE                               '
      TXTPLS(1,3)='PLASMA DENSITY (BULK PARTICLES)                  '
      TXTPLS(1,4)='PLASMA DENSITY (BULK PARTICLES)                  '
      TXTPLS(1,5)='DRIFT VELOCITY IN X-DIRECTION (BULK IONS)        '
      TXTPLS(1,6)='DRIFT VELOCITY IN Y-DIRECTION (BULK IONS)        '
      TXTPLS(1,7)='DRIFT VELOCITY IN Z-DIRECTION (BULK IONS)        '
      TXTPLS(1,8)='MAGN. FIELD UNIT VECTOR, X DIRECTION             '
      TXTPLS(1,9)='MAGN. FIELD UNIT VECTOR, Y DIRECTION             '
      TXTPLS(1,10)='MAGN. FIELD UNIT VECTOR, Z DIRECTION             '
      TXTPLS(1,11)='MAGN. FIELD STRENGTH                             '
cdr to be added here  TXTPLS(1,xx)='Magn. POTENTIAL, e.g. PSI fct.   '
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
      TXTPLS(1,23)='FLOW VELOCITY PARALLEL TO BFIELD                 '
      TXTPLS(1,24)='PARALLEL TO B MOMENTUM FLOW                      '

C  gradient tallies:  derivatives wrt. cartesian coordinates
      TXTPLS(1,25)='dTE/dX                                           '
      TXTPLS(1,26)='dTE/dY                                           '
      TXTPLS(1,27)='dTE/dZ                                           '
      TXTPLS(1,28)='dTI/dX                                           '
      TXTPLS(1,29)='dTI/dY                                           '
      TXTPLS(1,30)='dTI/dZ                                           '
      TXTPLS(1,31)='dNE/dX                                           '
      TXTPLS(1,32)='dNE/dY                                           '
      TXTPLS(1,33)='dNE/dZ                                           '
      TXTPLS(1,34)='dNI/dX                                           '
      TXTPLS(1,35)='dNI/dY                                           '
      TXTPLS(1,36)='dNI/dZ                                           '
      TXTPLS(1,37)='dVX/dX                                           '
      TXTPLS(1,38)='dVX/dY                                           '
      TXTPLS(1,39)='dVX/dZ                                           '
      TXTPLS(1,40)='dVY/dX                                           '
      TXTPLS(1,41)='dVY/dY                                           '
      TXTPLS(1,42)='dVY/dZ                                           '
      TXTPLS(1,43)='dVZ/dX                                           '
      TXTPLS(1,44)='dVZ/dY                                           '
      TXTPLS(1,45)='dVZ/dZ                                           '
      TXTPLS(1,46)='dBX/dX                                           '
      TXTPLS(1,47)='dBX/dY                                           '
      TXTPLS(1,48)='dBX/dZ                                           '
      TXTPLS(1,49)='dBY/dX                                           '
      TXTPLS(1,50)='dBY/dY                                           '
      TXTPLS(1,51)='dBY/dZ                                           '
      TXTPLS(1,52)='dBZ/dX                                           '
      TXTPLS(1,53)='dBZ/dY                                           '
      TXTPLS(1,54)='dBZ/dZ                                           '
      TXTPLS(1,55)='dBF/dX                                           '
      TXTPLS(1,56)='dBF/dY                                           '
      TXTPLS(1,57)='dBF/dZ                                           '
      TXTPLS(1,58)='dADIN/dX                                         '
cdr  psi funct. gradient to be added here ?
      TXTPLS(1,59)='dADIN/dY                                         '
      TXTPLS(1,60)='dADIN/dZ                                         '
      TXTPLS(1,61)='dEDRIFT/dX                                       '
      TXTPLS(1,62)='dEDRIFT/dY                                       '
      TXTPLS(1,63)='dEDRIFT/dZ                                       '
      TXTPLS(1,64)='dVOL/dX                                          '
      TXTPLS(1,65)='dVOL/dY                                          '
      TXTPLS(1,66)='dVOL/dZ                                          '
      TXTPLS(1,67)='dWGHT/dX                                         '
      TXTPLS(1,68)='dWGHT/dY                                         '
      TXTPLS(1,69)='dWGHT/dZ                                         '
      TXTPLS(1,70)='dBXPERP/dX                                       '
      TXTPLS(1,71)='dBXPERP/dY                                       '
      TXTPLS(1,72)='dBXPERP/dZ                                       '
      TXTPLS(1,73)='dBYPERP/dX                                       '
      TXTPLS(1,74)='dBYPERP/dY                                       '
      TXTPLS(1,75)='dBYPERP/dZ                                       '
      TXTPLS(1,76)='dEX/dX                                           '
      TXTPLS(1,77)='dEX/dY                                           '
      TXTPLS(1,78)='dEX/dZ                                           '
      TXTPLS(1,79)='dEY/dX                                           '
      TXTPLS(1,80)='dEY/dY                                           '
      TXTPLS(1,81)='dEY/dZ                                           '
      TXTPLS(1,82)='dEZ/dX                                           '
      TXTPLS(1,83)='dEZ/dY                                           '
      TXTPLS(1,84)='dEZ/dZ                                           '
      TXTPLS(1,85)='dEF/dX                                           '
      TXTPLS(1,86)='dEF/dY                                           '
      TXTPLS(1,87)='dEF/dZ                                           '
      TXTPLS(1,88)='dPOT/dX                                          '
      TXTPLS(1,89)='dPOT/dY                                          '
      TXTPLS(1,90)='dPOT/dZ                                          '
      TXTPLS(1,91)='dBVIN/dX                                          '
      TXTPLS(1,92)='dBVIN/dY                                          '
      TXTPLS(1,93)='dBVIN/dZ                                          '
      TXTPLS(1,94)='dPARMOM/dX                                        '
      TXTPLS(1,95)='dPARMOM/dY                                        '
      TXTPLS(1,96)='dPARMOM/dZ                                        '
C
c  currently: ntali=4*24=96
      DO J=1,NTALI
        IF (J.NE.12) THEN
          DO I=2,N1MX
            TEXT72=TXTPLS(1,J)
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
      TXTPUN(1,12)='ADDITIONAL TALLY UNITS  '
      TXTPUN(1,13)='EV                      '  ! EDRIFT  --> DERIVED QUANTITY
      TXTPUN(1,14)='CM**3                   '  ! VOL
      TXTPUN(1,15)=' ---                    '  ! WEIGHT WINDOW
      TXTPUN(1,16)=' ---                    '  ! BX_PERP --> DERIVED QUANTITY 
      TXTPUN(1,17)=' ---                    '  ! BY_PERP --> DERIVED QUANTITY
      TXTPUN(1,18)=' ---                    '
      TXTPUN(1,19)=' ---                    '
      TXTPUN(1,20)=' ---                    '
      TXTPUN(1,21)='V/CM                    '  ! EF
      TXTPUN(1,22)='V                       '  ! POT
      TXTPUN(1,23)='CM/S                    '  ! BVIN
      TXTPUN(1,24)='G*CM/S                  '  ! PARMOM

cdr derivaties in cart. coordinates  (gradient vector)  
      TXTPUN(1,25)='EV/CM                   '  ! grad(Te) 
      TXTPUN(1,26)='EV/CM                   '
      TXTPUN(1,27)='EV/CM                   '
      TXTPUN(1,28)='EV/CM                   '  ! grad(Ti)
      TXTPUN(1,29)='EV/CM                   '
      TXTPUN(1,30)='EV/CM                   '
      TXTPUN(1,31)='CM**-3/CM               '  ! grad(ne)
      TXTPUN(1,32)='CM**-3/CM               '
      TXTPUN(1,33)='CM**-3/CM               '
      TXTPUN(1,34)='CM**-3/CM               '  ! grad(ni)
      TXTPUN(1,35)='CM**-3/CM               '
      TXTPUN(1,36)='CM**-3/CM               '
      TXTPUN(1,37)='CM/SEC/CM               '
      TXTPUN(1,38)='CM/SEC/CM               '
      TXTPUN(1,39)='CM/SEC/CM               '
      TXTPUN(1,40)='CM/SEC/CM               '
      TXTPUN(1,41)='CM/SEC/CM               '
      TXTPUN(1,42)='CM/SEC/CM               '
      TXTPUN(1,43)='CM/SEC/CM               '
      TXTPUN(1,44)='CM/SEC/CM               '
      TXTPUN(1,45)='CM/SEC/CM               '
      TXTPUN(1,46)=' ---                    '
      TXTPUN(1,47)=' ---                    '
      TXTPUN(1,48)=' ---                    '
      TXTPUN(1,49)=' ---                    '
      TXTPUN(1,50)=' ---                    '
      TXTPUN(1,51)=' ---                    '
      TXTPUN(1,52)=' ---                    '
      TXTPUN(1,53)=' ---                    '
      TXTPUN(1,54)=' ---                    '
      TXTPUN(1,55)='TESLA/CM                '
      TXTPUN(1,56)='TESLA/CM                '
      TXTPUN(1,57)='TESLA/CM                '
C     TXTPUN(1,58)='TO BE READ, ADIN        '
C     TXTPUN(1,59)='TO BE READ, ADIN        '
C     TXTPUN(1,60)='TO BE READ, ADIN        '
      TXTPUN(1,61)='EV/CM                   '  ! EDRIFT  --> DERIVED QUANTITY
      TXTPUN(1,62)='EV/CM                   '  ! EDRIFT  --> DERIVED QUANTITY
      TXTPUN(1,63)='EV/CM                   '  ! EDRIFT  --> DERIVED QUANTITY
      TXTPUN(1,64)='CM**3/CM                '  ! VOL
      TXTPUN(1,65)='CM**3/CM                '  ! VOL
      TXTPUN(1,66)='CM**3/CM                '  ! VOL
      TXTPUN(1,67)=' ---                    '  ! WEIGHT WINDOW
      TXTPUN(1,68)=' ---                    '  ! WEIGHT WINDOW
      TXTPUN(1,69)=' ---                    '  ! WEIGHT WINDOW
      TXTPUN(1,70)=' ---                    '  ! BX_PERP --> DERIVED QUANTITY 
      TXTPUN(1,71)=' ---                    '  ! BX_PERP --> DERIVED QUANTITY 
      TXTPUN(1,72)=' ---                    '  ! BX_PERP --> DERIVED QUANTITY 
      TXTPUN(1,73)=' ---                    '  ! BY_PERP --> DERIVED QUANTITY
      TXTPUN(1,74)=' ---                    '  ! BY_PERP --> DERIVED QUANTITY
      TXTPUN(1,75)=' ---                    '  ! BY_PERP --> DERIVED QUANTITY
      TXTPUN(1,76)=' ---                    '
      TXTPUN(1,77)=' ---                    '
      TXTPUN(1,78)=' ---                    '
      TXTPUN(1,79)=' ---                    '
      TXTPUN(1,80)=' ---                    '
      TXTPUN(1,81)=' ---                    '
      TXTPUN(1,82)=' ---                    '
      TXTPUN(1,83)=' ---                    '
      TXTPUN(1,84)=' ---                    '
      TXTPUN(1,85)='V/CM/CM                 '  ! EF
      TXTPUN(1,86)='V/CM/CM                 '  ! EF
      TXTPUN(1,87)='V/CM/CM                 '  ! EF
      TXTPUN(1,88)='V/CM                    '  ! POT
      TXTPUN(1,89)='V/CM                    '  ! POT
      TXTPUN(1,90)='V/CM                    '  ! POT
      TXTPUN(1,91)='CM/S/CM                 '  ! BVIN
      TXTPUN(1,92)='CM/S/CM                 '  ! BVIN
      TXTPUN(1,93)='CM/S/CM                 '  ! BVIN
      TXTPUN(1,94)='G*CM/S/CM               '  ! PARMOM
      TXTPUN(1,95)='G*CM/S/CM               '  ! PARMOM
      TXTPUN(1,96)='G*CM/S/CM               '  ! PARMOM
C
      DO J=1,NTALI
        IF (J.NE.12) THEN
          DO I=2,N1MX
            TEXT24=TXTPUN(1,J)
            TXTPUN(I,J)=TEXT24
          ENDDO
        ENDIF
      ENDDO
C
!      ENTRY EIRENE_STTXT1_INTAL
C
      NFSTPI(1)=1
      NFSTPI(2)=NPLSTI
      NFSTPI(3)=1
      NFSTPI(4)=NPLSI
      NFSTPI(5)=NPLSV
      NFSTPI(6)=NPLSV
      NFSTPI(7)=NPLSV
      NFSTPI(8)=1
      NFSTPI(9)=1
      NFSTPI(10)=1
      NFSTPI(11)=1
      NFSTPI(12)=NAIN    ! use NAIN here, as NAINI is not yet known
      NFSTPI(13)=NPLSI
      NFSTPI(14)=1
      NFSTPI(15)=NATMI+NMOLI+NIONI
      NFSTPI(16)=1
      NFSTPI(17)=1
      NFSTPI(18)=1
      NFSTPI(19)=1
      NFSTPI(20)=1
      NFSTPI(21)=1
      NFSTPI(22)=1
      NFSTPI(23)=NPLSV
      NFSTPI(24)=NPLS

      NFSTPI(25)=1
      NFSTPI(26)=1
      NFSTPI(27)=1
      NFSTPI(28)=NPLSTI
      NFSTPI(29)=NPLSTI
      NFSTPI(30)=NPLSTI
      NFSTPI(31)=1
      NFSTPI(32)=1
      NFSTPI(33)=1
      NFSTPI(34)=NPLSI
      NFSTPI(35)=NPLSI
      NFSTPI(36)=NPLSI
      NFSTPI(37)=NPLSV
      NFSTPI(38)=NPLSV
      NFSTPI(39)=NPLSV
      NFSTPI(40)=NPLSV
      NFSTPI(41)=NPLSV
      NFSTPI(42)=NPLSV
      NFSTPI(43)=NPLSV
      NFSTPI(44)=NPLSV
      NFSTPI(45)=NPLSV
      NFSTPI(46)=1
      NFSTPI(47)=1
      NFSTPI(48)=1
      NFSTPI(49)=1
      NFSTPI(50)=1
      NFSTPI(51)=1
      NFSTPI(52)=1
      NFSTPI(53)=1
      NFSTPI(54)=1
      NFSTPI(55)=1
      NFSTPI(56)=1
      NFSTPI(57)=1
      NFSTPI(58)=NAIN    ! use NAIN here as NAINI is not yet known
      NFSTPI(59)=NAIN    ! use NAIN here as NAINI is not yet known
      NFSTPI(60)=NAIN    ! use NAIN here as NAINI is not yet known
      NFSTPI(61)=NPLSI
      NFSTPI(62)=NPLSI
      NFSTPI(63)=NPLSI
      NFSTPI(64)=1
      NFSTPI(65)=1
      NFSTPI(66)=1
      NFSTPI(67)=NATMI+NMOLI+NIONI
      NFSTPI(68)=NATMI+NMOLI+NIONI
      NFSTPI(69)=NATMI+NMOLI+NIONI
      NFSTPI(70)=1
      NFSTPI(71)=1
      NFSTPI(72)=1
      NFSTPI(73)=1
      NFSTPI(74)=1
      NFSTPI(75)=1
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
      NFSTPI(91)=NPLSV
      NFSTPI(92)=NPLSV
      NFSTPI(93)=NPLSV
      NFSTPI(94)=NPLS
      NFSTPI(95)=NPLS
      NFSTPI(96)=NPLS
C
      TXTPSP(1,1)='ELECTRONS               '
      TXTPSP(1,3)='ELECTRONS               '
      TXTPSP(1,8)=' ---                    '
      TXTPSP(1,9)=' ---                    '
      TXTPSP(1,10)=' ---                    '
      TXTPSP(1,11)=' ---                    '
      TXTPSP(1,14)=' ---                    '
      TXTPSP(1,16)=' ---                    '
      TXTPSP(1,17)=' ---                    '
      TXTPSP(1,18)=' ---                    '
      TXTPSP(1,19)=' ---                    '
      TXTPSP(1,20)=' ---                    '
      TXTPSP(1,21)=' ---                    '
      TXTPSP(1,22)=' ---                    '
      TXTPSP(1,25)='ELECTRONS               '
      TXTPSP(1,26)='ELECTRONS               '
      TXTPSP(1,27)='ELECTRONS               '
      TXTPSP(1,31)='ELECTRONS               '
      TXTPSP(1,32)='ELECTRONS               '
      TXTPSP(1,33)='ELECTRONS               '
      TXTPSP(1,46)=' ---                    '
      TXTPSP(1,47)=' ---                    '
      TXTPSP(1,48)=' ---                    '
      TXTPSP(1,49)=' ---                    '
      TXTPSP(1,50)=' ---                    '
      TXTPSP(1,51)=' ---                    '
      TXTPSP(1,52)=' ---                    '
      TXTPSP(1,53)=' ---                    '
      TXTPSP(1,54)=' ---                    '
      TXTPSP(1,55)=' ---                    '
      TXTPSP(1,56)=' ---                    '
      TXTPSP(1,57)=' ---                    '
      TXTPSP(1,64)=' ---                    '
      TXTPSP(1,65)=' ---                    '
      TXTPSP(1,66)=' ---                    '
      TXTPSP(1,70)=' ---                    '
      TXTPSP(1,71)=' ---                    '
      TXTPSP(1,72)=' ---                    '
      TXTPSP(1,73)=' ---                    '
      TXTPSP(1,74)=' ---                    '
      TXTPSP(1,75)=' ---                    '
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
C
C     TXTPSP(IAIN,12)='TO BE READ            '
C     TXTPSP(IAIN,58)='TO BE READ            '
C     TXTPSP(IAIN,59)='TO BE READ            '
C     TXTPSP(IAIN,60)='TO BE READ            '
C
      DO 50 ISPZ=1,NSPAMI
        TXTPSP(ISPZ,15)=TEXTS(ISPZ)
        TXTPSP(ISPZ,67)=TEXTS(ISPZ)
        TXTPSP(ISPZ,68)=TEXTS(ISPZ)
        TXTPSP(ISPZ,69)=TEXTS(ISPZ)
 50   CONTINUE
C
      DO 80 IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
        TXTPSP(IPLS,2)=TEXTS(ISPZ)
        TXTPSP(IPLS,4)=TEXTS(ISPZ)
        TXTPSP(IPLS,5)=TEXTS(ISPZ)
        TXTPSP(IPLS,6)=TEXTS(ISPZ)
        TXTPSP(IPLS,7)=TEXTS(ISPZ)
        TXTPSP(IPLS,13)=TEXTS(ISPZ)
        TXTPSP(IPLS,23)=TEXTS(ISPZ)
        TXTPSP(IPLS,24)=TEXTS(ISPZ)
        TXTPSP(IPLS,28)=TEXTS(ISPZ)
        TXTPSP(IPLS,29)=TEXTS(ISPZ)
        TXTPSP(IPLS,30)=TEXTS(ISPZ)
        TXTPSP(IPLS,34)=TEXTS(ISPZ)
        TXTPSP(IPLS,35)=TEXTS(ISPZ)
        TXTPSP(IPLS,36)=TEXTS(ISPZ)
        TXTPSP(IPLS,37)=TEXTS(ISPZ)
        TXTPSP(IPLS,38)=TEXTS(ISPZ)
        TXTPSP(IPLS,39)=TEXTS(ISPZ)
        TXTPSP(IPLS,40)=TEXTS(ISPZ)
        TXTPSP(IPLS,41)=TEXTS(ISPZ)
        TXTPSP(IPLS,42)=TEXTS(ISPZ)
        TXTPSP(IPLS,43)=TEXTS(ISPZ)
        TXTPSP(IPLS,44)=TEXTS(ISPZ)
        TXTPSP(IPLS,45)=TEXTS(ISPZ)
        TXTPSP(IPLS,61)=TEXTS(ISPZ)
        TXTPSP(IPLS,62)=TEXTS(ISPZ)
        TXTPSP(IPLS,63)=TEXTS(ISPZ)
        TXTPSP(IPLS,91)=TEXTS(ISPZ)
        TXTPSP(IPLS,92)=TEXTS(ISPZ)
        TXTPSP(IPLS,93)=TEXTS(ISPZ)
        TXTPSP(IPLS,94)=TEXTS(ISPZ)
        TXTPSP(IPLS,95)=TEXTS(ISPZ)
        TXTPSP(IPLS,96)=TEXTS(ISPZ)
 80     CONTINUE
C
      RETURN
      END
 
 
 
 
 
 
 
