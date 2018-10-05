      MODULE EIRMOD_CINIT
 
!  sep-05: specifications for databases added, ndbnames, dbhandle, dbfname
!  jul-06: database handle added for ADAS
!  20.06.07: deallocate tdmpar
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
 
      IMPLICIT NONE
 
      PRIVATE
 
      PUBLIC :: EIRENE_ALLOC_CINIT, EIRENE_DEALLOC_CINIT, 
     .          EIRENE_INIT_CINIT,
     .          TDENMODEL, TDENMODAR
 
      INTEGER, PUBLIC, SAVE ::
     I         NCINIT, MCINIT, LCINIT
 
      REAL(DP),  PUBLIC, TARGET, ALLOCATABLE, SAVE :: RCINIT(:)
      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE :: ICINIT(:)
      LOGICAL, PUBLIC, TARGET, SAVE :: LCNIT(4)
 
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R EP1IN,  EP1CH,  EP1OT,  EXEP1,
     R ELLIN,  ELLOT,  ELLCH,  EXELL,
     R TRIIN,  TRIOT,  TRICH,  EXTRI,
     R TE0,    TE1,    TE2,    TE3,    TE4,   TE5,
     R BP0,    BP1,    BP2,    BP3,    BP4,   BP5,
     R B0,     B1,     B2,     B3,     B4,    B5,
     R VL0,    VL1,    VL2,    VL3,    VL4,   VL5,
     R EF0,    EF1,    EF2,    EF3,    EF4,   EF5
     
 
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R TI0(:),TI1(:),TI2(:),TI3(:),TI4(:),TI5(:),
     R VX0(:),VX1(:),VX2(:),VX3(:),VX4(:),VX5(:),
     R VY0(:),VY1(:),VY2(:),VY3(:),VY4(:),VY5(:),
     R VZ0(:),VZ1(:),VZ2(:),VZ3(:),VZ4(:),VZ5(:),
     R DI0(:),DI1(:),DI2(:),DI3(:),DI4(:),DI5(:)
 
      INTEGER, PUBLIC, POINTER, SAVE ::
     I INDPRO(:),INDGRD(:),INDSRC(:)
 
      LOGICAL, PUBLIC, POINTER, SAVE ::
     L NLMACH,NLMLTI,NLMLV,NLPITCH
 
      CHARACTER(66), PUBLIC, SAVE :: CASENAME
      CHARACTER(10), PUBLIC, ALLOCATABLE, SAVE :: CDENMODEL(:)
 
      TYPE TDENMODEL
        REAL(DP) :: G_BOLTZ, DELTAE, A_CORONA, DVAL, TVAL,
     R              VXVAL, VYVAl, VZVAL, DFACTOR, TFACTOR,
     R              VFACTOR
        INTEGER :: NRE
        INTEGER, POINTER :: ISP(:), ITP(:), ISTR(:)
        CHARACTER(8), POINTER :: FNAME(:)
        CHARACTER(4), POINTER :: H123(:)
        CHARACTER(9), POINTER :: REACTION(:)
        CHARACTER(3), POINTER :: CR(:)
      END TYPE
 
      TYPE TDENMODAR
        TYPE(TDENMODEL), POINTER :: TDM
      END TYPE
 
      TYPE(TDENMODAR), PUBLIC, ALLOCATABLE, SAVE :: TDMPAR(:)
 
      INTEGER, PUBLIC, PARAMETER ::
     I NDBNAMES=17
 
      CHARACTER(400), PUBLIC, SAVE :: DBFNAME(NDBNAMES)
      CHARACTER(6), PUBLIC, SAVE :: DBHANDLE(NDBNAMES)
 
      CONTAINS
 
      SUBROUTINE EIRENE_ALLOC_CINIT
 
      IF (ALLOCATED(RCINIT)) RETURN
 
      NCINIT=42+30*NPLS
      MCINIT=12+3+NSTRA
      LCINIT=3
 
      ALLOCATE (RCINIT(NCINIT))
      ALLOCATE (ICINIT(MCINIT))
      ALLOCATE (CDENMODEL(NPLS))
      ALLOCATE (TDMPAR(NPLS))
 
      WRITE (55+IFOFF,'(A,T25,I15)')
     .      ' CINIT ',NCINIT*8 + MCINIT*4 +
     .                NPLS*LEN(CDENMODEL(1))
 
      EP1IN  => RCINIT(1)
      EP1CH  => RCINIT(2)
      EP1OT  => RCINIT(3)
      EXEP1  => RCINIT(4)
      ELLIN  => RCINIT(5)
      ELLOT  => RCINIT(6)
      ELLCH  => RCINIT(7)
      EXELL  => RCINIT(8)
      TRIIN  => RCINIT(9)
      TRIOT  => RCINIT(10)
      TRICH  => RCINIT(11)
      EXTRI  => RCINIT(12)
      TE0    => RCINIT(13)
      TE1    => RCINIT(14)
      TE2    => RCINIT(15)
      TE3    => RCINIT(16)
      TE4    => RCINIT(17)
      TE5    => RCINIT(18)
      BP0    => RCINIT(19)
      BP1    => RCINIT(20)
      BP2    => RCINIT(21)
      BP3    => RCINIT(22)
      BP4    => RCINIT(23)
      BP5    => RCINIT(24)
      B0     => RCINIT(25)
      B1     => RCINIT(26)
      B2     => RCINIT(27)
      B3     => RCINIT(28)
      B4     => RCINIT(29)
      B5     => RCINIT(30)
      VL0    => RCINIT(31)
      VL1    => RCINIT(32)
      VL2    => RCINIT(33)
      VL3    => RCINIT(34)
      VL4    => RCINIT(35)
      VL5    => RCINIT(36)
      EF0    => RCINIT(37)
      EF1    => RCINIT(38)
      EF2    => RCINIT(39)
      EF3    => RCINIT(40)
      EF4    => RCINIT(41)
      EF5    => RCINIT(42)
      TI0    => RCINIT(43         : 42+ 1*NPLS)
      TI1    => RCINIT(43+ 1*npls : 42+ 2*NPLS)
      TI2    => RCINIT(43+ 2*npls : 42+ 3*NPLS)
      TI3    => RCINIT(43+ 3*npls : 42+ 4*NPLS)
      TI4    => RCINIT(43+ 4*npls : 42+ 5*NPLS)
      TI5    => RCINIT(43+ 5*npls : 42+ 6*NPLS)
      VX0    => RCINIT(43+ 6*npls : 42+ 7*NPLS)
      VX1    => RCINIT(43+ 7*npls : 42+ 8*NPLS)
      VX2    => RCINIT(43+ 8*npls : 42+ 9*NPLS)
      VX3    => RCINIT(43+ 9*npls : 42+10*NPLS)
      VX4    => RCINIT(43+10*npls : 42+11*NPLS)
      VX5    => RCINIT(43+11*npls : 42+12*NPLS)
      VY0    => RCINIT(43+12*npls : 42+13*NPLS)
      VY1    => RCINIT(43+13*npls : 42+14*NPLS)
      VY2    => RCINIT(43+14*npls : 42+15*NPLS)
      VY3    => RCINIT(43+15*npls : 42+16*NPLS)
      VY4    => RCINIT(43+16*npls : 42+17*NPLS)
      VY5    => RCINIT(43+17*npls : 42+18*NPLS)
      VZ0    => RCINIT(43+18*npls : 42+19*NPLS)
      VZ1    => RCINIT(43+19*npls : 42+20*NPLS)
      VZ2    => RCINIT(43+20*npls : 42+21*NPLS)
      VZ3    => RCINIT(43+21*npls : 42+22*NPLS)
      VZ4    => RCINIT(43+22*npls : 42+23*NPLS)
      VZ5    => RCINIT(43+23*npls : 42+24*NPLS)
      DI0    => RCINIT(43+24*npls : 42+25*NPLS)
      DI1    => RCINIT(43+25*npls : 42+26*NPLS)
      DI2    => RCINIT(43+26*npls : 42+27*NPLS)
      DI3    => RCINIT(43+27*npls : 42+28*NPLS)
      DI4    => RCINIT(43+28*npls : 42+29*NPLS)
      DI5    => RCINIT(43+29*npls : 42+30*NPLS)
 
      INDPRO => ICINIT( 1 : 12)
      INDGRD => ICINIT(13 : 15)
      INDSRC => ICINIT(16 : 15+NSTRA)
 
      NLMACH => LCNIT(1)
      NLMLTI => LCNIT(2)
      NLMLV  => LCNIT(3)
      NLPITCH=> LCNIT(4)
 
      CALL EIRENE_INIT_CINIT
 
      RETURN
 
      END SUBROUTINE EIRENE_ALLOC_CINIT
 
 
      SUBROUTINE EIRENE_DEALLOC_CINIT
      INTEGER :: I
 
      IF (.NOT.ALLOCATED(RCINIT)) RETURN
 
      DEALLOCATE (RCINIT)
      DEALLOCATE (ICINIT)
      DEALLOCATE (CDENMODEL)
 
      DO I=1,NPLS
        IF (ASSOCIATED(TDMPAR(I)%TDM)) THEN
          DEALLOCATE (TDMPAR(I)%TDM%ISP)
          DEALLOCATE (TDMPAR(I)%TDM%ITP)
          DEALLOCATE (TDMPAR(I)%TDM%ISTR)
          DEALLOCATE (TDMPAR(I)%TDM%FNAME)
          DEALLOCATE (TDMPAR(I)%TDM%H123)
          DEALLOCATE (TDMPAR(I)%TDM%REACTION)
          DEALLOCATE (TDMPAR(I)%TDM%CR)
          DEALLOCATE (TDMPAR(I)%TDM)
        END IF
      END DO
      DEALLOCATE (TDMPAR)
 
      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CINIT
 
 
      SUBROUTINE EIRENE_INIT_CINIT

      INTEGER :: I
 
      RCINIT = 0._DP
      ICINIT = 0
      LCNIT  = .FALSE.
      CDENMODEL = REPEAT(' ',LEN(CDENMODEL))
      CASENAME = REPEAT(' ',LEN(CASENAME))

      DO I = 1, NPLS
         NULLIFY(TDMPAR(I)%TDM)
      END DO
C
C  INITIALIZE GEOMETRY DATA
C
      ELLIN=1._DP
      ELLOT=1._DP
      ELLCH=1._DP
      TRIIN=1._DP
      TRIOT=1._DP
      TRICH=1._DP
 
C
C  INITIALIZE DATABASE NAMES
C
      DBHANDLE(1) = 'AMJUEL'
      DBHANDLE(2) = 'METHAN'
      DBHANDLE(3) = 'HYDHEL'
      DBHANDLE(4) = 'H2VIBR'
      DBHANDLE(5) = 'SPECTR'
      DBHANDLE(6) = 'PHOTON'
      DBHANDLE(7) = 'PHTNEW'
      DBHANDLE(8) = 'SPUTER'
      DBHANDLE(9) = 'TRIM  '
      DBHANDLE(10) = 'POLARI'
      DBHANDLE(11) = 'gr_ext'
      DBHANDLE(12) = 'mo_ext'
      DBHANDLE(13) = 'ADAS  '
      DBHANDLE(14) = 'HYDRTC'
      DBHANDLE(15) = 'HYDCRS'
      DBHANDLE(16) = 'HYDREA'
      DBHANDLE(17) = 'TAB2D '
 
      DBFNAME = REPEAT(' ',LEN(DBFNAME(1)))
      DBFNAME(1) = 'AMJUEL'
      DBFNAME(2) = 'METHANE'
      DBFNAME(3) = 'HYDHEL'
      DBFNAME(4) = 'H2VIBR'
      DBFNAME(5) = 'SPECTR'
      DBFNAME(6) = 'PHOTON'
      DBFNAME(7) = 'PHTNEW'
      DBFNAME(8) = 'SPUTER'
      DBFNAME(9)(1:5) = 'fort.'
      if (21+ifoff >= 100) then
        write (DBFNAME(9)(6:8),'(I3)') 21+ifoff
      else
        write (DBFNAME(9)(6:7),'(I2)') 21+ifoff
      endif
      DBFNAME(10) = 'POLARI'
      DBFNAME(11) = 'graphite_ext.dat'
      DBFNAME(12) = 'mo_ext.dat'
      DBFNAME(13) = ' '
      DBFNAME(14) = 'HYDRTC'
      DBFNAME(15) = 'HYDCRS'
      DBFNAME(16) = ' '
      DBFNAME(17) = ' '
 
      RETURN
      END SUBROUTINE EIRENE_INIT_CINIT
 
      END MODULE EIRMOD_CINIT
