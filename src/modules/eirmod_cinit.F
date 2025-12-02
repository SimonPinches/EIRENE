      MODULE EIRMOD_CINIT

!  sep-05: specifications for databases added, ndbnames, dbhandle, dbfname
!  jul-06: database handle added for ADAS
!  20.06.07: deallocate tdmpar
cdr oct 18 : "density models" TDMPAR are now fully characterized
cdr          by reaction identifyer IRC and donor state
cdr          isp,itp,istr,
cdr          rather than by requiring a full reaction deck and executing SLREAC.f in here.
cdr dec 18 : remove obsolete photon-reflection data files graphite_ext and mo_ext
cdr aug 19 : add legal filename: AMMONX (ifile=12)
cdr nov.20 : remove CHELP, after cleanup of reaction data for density models.
cdr          Now only reaction number IRC read from block 4 is kept.
cdr nov 20 : remove obsolete POLARI file (was for high pressure mercury gas discharge lamps)
cdr dec 20 : remove unused PHTNEW file
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE
CDR FOR DIRECT TRANSFER OF INFO FROM FIND_PARAM.F INTO INPUT.F
      INTEGER, PUBLIC, SAVE :: INDPRO2_SAVE

      PUBLIC :: EIRENE_ALLOC_CINIT, EIRENE_DEALLOC_CINIT,
     .          EIRENE_INIT_CINIT, EIRENE_BROADCAST_CINIT,
     .          TDENMODEL, TDENMODAR

      INTEGER, PRIVATE, SAVE ::
     I         NCINIT, MCINIT, LCINIT

      REAL(DP),  PUBLIC, TARGET, ALLOCATABLE, SAVE :: RCINIT(:)
      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE :: ICINIT(:)
      LOGICAL, PUBLIC, TARGET, SAVE :: LCNIT(5)

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
     R DI0(:),DI1(:),DI2(:),DI3(:),DI4(:),DI5(:),
     R ZI0(:),ZI1(:),ZI2(:),ZI3(:),ZI4(:),ZI5(:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I INDPRO(:),INDGRD(:),INDSRC(:)

!pb   for NBMLT option, default .true.
      LOGICAL, PUBLIC, SAVE :: LCPYPRO(12)

      LOGICAL, PUBLIC, POINTER, SAVE ::
     L NLMACH, NLPITCH, NLKELV, NLMLTI, NLMLV

      CHARACTER(66), PUBLIC, SAVE :: CASENAME
      CHARACTER(10), PUBLIC, ALLOCATABLE, SAVE :: CDENMODEL(:)

      TYPE TDENMODEL
        REAL(DP) :: G_BOLTZ, DELTAE, A_CORONA, G_PLANCK,
     R              DVAL, TVAL,
     R              VXVAL, VYVAL, VZVAL,
     R              DFACTOR, TFACTOR, VFACTOR
cdr  NRE contributions to field particle, single (effective temperature)
cdr  currently NRE GT 1 available only for density model "COLRAD",
cdr  all others: NRE=1 (to be generalized still e.g. for density model "CORONA").
        INTEGER :: NRE
        INTEGER, POINTER :: ISP(:), ITP(:), ISTR(:), IRC(:)
      END TYPE

      TYPE TDENMODAR
        TYPE(TDENMODEL), POINTER :: TDM
      END TYPE

      TYPE(TDENMODAR), PUBLIC, ALLOCATABLE, SAVE :: TDMPAR(:)

      INTEGER, PUBLIC, PARAMETER :: NDBNAMES=17
      CHARACTER(400), PUBLIC, SAVE :: DBFNAME(NDBNAMES)
      CHARACTER(6), PUBLIC, SAVE :: DBHANDLE(NDBNAMES)
      LOGICAL, PUBLIC, SAVE :: LDBREAD(NDBNAMES)
cpg
      character(400), public, save :: MASTER_PATH
cpg

C> Default output file prefix
#ifndef HPUX
      CHARACTER(5), PUBLIC, PARAMETER :: FORT    = 'FORT.'
      CHARACTER(5), PUBLIC, PARAMETER :: FORT_LC = 'fort.'
#else
      CHARACTER(3), PUBLIC, PARAMETER :: FORT    = 'FTN'
      CHARACTER(3), PUBLIC, PARAMETER :: FORT_LC = 'ftn'
#endif

      CONTAINS

      SUBROUTINE EIRENE_ALLOC_CINIT
      INTEGER :: I

      IF (ALLOCATED(RCINIT)) RETURN

      NCINIT=42+36*NPLS
      MCINIT=12+3+NSTRA
      LCINIT=5    !dr: unused, so far.

      ALLOCATE (RCINIT(NCINIT))
      ALLOCATE (ICINIT(MCINIT))
      ALLOCATE (CDENMODEL(NPLS))
      ALLOCATE (TDMPAR(NPLS))
      DO I=1,NPLS
         NULLIFY(TDMPAR(I)%TDM)
      ENDDO
cpg
!pb   MASTER_PATH has already been set. Do not overwrite.
!pb   MASTER_PATH= ''
cpg

      WRITE (IUNMEM,'(A,T25,I15)')
     .      ' CINIT ',NCINIT*8 + MCINIT*4 +
     .                NPLS*LEN(CDENMODEL(1))
!dr  here should be added also LCINIT*.

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
      ZI0    => RCINIT(43+30*npls : 42+31*NPLS)
      ZI1    => RCINIT(43+31*npls : 42+32*NPLS)
      ZI2    => RCINIT(43+32*npls : 42+33*NPLS)
      ZI3    => RCINIT(43+33*npls : 42+34*NPLS)
      ZI4    => RCINIT(43+34*npls : 42+35*NPLS)
      ZI5    => RCINIT(43+35*npls : 42+36*NPLS)

      INDPRO => ICINIT( 1 : 12)
      INDGRD => ICINIT(13 : 15)
      INDSRC => ICINIT(16 : 15+NSTRA)

cdr units used for input of background profiles: V, B, Ti
      NLMACH => LCNIT(1)
      NLMLTI => LCNIT(2)
      NLMLV  => LCNIT(3)
      NLPITCH=> LCNIT(4)
      NLKELV => LCNIT(5)  ! not in use yet

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
          DEALLOCATE (TDMPAR(I)%TDM%IRC)
          DEALLOCATE (TDMPAR(I)%TDM)
        END IF
      END DO
      DEALLOCATE (TDMPAR)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CINIT


      SUBROUTINE EIRENE_INIT_CINIT
cdr called from ALLOC_CINIT,
cdr also called (EARLIER) from FIND_PARAM.
cdr In that case
cdr make sure then to not deal with unallocated variables
cdr i.e. not with RCINIT, ICINIT, CDENMODEL,
cdr but only with LCNIT, DBHANDLE, DBFNAME

      INTEGER :: I

      IF (ALLOCATED(RCINIT)) RCINIT = 0._DP
      IF (ALLOCATED(ICINIT)) ICINIT = 0
      LCNIT  = .FALSE.
      LCPYPRO = .TRUE.
      IF (ALLOCATED(CDENMODEL)) CDENMODEL = REPEAT(' ',LEN(CDENMODEL))

      CASENAME = REPEAT(' ',LEN(CASENAME))

      IF (ALLOCATED(TDMPAR)) THEN
        DO I = 1, NPLS
          NULLIFY(TDMPAR(I)%TDM)
        END DO
      END IF
C
C  INITIALIZE GEOMETRY DATA
C
      IF (ALLOCATED(RCINIT)) THEN
        ELLIN=1._DP
        ELLOT=1._DP
        ELLCH=1._DP
        TRIIN=1._DP
        TRIOT=1._DP
        TRICH=1._DP
      END IF

C
C  INITIALIZE DATABASE NAMES, KEY WORDS FOR PERMITTED NAMES
C
      DBHANDLE(1) = 'AMJUEL'
      DBHANDLE(2) = 'METHAN'
      DBHANDLE(3) = 'HYDHEL'
      DBHANDLE(4) = 'H2VIBR'
      DBHANDLE(5) = 'SPECTR'
      DBHANDLE(6) = 'PHOTON'
      DBHANDLE(7) = 'free'
      DBHANDLE(8) = 'SPUTER'
      DBHANDLE(9) = 'TRIM  '
      DBHANDLE(10) = 'free'
      DBHANDLE(11) = 'free'
      DBHANDLE(12) = 'AMMONX'
      DBHANDLE(13) = 'ADAS  '
      DBHANDLE(14) = 'free'
      DBHANDLE(15) = 'free'
      DBHANDLE(16) = 'free'
      DBHANDLE(17) = 'TAB2D '

C  FOR READING CFILE CARDS,
C  DEFAULT FOR DATABASE NAMES, CAN BE OVERWRITTEN BY CFILE CARDS
      DBFNAME = REPEAT(' ',LEN(DBFNAME(1)))
      DBFNAME(1) = 'AMJUEL'
      DBFNAME(2) = 'METHANE'
      DBFNAME(3) = 'HYDHEL'
      DBFNAME(4) = 'H2VIBR'
      DBFNAME(5) = 'SPECTR'  ! also out ?
      DBFNAME(6) = 'PHOTON'
      DBFNAME(7) = 'free'
      DBFNAME(8) = 'SPUTER'
#ifndef HPUX
      I = 5
#else
      I = 3
#endif
      DBFNAME(9)(1:I) = FORT_LC
CDR  TRIM USUALLY IS INCLUDED AS fort.21.
cdr  FOR SOME COMPILERS THIS IS ILLEGAL.
      if (21+ifoff >= 100) then
        write (DBFNAME(9)(I+1:I+3),'(I3)') 21+ifoff
      else
        write (DBFNAME(9)(I+1:I+2),'(I2)') 21+ifoff
      endif
      DBFNAME(10) = 'free'
      DBFNAME(11) = 'free'
      DBFNAME(12) = 'AMMONX'
      DBFNAME(13) = 'ADAS'
      DBFNAME(14) = 'free'
      DBFNAME(15) = 'free'
      DBFNAME(16) = 'free'
      DBFNAME(17) = 'TAB2D'

      LDBREAD = .FALSE.

      RETURN
      END SUBROUTINE EIRENE_INIT_CINIT


      SUBROUTINE EIRENE_BROADCAST_CINIT(ME)
      USE EIRMOD_MPI
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER, IPLS, NREF, IAN
      REAL(DP) :: RHELP(4)
      INTEGER, ALLOCATABLE :: IHELP(:)
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
      EXTERNAL :: MPI_BCAST
#endif

      IF (ME /= 0) CALL EIRENE_ALLOC_CINIT

      CALL MPI_BCAST (RCINIT,NCINIT,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICINIT,MCINIT,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LCNIT,LCINIT,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LCPYPRO,12,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

cdr  "density models", to set derived plasma background species IPLS
cdr  e.g. corona, colrad, const, saha, planck, .....
      CALL MPI_BCAST (CDENMODEL,10*NPLS,MPI_CHARACTER,0,
     .                MPI_COMM_WORLD,ier)

      DO IPLS=1,NPLS

        IF (VERIFY(CDENMODEL(IPLS),' ') == 0) CYCLE

        IF (ME == 0) NREF = TDMPAR(IPLS)%TDM%NRE
        CALL MPI_BCAST(NREF,1,MPI_INTEGER,0,
     .                 MPI_COMM_WORLD,ier)

cdr
        if(.not.allocated(ihelp)) then
          ALLOCATE (IHELP(4*NREF))
        endif

        IF (ME == 0) THEN
          RHELP(1) = TDMPAR(IPLS)%TDM%G_BOLTZ
          RHELP(2) = TDMPAR(IPLS)%TDM%DELTAE
          RHELP(3) = TDMPAR(IPLS)%TDM%A_CORONA
          RHELP(4) = TDMPAR(IPLS)%TDM%G_PLANCK

          IHELP(1:NREF) = TDMPAR(IPLS)%TDM%ISP(1:NREF)
          IHELP(1+NREF:2*NREF) = TDMPAR(IPLS)%TDM%ITP(1:NREF)
          IHELP(1+2*NREF:3*NREF) = TDMPAR(IPLS)%TDM%ISTR(1:NREF)
          IHELP(1+3*NREF:4*NREF) = TDMPAR(IPLS)%TDM%IRC(1:NREF)
          IAN = 0
        END IF

        CALL MPI_BCAST(RHELP,4,MPI_REAL8,0,
     .                MPI_COMM_WORLD,ier)
        CALL MPI_BCAST(IHELP,4*NREF,MPI_INTEGER,0,
     .                MPI_COMM_WORLD,ier)

        IF (ME > 0) THEN

          if(.not.associated(tdmpar(ipls)%tdm)) then

            ALLOCATE (TDMPAR(IPLS)%TDM)
            ALLOCATE (TDMPAR(IPLS)%TDM%ISP(NREF))
            ALLOCATE (TDMPAR(IPLS)%TDM%ITP(NREF))
            ALLOCATE (TDMPAR(IPLS)%TDM%ISTR(NREF))
            ALLOCATE (TDMPAR(IPLS)%TDM%IRC(NREF))

          endif

          TDMPAR(IPLS)%TDM%NRE = NREF
          TDMPAR(IPLS)%TDM%G_BOLTZ = RHELP(1)
          TDMPAR(IPLS)%TDM%DELTAE = RHELP(2)
          TDMPAR(IPLS)%TDM%A_CORONA = RHELP(3)
          TDMPAR(IPLS)%TDM%G_PLANCK = RHELP(4)

          TDMPAR(IPLS)%TDM%ISP(1:NREF)  = IHELP(1:NREF)
          TDMPAR(IPLS)%TDM%ITP(1:NREF)  = IHELP(1+NREF:2*NREF)
          TDMPAR(IPLS)%TDM%ISTR(1:NREF) = IHELP(1+2*NREF:3*NREF)
          TDMPAR(IPLS)%TDM%IRC(1:NREF)  = IHELP(1+3*NREF:4*NREF)
          IAN = 0
        END IF

        DEALLOCATE (IHELP)

      END DO

      CALL MPI_BCAST(DBFNAME,100*NDBNAMES,MPI_CHARACTER,0,
     .               MPI_COMM_WORLD,ier)

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END SUBROUTINE EIRENE_BROADCAST_CINIT

      END MODULE EIRMOD_CINIT
