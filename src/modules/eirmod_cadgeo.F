      MODULE EIRMOD_CADGEO

cdr Sept. 19
cdr Data for "additional surfaces": BREP formulation,
cdr      1st or 2nd order algebraic surfaces, and their boundaries
cdr Mostly read in input block 3B, (input.f), and some precomputed
cdr      coefficients set in TIMEA0 to speed up MC computation.

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_CADGEO, EIRENE_DEALLOC_CADGEO,
     P          EIRENE_INIT_CADGEO, EIRENE_BROADCAST_CADGEO

! PRIVATE VARIABLES USED FOR CALCULATION OF MEMORY USAGE
      INTEGER, PRIVATE, SAVE ::
     I NPLIM,  NADGEO, MADGEO

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        RADGEO(:,:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     I         IADGEO(:,:)

C NADGEO, REAL
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R A0LM(:),     A1LM(:),     A2LM(:),     A3LM(:),     A4LM(:),
     R A5LM(:),     A6LM(:),     A7LM(:),     A8LM(:),     A9LM(:),
     R ALIMS(:,:),  XLIMS(:,:),  YLIMS(:,:),  ZLIMS(:,:),
     R ALIMS0(:,:), XLIMS1(:,:), YLIMS1(:,:), ZLIMS1(:,:),
     R XLIMS2(:,:), YLIMS2(:,:), ZLIMS2(:,:),
     R XLIMS3(:,:), YLIMS3(:,:), ZLIMS3(:,:),
     R RLB(:),
     R P1(:,:),     P2(:,:),     P3(:,:),     P4(:,:),    P5(:,:),
     R P6(:,:),
     R ALM(:),      BLM(:),      CLM(:),
     R PS13(:,:),   PS23(:,:),
     R PS24(:,:),   PS34(:,:),   PS35(:,:),   PS45(:,:),
     R P1A(:),      P2A(:),      P1B(:),      P2B(:),
     R P1C(:),      P2C(:)

C MADGEO, INTEGER
      INTEGER, PUBLIC, POINTER, SAVE ::
     I ILIN(:),     ISCN(:)

      INTEGER, PUBLIC, SAVE :: NLIMI

      LOGICAL, PUBLIC, ALLOCATABLE, SAVE ::
     L RLBNOT(:)

! VARIABLES STORED FOR TRANSFER OF DATA TO JSON FILES
      INTEGER, PUBLIC, ALLOCATABLE, SAVE :: ITINI(:), ITEND(:)
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R XLCOR(:), YLCOR(:), ZLCOR(:),
     R XLREF(:), YLREF(:), ZLREF(:),
     R XLROT(:), YLROT(:), ZLROT(:), ALROT(:)

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_CADGEO

      IF (ALLOCATED(RADGEO)) RETURN

      NPLIM=10+14*9+10+12*3
      NADGEO=NLIM*NPLIM
      MADGEO=NLIM*2

      ALLOCATE (RADGEO(NPLIM,NLIM))
      ALLOCATE (IADGEO(2,NLIM))
      ALLOCATE (RLBNOT(NLIM))

      WRITE (IUNMEM,'(A,T25,I15)')
     .      ' CADGEO ',NLIM*(NPLIM+1)*8 + 2*NLIM*4

      A0LM => RADGEO(1,:)
      A1LM => RADGEO(2,:)
      A2LM => RADGEO(3,:)
      A3LM => RADGEO(4,:)
      A4LM => RADGEO(5,:)
      A5LM => RADGEO(6,:)
      A6LM => RADGEO(7,:)
      A7LM => RADGEO(8,:)
      A8LM => RADGEO(9,:)
      A9LM => RADGEO(10,:)
      ALIMS => RADGEO(11:19,:)
      XLIMS => RADGEO(20:28,:)
      YLIMS => RADGEO(29:37,:)
      ZLIMS => RADGEO(38:46,:)
      ALIMS0 => RADGEO(47:55,:)
      XLIMS1 => RADGEO(56:64,:)
      YLIMS1 => RADGEO(65:73,:)
      ZLIMS1 => RADGEO(74:82,:)
      XLIMS2 => RADGEO(83:91,:)
      YLIMS2 => RADGEO(92:100,:)
      ZLIMS2 => RADGEO(101:109,:)
      XLIMS3 => RADGEO(110:118,:)
      YLIMS3 => RADGEO(119:127,:)
      ZLIMS3 => RADGEO(128:136,:)
      RLB => RADGEO(137,:)
      P1 => RADGEO(138:140,:)
      P2 => RADGEO(141:143,:)
      P3 => RADGEO(144:146,:)
      P4 => RADGEO(147:149,:)
      P5 => RADGEO(150:152,:)
      P6 => RADGEO(153:155,:)
      ALM => RADGEO(156,:)
      BLM => RADGEO(157,:)
      CLM => RADGEO(158,:)
      PS13 => RADGEO(159:161,:)
      PS23 => RADGEO(162:164,:)
      PS24 => RADGEO(165:167,:)
      PS34 => RADGEO(168:170,:)
      PS35 => RADGEO(171:173,:)
      PS45 => RADGEO(174:176,:)
      P1A => RADGEO(177,:)
      P2A => RADGEO(178,:)
      P1B => RADGEO(179,:)
      P2B => RADGEO(180,:)
      P1C => RADGEO(181,:)
      P2C => RADGEO(182,:)

      ILIN => IADGEO(1,:)
      ISCN => IADGEO(2,:)

      ALLOCATE (ITINI(NLIM))
      ALLOCATE (ITEND(NLIM))
      ALLOCATE (XLCOR(NLIM))
      ALLOCATE (YLCOR(NLIM))
      ALLOCATE (ZLCOR(NLIM))
      ALLOCATE (XLREF(NLIM))
      ALLOCATE (YLREF(NLIM))
      ALLOCATE (ZLREF(NLIM))
      ALLOCATE (XLROT(NLIM))
      ALLOCATE (YLROT(NLIM))
      ALLOCATE (ZLROT(NLIM))
      ALLOCATE (ALROT(NLIM))

      CALL EIRENE_INIT_CADGEO

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CADGEO

      SUBROUTINE EIRENE_DEALLOC_CADGEO

      IF (.NOT.ALLOCATED(RADGEO)) RETURN

      DEALLOCATE (RADGEO)
      DEALLOCATE (IADGEO)
      DEALLOCATE (RLBNOT)

      DEALLOCATE (ITINI)
      DEALLOCATE (ITEND)
      DEALLOCATE (XLCOR)
      DEALLOCATE (YLCOR)
      DEALLOCATE (ZLCOR)
      DEALLOCATE (XLREF)
      DEALLOCATE (YLREF)
      DEALLOCATE (ZLREF)
      DEALLOCATE (XLROT)
      DEALLOCATE (YLROT)
      DEALLOCATE (ZLROT)
      DEALLOCATE (ALROT)
      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CADGEO


      SUBROUTINE EIRENE_INIT_CADGEO

      RADGEO = 0.D0
      IADGEO = 0
      RLBNOT = .FALSE.

      ITINI = 0
      ITEND = 0
      XLCOR = 0._DP
      YLCOR = 0._DP
      ZLCOR = 0._DP
      XLREF = 0._DP
      YLREF = 0._DP
      ZLREF = 0._DP
      XLROT = 0._DP
      YLROT = 0._DP
      ZLROT = 0._DP
      ALROT = 0._DP

      RETURN
      END SUBROUTINE EIRENE_INIT_CADGEO


      SUBROUTINE EIRENE_BROADCAST_CADGEO(ME)
      USE EIRMOD_MPI
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif

      IF (ME /= 0) CALL EIRENE_ALLOC_CADGEO

      CALL MPI_BCAST (RADGEO,NADGEO,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IADGEO,MADGEO,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NLIMI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END SUBROUTINE EIRENE_BROADCAST_CADGEO

      END MODULE EIRMOD_CADGEO
