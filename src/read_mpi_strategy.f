      subroutine eirene_read_mpi_strategy_fixed (line)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT
      USE EIRMOD_COMUSR, ONLY: NPRLL
      USE EIRMOD_CPES, ONLY: NPRS, NLIDENT,
     >    STRATEGY_UNDEFINED, STRATEGY_EMBARRASS,
     >    STRATEGY_ORIGINAL, STRATEGY_APCAS, STRATEGY_BALANCED,
     >    INPUT_DISTRIBUTION_STRATEGY

      implicit none
      character(80), intent(in) :: LINE
      character(80) :: zeile
C
C  READ MPI STRATEGY IF PRESENT
C
      IF (INDEX(LINE,'INFORMATION_FOR_MPI') /= 0) THEN
        CALL EIRENE_MASAGE('*** INFORMATION_FOR_MPI')
        DO
          READ (IUNIN,'(A80)',END=1596) ZEILE
          IF (INDEX(ZEILE,'ORIGINAL') /= 0) THEN
            NPRLL = 1
            input_distribution_strategy = STRATEGY_ORIGINAL
            CALL EIRENE_MASAGE
     .           ('APPLYING "ORIGINAL" PARALLELIZATION STRATEGY')
            goto 1597
          ELSE IF (INDEX(ZEILE,'APCAS') /= 0) THEN
            NPRLL = 2
            input_distribution_strategy = STRATEGY_APCAS
            CALL EIRENE_MASAGE
     .           ('APPLYING "APCAS" PARALLELIZATION STRATEGY')
            goto 1597
          ELSE IF (INDEX(ZEILE,'BALANCED') /= 0) THEN
            NPRLL = 3
            input_distribution_strategy = STRATEGY_BALANCED
            CALL EIRENE_MASAGE
     .           ('APPLYING "BALANCED" PARALLELIZATION STRATEGY')
            goto 1597
          ELSE IF (INDEX(ZEILE,'EMBARRASS') /= 0) THEN
            NPRLL = 0
            input_distribution_strategy = STRATEGY_EMBARRASS
            CALL EIRENE_MASAGE
     .           ('APPLYING "EMBARRASSINGLY PARALLEL" STRATEGY')
            goto 1597
          ELSE IF (INDEX(ZEILE,'AUTOMATIC') /=0) THEN
            CALL EIRENE_MASAGE
     .           ('APPLYING AUTOMATIC PARALLELIZATION STRATEGY')
#if MPI_VERSION < 3
            NPRLL = 1
            input_distribution_strategy = STRATEGY_ORIGINAL
            CALL EIRENE_MASAGE
     .           ('CURRENT AUTOMATIC ASSIGNMENT IS "ORIGINAL"')
#else
            NPRLL = 3
            input_distribution_strategy = STRATEGY_BALANCED
            CALL EIRENE_MASAGE
     .           ('CURRENT AUTOMATIC ASSIGNMENT IS "BALANCED"')
#endif
            goto 1597
          END IF
        END DO
 1596   CONTINUE
        if (input_distribution_strategy.eq.STRATEGY_UNDEFINED) THEN
          CALL EIRENE_MASAGE
     .         ('NO MPI PARALLELIZATION STRATEGY NAME RECOGNIZED')
#if MPI_VERSION < 3
          NPRLL = 1
          input_distribution_strategy = STRATEGY_ORIGINAL
          CALL EIRENE_MASAGE
     .         ('APPLYING DEFAULT "ORIGINAL" STRATEGY')
#else
          NPRLL = 3
          input_distribution_strategy = STRATEGY_BALANCED
          CALL EIRENE_MASAGE
     .         ('APPLYING DEFAULT "BALANCED" STRATEGY')
#endif
        end if
 1597   CONTINUE
        CALL EIRENE_LEER(1)
        GOTO 1599

      ELSE ! string 'INFORMATION_FOR_MPI' not found
        IF (NPRS > 1 .and.
     &     input_distribution_strategy.eq.STRATEGY_UNDEFINED) THEN
          CALL EIRENE_MASAGE
     .       ('NO MPI PARALLELIZATION STRATEGY PROVIDED')
#if MPI_VERSION < 3
          NPRLL = 1
          input_distribution_strategy = STRATEGY_ORIGINAL
          CALL EIRENE_MASAGE
     .       ('APPLYING DEFAULT "ORIGINAL" STRATEGY')
#else
          NPRLL = 3
          input_distribution_strategy = STRATEGY_BALANCED
          CALL EIRENE_MASAGE
     .       ('APPLYING DEFAULT "BALANCED" STRATEGY')
#endif
          CALL EIRENE_LEER(1)
        END IF
      END IF
 1599 CONTINUE
      
      return
      end subroutine eirene_read_mpi_strategy_fixed
