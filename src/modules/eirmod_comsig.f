cdr  aug. 17: added prspec, prargl;
cdr           Separate printout for energy (spectrally) resolved
cdr           from printout for spatially (along LOS) resolved data.
cdr           This was, so far, all mixed with TRCSIG (for debugging printout)
cdr  Jan. 2018  mod_addv added, as well as CNT data structure.
cdr  Oct 18   : rationalization of line emission specification,
cdr             now via reacdat(irc.....), block 4, exclusively.

      MODULE EIRMOD_COMSIG

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_COMSIG, EIRENE_DEALLOC_COMSIG,
     P          EIRENE_INIT_COMSIG, EIRENE_BROADCAST_COMSIG,
     p          TEMIS_MODEL, TCOMPO, TCONTRIB,
     P          ASSIGNMENT(=)


      INTEGER, PUBLIC, SAVE ::
     P         NCMSIG, MCMSIG

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE :: RCMSIG(:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R        FUFFER(:,:), ENERGY(:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R        XCHORD(:), YCHORD(:), ZCHORD(:),
     R        XPIVOT(:), YPIVOT(:), ZPIVOT(:),
     R        TILINE(:), TINP(:),   EMIN1(:),  EMAX1(:), ESHIFT(:)
      REAL(DP), PUBLIC, ALLOCATABLE :: RECADD(:,:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE :: ICMSIG(:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I         IPIVOT(:), ICHORD(:),
     I         NSPSTR(:), NSPSPZ(:),
     I         NSPBLC(:), NSPADD(:),
     I         NCHTAL(:), NSPINI(:), NSPEND(:),
     I         NSPSCL(:), NSPNEW(:)
      INTEGER, PUBLIC, ALLOCATABLE :: INTADD(:,:)
      INTEGER, PUBLIC :: NCTSIG

      INTEGER, PUBLIC, SAVE ::
     I         NCHORI,   NCHENI,   MOD_ADDV
      LOGICAL, PUBLIC, SAVE ::
     L         PRSPEC,   PRARGL

      LOGICAL, PUBLIC, ALLOCATABLE, SAVE :: NLSTCHR(:)

cdr data structures below are for line emission options,
cdr and generalize the older Ba_alpha, Ba_beta, ..., Ly_beta hard coded
cdr six hydrogenic line emission routines (6 lines), and, by coincidence, also
cdr six components per line: H, H+, H-, H2, H2+, H3+.

      CHARACTER(80), PUBLIC, ALLOCATABLE, SAVE :: CH_LINE_NAME(:)

      TYPE TCONTRIB
        INTEGER :: ISP, ITP, IRATIO, IRC_RAT(2), ISP_RAT(2), ITP_RAT(2) 
      END TYPE TCONTRIB

      TYPE TCOMPO
        CHARACTER(80) :: COMPO_NAME
        INTEGER :: NUM_CONTRIB, IADV, IRC
        TYPE(TCONTRIB), ALLOCATABLE :: CONTRIB(:)
      END TYPE TCOMPO

      TYPE TEMIS_MODEL
        CHARACTER(80) :: LINE_NAME
        INTEGER :: NUM_COMPO, IADV_TOTAL
        REAL(DP) :: EINSTEIN, TRANS_EN, ENERGY
        TYPE(TCOMPO), ALLOCATABLE :: COMPO(:)
      END TYPE TEMIS_MODEL

      TYPE(TEMIS_MODEL), PUBLIC, ALLOCATABLE, SAVE :: EMIS_LINES(:)

      INTERFACE ASSIGNMENT(=)  ! DEFINE ASSIGNMENT
        MODULE PROCEDURE EIRENE_CONTRIB_TO_CONTRIB
      END INTERFACE

      CONTAINS

      SUBROUTINE EIRENE_ALLOC_COMSIG

      IF (ALLOCATED(RCMSIG)) RETURN

      NCMSIG=11*NCHOR  ! 11 reals per line of sight, XCHORD,....
      MCMSIG=11*NCHOR  ! 11 integers per line of sight, IPIVOT,....

      ALLOCATE (RCMSIG(NCMSIG))
      ALLOCATE (FUFFER(NCHOR,NCHEN))
      ALLOCATE (ENERGY(NCHEN))
      ALLOCATE (ICMSIG(MCMSIG))
      ALLOCATE (NLSTCHR(NCHOR))

      ALLOCATE (CH_LINE_NAME(NCHOR))

      WRITE (IUNMEM,'(A,T25,I15)')
     .      ' COMSIG ',(NCMSIG+(NCHOR+1)*NCHEN)+8 + (MCMSIG+NCHOR)*4 +
     .                 NCHOR*80

      XCHORD => RCMSIG( 0*NCHOR+1 :  1*NCHOR)
      YCHORD => RCMSIG( 1*NCHOR+1 :  2*NCHOR)
      ZCHORD => RCMSIG( 2*NCHOR+1 :  3*NCHOR)
      XPIVOT => RCMSIG( 3*NCHOR+1 :  4*NCHOR)
      YPIVOT => RCMSIG( 4*NCHOR+1 :  5*NCHOR)
      ZPIVOT => RCMSIG( 5*NCHOR+1 :  6*NCHOR)
      TILINE => RCMSIG( 6*NCHOR+1 :  7*NCHOR)
      TINP   => RCMSIG( 7*NCHOR+1 :  8*NCHOR)
      EMIN1  => RCMSIG( 8*NCHOR+1 :  9*NCHOR)
      EMAX1  => RCMSIG( 9*NCHOR+1 : 10*NCHOR)
      ESHIFT => RCMSIG(10*NCHOR+1 : 11*NCHOR)

      IPIVOT => ICMSIG( 0*NCHOR+1 :  1*NCHOR)
      ICHORD => ICMSIG( 1*NCHOR+1 :  2*NCHOR)
      NSPSTR => ICMSIG( 2*NCHOR+1 :  3*NCHOR)
      NSPSPZ => ICMSIG( 3*NCHOR+1 :  4*NCHOR)
      NSPBLC => ICMSIG( 4*NCHOR+1 :  5*NCHOR)
      NSPADD => ICMSIG( 5*NCHOR+1 :  6*NCHOR)
      NCHTAL => ICMSIG( 6*NCHOR+1 :  7*NCHOR)
      NSPINI => ICMSIG( 7*NCHOR+1 :  8*NCHOR)
      NSPEND => ICMSIG( 8*NCHOR+1 :  9*NCHOR)
      NSPSCL => ICMSIG( 9*NCHOR+1 : 10*NCHOR)
      NSPNEW => ICMSIG(10*NCHOR+1 : 11*NCHOR)

      CH_LINE_NAME = REPEAT(' ',80)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_COMSIG


      SUBROUTINE EIRENE_DEALLOC_COMSIG

      INTEGER :: I, J

      IF (.NOT.ALLOCATED(RCMSIG)) RETURN

      DEALLOCATE (RCMSIG)
      DEALLOCATE (FUFFER)
      DEALLOCATE (ENERGY)
      DEALLOCATE (ICMSIG)
      DEALLOCATE (NLSTCHR)

      DEALLOCATE (CH_LINE_NAME)

      IF (ALLOCATED(EMIS_LINES) .AND. (NUM_LINES > 0)) THEN

         DO I = 1, NUM_LINES

           IF (EMIS_LINES(I)%NUM_COMPO > 0) THEN

             DO J=1, EMIS_LINES(I)%NUM_COMPO
               DEALLOCATE (EMIS_LINES(I)%COMPO(J)%CONTRIB)
             END DO

             DEALLOCATE (EMIS_LINES(I)%COMPO)

           END IF

        ENDDO

        DEALLOCATE (EMIS_LINES)
      END IF

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_COMSIG


      SUBROUTINE EIRENE_INIT_COMSIG

      RCMSIG = 0._DP
      FUFFER = 0._DP
      ENERGY = 0._DP
      ICMSIG = 0
C
C  SET DEFAULTS FOR LINE OF SIGHT INTEGRATION (BLOCK 12)
C
      NSPBLC = 1
      NSPADD = 0
      NSPNEW = 0

      NLSTCHR = .FALSE.

      RETURN
      END SUBROUTINE EIRENE_INIT_COMSIG


      SUBROUTINE EIRENE_CONTRIB_TO_CONTRIB (CONA, CONB)

      TYPE(TCONTRIB), INTENT(OUT) :: CONA
      TYPE(TCONTRIB), INTENT(IN) :: CONB

      CONA%ISP          = CONB%ISP
      CONA%ITP          = CONB%ITP
      CONA%IRATIO       = CONB%IRATIO
      CONA%IRC_RAT      = CONB%IRC_RAT
      CONA%ISP_RAT = CONB%ISP_RAT
      CONA%ITP_RAT = CONB%ITP_RAT

      RETURN
      END SUBROUTINE EIRENE_CONTRIB_TO_CONTRIB


      SUBROUTINE EIRENE_BROADCAST_COMSIG(ME)
      USE EIRMOD_MPI
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER

      IF ((ME /= 0) .AND. (NCHOR > 0)) CALL EIRENE_ALLOC_COMSIG

      IF (ALLOCATED(RCMSIG)) THEN
        CALL MPI_BCAST (RCMSIG,NCMSIG,MPI_REAL8,0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (FUFFER,NCHOR*NCHEN,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (ENERGY,NCHEN,MPI_REAL8,0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (ICMSIG,MCMSIG,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (NLSTCHR,NCHOR,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (CH_LINE_NAME,80*NCHOR,MPI_CHARACTER,
     .                   0,MPI_COMM_WORLD,ier)
      END IF
      CALL MPI_BCAST (NCHORI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NCHENI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
cdr  additional output tallies added by code itself (rather than via input block 14).
      CALL MPI_BCAST (MOD_ADDV,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      IF (NUM_LINES > 0) THEN
        CALL EIRENE_BROAD_EMIS_LINES(ME)
      END IF
      
      END SUBROUTINE EIRENE_BROADCAST_COMSIG


      SUBROUTINE EIRENE_BROAD_EMIS_LINES(ME)
      USE EIRMOD_MPI
      INTEGER, INTENT(IN) :: ME
      INTEGER :: I, J, K, NUM_COMPO, NUM_CONTRIB, IER
      TYPE(TCONTRIB) :: CNT

      IF (.NOT.ALLOCATED(EMIS_LINES)) THEN
        ALLOCATE (EMIS_LINES(NUM_LINES))
        EMIS_LINES%LINE_NAME = REPEAT(' ',80)
        EMIS_LINES%NUM_COMPO = 0
      END IF

      DO I = 1, NUM_LINES

        CALL MPI_BCAST (EMIS_LINES(I)%LINE_NAME,80,MPI_CHARACTER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (EMIS_LINES(I)%NUM_COMPO,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (EMIS_LINES(I)%IADV_TOTAL,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (EMIS_LINES(I)%EINSTEIN,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (EMIS_LINES(I)%TRANS_EN,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)

        NUM_COMPO = EMIS_LINES(I)%NUM_COMPO

        IF (ME /= 0) THEN
          ALLOCATE (EMIS_LINES(I)%COMPO(NUM_COMPO))
        END IF

        DO J = 1, NUM_COMPO
          CALL MPI_BCAST (EMIS_LINES(I)%COMPO(J)%COMPO_NAME,80,
     .                    MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
          CALL MPI_BCAST (EMIS_LINES(I)%COMPO(J)%NUM_CONTRIB,1,
     .                    MPI_INTEGER,0,MPI_COMM_WORLD,ier)
          CALL MPI_BCAST (EMIS_LINES(I)%COMPO(J)%IADV,1,
     .                    MPI_INTEGER,0,MPI_COMM_WORLD,ier)
          CALL MPI_BCAST (EMIS_LINES(I)%COMPO(J)%IRC,1,
     .                    MPI_INTEGER,0,MPI_COMM_WORLD,ier)

          NUM_CONTRIB = EMIS_LINES(I)%COMPO(J)%NUM_CONTRIB

          IF (ME /= 0) THEN
            ALLOCATE (EMIS_LINES(I)%COMPO(J)%CONTRIB(NUM_CONTRIB))
          END IF

          DO K = 1, NUM_CONTRIB

            IF (ME == 0) CNT = EMIS_LINES(I)%COMPO(J)%CONTRIB(K)

            CALL MPI_BCAST (CNT%ISP,3,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
            CALL MPI_BCAST (CNT%ITP,3,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
            CALL MPI_BCAST (CNT%IRATIO,1,MPI_INTEGER,
     .                      0,MPI_COMM_WORLD,ier)
cdr donor state densities may be converted to other species/isotopes densities,
cdr by multiplication with a density "ratio": e.g. nH2+ = nH2 * ratio, 
cdr with ratio = CR equilibrium ratio:(nH2+/nH2) 

c  density ratio factor, for this donor state 
c  plus,(for ratio2) two further density factors
            CALL MPI_BCAST (CNT%ISP_RAT,2,MPI_INTEGER,
     .                      0,MPI_COMM_WORLD,ier)
            CALL MPI_BCAST (CNT%ITP_RAT,2,MPI_INTEGER,
     .                      0,MPI_COMM_WORLD,ier)
            CALL MPI_BCAST (CNT%IRC_RAT,2,MPI_INTEGER,
     .                      0,MPI_COMM_WORLD,ier)

            IF (ME /= 0) EMIS_LINES(I)%COMPO(J)%CONTRIB(K) = CNT

          END DO
        END DO


      END DO

      RETURN
      END SUBROUTINE EIRENE_BROAD_EMIS_LINES

      END MODULE EIRMOD_COMSIG
