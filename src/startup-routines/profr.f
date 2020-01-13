C
C
      SUBROUTINE EIRENE_PROFR (PRO,IINDEX,NSPZI,NSPZ1,NDAT)
C
C  READ ENTIRE PROFILE FROM WORK ARRAY PLASMA_BCKGROUND   (EIRMOD_cspei)
C  NSPZ1:  first dimension of PRO array as in calling program
C  NSPZI:  fill the first NSPZI fields. NZPZI LE NSPZ1 necessarily

C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CSPEI
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IINDEX, NSPZI, NSPZ1, NDAT
      REAL(DP), INTENT(OUT) :: PRO(NSPZ1,*)
      if (nspzi.gt.nspz1 .or. nspz1.le.0) then
        write (iunout,*) 'error in PROFR' 
        write (iunout,*) 'PRO: incorrect dimension in calling program'
        write (iunout,*) 'nspz1, nspzi= ',NSPZ1, NSPZI
      endif

      PRO(1:NSPZI,1:NDAT) = PLASMA_BCKGRND(IINDEX+1:IINDEX+NSPZI,1:NDAT)
      RETURN
      END
