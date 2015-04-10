C
C
      SUBROUTINE EIRENE_OUTUSR
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CPLOT
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_COUTAU
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CSPEZ
      IMPLICIT NONE
      real(dp) :: dummy(nrtal)
      integer :: iadv

      if (nadv <= nadvi+7) return

      ADDV(1,1:NSBOX_TAL) = ADDV(NADVI+1,1:NSBOX_TAL)
      ADDV(2,1:NSBOX_TAL) = ADDV(NADVI+2,1:NSBOX_TAL)
      ADDV(3,1:NSBOX_TAL) = ADDV(NADVI+3,1:NSBOX_TAL)
      ADDV(4,1:NSBOX_TAL) = ADDV(NADVI+4,1:NSBOX_TAL)
      ADDV(5,1:NSBOX_TAL) = ADDV(NADVI+5,1:NSBOX_TAL)
      ADDV(6,1:NSBOX_TAL) = ADDV(NADVI+6,1:NSBOX_TAL)
      ADDV(7,1:NSBOX_TAL) = ADDV(NADVI+7,1:NSBOX_TAL)

      DO 455 IADV=1,NADVI
        DUMMY(1:NSBOX_TAL) = ADDV(IADV,1:NSBOX_TAL)
        CALL EIRENE_INTTAL (DUMMY,VOLTAL,1,1,NSBOX_TAL,
     .                      ADDVI(IADV,0),
     .                      nr1tal, np2tal, nt3tal,nbmlt)
        ADDV(IADV,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)
455   end do

      iraps = 0
      nraps = 60

      call eirene_plteir(0)
      call eirene_rpsout

      RETURN
      END
