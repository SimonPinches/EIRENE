      SUBROUTINE EIRENE_REINITIALIZATION_OF_EIRENE

      use EIRMOD_PRECISION
      use EIRMOD_PARMMOD
      use EIRMOD_COMUSR
      use EIRMOD_COUTAU
      use EIRMOD_CRECH
      use EIRMOD_LININT, ONLY: EIRENE_LININT_REINIT
      use EIRMOD_SIGLINE, ONLY: EIRENE_SIGLINE_REINIT
      use EIRMOD_COLRAD, ONLY: EIRENE_COLRAD_REINIT
      use EIRMOD_SAMVOL, ONLY: EIRENE_SAMVOL_REINIT
      use EIRMOD_UPTBGK, ONLY: EIRENE_UPTBGK_REINIT
      use EIRMOD_STATIS_BGK, ONLY: EIRENE_STATIS_BGK_REINIT
      use EIRMOD_RANF, ONLY: RANF_EIRENE_REINIT
      use EIRMOD_RANSET, ONLY: RANSET_EIRENE_REINIT
      use EIRMOD_VELOPI, ONLY: EIRENE_VELOPI_REINIT
      use EIRMOD_VELOEL, ONLY: EIRENE_VELOEL_REINIT
      use EIRMOD_VELOCX, ONLY: EIRENE_VELOCX_REINIT
      use EIRMOD_SWITCH_PARTINFO, ONLY: EIRENE_REINIT_PARTINFO

      implicit none

      REAL(DP) :: dummy
ctk     , ranf_eirene_reinit
C      REAL(DP) :: H1RN_REINIT
      integer :: idummy
ctk      INTEGER, EXTERNAL :: ranset_eirene_reinit

C     reinitialization start
      call EIRENE_EIRENE_REINIT
      call EIRENE_SIGLINE_REINIT

      dummy = ranf_eirene_reinit()
      idummy = ranset_eirene_reinit()

!pb   call INIT_COUTAU_REINIT
      call EIRENE_CRECH_REINIT
      call EIRENE_COMUSR_REINIT

      call EIRENE_STCOOR_REINIT
      call EIRENE_SAMVOL_REINIT
      call EIRENE_VELOPI_REINIT
      call EIRENE_VELOCX_REINIT
      call EIRENE_VELOEL_REINIT

      call EIRENE_PL3D_REINIT
      call EIRENE_PLT2D_REINIT
      call EIRENE_PLTEIR_REINIT

      call EIRENE_REFLEC_REINIT

!pb   call UPTCOP_REINIT

      call EIRENE_STATIS_BGK_REINIT

csw 18apr07
      call EIRENE_SPUTER_REINIT

cdr   call EIRENE_BA_ALPHA_REINIT
!pb   call EIRENE_BA_GAMMA_REINIT
!pb   call EIRENE_LY_BETA_REINIT
      call EIRENE_LININT_REINIT  !cdr, july 17, added

      call EIRENE_UPTBGK_REINIT
!out  call EIRENE_MKCENS_REINIT
      call EIRENE_update_reinit
      call EIRENE_update_spectrum_reinit
!pb 03aug17
      call eirene_colrad_reinit

      call eirene_reinit_partinfo

C     reinitialization end


      end
