      SUBROUTINE EIRENE_REINITIALIZATION_OF_EIRENE
 
      use EIRMOD_PRECISION
      use EIRMOD_PARMMOD
      use EIRMOD_COUTAU
      use EIRMOD_CRECH
 
      implicit none
 
      REAL(DP) :: dummy, ranf_eirene_reinit
C      REAL(DP) :: H1RN_REINIT
      integer :: idummy
      INTEGER, EXTERNAL :: ranset_eirene_reinit
 
C     reinitialization start
      call EIRENE_EIRENE_REINIT
      call EIRENE_SIGHA_REINIT
 
      dummy = ranf_eirene_reinit()
      idummy = ranset_eirene_reinit()
 
!pb   call INIT_COUTAU_REINIT
      call EIRENE_CRECH_REINIT
 
      call EIRENE_STCOOR_REINIT
      call EIRENE_SAMVOL_REINIT
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
      call EIRENE_BA_ALPHA_REINIT
      call EIRENE_EMIS_PROFILES_REINIT
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
 
