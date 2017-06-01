      SUBROUTINE EIRENE_REINITIALIZATION_OF_EIRENE
 
      use EIRMOD_PRECISION
      use EIRMOD_PARMMOD
      use EIRMOD_COUTAU
      use EIRMOD_CRECH
 
      implicit none
 
      REAL(DP) :: dummy, H1RN_REINIT,  ranf_eirene_reinit,
     &     sheath_reinit
      INTEGER, EXTERNAL :: RANSET_EIRENE_REINIT
 
C     reinitialization start
      call EIRENE_EIRENE_REINIT
      call EIRENE_SIGHA_REINIT
 
!  OUT Of USE
!      dummy = H1RN_REINIT
!      call H1RNV_REINIT
 
      dummy = ranf_eirene_reinit()
      dummy = ranset_eirene_reinit()
 
!pb      call INIT_COUTAU_REINIT
      call EIRENE_CRECH_REINIT
 
      call EIRENE_STCOOR_REINIT
      call EIRENE_SAMVOL_REINIT
      call EIRENE_VELOCX_REINIT
      call EIRENE_VELOEL_REINIT
 
      call EIRENE_PL3D_REINIT
      call EIRENE_PLT2D_REINIT
      call EIRENE_PLTEIR_REINIT
 
      call EIRENE_REFLEC_REINIT
 
!pb      call UPTCOP_REINIT
 
      call EIRENE_STATIS_BGK_REINIT

csw 18apr07
      call EIRENE_SPUTER_REINIT
      call EIRENE_BA_ALPHA_REINIT
!pb      call EIRENE_BA_GAMMA_REINIT
!pb      call EIRENE_LY_BETA_REINIT
      call EIRENE_UPTBGK_REINIT
      call EIRENE_STORE_REINIT
!out      call EIRENE_MKCENS_REINIT
      call EIRENE_update_reinit
      call EIRENE_update_spectrum_reinit
csw
 
C     reinitialization end
 
 
      end
 
