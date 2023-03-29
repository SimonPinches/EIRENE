      subroutine eirene_defaults_usr

      use eirmod_parmmod 
cpg    
      use eirmod_cinit , only : MASTER_PATH
      use eirmod_cpes , only : STRATEGY_DEFAULT, NPRLL_DEFAULT,
     .                         STRATEGY_BALANCED
cpg
      implicit none 
cpg
      character*256 :: get_solpstop       
      external get_solpstop 
cpg
c     set file units that allow direct control
      IUNMEM = 85
      IUNRAPS = 80
      IUNRAPSVEC = 70
      LOUTAPP = .FALSE. 
      LIF3COP_FROM_LOOP = .TRUE.
      LPE0_TO_STDOUT = .TRUE.
cpg
      MASTER_PATH = get_solpstop()
cpg
      STRATEGY_DEFAULT = STRATEGY_BALANCED
      NPRLL_DEFAULT = 3
      
      return
      end subroutine eirene_defaults_usr
