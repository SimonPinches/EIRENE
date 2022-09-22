      subroutine eirene_defaults_usr

      use eirmod_parmmod 
cpg    
      use eirmod_cinit , only : MASTER_PATH
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
cpg
      MASTER_PATH = get_solpstop()
cpg
      return
      end subroutine eirene_defaults_usr
