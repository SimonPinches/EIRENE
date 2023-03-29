cpg called from find_param.f 
 
      subroutine eirene_couple_param_consistency(nlimi,nstsi,textal,ntx)
       use eirmod_parmmod
       use eirmod_comusr, only : natmi,nmoli,nioni
       
       IMPLICIT NONE
       integer, intent(in) :: nlimi,nstsi,ntx
       character(8), intent(in) :: textal(ntx)
     
      return
      end subroutine eirene_couple_param_consistency
