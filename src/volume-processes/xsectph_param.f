C  28.6.05
C nnrot --> nrPH, for subr. find_param, setamd, after removing phv_nrota,
C                 and phv_nrotph
C
      SUBROUTINE EIRENE_XSECTPH_PARAM

      USE EIRMOD_PARMMOD
      USE EIRMOD_COMXS
      USE EIRMOD_COMUSR
      USE EIRMOD_PHOTON
      IMPLICIT NONE
csw
csw  PHOTON COLLISION, PH - type
csw
      integer :: iphot,nrc,kk

      do iphot=1,nphoti
         if(nrcph(iphot) > 0) then
            do nrc=1,nrcph(iphot)
               kk=ireacph(iphot,nrc)
               if(iswr(kk) == 7) then
                 NRPH=NRPH+1
               endif
            enddo
         endif
      enddo
cdr
c  this call is still necessary because not all PH-process data
c  have already been moved to module COMXS.
C  Still some clean-up work to be done
cdr
      call EIRENE_PH_ALLOC_XSECTPH(nrph)

      RETURN

      END SUBROUTINE EIRENE_XSECTPH_PARAM
