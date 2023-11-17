cpg called from find_param.f

      subroutine eirene_couple_param_consistency(nlimi,nstsi,textal,ntx)
      use eirmod_parmmod
      use eirmod_extrab25
      use eirmod_comusr, only : natmi,nmoli,nioni

      IMPLICIT NONE
      integer, intent(in) :: nlimi,nstsi,ntx
      character(8), intent(in) :: textal(ntx)


#ifndef DEF_ISOEXTRA
#define DEF_ISOEXTRA 0
#endif

#ifndef DEF_NPHID
#define DEF_NPHID 1
      N3RD=DEF_NPHID
#endif

#ifndef DEF_NGSTAL
#define DEF_NGSTAL 0
      NGSTAL=DEF_NGSTAL
#endif

#ifdef B25_EIRENE
#include <DIMENSIONS.F>

      IF (NGSTAL.GT.DEF_NGSTAL) THEN
        WRITE(iunout,*)
     .   'NGSTAL from KOPPLDIM.F (',DEF_NGSTAL,') is too small'
        WRITE(iunout,*)
     .   'compared to value from input file (',NGSTAL,')'
        WRITE(iunout,*) 'Expect trouble at coupling time !'
        WRITE(iunout,*)
     .   'Increase DEF_NGSTAL in DIMENSIONS.F and recompile.'
        WRITE(0,*)
     .   'NGSTAL from KOPPLDIM.F (',DEF_NGSTAL,') is too small'
        WRITE(0,*)
     .   'compared to value from input file (',NGSTAL,')'
        WRITE(0,*) 'Expect trouble at coupling time !'
        WRITE(0,*)
     .   'Increase DEF_NGSTAL in DIMENSIONS.F and recompile.'
      ENDIF

#ifndef ALLOCATE_AND_NAMELIST
      IF (NR1ST.GT.DEF_NYD+1) THEN
        WRITE(iunout,*)
     .   'N1ST from KOPPLDIM.F is too small (',DEF_NYD+1,')'
        WRITE(iunout,*)
     .   'compared to NR1ST from input file (',NR1ST,')'
        WRITE(iunout,*) 'Expect trouble at coupling time !'
        WRITE(iunout,*)
     .   'Increase DEF_NYD in DIMENSIONS.F and recompile.'
        WRITE(0,*)
     .   'N1ST from KOPPLDIM.F is too small (',DEF_NYD+1,')'
        WRITE(0,*)
     .   'compared to NR1ST from input file (',NR1ST,')'
        WRITE(0,*) 'Expect trouble at coupling time !'
        WRITE(0,*)
     .   'Increase DEF_NYD in DIMENSIONS.F and recompile.'
      ENDIF
#endif

#ifndef ALLOCATE_AND_NAMELIST
      IF (NLPOL.AND.NP2ND.GT.
     . DEF_NXD+1+DEF_NCUT+(DEF_NCUT/2-1)*(1+DEF_ISOEXTRA)) THEN
        WRITE(iunout,*)
     .   'N2ND from KOPPLDIM.F is too small (',
     .    DEF_NXD+1+DEF_NCUT+(DEF_NCUT/2-1)*(1+DEF_ISOEXTRA),')'
        WRITE(iunout,*) 'compared to NP2ND from input file (',NP2ND,')'
        WRITE(iunout,*) 'Expect trouble at coupling time !'
        WRITE(iunout,*)
     .   'Increase DEF_NXD in DIMENSIONS.F and recompile.'
        WRITE(0,*)
     .   'N2ND from KOPPLDIM.F is too small (',
     .    DEF_NXD+1+DEF_NCUT+(DEF_NCUT/2-1)*(1+DEF_ISOEXTRA),')'
        WRITE(0,*) 'compared to NP2ND from input file (',NP2ND,')'
        WRITE(0,*) 'Expect trouble at coupling time !'
        WRITE(0,*)
     .   'Increase DEF_NXD in DIMENSIONS.F and recompile.'
      ENDIF
#endif

#ifndef ALLOCATE_AND_NAMELIST
      IF (NSTSI.GT.DEF_NSTS) THEN
        WRITE(iunout,*)
     .   'NSTS from KOPPLDIM.F is too small (',DEF_NSTS,')'
        WRITE(iunout,*)
     .   'compared to NSTSI from input file (',NSTSI,')'
        WRITE(iunout,*) 'Expect trouble at coupling time !'
        WRITE(iunout,*)
     .   'Increase DEF_NSTS in DIMENSIONS.F and recompile.'
        WRITE(0,*)
     .   'NSTS from KOPPLDIM.F is too small (',DEF_NSTS,')'
        WRITE(0,*) 'compared to NSTSI from input file (',NSTSI,')'
        WRITE(0,*) 'Expect trouble at coupling time !'
        WRITE(0,*)
     .   'Increase DEF_NSTS in DIMENSIONS.F and recompile.'
      ENDIF

#endif 

#ifndef ALLOCATE_AND_NAMELIST
      IF (NLIMI.GT.DEF_NLIM) THEN
        WRITE(iunout,*)
     .   'NLIM from KOPPLDIM.F is too small (',DEF_NLIM,')'
        WRITE(iunout,*) 'compared to NLIMI from input file (',NLIMI,')'
        WRITE(iunout,*) 'Expect trouble at coupling time !'
        WRITE(iunout,*)
     .   'Increase DEF_NLIM in DIMENSIONS.F and recompile.'
        WRITE(0,*)
     .   'NLIM from KOPPLDIM.F is too small (',DEF_NLIM,')'
        WRITE(0,*) 'compared to NLIMI from input file (',NLIMI,')'
        WRITE(0,*) 'Expect trouble at coupling time !'
        WRITE(0,*)
     .   'Increase DEF_NLIM in DIMENSIONS.F and recompile.'
      ENDIF
#endif

#ifndef ALLOCATE_AND_NAMELIST
      IF (NATMI.GT.DEF_NATM) THEN
        WRITE(iunout,*)
     .   'NATM from KOPPLDIM.F is too small (',DEF_NATM,')'
        WRITE(iunout,*) 'compared to NATMI from input file (',NATMI,')'
        WRITE(iunout,*) 'Expect trouble at coupling time !'
        WRITE(iunout,*)
     .   'Increase DEF_NATM in DIMENSIONS.F and recompile.'
        WRITE(0,*)
     .   'NATM from KOPPLDIM.F is too small (',DEF_NATM,')'
        WRITE(0,*) 'compared to NATMI from input file (',NATMI,')'
        WRITE(0,*) 'Expect trouble at coupling time !'
        WRITE(0,*)
     .   'Increase DEF_NATM in DIMENSIONS.F and recompile.'
      ENDIF
#endif
      if(.not.allocated(TEXTA)) allocate(TEXTA(NATM))

#ifndef ALLOCATE_AND_NAMELIST
      IF (NMOLI.GT.DEF_NMOL) THEN
        WRITE(iunout,*)
     .   'NMOL from KOPPLDIM.F is too small (',DEF_NMOL,')'
        WRITE(iunout,*) 'compared to NMOLI from input file (',NMOLI,')'
        WRITE(iunout,*) 'Expect trouble at coupling time !'
        WRITE(iunout,*)
     .   'Increase DEF_NMOL in DIMENSIONS.F and recompile.'
        WRITE(0,*)
     .   'NMOL from KOPPLDIM.F is too small (',DEF_NMOL,')'
        WRITE(0,*) 'compared to NMOLI from input file (',NMOLI,')'
        WRITE(0,*) 'Expect trouble at coupling time !'
        WRITE(0,*)
     .   'Increase DEF_NMOL in DIMENSIONS.F and recompile.'
      ENDIF

#endif

#ifndef ALLOCATE_AND_NAMELIST
      IF (NIONI.GT.DEF_NION) THEN
        WRITE(iunout,*)
     .   'NION from KOPPLDIM.F is too small (',DEF_NION,')'
        WRITE(iunout,*) 'compared to NIONI from input file (',NIONI,')'
        WRITE(iunout,*) 'Expect trouble at coupling time !'
        WRITE(iunout,*)
     .   'Increase DEF_NION in DIMENSIONS.F and recompile.'
        WRITE(0,*)
     .   'NION from KOPPLDIM.F is too small (',DEF_NION,')'
        WRITE(0,*) 'compared to NIONI from input file (',NIONI,')'
        WRITE(0,*) 'Expect trouble at coupling time !'
        WRITE(0,*)
     .   'Increase DEF_NION in DIMENSIONS.F and recompile.'
      ENDIF
#endif

#ifndef ALLOCATE_AND_NAMELIST
      IF (NSTRAI.GT.DEF_NSTRA) THEN
        WRITE(iunout,*)
     .   'NSTRA from KOPPLDIM.F is too small (',DEF_NSTRA,')'
        WRITE(iunout,*)
     .   'compared to NSTRAI from input file (',NSTRAI,')'
        WRITE(iunout,*) 'Expect trouble at coupling time !'
        WRITE(iunout,*)
     .   'Increase DEF_NSTRA in DIMENSIONS.F and recompile.'
        WRITE(0,*)
     .   'NSTRA from KOPPLDIM.F is too small (',DEF_NSTRA,')'
        WRITE(0,*)
     .   'compared to NSTRAI from input file (',NSTRAI,')'
        WRITE(0,*) 'Expect trouble at coupling time !'
        WRITE(0,*)
     .   'Increase DEF_NSTRA in DIMENSIONS.F and recompile.'
      ENDIF
#endif

#endif

cxpb Define some numbers needed by eirmod_extrab25
!cank 960623
      nnatmi=natmi
      nnmoli=nmoli
      nnioni=nioni
!cank 960513
      nnlimi=nlimi
      nnstsi=nstsi

!cym get texta from local variable in find_param
      allocate(texta(size(textal)))
      texta=textal

      return
      end subroutine eirene_couple_param_consistency
