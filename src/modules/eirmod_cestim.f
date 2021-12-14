c   march 19, 2006:  corrected pointer for spttot in "associate_cestim"
!   20.06.07:        deallocate ESTIML and SMESTL
cdr 14.10.14:        naming of arrays in smestl adapted to those of other eirene std. dev. tallies
cdr                  two further tallies introduced (gg, stv) for stand. dev.
cdr                  of sum over strata
cdr dec 15:  species index added for eapl,empl,eipl,ephpl,eppl
cdr mar 17:  comments added
cpb Dec. 17: remove type SPECT_ARRAY, not needed in Fortran 2003

      MODULE EIRMOD_CESTIM

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_CESTIM, EIRENE_DEALLOC_CESTIM,
     P          EIRENE_ASSOCIATE_CESTIM,
     P          EIRENE_INIT_CESTIM, EIRENE_BROADCAST_CESTIM


      TYPE(EIRENE_SPECTRUM), PUBLIC, ALLOCATABLE, TARGET, SAVE ::
     .        ESTIML(:)
      TYPE(EIRENE_SPECTRUM), PUBLIC, ALLOCATABLE, TARGET, SAVE ::
     .        SMESTL(:)

      INTEGER, PUBLIC, SAVE ::
     I NESTM1, NESTM2, NESTIM

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        ESTIMV(:,:), ESTIMS(:,:)

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R          CEMETERYV(:,:), CEMETERYS(:,:)

C  NESTM1, REAL, VOLUME-AVERAGED TALLIES
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R PDENA(:,:), PDENM(:,:), PDENI(:,:), PDENPH(:,:),
     R EDENA(:,:), EDENM(:,:), EDENI(:,:), EDENPH(:,:),
     R PAEL(:),    PAAT(:,:),  PAML(:,:),  PAIO(:,:),  PAPHT(:,:),
     R PAPL(:,:),
     R PMEL(:),    PMAT(:,:),  PMML(:,:),  PMIO(:,:),  PMPHT(:,:),
     R PMPL(:,:),
     R PIEL(:),    PIAT(:,:),  PIML(:,:),  PIIO(:,:),  PIPHT(:,:),
     R PIPL(:,:),
     R PPHEL(:),   PPHAT(:,:), PPHML(:,:), PPHIO(:,:), PPHPHT(:,:),
     R PPHPL(:,:),
     R EAEL(:),  EAAT(:),  EAML(:),  EAIO(:),  EAPHT(:),  EAPL(:,:),
     R EMEL(:),  EMAT(:),  EMML(:),  EMIO(:),  EMPHT(:),  EMPL(:,:),
     R EIEL(:),  EIAT(:),  EIML(:),  EIIO(:),  EIPHT(:),  EIPL(:,:),
     R EPHEL(:), EPHAT(:), EPHML(:), EPHIO(:), EPHPHT(:), EPHPL(:,:),
c  additional (non-default) tallies: range: 57 --62. special treatment re scaling
     R ADDV(:,:),  COLV(:,:),  SNAPV(:,:),
     R COPV(:,:),  BGKV(:,:),  ALGV(:,:),
c  more recent tallies  63 --100
     R PGENA(:,:), PGENM(:,:), PGENI(:,:), PGENPH(:,:),
     R EGENA(:,:), EGENM(:,:), EGENI(:,:), EGENPH(:,:),
     R VGENA(:,:), VGENM(:,:), VGENI(:,:), VGENPH(:,:),
     R PPAT(:,:),  PPML(:,:),  PPIO(:,:),  PPPHT(:,:), PPPL(:,:),
     R EPAT(:),    EPML(:),    EPIO(:),    EPPHT(:),   EPPL(:,:),
     R VXDENA(:,:), VXDENM(:,:), VXDENI(:,:), VXDENPH(:,:),
     R VYDENA(:,:), VYDENM(:,:), VYDENI(:,:), VYDENPH(:,:),
     R VZDENA(:,:), VZDENM(:,:), VZDENI(:,:), VZDENPH(:,:),
     R MAPL(:,:), MMPL(:,:), MIPL(:,:), MPHPL(:,:)

c  POINTER FOR "A,M,I,PH"-UNIFIED SUBROUTINES
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R PDENX(:),  EDENX(:),
     R PXEL(:),   PXAT(:,:),  PXML(:,:),  PXIO(:,:), PXPL(:,:),
     R EXEL(:),   EXAT(:),    EXML(:),    EXIO(:),   EXPL(:,:),
     R VXDENX(:), VYDENX(:),  VZDENX(:),
     R MXPL(:,:), PXX(:,:),   EXX(:)

C  NESTM2, REAL, SURFACE-AVERAGED TALLIES
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R POTAT(:,:),
     R PRFAAT(:,:), PRFMAT(:,:), PRFIAT(:,:), PRFPHAT(:,:),
     R PRFPAT(:,:),
C
     R POTML(:,:),
     R PRFAML(:,:), PRFMML(:,:), PRFIML(:,:), PRFPHML(:,:),
     R PRFPML(:,:),
C
     R POTIO(:,:),
     R PRFAIO(:,:), PRFMIO(:,:), PRFIIO(:,:), PRFPHIO(:,:),
     R PRFPIO(:,:),
C
     R POTPHT(:,:),
     R PRFAPHT(:,:), PRFMPHT(:,:), PRFIPHT(:,:), PRFPHPHT(:,:),
     R PRFPPHT(:,:),
C
     R POTPL(:,:)
C
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R EOTAT(:,:),
     R ERFAAT(:,:), ERFMAT(:,:), ERFIAT(:,:), ERFPHAT(:,:),
     R ERFPAT(:,:),
C
     R EOTML(:,:),
     R ERFAML(:,:), ERFMML(:,:), ERFIML(:,:), ERFPHML(:,:),
     R ERFPML(:,:),
C
     R EOTIO(:,:),
     R ERFAIO(:,:), ERFMIO(:,:), ERFIIO(:,:), ERFPHIO(:,:),
     R ERFPIO(:,:),
C
     R EOTPHT(:,:),
     R ERFAPHT(:,:), ERFMPHT(:,:), ERFIPHT(:,:), ERFPHPHT(:,:),
     R ERFPPHT(:,:),
C
     R EOTPL(:,:),
C  FULL MATRIX: SPUTTERED FLUXES RESOLVED BY INCIDENT TYPE 
C               AND EMITTED TYPE AND SPECIES
     R SPTAAT(:,:), SPTMAT(:,:), SPTIAT(:,:), SPTPHAT(:,:),
     R SPTPAT(:,:),
     R SPTAML(:,:), SPTMML(:,:), SPTIML(:,:), SPTPHML(:,:),
     R SPTPML(:,:),
     R SPTAIO(:,:), SPTMIO(:,:), SPTIIO(:,:), SPTPHIO(:,:),
     R SPTPIO(:,:),
     R SPTAPHT(:,:), SPTMPHT(:,:), SPTIPHT(:,:), SPTPHPHT(:,:),
     R SPTPPHT(:,:),
     R SPTAPL(:,:), SPTMPL(:,:), SPTIPL(:,:), SPTPHPL(:,:),
     R SPTPPL(:,:),
! next: incident type: atoms, molecs., test ions, photons, bulk ions,
! but no emitted species index, only surface index
! analog to spttot, but: a,m,i,pl,ph, in tally name, incident type-resolved
     R sptatot(:), sptmtot(:), sptitot(:), sptphtot(:), sptpltot(:),
     R SPTTOT(:),
C
     R ADDS(:,:),  ALGS(:,:),
     R SPUMP(:,:)

C  FROM HERE: NO POINTERS ?
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NFIRST(:), NADDV(:),
     I IRESC1(:), IRESC2(:),
     I NFRSTW(:), NADDW(:)

      LOGICAL, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     L LIVTALV(:), LIVTALS(:)
      LOGICAL, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     L LMISTALV(:), LMISTALS(:)
c  logical, for each volume-averaged tally.
c  either active tally (if true) or deactivated tally, no storage (if false)
      LOGICAL, PUBLIC, POINTER, SAVE ::
     L LPDENA, LPDENM, LPDENI, LPDENPH,
     L LEDENA, LEDENM, LEDENI, LEDENPH,
     L LPAEL,  LPAAT,  LPAML,  LPAIO,   LPAPHT,  LPAPL,
     L LPMEL,  LPMAT,  LPMML,  LPMIO,   LPMPHT,  LPMPL,
     L LPIEL,  LPIAT,  LPIML,  LPIIO,   LPIPHT,  LPIPL,
     L LPPHEL, LPPHAT, LPPHML, LPPHIO,  LPPHPHT, LPPHPL,
     L LEAEL,  LEAAT,  LEAML,  LEAIO,   LEAPHT,  LEAPL,
     L LEMEL,  LEMAT,  LEMML,  LEMIO,   LEMPHT,  LEMPL,
     L LEIEL,  LEIAT,  LEIML,  LEIIO,   LEIPHT,  LEIPL,
     L LEPHEL, LEPHAT, LEPHML, LEPHIO,  LEPHPHT, LEPHPL,
c  additional tallies
     L LADDV,  LCOLV,  LSNAPV,
     L LCOPV,  LBGKV,  LALGV,
c  generation (fluid) limit tallies
     L LPGENA, LPGENM, LPGENI, LPGENPH,
     L LEGENA, LEGENM, LEGENI, LEGENPH,
     L LVGENA, LVGENM, LVGENI, LVGENPH,
c  volumetric primary source (field particle) tallies
     L LPPAT,  LPPML,  LPPIO,  LPPPHT,  LPPPL,
     L LEPAT,  LEPML,  LEPIO,  LEPPHT,  LEPPL,
c  test particle flow velocity densities
     L LVXDENA, LVXDENM, LVXDENI, LVXDENPH,
     L LVYDENA, LVYDENM, LVYDENI, LVYDENPH,
     L LVZDENA, LVZDENM, LVZDENI, LVZDENPH,
c  parallel (to B field) momentum source tallies
     L LMAPL,  LMMPL,  LMIPL,  LMPHPL

c  POINTER FOR "A,M,I,PH"-UNIFIED SUBROUTINES
      LOGICAL, PUBLIC, POINTER, SAVE ::
     L LPDENX,  LEDENX,
     L LPXEL,   LPXAT,   LPXML,   LPXIO, LPXPL,
     L LEXEL,   LEXAT,   LEXML,   LEXIO, LEXPL,
     L LVXDENX, LVYDENX, LVZDENX,
     L LMXPL,   LPXX,    LEXX

      LOGICAL, PUBLIC, POINTER, SAVE ::
     L LMSPDENA, LMSPDENM, LMSPDENI, LMSPDENPH,
     L LMSEDENA, LMSEDENM, LMSEDENI, LMSEDENPH,
     L LMSPAEL,  LMSPAAT,  LMSPAML,  LMSPAIO,   LMSPAPHT,  LMSPAPL,
     L LMSPMEL,  LMSPMAT,  LMSPMML,  LMSPMIO,   LMSPMPHT,  LMSPMPL,
     L LMSPIEL,  LMSPIAT,  LMSPIML,  LMSPIIO,   LMSPIPHT,  LMSPIPL,
     L LMSPPHEL, LMSPPHAT, LMSPPHML, LMSPPHIO,  LMSPPHPHT, LMSPPHPL,
     L LMSEAEL,  LMSEAAT,  LMSEAML,  LMSEAIO,   LMSEAPHT,  LMSEAPL,
     L LMSEMEL,  LMSEMAT,  LMSEMML,  LMSEMIO,   LMSEMPHT,  LMSEMPL,
     L LMSEIEL,  LMSEIAT,  LMSEIML,  LMSEIIO,   LMSEIPHT,  LMSEIPL,
     L LMSEPHEL, LMSEPHAT, LMSEPHML, LMSEPHIO,  LMSEPHPHT, LMSEPHPL,
     L LMSADDV,  LMSCOLV,  LMSSNAPV,
     L LMSCOPV,  LMSBGKV,  LMSALGV,
     L LMSPGENA, LMSPGENM, LMSPGENI, LMSPGENPH,
     L LMSEGENA, LMSEGENM, LMSEGENI, LMSEGENPH,
     L LMSVGENA, LMSVGENM, LMSVGENI, LMSVGENPH,
     L LMSPPAT,  LMSPPML,  LMSPPIO,  LMSPPPHT,  LMSPPPL,
     L LMSEPAT,  LMSEPML,  LMSEPIO,  LMSEPPHT,  LMSEPPL,
     L LMSVXDENA, LMSVXDENM, LMSVXDENI, LMSVXDENPH,
     L LMSVYDENA, LMSVYDENM, LMSVYDENI, LMSVYDENPH,
     L LMSVZDENA, LMSVZDENM, LMSVZDENI, LMSVZDENPH,
     L LMSMAPL,  LMSMMPL,  LMSMIPL,  LMSMPHPL

c  logical, for each surface-averaged tally, particle flux.
c  either active tally (if true) or deactivated tally, no storage (if false)
      LOGICAL, PUBLIC, POINTER, SAVE ::
     L LPOTAT,
     L LPRFAAT, LPRFMAT, LPRFIAT, LPRFPHAT,
     L LPRFPAT,
C
     L LPOTML,
     L LPRFAML, LPRFMML, LPRFIML, LPRFPHML,
     L LPRFPML,
C
     L LPOTIO,
     L LPRFAIO, LPRFMIO, LPRFIIO, LPRFPHIO,
     L LPRFPIO,
C
     L LPOTPHT,
     L LPRFAPHT, LPRFMPHT, LPRFIPHT, LPRFPHPHT,
     L LPRFPPHT,
C
     L LPOTPL

      LOGICAL, PUBLIC, POINTER, SAVE ::
     L LMSPOTAT,
     L LMSPRFAAT, LMSPRFMAT, LMSPRFIAT, LMSPRFPHAT,
     L LMSPRFPAT,
C
     L LMSPOTML,
     L LMSPRFAML, LMSPRFMML, LMSPRFIML, LMSPRFPHML,
     L LMSPRFPML,
C
     L LMSPOTIO,
     L LMSPRFAIO, LMSPRFMIO, LMSPRFIIO, LMSPRFPHIO,
     L LMSPRFPIO,
C
     L LMSPOTPHT,
     L LMSPRFAPHT, LMSPRFMPHT, LMSPRFIPHT, LMSPRFPHPHT,
     L LMSPRFPPHT,
C
     L LMSPOTPL
C
c  logical, for each surface-averaged tally, energy flux.
c  either active tally (if true) or deactivated tally, no storage (if false)
      LOGICAL, PUBLIC, POINTER, SAVE ::
     L LEOTAT,
     L LERFAAT, LERFMAT, LERFIAT, LERFPHAT, LERFPAT,
C
     L LEOTML,
     L LERFAML, LERFMML, LERFIML, LERFPHML, LERFPML,
C
     L LEOTIO,
     L LERFAIO, LERFMIO, LERFIIO, LERFPHIO, LERFPIO,
C
     L LEOTPHT,
     L LERFAPHT, LERFMPHT, LERFIPHT, LERFPHPHT, LERFPPHT,
C
     L LEOTPL,
C
     L LSPTAAT,  LSPTMAT,  LSPTIAT,  LSPTPHAT,  LSPTPAT,
     L LSPTAML,  LSPTMML,  LSPTIML,  LSPTPHML,  LSPTPML,
     L LSPTAIO,  LSPTMIO,  LSPTIIO,  LSPTPHIO,  LSPTPIO,
     L LSPTAPHT, LSPTMPHT, LSPTIPHT, LSPTPHPHT, LSPTPPHT,
     L LSPTAPL,  LSPTMPL,  LSPTIPL,  LSPTPHPL,  LSPTPPL,
     L Lsptatot, Lsptmtot, Lsptitot, Lsptpltot, Lsptphtot,
     L LSPTTOT,

     L LADDS,  LALGS,
     L LSPUMP
C
      LOGICAL, PUBLIC, POINTER, SAVE ::
     L LMSEOTAT,
     L LMSERFAAT, LMSERFMAT, LMSERFIAT, LMSERFPHAT, LMSERFPAT,
C
     L LMSEOTML,
     L LMSERFAML, LMSERFMML, LMSERFIML, LMSERFPHML, LMSERFPML,
C
     L LMSEOTIO,
     L LMSERFAIO, LMSERFMIO, LMSERFIIO, LMSERFPHIO, LMSERFPIO,
C
     L LMSEOTPHT,
     L LMSERFAPHT, LMSERFMPHT, LMSERFIPHT, LMSERFPHPHT, LMSERFPPHT,
C
     L LMSEOTPL,
C
     L LMSSPTAAT,  LMSSPTMAT,  LMSSPTIAT,  LMSSPTPHAT,  LMSSPTPAT,
     L LMSSPTAML,  LMSSPTMML,  LMSSPTIML,  LMSSPTPHML,  LMSSPTPML,
     L LMSSPTAIO,  LMSSPTMIO,  LMSSPTIIO,  LMSSPTPHIO,  LMSSPTPIO,
     L LMSSPTAPHT, LMSSPTMPHT, LMSSPTIPHT, LMSSPTPHPHT, LMSSPTPPHT,
     L LMSSPTAPL,  LMSSPTMPL,  LMSSPTIPL,  LMSSPTPHPL,  LMSSPTPPL,
     L LMSSPTTOT,
     L Lmssptatot, Lmssptmtot, Lmssptitot, Lmssptpltot, Lmssptphtot,
     L LMSADDS,  LMSALGS,
     L LMSSPUMP

!  DECLARATION AS TARGET ARRAYS FOR POINTERS USED BY UNIFIED SUBROUTINES
      LOGICAL, PUBLIC, TARGET, SAVE :: LEA, LEM, LEIO, LEPH

!  POINTER FOR "A,M,I,PH"-UNIFIED SUBROUTINES: ENERGY RATE TALLIES
      LOGICAL, PUBLIC, POINTER, SAVE :: LEX
C

C

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_CESTIM(ICAL)

      INTEGER, INTENT(IN) :: ICAL

      IF (ICAL == 1) THEN

        IF (ALLOCATED(LIVTALV)) RETURN

        ALLOCATE (LIVTALV(NTALV))
        ALLOCATE (LIVTALS(NTALS))
        ALLOCATE (LMISTALV(NTALV))
        ALLOCATE (LMISTALS(NTALS))

        ALLOCATE (NFIRST(NTALV))
        ALLOCATE (NADDV(NTALV))
        ALLOCATE (IRESC1(NTALV))
        ALLOCATE (IRESC2(NTALV))
        ALLOCATE (NFRSTW(NTALS))
        ALLOCATE (NADDW(NTALS))

        WRITE (IUNMEM,'(A,T25,I15)')
     .      ' CESTIM(1) ',4*(NTALV+NTALS) + (4*NTALV+2*NTALS)*4

      ELSE IF (ICAL == 2) THEN

        IF (ALLOCATED(ESTIMV)) RETURN

        NESTM1=NVOLTL*NRTAL
        NESTM2=NSRFTL*NLMPGS
        NESTIM=NESTM1+NESTM2

        ALLOCATE (ESTIMV(NVOLTL,NRTAL))
        ALLOCATE (ESTIMS(NSRFTL,NLMPGS))

        ALLOCATE (CEMETERYV(0:0,NRTAL))
        ALLOCATE (CEMETERYS(0:0,NLMPGS))

        WRITE (IUNMEM,'(A,T25,I15)')
     .      ' CESTIM(2) ',(NESTIM+NRTAL+NLMPGS)*8

      END IF

      CALL EIRENE_INIT_CESTIM(ICAL)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CESTIM


      SUBROUTINE EIRENE_ASSOCIATE_CESTIM

C  VOLUME-AVERAGED TALLIES:
C     if tally is active in this run     : Pointer to allocatable array ESTIMV
C     if tally is deactivated in this run: Pointer to CEMETERYV

      IF (LPDENA) THEN
        PDENA => ESTIMV(NADDV(1)+1:NADDV(2),:)
      ELSE
        PDENA => CEMETERYV(0:0,:)
      END IF
      IF (LPDENM) THEN
        PDENM => ESTIMV(NADDV(2)+1:NADDV(3),:)
      ELSE
        PDENM => CEMETERYV(0:0,:)
      END IF
      IF (LPDENI) THEN
        PDENI => ESTIMV(NADDV(3)+1:NADDV(4),:)
      ELSE
        PDENI => CEMETERYV(0:0,:)
      END IF
      IF (LPDENPH) THEN
        PDENPH => ESTIMV(NADDV(4)+1:NADDV(5),:)
      ELSE
        PDENPH => CEMETERYV(0:0,:)
      END IF

      IF (LEDENA) THEN
        EDENA => ESTIMV(NADDV(5)+1:NADDV(6),:)
      ELSE
        EDENA => CEMETERYV(0:0,:)
      END IF
      IF (LEDENM) THEN
        EDENM => ESTIMV(NADDV(6)+1:NADDV(7),:)
      ELSE
        EDENM => CEMETERYV(0:0,:)
      END IF
      IF (LEDENI) THEN
        EDENI => ESTIMV(NADDV(7)+1:NADDV(8),:)
      ELSE
        EDENI => CEMETERYV(0:0,:)
      END IF
      IF (LEDENPH) THEN
        EDENPH => ESTIMV(NADDV(8)+1:NADDV(9),:)
      ELSE
        EDENPH => CEMETERYV(0:0,:)
      END IF

      IF (LPAEL) THEN
        PAEL => ESTIMV(NADDV(10),:)
      ELSE
        PAEL => CEMETERYV(0,:)
      END IF
      IF (LPAAT) THEN
        PAAT => ESTIMV(NADDV(10)+1:NADDV(11),:)
      ELSE
        PAAT => CEMETERYV(0:0,:)
      END IF
      IF (LPAML) THEN
        PAML => ESTIMV(NADDV(11)+1:NADDV(12),:)
      ELSE
        PAML => CEMETERYV(0:0,:)
      END IF
      IF (LPAIO) THEN
        PAIO => ESTIMV(NADDV(12)+1:NADDV(13),:)
      ELSE
        PAIO => CEMETERYV(0:0,:)
      END IF
      IF (LPAPHT) THEN
        PAPHT => ESTIMV(NADDV(13)+1:NADDV(14),:)
      ELSE
        PAPHT => CEMETERYV(0:0,:)
      END IF
      IF (LPAPL) THEN
        PAPL => ESTIMV(NADDV(14)+1:NADDV(15),:)
      ELSE
        PAPL => CEMETERYV(0:0,:)
      END IF

      IF (LPMEL) THEN
        PMEL => ESTIMV(NADDV(16),:)
      ELSE
        PMEL => CEMETERYV(0,:)
      END IF
      IF (LPMAT) THEN
        PMAT => ESTIMV(NADDV(16)+1:NADDV(17),:)
      ELSE
        PMAT => CEMETERYV(0:0,:)
      END IF
      IF (LPMML) THEN
        PMML => ESTIMV(NADDV(17)+1:NADDV(18),:)
      ELSE
        PMML => CEMETERYV(0:0,:)
      END IF
      IF (LPMIO) THEN
        PMIO => ESTIMV(NADDV(18)+1:NADDV(19),:)
      ELSE
        PMIO => CEMETERYV(0:0,:)
      END IF
      IF (LPMPHT) THEN
        PMPHT => ESTIMV(NADDV(19)+1:NADDV(20),:)
      ELSE
        PMPHT => CEMETERYV(0:0,:)
      END IF
      IF (LPMPL) THEN
        PMPL => ESTIMV(NADDV(20)+1:NADDV(21),:)
      ELSE
        PMPL => CEMETERYV(0:0,:)
      END IF

      IF (LPIEL) THEN
        PIEL => ESTIMV(NADDV(22),:)
      ELSE
        PIEL => CEMETERYV(0,:)
      END IF
      IF (LPIAT) THEN
        PIAT => ESTIMV(NADDV(22)+1:NADDV(23),:)
      ELSE
        PIAT => CEMETERYV(0:0,:)
      END IF
      IF (LPIML) THEN
        PIML => ESTIMV(NADDV(23)+1:NADDV(24),:)
      ELSE
        PIML => CEMETERYV(0:0,:)
      END IF
      IF (LPIIO) THEN
        PIIO => ESTIMV(NADDV(24)+1:NADDV(25),:)
      ELSE
        PIIO => CEMETERYV(0:0,:)
      END IF
      IF (LPIPHT) THEN
        PIPHT => ESTIMV(NADDV(25)+1:NADDV(26),:)
      ELSE
        PIPHT => CEMETERYV(0:0,:)
      END IF
      IF (LPIPL) THEN
        PIPL => ESTIMV(NADDV(26)+1:NADDV(27),:)
      ELSE
        PIPL => CEMETERYV(0:0,:)
      END IF

      IF (LPPHEL) THEN
        PPHEL => ESTIMV(NADDV(28),:)
      ELSE
        PPHEL => CEMETERYV(0,:)
      END IF
      IF (LPPHAT) THEN
        PPHAT => ESTIMV(NADDV(28)+1:NADDV(29),:)
      ELSE
        PPHAT => CEMETERYV(0:0,:)
      END IF
      IF (LPPHML) THEN
        PPHML => ESTIMV(NADDV(29)+1:NADDV(30),:)
      ELSE
        PPHML => CEMETERYV(0:0,:)
      END IF
      IF (LPPHIO) THEN
        PPHIO => ESTIMV(NADDV(30)+1:NADDV(31),:)
      ELSE
        PPHIO => CEMETERYV(0:0,:)
      END IF
      IF (LPPHPHT) THEN
        PPHPHT => ESTIMV(NADDV(31)+1:NADDV(32),:)
      ELSE
        PPHPHT => CEMETERYV(0:0,:)
      END IF
      IF (LPPHPL) THEN
        PPHPL => ESTIMV(NADDV(32)+1:NADDV(33),:)
      ELSE
        PPHPL => CEMETERYV(0:0,:)
      END IF

      IF (LEAEL) THEN
        EAEL => ESTIMV(NADDV(34),:)
      ELSE
        EAEL => CEMETERYV(0,:)
      END IF
      IF (LEAAT) THEN
        EAAT => ESTIMV(NADDV(35),:)
      ELSE
        EAAT => CEMETERYV(0,:)
      END IF
      IF (LEAML) THEN
        EAML => ESTIMV(NADDV(36),:)
      ELSE
        EAML => CEMETERYV(0,:)
      END IF
      IF (LEAIO) THEN
        EAIO => ESTIMV(NADDV(37),:)
      ELSE
        EAIO => CEMETERYV(0,:)
      END IF
      IF (LEAPHT) THEN
        EAPHT => ESTIMV(NADDV(38),:)
      ELSE
        EAPHT => CEMETERYV(0,:)
      END IF
      IF (LEAPL) THEN
        EAPL => ESTIMV(NADDV(38)+1:NADDV(39),:)
      ELSE
        EAPL => CEMETERYV(0:0,:)
      END IF

      IF (LEMEL) THEN
        EMEL => ESTIMV(NADDV(40),:)
      ELSE
        EMEL => CEMETERYV(0,:)
      END IF
      IF (LEMAT) THEN
        EMAT => ESTIMV(NADDV(41),:)
      ELSE
        EMAT => CEMETERYV(0,:)
      END IF
      IF (LEMML) THEN
        EMML => ESTIMV(NADDV(42),:)
      ELSE
        EMML => CEMETERYV(0,:)
      END IF
      IF (LEMIO) THEN
        EMIO => ESTIMV(NADDV(43),:)
      ELSE
        EMIO => CEMETERYV(0,:)
      END IF
      IF (LEMPHT) THEN
        EMPHT => ESTIMV(NADDV(44),:)
      ELSE
        EMPHT => CEMETERYV(0,:)
      END IF
      IF (LEMPL) THEN
        EMPL => ESTIMV(NADDV(44)+1:NADDV(45),:)
      ELSE
        EMPL => CEMETERYV(0:0,:)
      END IF

      IF (LEIEL) THEN
        EIEL => ESTIMV(NADDV(46),:)
      ELSE
        EIEL => CEMETERYV(0,:)
      END IF
      IF (LEIAT) THEN
        EIAT => ESTIMV(NADDV(47),:)
      ELSE
        EIAT => CEMETERYV(0,:)
      END IF
      IF (LEIML) THEN
        EIML => ESTIMV(NADDV(48),:)
      ELSE
        EIML => CEMETERYV(0,:)
      END IF
      IF (LEIIO) THEN
        EIIO => ESTIMV(NADDV(49),:)
      ELSE
        EIIO => CEMETERYV(0,:)
      END IF
      IF (LEIPHT) THEN
        EIPHT => ESTIMV(NADDV(50),:)
      ELSE
        EIPHT => CEMETERYV(0,:)
      END IF
      IF (LEIPL) THEN
        EIPL => ESTIMV(NADDV(50)+1:NADDV(51),:)
      ELSE
        EIPL => CEMETERYV(0:0,:)
      END IF

      IF (LEPHEL) THEN
        EPHEL => ESTIMV(NADDV(52),:)
      ELSE
        EPHEL => CEMETERYV(0,:)
      END IF
      IF (LEPHAT) THEN
        EPHAT => ESTIMV(NADDV(53),:)
      ELSE
        EPHAT => CEMETERYV(0,:)
      END IF
      IF (LEPHML) THEN
        EPHML => ESTIMV(NADDV(54),:)
      ELSE
        EPHML => CEMETERYV(0,:)
      END IF
      IF (LEPHIO) THEN
        EPHIO => ESTIMV(NADDV(55),:)
      ELSE
        EPHIO => CEMETERYV(0,:)
      END IF
      IF (LEPHPHT) THEN
        EPHPHT => ESTIMV(NADDV(56),:)
      ELSE
        EPHPHT => CEMETERYV(0,:)
      END IF
      IF (LEPHPL) THEN
        EPHPL => ESTIMV(NADDV(56)+1:NADDV(57),:)
      ELSE
        EPHPL => CEMETERYV(0:0,:)
      END IF

      IF (LADDV) THEN
c  ntala =57
        ADDV => ESTIMV(NADDV(NTALA)+1:NADDV(NTALA+1),:)
      ELSE
        ADDV => CEMETERYV(0:0,:)
      END IF
      IF (LCOLV) THEN
c  ntalc =58
        COLV => ESTIMV(NADDV(NTALC)+1:NADDV(NTALC+1),:)
      ELSE
        COLV => CEMETERYV(0:0,:)
      END IF
      IF (LSNAPV) THEN
c  ntalt =59
        SNAPV => ESTIMV(NADDV(NTALT)+1:NADDV(NTALT+1),:)
      ELSE
        SNAPV => CEMETERYV(0:0,:)
      END IF
      IF (LCOPV) THEN
c  ntalm =60
        COPV => ESTIMV(NADDV(NTALM)+1:NADDV(NTALM+1),:)
      ELSE
        COPV => CEMETERYV(0:0,:)
      END IF
      IF (LBGKV) THEN
c  ntalb =61
        BGKV => ESTIMV(NADDV(NTALB)+1:NADDV(NTALB+1),:)
      ELSE
        BGKV => CEMETERYV(0:0,:)
      END IF
      IF (LALGV) THEN
c  ntalr =62
        ALGV => ESTIMV(NADDV(NTALR)+1:NADDV(NTALR+1),:)
      ELSE
        ALGV => CEMETERYV(0:0,:)
      END IF

      IF (LPGENA) THEN
        PGENA => ESTIMV(NADDV(63)+1:NADDV(64),:)
      ELSE
        PGENA => CEMETERYV(0:0,:)
      END IF
      IF (LPGENM) THEN
        PGENM => ESTIMV(NADDV(64)+1:NADDV(65),:)
      ELSE
        PGENM => CEMETERYV(0:0,:)
      END IF
      IF (LPGENI) THEN
        PGENI => ESTIMV(NADDV(65)+1:NADDV(66),:)
      ELSE
        PGENI => CEMETERYV(0:0,:)
      END IF
      IF (LPGENPH) THEN
        PGENPH => ESTIMV(NADDV(66)+1:NADDV(67),:)
      ELSE
        PGENPH => CEMETERYV(0:0,:)
      END IF
      IF (LEGENA) THEN
        EGENA => ESTIMV(NADDV(67)+1:NADDV(68),:)
      ELSE
        EGENA => CEMETERYV(0:0,:)
      END IF
      IF (LEGENM) THEN
        EGENM => ESTIMV(NADDV(68)+1:NADDV(69),:)
      ELSE
        EGENM => CEMETERYV(0:0,:)
      END IF
      IF (LEGENI) THEN
        EGENI => ESTIMV(NADDV(69)+1:NADDV(70),:)
      ELSE
        EGENI => CEMETERYV(0:0,:)
      END IF
      IF (LEGENPH) THEN
        EGENPH => ESTIMV(NADDV(70)+1:NADDV(71),:)
      ELSE
        EGENPH => CEMETERYV(0:0,:)
      END IF
      IF (LVGENA) THEN
        VGENA => ESTIMV(NADDV(71)+1:NADDV(72),:)
      ELSE
        VGENA => CEMETERYV(0:0,:)
      END IF
      IF (LVGENM) THEN
        VGENM => ESTIMV(NADDV(72)+1:NADDV(73),:)
      ELSE
        VGENM => CEMETERYV(0:0,:)
      END IF
      IF (LVGENI) THEN
        VGENI => ESTIMV(NADDV(73)+1:NADDV(74),:)
      ELSE
        VGENI => CEMETERYV(0:0,:)
      END IF
      IF (LVGENPH) THEN
        VGENPH => ESTIMV(NADDV(74)+1:NADDV(75),:)
      ELSE
        VGENPH => CEMETERYV(0:0,:)
      END IF

      IF (LPPAT) THEN
        PPAT => ESTIMV(NADDV(75)+1:NADDV(76),:)
      ELSE
        PPAT => CEMETERYV(0:0,:)
      END IF
      IF (LPPML) THEN
        PPML => ESTIMV(NADDV(76)+1:NADDV(77),:)
      ELSE
        PPML => CEMETERYV(0:0,:)
      END IF
      IF (LPPIO) THEN
        PPIO => ESTIMV(NADDV(77)+1:NADDV(78),:)
      ELSE
        PPIO => CEMETERYV(0:0,:)
      END IF
      IF (LPPPHT) THEN
        PPPHT => ESTIMV(NADDV(78)+1:NADDV(79),:)
      ELSE
        PPPHT => CEMETERYV(0:0,:)
      END IF
      IF (LPPPL) THEN
        PPPL => ESTIMV(NADDV(79)+1:NADDV(80),:)
      ELSE
        PPPL => CEMETERYV(0:0,:)
      END IF

      IF (LEPAT) THEN
        EPAT => ESTIMV(NADDV(81),:)
      ELSE
        EPAT => CEMETERYV(0,:)
      END IF
      IF (LEPML) THEN
        EPML => ESTIMV(NADDV(82),:)
      ELSE
        EPML => CEMETERYV(0,:)
      END IF
      IF (LEPIO) THEN
        EPIO => ESTIMV(NADDV(83),:)
      ELSE
        EPIO => CEMETERYV(0,:)
      END IF
      IF (LEPPHT) THEN
        EPPHT => ESTIMV(NADDV(84),:)
      ELSE
        EPPHT => CEMETERYV(0,:)
      END IF
      IF (LEPPL) THEN
        EPPL => ESTIMV(NADDV(84)+1:NADDV(85),:)
      ELSE
        EPPL => CEMETERYV(0:0,:)
      END IF
      IF (LVXDENA) THEN
        VXDENA => ESTIMV(NADDV(85)+1:NADDV(86),:)
      ELSE
        VXDENA => CEMETERYV(0:0,:)
      END IF
      IF (LVXDENM) THEN
        VXDENM => ESTIMV(NADDV(86)+1:NADDV(87),:)
      ELSE
        VXDENM => CEMETERYV(0:0,:)
      END IF
      IF (LVXDENI) THEN
        VXDENI => ESTIMV(NADDV(87)+1:NADDV(88),:)
      ELSE
        VXDENI => CEMETERYV(0:0,:)
      END IF
      IF (LVXDENPH) THEN
        VXDENPH => ESTIMV(NADDV(88)+1:NADDV(89),:)
      ELSE
        VXDENPH => CEMETERYV(0:0,:)
      END IF
      IF (LVYDENA) THEN
        VYDENA => ESTIMV(NADDV(89)+1:NADDV(90),:)
      ELSE
        VYDENA => CEMETERYV(0:0,:)
      END IF
      IF (LVYDENM) THEN
        VYDENM => ESTIMV(NADDV(90)+1:NADDV(91),:)
      ELSE
        VYDENM => CEMETERYV(0:0,:)
      END IF
      IF (LVYDENI) THEN
        VYDENI => ESTIMV(NADDV(91)+1:NADDV(92),:)
      ELSE
        VYDENI => CEMETERYV(0:0,:)
      END IF
      IF (LVYDENPH) THEN
        VYDENPH => ESTIMV(NADDV(92)+1:NADDV(93),:)
      ELSE
        VYDENPH => CEMETERYV(0:0,:)
      END IF
      IF (LVZDENA) THEN
        VZDENA => ESTIMV(NADDV(93)+1:NADDV(94),:)
      ELSE
        VZDENA => CEMETERYV(0:0,:)
      END IF
      IF (LVZDENM) THEN
        VZDENM => ESTIMV(NADDV(94)+1:NADDV(95),:)
      ELSE
        VZDENM => CEMETERYV(0:0,:)
      END IF
      IF (LVZDENI) THEN
        VZDENI => ESTIMV(NADDV(95)+1:NADDV(96),:)
      ELSE
        VZDENI => CEMETERYV(0:0,:)
      END IF
      IF (LVZDENPH) THEN
        VZDENPH => ESTIMV(NADDV(96)+1:NADDV(97),:)
      ELSE
        VZDENPH => CEMETERYV(0:0,:)
      END IF
      IF (LMAPL) THEN
        MAPL => ESTIMV(NADDV(97)+1:NADDV(98),:)
      ELSE
        MAPL => CEMETERYV(0:0,:)
      END IF
      IF (LMMPL) THEN
        MMPL => ESTIMV(NADDV(98)+1:NADDV(99),:)
      ELSE
        MMPL => CEMETERYV(0:0,:)
      END IF
      IF (LMIPL) THEN
        MIPL => ESTIMV(NADDV(99)+1:NADDV(100),:)
      ELSE
        MIPL => CEMETERYV(0:0,:)
      END IF
      IF (LMPHPL) THEN
!pb   I would have expected the compiler to find the length of
!pb   array ESTIMV automatically but ifort version 12.0.4 does not
!pb 20.06.2012        MPHPL => ESTIMV(NADDV(100)+1: ,:)
        MPHPL => ESTIMV(NADDV(100)+1:NVOLTL ,:)
      ELSE
        MPHPL => CEMETERYV(0:0,:)
      END IF

C  SURFACE-AVERAGED TALLIES:
C     if tally is active in this run     : Pointer to allocatable array ESTIMS
C     if tally is deactivated in this run: Pointer to CEMETERYS
      IF (LPOTAT) THEN
        POTAT => ESTIMS(1:NADDW(2),:)
      ELSE
        POTAT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFAAT) THEN
        PRFAAT => ESTIMS(NADDW(2)+1:NADDW(3),:)
      ELSE
        PRFAAT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFMAT) THEN
        PRFMAT => ESTIMS(NADDW(3)+1:NADDW(4),:)
      ELSE
        PRFMAT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFIAT) THEN
        PRFIAT => ESTIMS(NADDW(4)+1:NADDW(5),:)
      ELSE
        PRFIAT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPHAT) THEN
        PRFPHAT => ESTIMS(NADDW(5)+1:NADDW(6),:)
      ELSE
        PRFPHAT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPAT) THEN
        PRFPAT => ESTIMS(NADDW(6)+1:NADDW(7),:)
      ELSE
        PRFPAT => CEMETERYS(0:0,:)
      END IF
C
      IF (LPOTML) THEN
        POTML => ESTIMS(NADDW(7)+1:NADDW(8),:)
      ELSE
        POTML => CEMETERYS(0:0,:)
      END IF
      IF (LPRFAML) THEN
        PRFAML => ESTIMS(NADDW(8)+1:NADDW(9),:)
      ELSE
        PRFAML => CEMETERYS(0:0,:)
      END IF
      IF (LPRFMML) THEN
        PRFMML => ESTIMS(NADDW(9)+1:NADDW(10),:)
      ELSE
        PRFMML => CEMETERYS(0:0,:)
      END IF
      IF (LPRFIML) THEN
        PRFIML => ESTIMS(NADDW(10)+1:NADDW(11),:)
      ELSE
        PRFIML => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPHML) THEN
        PRFPHML => ESTIMS(NADDW(11)+1:NADDW(12),:)
      ELSE
        PRFPHML => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPML) THEN
        PRFPML => ESTIMS(NADDW(12)+1:NADDW(13),:)
      ELSE
        PRFPML => CEMETERYS(0:0,:)
      END IF
C
      IF (LPOTIO) THEN
        POTIO => ESTIMS(NADDW(13)+1:NADDW(14),:)
      ELSE
        POTIO => CEMETERYS(0:0,:)
      END IF
      IF (LPRFAIO) THEN
        PRFAIO => ESTIMS(NADDW(14)+1:NADDW(15),:)
      ELSE
        PRFAIO => CEMETERYS(0:0,:)
      END IF
      IF (LPRFMIO) THEN
        PRFMIO => ESTIMS(NADDW(15)+1:NADDW(16),:)
      ELSE
        PRFMIO => CEMETERYS(0:0,:)
      END IF
      IF (LPRFIIO) THEN
        PRFIIO => ESTIMS(NADDW(16)+1:NADDW(17),:)
      ELSE
        PRFIIO => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPHIO) THEN
        PRFPHIO => ESTIMS(NADDW(17)+1:NADDW(18),:)
      ELSE
        PRFPHIO => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPIO) THEN
        PRFPIO => ESTIMS(NADDW(18)+1:NADDW(19),:)
      ELSE
        PRFPIO => CEMETERYS(0:0,:)
      END IF
C
      IF (LPOTPHT) THEN
        POTPHT => ESTIMS(NADDW(19)+1:NADDW(20),:)
      ELSE
        POTPHT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFAPHT) THEN
        PRFAPHT => ESTIMS(NADDW(20)+1:NADDW(21),:)
      ELSE
        PRFAPHT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFMPHT) THEN
        PRFMPHT => ESTIMS(NADDW(21)+1:NADDW(22),:)
      ELSE
        PRFMPHT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFIPHT) THEN
        PRFIPHT => ESTIMS(NADDW(22)+1:NADDW(23),:)
      ELSE
        PRFIPHT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPHPHT) THEN
        PRFPHPHT =>ESTIMS(NADDW(23)+1:NADDW(24),:)
      ELSE
        PRFPHPHT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPPHT) THEN
        PRFPPHT => ESTIMS(NADDW(24)+1:NADDW(25),:)
      ELSE
        PRFPPHT => CEMETERYS(0:0,:)
      END IF
C
      IF (LPOTPL) THEN
        POTPL => ESTIMS(NADDW(25)+1:NADDW(26),:)
      ELSE
        POTPL => CEMETERYS(0:0,:)
      END IF
C
      IF (LEOTAT) THEN
        EOTAT => ESTIMS(NADDW(26)+1:NADDW(27),:)
      ELSE
        EOTAT => CEMETERYS(0:0,:)
      END IF
      IF (LERFAAT) THEN
        ERFAAT => ESTIMS(NADDW(27)+1:NADDW(28),:)
      ELSE
        ERFAAT => CEMETERYS(0:0,:)
      END IF
      IF (LERFMAT) THEN
        ERFMAT => ESTIMS(NADDW(28)+1:NADDW(29),:)
      ELSE
        ERFMAT => CEMETERYS(0:0,:)
      END IF
      IF (LERFIAT) THEN
        ERFIAT => ESTIMS(NADDW(29)+1:NADDW(30),:)
      ELSE
        ERFIAT => CEMETERYS(0:0,:)
      END IF
      IF (LERFPHAT) THEN
        ERFPHAT => ESTIMS(NADDW(30)+1:NADDW(31),:)
      ELSE
        ERFPHAT => CEMETERYS(0:0,:)
      END IF
      IF (LERFPAT) THEN
        ERFPAT => ESTIMS(NADDW(31)+1:NADDW(32),:)
      ELSE
        ERFPAT => CEMETERYS(0:0,:)
      END IF
C
      IF (LEOTML) THEN
        EOTML => ESTIMS(NADDW(32)+1:NADDW(33),:)
      ELSE
        EOTML => CEMETERYS(0:0,:)
      END IF
      IF (LERFAML) THEN
        ERFAML => ESTIMS(NADDW(33)+1:NADDW(34),:)
      ELSE
        ERFAML => CEMETERYS(0:0,:)
      END IF
      IF (LERFMML) THEN
        ERFMML => ESTIMS(NADDW(34)+1:NADDW(35),:)
      ELSE
        ERFMML => CEMETERYS(0:0,:)
      END IF
      IF (LERFIML) THEN
        ERFIML => ESTIMS(NADDW(35)+1:NADDW(36),:)
      ELSE
        ERFIML => CEMETERYS(0:0,:)
      END IF
      IF (LERFPHML) THEN
        ERFPHML => ESTIMS(NADDW(36)+1:NADDW(37),:)
      ELSE
        ERFPHML => CEMETERYS(0:0,:)
      END IF
      IF (LERFPML) THEN
        ERFPML => ESTIMS(NADDW(37)+1:NADDW(38),:)
      ELSE
        ERFPML => CEMETERYS(0:0,:)
      END IF
C
      IF (LEOTIO) THEN
        EOTIO => ESTIMS(NADDW(38)+1:NADDW(39),:)
      ELSE
        EOTIO => CEMETERYS(0:0,:)
      END IF
      IF (LERFAIO) THEN
        ERFAIO => ESTIMS(NADDW(39)+1:NADDW(40),:)
      ELSE
        ERFAIO => CEMETERYS(0:0,:)
      END IF
      IF (LERFMIO) THEN
        ERFMIO => ESTIMS(NADDW(40)+1:NADDW(41),:)
      ELSE
        ERFMIO => CEMETERYS(0:0,:)
      END IF
      IF (LERFIIO) THEN
        ERFIIO => ESTIMS(NADDW(41)+1:NADDW(42),:)
      ELSE
        ERFIIO => CEMETERYS(0:0,:)
      END IF
      IF (LERFPHIO) THEN
        ERFPHIO => ESTIMS(NADDW(42)+1:NADDW(43),:)
      ELSE
        ERFPHIO => CEMETERYS(0:0,:)
      END IF
      IF (LERFPIO) THEN
        ERFPIO => ESTIMS(NADDW(43)+1:NADDW(44),:)
      ELSE
        ERFPIO => CEMETERYS(0:0,:)
      END IF
C
      IF (LEOTPHT) THEN
        EOTPHT => ESTIMS(NADDW(44)+1:NADDW(45),:)
      ELSE
        EOTPHT => CEMETERYS(0:0,:)
      END IF
      IF (LERFAPHT) THEN
        ERFAPHT => ESTIMS(NADDW(45)+1:NADDW(46),:)
      ELSE
        ERFAPHT => CEMETERYS(0:0,:)
      END IF
      IF (LERFMPHT) THEN
        ERFMPHT => ESTIMS(NADDW(46)+1:NADDW(47),:)
      ELSE
        ERFMPHT => CEMETERYS(0:0,:)
      END IF
      IF (LERFIPHT) THEN
        ERFIPHT => ESTIMS(NADDW(47)+1:NADDW(48),:)
      ELSE
        ERFIPHT => CEMETERYS(0:0,:)
      END IF
      IF (LERFPHPHT) THEN
        ERFPHPHT =>ESTIMS(NADDW(48)+1:NADDW(49),:)
      ELSE
        ERFPHPHT => CEMETERYS(0:0,:)
      END IF
      IF (LERFPPHT) THEN
        ERFPPHT => ESTIMS(NADDW(49)+1:NADDW(50),:)
      ELSE
        ERFPPHT => CEMETERYS(0:0,:)
      END IF
C
      IF (LEOTPL) THEN
        EOTPL => ESTIMS(NADDW(50)+1:NADDW(51),:)
      ELSE
        EOTPL => CEMETERYS(0:0,:)
      END IF
C
C  SPUTTER TALLIES  
C
C  EMITTED TYPE: ATOMS
      IF (LSPTAAT) THEN
        SPTAAT => ESTIMS(NADDW(51)+1:NADDW(52),:)
      ELSE
        SPTAAT => CEMETERYS(0:0,:)
      END IF
      IF (LSPTMAT) THEN
        SPTMAT => ESTIMS(NADDW(52)+1:NADDW(53),:)
      ELSE
        SPTMAT => CEMETERYS(0:0,:)
      END IF
      IF (LSPTIAT) THEN
        SPTIAT => ESTIMS(NADDW(53)+1:NADDW(54),:)
      ELSE
        SPTIAT => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPHAT) THEN
        SPTPHAT => ESTIMS(NADDW(54)+1:NADDW(55),:)
      ELSE
        SPTPHAT => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPAT) THEN
        SPTPAT => ESTIMS(NADDW(55)+1:NADDW(56),:)
      ELSE
        SPTPAT => CEMETERYS(0:0,:)
      END IF

C  EMITTED TYPE: MOLECULES
      IF (LSPTAML) THEN
        SPTAML => ESTIMS(NADDW(56)+1:NADDW(57),:)
      ELSE
        SPTAML => CEMETERYS(0:0,:)
      END IF
      IF (LSPTMML) THEN
        SPTMML => ESTIMS(NADDW(57)+1:NADDW(58),:)
      ELSE
        SPTMML => CEMETERYS(0:0,:)
      END IF
      IF (LSPTIML) THEN
        SPTIML => ESTIMS(NADDW(58)+1:NADDW(59),:)
      ELSE
        SPTIML => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPHML) THEN
        SPTPHML => ESTIMS(NADDW(59)+1:NADDW(60),:)
      ELSE
        SPTPHML => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPML) THEN
        SPTPML => ESTIMS(NADDW(60)+1:NADDW(61),:)
      ELSE
        SPTPML => CEMETERYS(0:0,:)
      END IF

C  EMITTED TYPE: TEST IONS
      IF (LSPTAIO) THEN
        SPTAIO => ESTIMS(NADDW(61)+1:NADDW(62),:)
      ELSE
        SPTAIO => CEMETERYS(0:0,:)
      END IF
      IF (LSPTMIO) THEN
        SPTMIO => ESTIMS(NADDW(62)+1:NADDW(63),:)
      ELSE
        SPTMIO => CEMETERYS(0:0,:)
      END IF
      IF (LSPTIIO) THEN
        SPTIIO => ESTIMS(NADDW(63)+1:NADDW(64),:)
      ELSE
        SPTIIO => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPHIO) THEN
        SPTPHIO => ESTIMS(NADDW(64)+1:NADDW(65),:)
      ELSE
        SPTPHIO => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPIO) THEN
        SPTPIO => ESTIMS(NADDW(65)+1:NADDW(66),:)
      ELSE
        SPTPIO => CEMETERYS(0:0,:)
      END IF

C  EMITTED TYPE: PHOTONS
      IF (LSPTAPHT) THEN
        SPTAPHT => ESTIMS(NADDW(66)+1:NADDW(67),:)
      ELSE
        SPTAPHT => CEMETERYS(0:0,:)
      END IF
      IF (LSPTMPHT) THEN
        SPTMPHT => ESTIMS(NADDW(67)+1:NADDW(68),:)
      ELSE
        SPTMPHT => CEMETERYS(0:0,:)
      END IF
      IF (LSPTIPHT) THEN
        SPTIPHT => ESTIMS(NADDW(68)+1:NADDW(69),:)
      ELSE
        SPTIPHT => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPHPHT) THEN
        SPTPHPHT => ESTIMS(NADDW(69)+1:NADDW(70),:)
      ELSE
        SPTPHPHT => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPPHT) THEN
        SPTPPHT => ESTIMS(NADDW(70)+1:NADDW(71),:)
      ELSE
        SPTPPHT => CEMETERYS(0:0,:)
      END IF

C  EMITTED TYPE: BULK IONS
      IF (LSPTAPL) THEN
        SPTAPL => ESTIMS(NADDW(71)+1:NADDW(72),:)
      ELSE
        SPTAPL => CEMETERYS(0:0,:)
      END IF
      IF (LSPTMPL) THEN
        SPTMPL => ESTIMS(NADDW(72)+1:NADDW(73),:)
      ELSE
        SPTMPL => CEMETERYS(0:0,:)
      END IF
      IF (LSPTIPL) THEN
        SPTIPL => ESTIMS(NADDW(73)+1:NADDW(74),:)
      ELSE
        SPTIPL => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPHPL) THEN
        SPTPHPL => ESTIMS(NADDW(74)+1:NADDW(75),:)
      ELSE
        SPTPHPL => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPPL) THEN
        SPTPPL => ESTIMS(NADDW(75)+1:NADDW(76),:)
      ELSE
        SPTPPL => CEMETERYS(0:0,:)
      END IF

C  TOTALS
      IF (LSPTATOT) THEN
        SPTATOT => ESTIMS(NADDW(76)+1,:)
      ELSE
        SPTATOT => CEMETERYS(0,:)
      END IF
      IF (LSPTMTOT) THEN
        SPTMTOT => ESTIMS(NADDW(77)+1,:)
      ELSE
        SPTMTOT => CEMETERYS(0,:)
      END IF
      IF (LSPTITOT) THEN
        SPTITOT => ESTIMS(NADDW(78)+1,:)
      ELSE
        SPTITOT => CEMETERYS(0,:)
      END IF
      IF (LSPTPHTOT) THEN
        SPTPHTOT => ESTIMS(NADDW(79)+1,:)
      ELSE
        SPTPHTOT => CEMETERYS(0,:)
      END IF
      IF (LSPTPLTOT) THEN
        SPTPLTOT => ESTIMS(NADDW(80)+1,:)
      ELSE
        SPTPLTOT => CEMETERYS(0,:)
      END IF
      IF (LSPTTOT) THEN
        SPTTOT => ESTIMS(NADDW(81)+1,:)
      ELSE
        SPTTOT => CEMETERYS(0,:)
      END IF

C  ADDITIONAL SURFACE TALLIES
      IF (LADDS) THEN
        ADDS => ESTIMS(NADDW(82)+1:NADDW(83),:)
      ELSE
        ADDS => CEMETERYS(0:0,:)
      END IF

C  ALGEBRAIC SURFACE TALLIES
      IF (LALGS) THEN
        ALGS => ESTIMS(NADDW(83)+1:NADDW(84),:)
      ELSE
        ALGS => CEMETERYS(0:0,:)
      END IF

C  PUMPED FLUXES
      IF (LSPUMP) THEN
cdr this coding now seems inconsistent with corresponding coding for last vol. av. tally.
        SPUMP => ESTIMS(NADDW(84)+1:,:)
cdr     SPUMP => ESTIMS(NADDW(84)+1:NADDW(85),:)   ! clearer, safer ? see above, for vol. av. tallies
      ELSE
        SPUMP => CEMETERYS(0:0,:)
      END IF

      RETURN
      END SUBROUTINE EIRENE_ASSOCIATE_CESTIM

      SUBROUTINE EIRENE_DEALLOC_CESTIM
      INTEGER :: I
C
      IF (ALLOCATED(ESTIMV)) THEN
         DEALLOCATE (ESTIMV)
         DEALLOCATE (ESTIMS)
         IF (NADSPC > 0) THEN
c  spectra tallies: standard deviation
           DO I=1,NADSPC
             DEALLOCATE(ESTIML(I)%SPC)
             IF (ASSOCIATED(ESTIML(I)%SDV)) THEN
               DEALLOCATE(ESTIML(I)%SDV)
               DEALLOCATE(ESTIML(I)%SGM)
               DEALLOCATE(ESTIML(I)%GG)
               DEALLOCATE(ESTIML(I)%STV)
             END IF
c  spectra tallies: standard deviation for sum over strata, intermediate storage
             IF (NSMSTRA > 0) THEN
               DEALLOCATE(SMESTL(I)%SPC)
               IF (ASSOCIATED(SMESTL(I)%SDV)) THEN
                 DEALLOCATE(SMESTL(I)%SDV)
                 DEALLOCATE(SMESTL(I)%SGM)
                 DEALLOCATE(SMESTL(I)%GG)
                 DEALLOCATE(SMESTL(I)%STV)
               END IF
             END IF
           END DO
           DEALLOCATE (ESTIML)
cdr  heinke frerichs, juli 2016
           IF (NSMSTRA > 0) DEALLOCATE (SMESTL)
         ELSE
           IF (ALLOCATED(ESTIML)) DEALLOCATE (ESTIML)
           IF (ALLOCATED(SMESTL)) DEALLOCATE (SMESTL)
         END IF
         DEALLOCATE (NFIRST)
         DEALLOCATE (NADDV)
         DEALLOCATE (IRESC1)
         DEALLOCATE (IRESC2)
         DEALLOCATE (NFRSTW)
         DEALLOCATE (NADDW)

         DEALLOCATE (CEMETERYV)
         DEALLOCATE (CEMETERYS)
      END IF

      IF (ALLOCATED(LIVTALV)) THEN
         DEALLOCATE (LIVTALV)
         DEALLOCATE (LIVTALS)
         DEALLOCATE (LMISTALV)
         DEALLOCATE (LMISTALS)
      END IF

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CESTIM


      SUBROUTINE EIRENE_INIT_CESTIM(ICAL)

      INTEGER, INTENT(IN) :: ICAL
C
      IF (ICAL == 1) THEN

        LIVTALV = .TRUE.
        LIVTALS = .TRUE.
        LMISTALV = .FALSE.
        LMISTALS = .FALSE.

! volume-averaged tallies

        LPDENA   => LIVTALV(1)
        LPDENM   => LIVTALV(2)
        LPDENI   => LIVTALV(3)
        LPDENPH  => LIVTALV(4)
        LEDENA   => LIVTALV(5)
        LEDENM   => LIVTALV(6)
        LEDENI   => LIVTALV(7)
        LEDENPH  => LIVTALV(8)
        LPAEL    => LIVTALV(9)
        LPAAT    => LIVTALV(10)
        LPAML    => LIVTALV(11)
        LPAIO    => LIVTALV(12)
        LPAPHT   => LIVTALV(13)
        LPAPL    => LIVTALV(14)
        LPMEL    => LIVTALV(15)
        LPMAT    => LIVTALV(16)
        LPMML    => LIVTALV(17)
        LPMIO    => LIVTALV(18)
        LPMPHT   => LIVTALV(19)
        LPMPL    => LIVTALV(20)
        LPIEL    => LIVTALV(21)
        LPIAT    => LIVTALV(22)
        LPIML    => LIVTALV(23)
        LPIIO    => LIVTALV(24)
        LPIPHT   => LIVTALV(25)
        LPIPL    => LIVTALV(26)
        LPPHEL   => LIVTALV(27)
        LPPHAT   => LIVTALV(28)
        LPPHML   => LIVTALV(29)
        LPPHIO   => LIVTALV(30)
        LPPHPHT  => LIVTALV(31)
        LPPHPL   => LIVTALV(32)
        LEAEL    => LIVTALV(33)
        LEAAT    => LIVTALV(34)
        LEAML    => LIVTALV(35)
        LEAIO    => LIVTALV(36)
        LEAPHT   => LIVTALV(37)
        LEAPL    => LIVTALV(38)
        LEMEL    => LIVTALV(39)
        LEMAT    => LIVTALV(40)
        LEMML    => LIVTALV(41)
        LEMIO    => LIVTALV(42)
        LEMPHT   => LIVTALV(43)
        LEMPL    => LIVTALV(44)
        LEIEL    => LIVTALV(45)
        LEIAT    => LIVTALV(46)
        LEIML    => LIVTALV(47)
        LEIIO    => LIVTALV(48)
        LEIPHT   => LIVTALV(49)
        LEIPL    => LIVTALV(50)
        LEPHEL   => LIVTALV(51)
        LEPHAT   => LIVTALV(52)
        LEPHML   => LIVTALV(53)
        LEPHIO   => LIVTALV(54)
        LEPHPHT  => LIVTALV(55)
        LEPHPL   => LIVTALV(56)
        LADDV    => LIVTALV(57)
        LCOLV    => LIVTALV(58)
        LSNAPV   => LIVTALV(59)
        LCOPV    => LIVTALV(60)
        LBGKV    => LIVTALV(61)
        LALGV    => LIVTALV(62)
        LPGENA   => LIVTALV(63)
        LPGENM   => LIVTALV(64)
        LPGENI   => LIVTALV(65)
        LPGENPH  => LIVTALV(66)
        LEGENA   => LIVTALV(67)
        LEGENM   => LIVTALV(68)
        LEGENI   => LIVTALV(69)
        LEGENPH  => LIVTALV(70)
        LVGENA   => LIVTALV(71)
        LVGENM   => LIVTALV(72)
        LVGENI   => LIVTALV(73)
        LVGENPH  => LIVTALV(74)
        LPPAT    => LIVTALV(75)
        LPPML    => LIVTALV(76)
        LPPIO    => LIVTALV(77)
        LPPPHT   => LIVTALV(78)
        LPPPL    => LIVTALV(79)
        LEPAT    => LIVTALV(80)
        LEPML    => LIVTALV(81)
        LEPIO    => LIVTALV(82)
        LEPPHT   => LIVTALV(83)
        LEPPL    => LIVTALV(84)
        LVXDENA  => LIVTALV(85)
        LVXDENM  => LIVTALV(86)
        LVXDENI  => LIVTALV(87)
        LVXDENPH => LIVTALV(88)
        LVYDENA  => LIVTALV(89)
        LVYDENM  => LIVTALV(90)
        LVYDENI  => LIVTALV(91)
        LVYDENPH => LIVTALV(92)
        LVZDENA  => LIVTALV(93)
        LVZDENM  => LIVTALV(94)
        LVZDENI  => LIVTALV(95)
        LVZDENPH => LIVTALV(96)
        LMAPL    => LIVTALV(97)
        LMMPL    => LIVTALV(98)
        LMIPL    => LIVTALV(99)
        LMPHPL   => LIVTALV(100)

        LMSPDENA   => LMISTALV(1)
        LMSPDENM   => LMISTALV(2)
        LMSPDENI   => LMISTALV(3)
        LMSPDENPH  => LMISTALV(4)
        LMSEDENA   => LMISTALV(5)
        LMSEDENM   => LMISTALV(6)
        LMSEDENI   => LMISTALV(7)
        LMSEDENPH  => LMISTALV(8)
        LMSPAEL    => LMISTALV(9)
        LMSPAAT    => LMISTALV(10)
        LMSPAML    => LMISTALV(11)
        LMSPAIO    => LMISTALV(12)
        LMSPAPHT   => LMISTALV(13)
        LMSPAPL    => LMISTALV(14)
        LMSPMEL    => LMISTALV(15)
        LMSPMAT    => LMISTALV(16)
        LMSPMML    => LMISTALV(17)
        LMSPMIO    => LMISTALV(18)
        LMSPMPHT   => LMISTALV(19)
        LMSPMPL    => LMISTALV(20)
        LMSPIEL    => LMISTALV(21)
        LMSPIAT    => LMISTALV(22)
        LMSPIML    => LMISTALV(23)
        LMSPIIO    => LMISTALV(24)
        LMSPIPHT   => LMISTALV(25)
        LMSPIPL    => LMISTALV(26)
        LMSPPHEL   => LMISTALV(27)
        LMSPPHAT   => LMISTALV(28)
        LMSPPHML   => LMISTALV(29)
        LMSPPHIO   => LMISTALV(30)
        LMSPPHPHT  => LMISTALV(31)
        LMSPPHPL   => LMISTALV(32)
        LMSEAEL    => LMISTALV(33)
        LMSEAAT    => LMISTALV(34)
        LMSEAML    => LMISTALV(35)
        LMSEAIO    => LMISTALV(36)
        LMSEAPHT   => LMISTALV(37)
        LMSEAPL    => LMISTALV(38)
        LMSEMEL    => LMISTALV(39)
        LMSEMAT    => LMISTALV(40)
        LMSEMML    => LMISTALV(41)
        LMSEMIO    => LMISTALV(42)
        LMSEMPHT   => LMISTALV(43)
        LMSEMPL    => LMISTALV(44)
        LMSEIEL    => LMISTALV(45)
        LMSEIAT    => LMISTALV(46)
        LMSEIML    => LMISTALV(47)
        LMSEIIO    => LMISTALV(48)
        LMSEIPHT   => LMISTALV(49)
        LMSEIPL    => LMISTALV(50)
        LMSEPHEL   => LMISTALV(51)
        LMSEPHAT   => LMISTALV(52)
        LMSEPHML   => LMISTALV(53)
        LMSEPHIO   => LMISTALV(54)
        LMSEPHPHT  => LMISTALV(55)
        LMSEPHPL   => LMISTALV(56)
        LMSADDV    => LMISTALV(57)
        LMSCOLV    => LMISTALV(58)
        LMSSNAPV   => LMISTALV(59)
        LMSCOPV    => LMISTALV(60)
        LMSBGKV    => LMISTALV(61)
        LMSALGV    => LMISTALV(62)
        LMSPGENA   => LMISTALV(63)
        LMSPGENM   => LMISTALV(64)
        LMSPGENI   => LMISTALV(65)
        LMSPGENPH  => LMISTALV(66)
        LMSEGENA   => LMISTALV(67)
        LMSEGENM   => LMISTALV(68)
        LMSEGENI   => LMISTALV(69)
        LMSEGENPH  => LMISTALV(70)
        LMSVGENA   => LMISTALV(71)
        LMSVGENM   => LMISTALV(72)
        LMSVGENI   => LMISTALV(73)
        LMSVGENPH  => LMISTALV(74)
        LMSPPAT    => LMISTALV(75)
        LMSPPML    => LMISTALV(76)
        LMSPPIO    => LMISTALV(77)
        LMSPPPHT   => LMISTALV(78)
        LMSPPPL    => LMISTALV(79)
        LMSEPAT    => LMISTALV(80)
        LMSEPML    => LMISTALV(81)
        LMSEPIO    => LMISTALV(82)
        LMSEPPHT   => LMISTALV(83)
        LMSEPPL    => LMISTALV(84)
        LMSVXDENA  => LMISTALV(85)
        LMSVXDENM  => LMISTALV(86)
        LMSVXDENI  => LMISTALV(87)
        LMSVXDENPH => LMISTALV(88)
        LMSVYDENA  => LMISTALV(89)
        LMSVYDENM  => LMISTALV(90)
        LMSVYDENI  => LMISTALV(91)
        LMSVYDENPH => LMISTALV(92)
        LMSVZDENA  => LMISTALV(93)
        LMSVZDENM  => LMISTALV(94)
        LMSVZDENI  => LMISTALV(95)
        LMSVZDENPH => LMISTALV(96)
        LMSMAPL    => LMISTALV(97)
        LMSMMPL    => LMISTALV(98)
        LMSMIPL    => LMISTALV(99)
        LMSMPHPL   => LMISTALV(100)

! surface-averaged tallies

        LPOTAT    => LIVTALS(1)
        LPRFAAT   => LIVTALS(2)
        LPRFMAT   => LIVTALS(3)
        LPRFIAT   => LIVTALS(4)
        LPRFPHAT  => LIVTALS(5)
        LPRFPAT   => LIVTALS(6)
        LPOTML    => LIVTALS(7)
        LPRFAML   => LIVTALS(8)
        LPRFMML   => LIVTALS(9)
        LPRFIML   => LIVTALS(10)
        LPRFPHML  => LIVTALS(11)
        LPRFPML   => LIVTALS(12)
        LPOTIO    => LIVTALS(13)
        LPRFAIO   => LIVTALS(14)
        LPRFMIO   => LIVTALS(15)
        LPRFIIO   => LIVTALS(16)
        LPRFPHIO  => LIVTALS(17)
        LPRFPIO   => LIVTALS(18)
        LPOTPHT   => LIVTALS(19)
        LPRFAPHT  => LIVTALS(20)
        LPRFMPHT  => LIVTALS(21)
        LPRFIPHT  => LIVTALS(22)
        LPRFPHPHT => LIVTALS(23)
        LPRFPPHT  => LIVTALS(24)
        LPOTPL    => LIVTALS(25)
        LEOTAT    => LIVTALS(26)
        LERFAAT   => LIVTALS(27)
        LERFMAT   => LIVTALS(28)
        LERFIAT   => LIVTALS(29)
        LERFPHAT  => LIVTALS(30)
        LERFPAT   => LIVTALS(31)
        LEOTML    => LIVTALS(32)
        LERFAML   => LIVTALS(33)
        LERFMML   => LIVTALS(34)
        LERFIML   => LIVTALS(35)
        LERFPHML  => LIVTALS(36)
        LERFPML   => LIVTALS(37)
        LEOTIO    => LIVTALS(38)
        LERFAIO   => LIVTALS(39)
        LERFMIO   => LIVTALS(40)
        LERFIIO   => LIVTALS(41)
        LERFPHIO  => LIVTALS(42)
        LERFPIO   => LIVTALS(43)
        LEOTPHT   => LIVTALS(44)
        LERFAPHT  => LIVTALS(45)
        LERFMPHT  => LIVTALS(46)
        LERFIPHT  => LIVTALS(47)
        LERFPHPHT => LIVTALS(48)
        LERFPPHT  => LIVTALS(49)
        LEOTPL    => LIVTALS(50)
C  SPUTTER TALLIES:  51 -- 75
        LSPTAAT   => LIVTALS(51)
        LSPTAML   => LIVTALS(52)
        LSPTAIO   => LIVTALS(53)
        LSPTAPHT  => LIVTALS(54)
        LSPTAPL   => LIVTALS(55)
        LSPTMAT   => LIVTALS(56)
        LSPTMML   => LIVTALS(57)
        LSPTMIO   => LIVTALS(58)
        LSPTMPHT  => LIVTALS(59)
        LSPTMPL   => LIVTALS(60)
        LSPTIAT   => LIVTALS(61)
        LSPTIML   => LIVTALS(62)
        LSPTIIO   => LIVTALS(63)
        LSPTIPHT  => LIVTALS(64)
        LSPTIPL   => LIVTALS(65)
        LSPTPHAT  => LIVTALS(66)
        LSPTPHML  => LIVTALS(67)
        LSPTPHIO  => LIVTALS(68)
        LSPTPHPHT => LIVTALS(69)
        LSPTPHPL  => LIVTALS(70)
        LSPTPAT   => LIVTALS(71)
        LSPTPML   => LIVTALS(72)
        LSPTPIO   => LIVTALS(73)
        LSPTPPHT  => LIVTALS(74)
        LSPTPPL   => LIVTALS(75)
        LSPTATOT  => LIVTALS(76)
        LSPTMTOT  => LIVTALS(77)
        LSPTITOT  => LIVTALS(78)
        LSPTPHTOT => LIVTALS(79)
        LSPTPLTOT => LIVTALS(80)
C
        LSPTTOT   => LIVTALS(81)
C
        LADDS     => LIVTALS(82)
        LALGS     => LIVTALS(83)
        LSPUMP    => LIVTALS(84)

        LMSPOTAT    => LMISTALS(1)
        LMSPRFAAT   => LMISTALS(2)
        LMSPRFMAT   => LMISTALS(3)
        LMSPRFIAT   => LMISTALS(4)
        LMSPRFPHAT  => LMISTALS(5)
        LMSPRFPAT   => LMISTALS(6)
        LMSPOTML    => LMISTALS(7)
        LMSPRFAML   => LMISTALS(8)
        LMSPRFMML   => LMISTALS(9)
        LMSPRFIML   => LMISTALS(10)
        LMSPRFPHML  => LMISTALS(11)
        LMSPRFPML   => LMISTALS(12)
        LMSPOTIO    => LMISTALS(13)
        LMSPRFAIO   => LMISTALS(14)
        LMSPRFMIO   => LMISTALS(15)
        LMSPRFIIO   => LMISTALS(16)
        LMSPRFPHIO  => LMISTALS(17)
        LMSPRFPIO   => LMISTALS(18)
        LMSPOTPHT   => LMISTALS(19)
        LMSPRFAPHT  => LMISTALS(20)
        LMSPRFMPHT  => LMISTALS(21)
        LMSPRFIPHT  => LMISTALS(22)
        LMSPRFPHPHT => LMISTALS(23)
        LMSPRFPPHT  => LMISTALS(24)
        LMSPOTPL    => LMISTALS(25)
        LMSEOTAT    => LMISTALS(26)
        LMSERFAAT   => LMISTALS(27)
        LMSERFMAT   => LMISTALS(28)
        LMSERFIAT   => LMISTALS(29)
        LMSERFPHAT  => LMISTALS(30)
        LMSERFPAT   => LMISTALS(31)
        LMSEOTML    => LMISTALS(32)
        LMSERFAML   => LMISTALS(33)
        LMSERFMML   => LMISTALS(34)
        LMSERFIML   => LMISTALS(35)
        LMSERFPHML  => LMISTALS(36)
        LMSERFPML   => LMISTALS(37)
        LMSEOTIO    => LMISTALS(38)
        LMSERFAIO   => LMISTALS(39)
        LMSERFMIO   => LMISTALS(40)
        LMSERFIIO   => LMISTALS(41)
        LMSERFPHIO  => LMISTALS(42)
        LMSERFPIO   => LMISTALS(43)
        LMSEOTPHT   => LMISTALS(44)
        LMSERFAPHT  => LMISTALS(45)
        LMSERFMPHT  => LMISTALS(46)
        LMSERFIPHT  => LMISTALS(47)
        LMSERFPHPHT => LMISTALS(48)
        LMSERFPPHT  => LMISTALS(49)
        LMSEOTPL    => LMISTALS(50)
C  SPUTTER TALLIES: 51 -- 75
        LMSSPTAAT   => LMISTALS(51)
        LMSSPTAML   => LMISTALS(52)
        LMSSPTAIO   => LMISTALS(53)
        LMSSPTAPHT  => LMISTALS(54)
        LMSSPTAPL   => LMISTALS(55)
        LMSSPTMAT   => LMISTALS(56)
        LMSSPTMML   => LMISTALS(57)
        LMSSPTMIO   => LMISTALS(58)
        LMSSPTMPHT  => LMISTALS(59)
        LMSSPTMPL   => LMISTALS(60)
        LMSSPTIAT   => LMISTALS(61)
        LMSSPTIML   => LMISTALS(62)
        LMSSPTIIO   => LMISTALS(63)
        LMSSPTIPHT  => LMISTALS(64)
        LMSSPTIPL   => LMISTALS(65)
        LMSSPTPHAT  => LMISTALS(66)
        LMSSPTPHML  => LMISTALS(67)
        LMSSPTPHIO  => LMISTALS(68)
        LMSSPTPHPHT => LMISTALS(69)
        LMSSPTPHPL  => LMISTALS(70)
        LMSSPTPAT   => LMISTALS(71)
        LMSSPTPML   => LMISTALS(72)
        LMSSPTPIO   => LMISTALS(73)
        LMSSPTPPHT  => LMISTALS(74)
        LMSSPTPPL   => LMISTALS(75)
C
        LMSSPTATOT  => LMISTALS(76)
        LMSSPTMTOT  => LMISTALS(77)
        LMSSPTITOT  => LMISTALS(78)
        LMSSPTPHTOT => LMISTALS(79)
        LMSSPTPLTOT => LMISTALS(80)
        LMSSPTTOT   => LMISTALS(81)
        LMSADDS     => LMISTALS(82)
        LMSALGS     => LMISTALS(83)
        LMSSPUMP    => LMISTALS(84)

        NFIRST = 0
        NADDV  = 0
        IRESC1 = 0
        IRESC2 = 0
        NFRSTW = 0
        NADDW  = 0

      ELSE IF (ICAL == 2) THEN

        ESTIMV = 0._DP
        ESTIMS = 0._DP

        CEMETERYV = 0._DP
        CEMETERYS = 0._DP

      END IF

      RETURN
      END SUBROUTINE EIRENE_INIT_CESTIM


      SUBROUTINE EIRENE_BROADCAST_CESTIM(ME)
      USE EIRMOD_MPI
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER, I, NSPS

      IF (ME /= 0) THEN
        CALL EIRENE_ALLOC_CESTIM(1)
      END IF

      CALL MPI_BCAST (NFIRST,NTALV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NADDV,NTALV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IRESC1,NTALV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IRESC2,NTALV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFRSTW,NTALS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NADDW,NTALS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

c  active and inactive tallies:
C  OUTPUT:
      CALL MPI_BCAST (LIVTALV,NTALV,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LIVTALS,NTALS,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LMISTALV,NTALV,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LMISTALS,NTALS,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LEA,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LEM,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LEIO,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LEPH,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

      IF (ME .NE. 0) THEN
        CALL EIRENE_ALLOC_CESTIM(2)
        CALL EIRENE_ASSOCIATE_CESTIM
      END IF

      IF (NADSPC > 0) THEN
         IF (ME .NE. 0) THEN
            IF (.NOT.ALLOCATED(ESTIML)) THEN
              ALLOCATE(ESTIML(NADSPC))
              ALLOCATE(SMESTL(NADSPC))
            END IF
         END IF
         DO I=1,NADSPC
           CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%SPCMIN,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%SPCMAX,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%SPCDEL,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%SPCDELI,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%ESP_MIN,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%ESP_MAX,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%ESP_00,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%SPC_XPLT,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%SPC_YPLT,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%SPC_SAME,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%SPCVX,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%SPCVY,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%SPCVZ,1,MPI_REAL8,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%NSPC,1,MPI_INTEGER,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%ISPCTYP,1,MPI_INTEGER,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%ISPCSRF,1,MPI_INTEGER,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%IPRTYP,1,MPI_INTEGER,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%IPRSP,1,MPI_INTEGER,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%IMETSP,1,MPI_INTEGER,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%ISRFCLL,1,MPI_INTEGER,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%IDIREC,1,MPI_INTEGER,0,
     .                     MPI_COMM_WORLD,ier)
           CALL MPI_BCAST (ESTIML(I)%LOG,1,MPI_LOGICAL,0,
     .                     MPI_COMM_WORLD,ier)

           IF (ME .NE. 0) THEN
             NSPS = ESTIML(I)%NSPC

             IF (.NOT.ASSOCIATED(ESTIML(I)%SPC)) THEN
               ALLOCATE(ESTIML(I)%SPC(0:NSPS+1))
               ALLOCATE(ESTIML(I)%SDV(0:NSPS+1))
               ALLOCATE(ESTIML(I)%SGM(0:NSPS+1))
               ALLOCATE(ESTIML(I)%STV(0:NSPS+1))
               ALLOCATE(ESTIML(I)%GG(0:NSPS+1))
C  variances for sum over strata
               ALLOCATE(SMESTL(I)%SPC(0:NSPS+1))
               ALLOCATE(SMESTL(I)%SDV(0:NSPS+1))
               ALLOCATE(SMESTL(I)%SGM(0:NSPS+1))
               ALLOCATE(SMESTL(I)%STV(0:NSPS+1))
               ALLOCATE(SMESTL(I)%GG(0:NSPS+1))
             END IF
             ESTIML(I)%SPC = 0._DP
             SMESTL(I) = ESTIML(I)
           END IF
         END DO
      ELSE
        IF (ME .NE. 0) THEN
          IF (.NOT.ALLOCATED(ESTIML)) THEN
            ALLOCATE(ESTIML(1))
          END IF
        END IF
      END IF
      
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      END SUBROUTINE EIRENE_BROADCAST_CESTIM

      END MODULE EIRMOD_CESTIM
