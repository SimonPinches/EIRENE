cdr  may 2017:  preparing for storage reduction by elimination of unnecessary input tallies:
cdr             commenting, 
cdr             lusr, musr, nusr, nplpr1, nplpr2, nsfprm made local, 
cdr             rather than public
cpb  Dec. 2017: remove type SPECT_ARRAY, not needed in Fortran 2003
cpb  jan 2018:  remove unused arrays TEDTEDX, TEDTEDY, TEDTEDZ
cdr             added: nstpi (formerly: coutou), naddcor
cdr  oct 2018:  bvin moved into LBSMO condition
cdr             POT  moved into LESMO condition
cdr             tbd:  BXPERP, BYPERP:  move into LBSMO condition
cdr             missing:  dealloc_corners  ??

      MODULE EIRMOD_COMUSR
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
 
      IMPLICIT NONE
 
      PRIVATE
 
      PUBLIC :: EIRENE_ALLOC_COMUSR, EIRENE_DEALLOC_COMUSR,
     P          EIRENE_INIT_COMUSR, 
     P          EIRENE_ALLOC_CORNERS,
     P          EIRENE_ASSOCIATE_COMUSR,
     P          EIRENE_COMUSR_REINIT
 
      INTEGER, SAVE ::
     P NPLPR1, NSFPRM, NPLPR2  ! internal, not public. former storage tests in setprm are abandoned
      INTEGER, PUBLIC, SAVE ::
     P NPLPRM  ! nplprm, is also used in setprm, for a storage test.
c 
      INTEGER, SAVE ::               
     P NUSR,   MUSR,   LUSR             ! also only local in this module, apparently
      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R         PLSTLS(:,:)
 
      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R          CEMETERYP(:,:)
 
      REAL(DP), POINTER, PUBLIC, SAVE ::
C  NPLPRM, REAL.
C  THE FIRST NPLPR1 DATA ARE PRIMARY INPUT PROFILES, SET IN SUBROUTINE PLASMA
     R        TEIN(:),        TIIN(:,:),      DEIN(:),     DIIN(:,:),
     R        VXIN(:,:),      VYIN(:,:),      VZIN(:,:),
     R        BXIN(:),        BYIN(:),        BZIN(:),     BFIN(:),
     R        ADIN(:,:),      VOL(:),         WGHT(:,:),
     R        EXIN(:),        EYIN(:),        EZIN(:),     EFIN(:),
     R        POT(:),
c  derived from primary input profils, in subr. PLASMA_DERIV  
c  (strictly: DEIN is also a derived tally) :
     R        BXPERP(:),      BYPERP(:),
     R        BVIN(:,:),      PARMOM(:,:),    EDRIFT(:,:),

c  optional: gradients of all scalar input tallies.
c  For the vectorial input tallies (V_IN, B_IN, E_IN) these gradients
c  then provide the full dyadic (all nine components).

     R        DTEDX(:),       DTEDY(:),       DTEDZ(:),
     R        DTIDX(:,:),     DTIDY(:,:),     DTIDZ(:,:),
     R        DDEDX(:),       DDEDY(:),       DDEDZ(:),
     R        DDIDX(:,:),     DDIDY(:,:),     DDIDZ(:,:),
     R        DVXDX(:,:),     DVXDY(:,:),     DVXDZ(:,:),
     R        DVYDX(:,:),     DVYDY(:,:),     DVYDZ(:,:),
     R        DVZDX(:,:),     DVZDY(:,:),     DVZDZ(:,:),
     R        DBXDX(:),       DBXDY(:),       DBXDZ(:),
     R        DBYDX(:),       DBYDY(:),       DBYDZ(:),
     R        DBZDX(:),       DBZDY(:),       DBZDZ(:),
     R        DBFDX(:),       DBFDY(:),       DBFDZ(:),
     R        DADINDX(:,:),   DADINDY(:,:),   DADINDZ(:,:),
     R        DVOLDX(:),      DVOLDY(:),      DVOLDZ(:),
     R        DWGHTDX(:,:),   DWGHTDY(:,:),   DWGHTDZ(:,:),
     R        DEXDX(:),       DEXDY(:),       DEXDZ(:),
     R        DEYDX(:),       DEYDY(:),       DEYDZ(:),
     R        DEZDX(:),       DEZDY(:),       DEZDZ(:),
     R        DEFDX(:),       DEFDY(:),       DEFDZ(:),
     R        DPOTDX(:),      DPOTDY(:),      DPOTDZ(:),
c 
     R        DBXPERPDX(:),   DBXPERPDY(:),   DBXPERPDZ(:),
     R        DBYPERPDX(:),   DBYPERPDY(:),   DBYPERPDZ(:),
     R        DBVINDX(:,:),   DBVINDY(:,:),   DBVINDZ(:,:),
     R        DPARMOMDX(:,:), DPARMOMDY(:,:), DPARMOMDZ(:,:),
     R        DEDRIFTDX(:,:), DEDRIFTDY(:,:), DEDRIFTDZ(:,:)

      REAL(DP), ALLOCATABLE, PUBLIC, SAVE ::
C  NPLPR2, REAL.
C  THIS SECOND SET OF DATA ARE DERIVED INPUT PROFILES, SET IN SUBROUTINE PLASMA_DERIV
C  (STRICTLY ALSO DEIN (ELECTRON DENSITY) FROM THE NPLPR1 BLOCK ABOVE
C   IS SUCH A DERIVED QUANTITY)  
     R        TEINL(:),  TIINL(:,:),  DEINL(:),  DIINL(:,:),
C  NSFPRM
     R        FLXOUT(:), SAREA(:),
C
     R        DIOD(:),   DATD(:),     DMLD(:),   DPLD(:),    DPHD(:),
     R        DION(:),   DATM(:),     DMOL(:),   DPLS(:),    DPHOT(:)
 
!  DECLARATION AS TARGET ARRAYS FOR POINTERS USED BY UNIFIED SUBROUTINES
      REAL(DP), TARGET, ALLOCATABLE, PUBLIC, SAVE ::
     R        RMASSI(:), RMASSA(:),   RMASSM(:), RMASSPH(:), RMASSP(:)

!     POINTER FOR UNIFIED SUBROUTINES
      REAL(DP), POINTER, PUBLIC, SAVE :: RMASSX

C     PLASMA PROFILES ON CELL VERTICES
      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        CORNER_PROFILES(:,:)

c  storage for setting tallies at cell vertices, rather than cell centres,
c  for interpolations
      REAL(DP), POINTER, PUBLIC, SAVE ::
     .        TEINCORNER(:),   TIINCORNER(:,:), DEINCORNER(:),
     .        DIINCORNER(:,:),
     .        VXINCORNER(:,:), VYINCORNER(:,:), VZINCORNER(:,:),
     .        BXINCORNER(:),   BYINCORNER(:),   BZINCORNER(:),
     .        BFINCORNER(:),   
     .        EXCORNER(:),     EYCORNER(:),     EZCORNER(:),
     .        EFCORNER(:),     POTCORNER(:),
     .        ADCORNER(:,:),   VOLCORNER(:),    WGHTCORNER(:,:),

     .        BXPERPCORNER(:), BYPERPCORNER(:), BVINCORNER(:,:),
     .        PARMOMCORNER(:,:), EDRIFTCORNER(:,:)
 
      REAL(DP), PUBLIC, SAVE :: TVAC, DVAC, VVAC, ALLOC
 
      CHARACTER(8), ALLOCATABLE, PUBLIC, SAVE :: TEXTS(:)
 
C  MUSR, INTEGER
      INTEGER, PUBLIC, SAVE ::
     I         NSPH  , NPHOTI, NPHOTIM, NPHOTI_IN,
     I         NSPA  , NATMI,  NATMIM,  NATMI_IN,
     I         NSPAM , NMOLI,  NMOLIM,  NMOLI_IN,
     I         NSPAMI, NIONI,  NIONIM,  NIONI_IN,
     I         NSPTOT, NPLSI,  NPLSIM,  NPLSI_IN,
     I         NSNVI,  NCPVI,  NADVI,   NBGVI,
     I         NALVI,  NCLVI,  NADSI,   NALSI, NAINI, NBITS
      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         NMASSA(:), NCHARA(:), NFOLA(:),  NGENA(:),
     I         NMASSM(:), NCHARM(:), NFOLM(:),  NGENM(:),
     I         NMASSI(:), NCHARI(:), NCHRGI(:), NFOLI(:), NGENI(:),
     I         NMASSP(:), NCHARP(:), NCHRGP(:),
     I         NFOLPH(:), NGENPH(:),
     I         NPRT(:),   ISPEZ(:,:,:,:,:,:),     ISPEZI(:,:),
     I         MPLSTI(:), MPLSV(:)
      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         ISPZ_BACK(:,:)
 
C  LUSR, LOGICAL
      LOGICAL, ALLOCATABLE, PUBLIC, SAVE ::
     L         LGVAC(:,:), LGDFT(:)

      LOGICAL, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     L         LIVTALI(:)
c
C  SMOOTHED INPUT TALLIES (interpolation from cell vertices into cell)
      LOGICAL, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     L         LSMOPRO(:)

      LOGICAL, PUBLIC, POINTER, SAVE ::
     L         LTESMO,     LTISMO,     LDESMO,    LDISMO,
     L         LVXSMO,     LVYSMO,     LVZSMO,
     L         LBXSMO,     LBYSMO,     LBZSMO,    LBFSMO,
     L         LADSMO,     LVOLSMO,    LWGHTSMO,
     L         LEXSMO,     LEYSMO,     LEZSMO,    LEFSMO,
     L         LPOTSMO,
C
     L         LBXPSMO,    LBYPSMO,
     L         LBVSMO,     LPARMOMSMO, LEDRIFTSMO

      LOGICAL, PUBLIC, SAVE ::
     L         LDSMO, LVSMO,  LBSMO,  LESMO

c  pointer to LIVTALI: active or inactive input tallies
      LOGICAL, PUBLIC, POINTER, SAVE ::
c  background, drifting maxwellian parameters
     L         LTEIN,      LTIIN,      LDEIN,     LDIIN,
     L         LVXIN,      LVYIN,      LVZIN,
c  magn. field
     L         LBXIN,      LBYIN,      LBZIN,     LBFIN,
c  additional stuff....
     L         LADIN,      LVOL,       LWGHT,
c  electr. field
     L         LEXIN,      LEYIN,      LEZIN,     LEFIN,
     L         LPOT,
c  derived tallies
     L         LBXPERP,    LBYPERP,
     L         LBVIN,      LPARMOM,    LEDRIFT,
c  gradient tallies
     L         LDTEDX,     LDTEDY,     LDTEDZ,
     L         LDTIDX,     LDTIDY,     LDTIDZ,
     L         LDDEDX,     LDDEDY,     LDDEDZ,
     L         LDDIDX,     LDDIDY,     LDDIDZ,
     L         LDVXDX,     LDVXDY,     LDVXDZ,
     L         LDVYDX,     LDVYDY,     LDVYDZ,
     L         LDVZDX,     LDVZDY,     LDVZDZ,
     L         LDBXDX,     LDBXDY,     LDBXDZ,
     L         LDBYDX,     LDBYDY,     LDBYDZ,
     L         LDBZDX,     LDBZDY,     LDBZDZ,
     L         LDBFDX,     LDBFDY,     LDBFDZ,
     L         LDADINDX,   LDADINDY,   LDADINDZ,
     L         LDVOLDX,    LDVOLDY,    LDVOLDZ,
     L         LDWGHTDX,   LDWGHTDY,   LDWGHTDZ,
     L         LDEXDX,     LDEXDY,     LDEXDZ,
     L         LDEYDX,     LDEYDY,     LDEYDZ,
     L         LDEZDX,     LDEZDY,     LDEZDZ,
     L         LDEFDX,     LDEFDY,     LDEFDZ,
     L         LDPOTDX,    LDPOTDY,    LDPOTDZ,
C  gradients of derived tallies
     L         LDBXPERPDX, LDBXPERPDY, LDBXPERPDZ,
     L         LDBYPERPDX, LDBYPERPDY, LDBYPERPDZ,
     L         LDBVINDX,   LDBVINDY,   LDBVINDZ,
     L         LDPARMOMDX, LDPARMOMDY, LDPARMOMDZ,
     L         LDEDRIFTDX, LDEDRIFTDY, LDEDRIFTDZ
 
      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         IADVE(:),  IADVS(:), IADVT(:),  IADRC(:),
     I         ICLVE(:),  ICLVS(:), ICLVT(:),  ICLRC(:),
     I         ISNVE(:),  ISNVS(:), ISNVT(:),  ISNRC(:),
     I         ICPVE(:),  ICPVS(:), ICPVT(:),  ICPRC(:),
     I         IBGVE(:),  IBGVS(:), IBGVT(:),  IBGRC(:),
     I         IADSE(:),  IADSS(:), IADST(:),  IADSC(:),
     I         NFRSTP(:), NADDP(:), NFSTPI(:), NADDCOR(:),
     I         NSPAN(:),  NSPEN(:),
     I         NSPANW(:), NSPENW(:)

      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         INTLOPTS(:)
 
      INTEGER, PUBLIC, SAVE ::
     I         NPRLL, NMODE,  NTCPU,
     I         NFILE, NFILEN, NFILEM, NFILEL, NFILEK, NFILEJ,
     I         NITER, IITER,  NTIME,  ITIMV
 
!      TYPE(SPECT_ARRAY), PUBLIC, ALLOCATABLE, SAVE :: BACK_SPEC(:)
      TYPE(EIRENE_SPECTRUM), PUBLIC, ALLOCATABLE, SAVE :: BACK_SPEC(:)
      LOGICAL, PUBLIC, ALLOCATABLE, SAVE :: LSPCCLL(:)
 
 
      CONTAINS
 
 
      SUBROUTINE EIRENE_ALLOC_COMUSR (ICAL)
 
      INTEGER, INTENT(IN) :: ICAL
 

      IF (.NOT.ALLOCATED(LSMOPRO))  ALLOCATE (LSMOPRO(NTALG))

      IF (ICAL == 1) THEN
 
        IF (ALLOCATED(RMASSI)) RETURN

        NPLPR2= 3*(NATM+NMOL+NION+NPLS)+4+NSPZ+2*NPHOT ! species (test particle and background) related data
 
        MUSR=4*NATM+4*NMOL+5*NION+3*NPLS+30+NSPZ+2*NPHOT+
     .       6*(1+NPHOTP)*(1+NATMP)*(1+NMOLP)*(1+NIONP)*(1+NPLSP)+NSPZ*6
     .       +2*NPLS+NSPZ*NPLS
     .       +4*NADV+4*NCLV+4*NSNV+4*NADS+4*NTALI+NTALG+2*NTALV+2*NTALS

        LUSR=NRAD*(NPLS+2)+2*NRAD+NTALI

        
        ALLOCATE (RMASSA(MAX(1,NATM)))
        ALLOCATE (RMASSM(MAX(1,NMOL)))
        ALLOCATE (RMASSI(MAX(1,NION)))
        ALLOCATE (RMASSPH(MAX(1,NPHOT)))
        ALLOCATE (RMASSP(MAX(1,NPLS)))

        ALLOCATE (DIOD(MAX(1,NION)))
        ALLOCATE (DATD(MAX(1,NATM)))
        ALLOCATE (DMLD(MAX(1,NMOL)))
        ALLOCATE (DPLD(MAX(1,NPLS)))
        ALLOCATE (DPHD(MAX(1,NPHOT)))
        ALLOCATE (DION(MAX(1,NION)))
        ALLOCATE (DATM(MAX(1,NATM)))
        ALLOCATE (DMOL(MAX(1,NMOL)))
        ALLOCATE (DPLS(MAX(1,NPLS)))
        ALLOCATE (DPHOT(MAX(1,NPHOT)))
 
        ALLOCATE (TEXTS(NSPZ))
 
c  integer  species and background tally data

        ALLOCATE (NMASSA(MAX(1,NATM)))
        ALLOCATE (NCHARA(MAX(1,NATM)))
        ALLOCATE (NFOLA(MAX(1,NATM)))
        ALLOCATE (NGENA(MAX(1,NATM)))
        ALLOCATE (NMASSM(MAX(1,NMOL)))
        ALLOCATE (NCHARM(MAX(1,NMOL)))
        ALLOCATE (NFOLM(MAX(1,NMOL)))
        ALLOCATE (NGENM(MAX(1,NMOL)))
        ALLOCATE (NMASSP(MAX(1,NPLS)))
        ALLOCATE (NCHARP(MAX(1,NPLS)))
        ALLOCATE (NCHRGP(MAX(1,NPLS)))
        ALLOCATE (NMASSI(MAX(1,NION)))
        ALLOCATE (NCHARI(MAX(1,NION)))
        ALLOCATE (NCHRGI(MAX(1,NION)))
        ALLOCATE (NFOLI(MAX(1,NION)))
        ALLOCATE (NGENI(MAX(1,NION)))
        ALLOCATE (NFOLPH(MAX(1,NPHOT)))
        ALLOCATE (NGENPH(MAX(1,NPHOT)))
        ALLOCATE (NPRT(MAX(1,NSPZ)))
        ALLOCATE (ISPEZ(-1:4,0:NPHOTP,0:NATMP,0:NMOLP,0:NIONP,0:NPLSP))
        ALLOCATE (ISPEZI(NSPZ,-1:4))

        ALLOCATE (MPLSTI(MAX(1,NPLS)))
        ALLOCATE (MPLSV(MAX(1,NPLS)))
        ALLOCATE (ISPZ_BACK(NSPZ,NPLS))
        ALLOCATE (IADVE(MAX(1,NADV)))
        ALLOCATE (IADVS(MAX(1,NADV)))
        ALLOCATE (IADVT(MAX(1,NADV)))
        ALLOCATE (IADRC(MAX(1,NADV)))
        ALLOCATE (ICLVE(MAX(1,NCLV)))
        ALLOCATE (ICLVS(MAX(1,NCLV)))
        ALLOCATE (ICLVT(MAX(1,NCLV)))
        ALLOCATE (ICLRC(MAX(1,NCLV)))
        ALLOCATE (ISNVE(MAX(1,NSNV)))
        ALLOCATE (ISNVS(MAX(1,NSNV)))
        ALLOCATE (ISNVT(MAX(1,NSNV)))
        ALLOCATE (ISNRC(MAX(1,NSNV)))
        ALLOCATE (IADSE(MAX(1,NADS)))
        ALLOCATE (IADSS(MAX(1,NADS)))
        ALLOCATE (IADST(MAX(1,NADS)))
        ALLOCATE (IADSC(MAX(1,NADS)))
        ALLOCATE (NFRSTP(NTALI))
        ALLOCATE (NFSTPI(NTALI))
        ALLOCATE (NADDP(NTALI))
        ALLOCATE (NADDCOR(NTALG))
        ALLOCATE (NSPAN(NTALV))
        ALLOCATE (NSPEN(NTALV))
        ALLOCATE (NSPANW(NTALS))
        ALLOCATE (NSPENW(NTALS))

        ALLOCATE (INTLOPTS(NTALI))

c  logicals
        ALLOCATE (LGVAC(NRAD,0:NPLS+1))
        ALLOCATE (LGDFT(NRAD))
        ALLOCATE (LSPCCLL(NRAD))
        ALLOCATE (LIVTALI(NTALI))
 
        WRITE (55+IFOFF,'(A,T25,I15)')
     .        ' COMUSR(1) ', NPLPR2*8 + MUSR*4 + LUSR*4 
 
      ELSE IF (ICAL == 2) THEN
 
        IF (ALLOCATED(PLSTLS)) RETURN

c
!        NPLPR1=(12+1*NPLS+NPLSTI+3*NPLSV)*NRAD  ! background data, set in plasma.f, 17 arrays
!        NPLPRM=NPLPR1+(NAIN+NSPZMC)*NRAD        !  adin, wght,...??? adin is allocated in call with ICAL == 2
cdr BVIN: add nplsv to nplpr2 and remove npls from nplprm. tbd:  check correct dimension of bvin !


        NPLPR1 = NINPTL*NRAD            ! storage only for active input tallies

        ALLOCATE (PLSTLS(NINPTL,NRAD))
 
        ALLOCATE (CEMETERYP(0:0,NRAD))  ! storage for in-active input tallies


        ALLOCATE (TEINL(NRAD))
        ALLOCATE (TIINL(NPLSTI,NRAD))
        ALLOCATE (DEINL(NRAD))
        ALLOCATE (DIINL(NPLS,NRAD))

c  NCPV, NBGV are now set
        ALLOCATE (ICPVE(NCPV))
        ALLOCATE (ICPVS(NCPV))
        ALLOCATE (ICPVT(NCPV))
        ALLOCATE (ICPRC(NCPV))
        ALLOCATE (IBGVE(NBGV))
        ALLOCATE (IBGVS(NBGV))
        ALLOCATE (IBGVT(NBGV))
        ALLOCATE (IBGRC(NBGV))

        WRITE (55+IFOFF,'(A,T25,I15)')
     .        ' COMUSR(2) ',(NPLPR1+(3+NPLSTI+NPLS)*NRAD)*8 + 
     .                      4*(NCPV+NBGV)*4 

      ELSE IF (ICAL == 3) THEN
 
        IF (ALLOCATED(FLXOUT)) RETURN
 
        NSFPRM=2*NLMPGS
        ALLOCATE (FLXOUT(NLMPGS))
        ALLOCATE (SAREA(NLMPGS))
 
        WRITE (55+IFOFF,'(A,T25,I15)')
     .         ' COMUSR(3) ',NSFPRM*8
 
      END IF
 
      CALL EIRENE_INIT_COMUSR(ICAL)
 
      RETURN
      END SUBROUTINE EIRENE_ALLOC_COMUSR


 
      SUBROUTINE EIRENE_ASSOCIATE_COMUSR
cdr special treatment of Ti:  intlopts.....

      IF (LTEIN) THEN
        TEIN => PLSTLS(NADDP(1)+1,:)
      ELSE 
        TEIN => CEMETERYP(0,:)
      END IF
      IF (LTIIN) THEN
        IF (INTLOPTS(2) >= 0) THEN
          TIIN => PLSTLS(NADDP(2)+1:NADDP(3),:)
        ELSE 
          TIIN => PLSTLS(NADDP(1)+1:NADDP(1)+1,:)    ! Ti = Te, no own storage for Ti
        END IF
      END IF
      IF (LDEIN) THEN
        DEIN => PLSTLS(NADDP(3)+1,:)
      ELSE 
        DEIN => CEMETERYP(0,:)
      END IF
      IF (LDIIN) THEN
        DIIN => PLSTLS(NADDP(4)+1:NADDP(5),:)
      ELSE 
        DIIN => CEMETERYP(0:0,:)
      END IF
      IF (LVXIN) THEN
        VXIN => PLSTLS(NADDP(5)+1:NADDP(6),:)
      ELSE 
        VXIN => CEMETERYP(0:0,:)
      END IF
      IF (LVYIN) THEN
        VYIN => PLSTLS(NADDP(6)+1:NADDP(7),:)
      ELSE 
        VYIN => CEMETERYP(0:0,:)
      END IF
      IF (LVZIN) THEN
        VZIN => PLSTLS(NADDP(7)+1:NADDP(8),:)
      ELSE 
        VZIN => CEMETERYP(0:0,:)
      END IF
      IF (LBXIN) THEN
        BXIN => PLSTLS(NADDP(8)+1,:)
      ELSE 
        BXIN => CEMETERYP(0,:)
      END IF
      IF (LBYIN) THEN
        BYIN => PLSTLS(NADDP(9)+1,:)
      ELSE 
        BYIN => CEMETERYP(0,:)
      END IF
      IF (LBZIN) THEN
        BZIN => PLSTLS(NADDP(10)+1,:)
      ELSE 
        BZIN => CEMETERYP(0,:)
      END IF
      IF (LBFIN) THEN
        BFIN => PLSTLS(NADDP(11)+1,:)
      ELSE 
        BFIN => CEMETERYP(0,:)
      END IF
      IF (LADIN) THEN
        ADIN => PLSTLS(NADDP(12)+1:NADDP(13),:)
      ELSE 
        ADIN => CEMETERYP(0:0,:)
      END IF
      IF (LEDRIFT) THEN
        EDRIFT => PLSTLS(NADDP(13)+1:NADDP(14),:)
      ELSE 
        EDRIFT => CEMETERYP(0:0,:)
      END IF
      IF (LVOL) THEN
        VOL => PLSTLS(NADDP(14)+1,:)
      ELSE 
        VOL => CEMETERYP(0,:)
      END IF
      IF (LWGHT) THEN
        WGHT => PLSTLS(NADDP(15)+1:NADDP(16),:)
      ELSE 
        WGHT => CEMETERYP(0:0,:)
      END IF
      IF (LBXPERP) THEN
        BXPERP => PLSTLS(NADDP(16)+1,:)
      ELSE 
        BXPERP => CEMETERYP(0,:)
      END IF
       IF (LBYPERP) THEN
        BYPERP => PLSTLS(NADDP(17)+1,:)
      ELSE 
        BYPERP => CEMETERYP(0,:)
      END IF
      IF (LEXIN) THEN
        EXIN => PLSTLS(NADDP(18)+1,:)
      ELSE 
        EXIN => CEMETERYP(0,:)
      END IF
      IF (LEYIN) THEN
        EYIN => PLSTLS(NADDP(19)+1,:)
      ELSE 
        EYIN => CEMETERYP(0,:)
      END IF
      IF (LEZIN) THEN
        EZIN => PLSTLS(NADDP(20)+1,:)
      ELSE 
        EZIN => CEMETERYP(0,:)
      END IF
      IF (LEFIN) THEN
        EFIN => PLSTLS(NADDP(21)+1,:)
      ELSE 
        EFIN => CEMETERYP(0,:)
      END IF
      IF (LPOT) THEN
        POT => PLSTLS(NADDP(22)+1,:)
      ELSE 
        POT => CEMETERYP(0,:)
      END IF
      IF (LBVIN) THEN
        BVIN => PLSTLS(NADDP(23)+1:NADDP(24),:)
      ELSE 
        BVIN => CEMETERYP(0:0,:)
      END IF
      IF (LPARMOM) THEN
        PARMOM => PLSTLS(NADDP(24)+1:NADDP(25),:)
      ELSE 
        PARMOM => CEMETERYP(0:0,:)
      END IF
     
      IF (LDTEDX) THEN
        DTEDX => PLSTLS(NADDP(25)+1,:)
      ELSE 
        DTEDX => CEMETERYP(0,:)
      END IF
      IF (LDTEDY) THEN
        DTEDY => PLSTLS(NADDP(26)+1,:)
      ELSE 
        DTEDY => CEMETERYP(0,:)
      END IF
      IF (LDTEDZ) THEN
        DTEDZ => PLSTLS(NADDP(27)+1,:)
      ELSE 
        DTEDZ => CEMETERYP(0,:)
      END IF
      IF (LDTIDX) THEN
        DTIDX => PLSTLS(NADDP(28)+1:NADDP(29),:)
      ELSE 
        DTIDX => CEMETERYP(0:0,:)
      END IF
      IF (LDTIDY) THEN
        DTIDY => PLSTLS(NADDP(29)+1:NADDP(30),:)
      ELSE 
        DTIDY => CEMETERYP(0:0,:)
      END IF
      IF (LDTIDZ) THEN
        DTIDZ => PLSTLS(NADDP(30)+1:NADDP(31),:)
      ELSE 
        DTIDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDDEDX) THEN
        DDEDX => PLSTLS(NADDP(31)+1,:)
      ELSE 
        DDEDX => CEMETERYP(0,:)
      END IF
      IF (LDDEDY) THEN
        DDEDY => PLSTLS(NADDP(32)+1,:)
      ELSE 
        DDEDY => CEMETERYP(0,:)
      END IF
      IF (LDDEDZ) THEN
        DDEDZ => PLSTLS(NADDP(33)+1,:)
      ELSE 
        DDEDZ => CEMETERYP(0,:)
      END IF
      IF (LDDIDX) THEN
        DDIDX => PLSTLS(NADDP(34)+1:NADDP(35),:)
      ELSE 
        DDIDX => CEMETERYP(0:0,:)
      END IF
      IF (LDDIDY) THEN
        DDIDY => PLSTLS(NADDP(35)+1:NADDP(36),:)
      ELSE 
        DDIDY => CEMETERYP(0:0,:)
      END IF
      IF (LDDIDZ) THEN
        DDIDZ => PLSTLS(NADDP(36)+1:NADDP(37),:)
      ELSE 
        DDIDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDVXDX) THEN
        DVXDX => PLSTLS(NADDP(37)+1:NADDP(38),:)
      ELSE 
        DVXDX => CEMETERYP(0:0,:)
      END IF
      IF (LDVXDY) THEN
        DVXDY => PLSTLS(NADDP(38)+1:NADDP(39),:)
      ELSE 
        DVXDY => CEMETERYP(0:0,:)
      END IF
      IF (LDVXDZ) THEN
        DVXDZ => PLSTLS(NADDP(39)+1:NADDP(40),:)
      ELSE 
        DVXDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDVYDX) THEN
        DVYDX => PLSTLS(NADDP(40)+1:NADDP(41),:)
      ELSE 
        DVYDX => CEMETERYP(0:0,:)
      END IF
      IF (LDVYDY) THEN
        DVYDY => PLSTLS(NADDP(41)+1:NADDP(42),:)
      ELSE 
        DVYDY => CEMETERYP(0:0,:)
      END IF
      IF (LDVYDZ) THEN
        DVYDZ => PLSTLS(NADDP(42)+1:NADDP(43),:)
      ELSE 
        DVYDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDVZDX) THEN
        DVZDX => PLSTLS(NADDP(43)+1:NADDP(44),:)
      ELSE 
        DVZDX => CEMETERYP(0:0,:)
      END IF
      IF (LDVZDY) THEN
        DVZDY => PLSTLS(NADDP(44)+1:NADDP(45),:)
      ELSE 
        DVZDY => CEMETERYP(0:0,:)
      END IF
      IF (LDVZDZ) THEN
        DVZDZ => PLSTLS(NADDP(45)+1:NADDP(46),:)
      ELSE 
        DVZDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDBXDX) THEN
        DBXDX => PLSTLS(NADDP(46)+1,:)
      ELSE 
        DBXDX => CEMETERYP(0,:)
      END IF
      IF (LDBXDY) THEN
        DBXDY => PLSTLS(NADDP(47)+1,:)
      ELSE 
        DBXDY => CEMETERYP(0,:)
      END IF
      IF (LDBXDZ) THEN
        DBXDZ => PLSTLS(NADDP(48)+1,:)
      ELSE 
        DBXDZ => CEMETERYP(0,:)
      END IF
      IF (LDBYDX) THEN
        DBYDX => PLSTLS(NADDP(49)+1,:)
      ELSE 
        DBYDX => CEMETERYP(0,:)
      END IF
      IF (LDBYDY) THEN
        DBYDY => PLSTLS(NADDP(50)+1,:)
      ELSE 
        DBYDY => CEMETERYP(0,:)
      END IF
      IF (LDBYDZ) THEN
        DBYDZ => PLSTLS(NADDP(51)+1,:)
      ELSE 
        DBYDZ => CEMETERYP(0,:)
      END IF
      IF (LDBZDX) THEN
        DBZDX => PLSTLS(NADDP(52)+1,:)
      ELSE 
        DBZDX => CEMETERYP(0,:)
      END IF
      IF (LDBZDY) THEN
        DBZDY => PLSTLS(NADDP(53)+1,:)
      ELSE 
        DBZDY => CEMETERYP(0,:)
      END IF
      IF (LDBZDZ) THEN
        DBZDZ => PLSTLS(NADDP(54)+1,:)
      ELSE 
        DBZDZ => CEMETERYP(0,:)
      END IF
      IF (LDBFDX) THEN
        DBFDX => PLSTLS(NADDP(55)+1,:)
      ELSE 
        DBFDX => CEMETERYP(0,:)
      END IF
      IF (LDBFDY) THEN
        DBFDY => PLSTLS(NADDP(56)+1,:)
      ELSE 
        DBFDY => CEMETERYP(0,:)
      END IF
      IF (LDBFDZ) THEN
        DBFDZ => PLSTLS(NADDP(57)+1,:)
      ELSE 
        DBFDZ => CEMETERYP(0,:)
      END IF
      IF (LDADINDX) THEN
        DADINDX => PLSTLS(NADDP(58)+1:NADDP(59),:)
      ELSE 
        DADINDX => CEMETERYP(0:0,:)
      END IF
      IF (LDADINDY) THEN
        DADINDY => PLSTLS(NADDP(59)+1:NADDP(60),:)
      ELSE 
        DADINDY => CEMETERYP(0:0,:)
      END IF
      IF (LDADINDZ) THEN
        DADINDZ => PLSTLS(NADDP(60)+1:NADDP(61),:)
      ELSE 
        DADINDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDEDRIFTDX) THEN
        DEDRIFTDX => PLSTLS(NADDP(61)+1:NADDP(62),:)
      ELSE 
        DEDRIFTDX => CEMETERYP(0:0,:)
      END IF
      IF (LDEDRIFTDY) THEN
        DEDRIFTDY => PLSTLS(NADDP(62)+1:NADDP(63),:)
      ELSE 
        DEDRIFTDY => CEMETERYP(0:0,:)
      END IF
      IF (LDEDRIFTDZ) THEN
        DEDRIFTDZ => PLSTLS(NADDP(63)+1:NADDP(64),:)
      ELSE 
        DEDRIFTDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDVOLDX) THEN
        DVOLDX => PLSTLS(NADDP(64)+1,:)
      ELSE 
        DVOLDX => CEMETERYP(0,:)
      END IF
      IF (LDVOLDY) THEN
        DVOLDY => PLSTLS(NADDP(65)+1,:)
      ELSE 
        DVOLDY => CEMETERYP(0,:)
      END IF
      IF (LDVOLDZ) THEN
        DVOLDZ => PLSTLS(NADDP(66)+1,:)
      ELSE 
        DVOLDZ => CEMETERYP(0,:)
      END IF
      IF (LDWGHTDX) THEN
        DWGHTDX => PLSTLS(NADDP(67)+1:NADDP(68),:)
      ELSE 
        DWGHTDX => CEMETERYP(0:0,:)
      END IF
      IF (LDWGHTDY) THEN
        DWGHTDY => PLSTLS(NADDP(68)+1:NADDP(69),:)
      ELSE 
        DWGHTDY => CEMETERYP(0:0,:)
      END IF
      IF (LDWGHTDZ) THEN
        DWGHTDZ => PLSTLS(NADDP(69)+1:NADDP(70),:)
      ELSE 
        DWGHTDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDBXPERPDX) THEN
        DBXPERPDX => PLSTLS(NADDP(70)+1,:)
      ELSE 
        DBXPERPDX => CEMETERYP(0,:)
      END IF
      IF (LDBXPERPDY) THEN
        DBXPERPDY => PLSTLS(NADDP(71)+1,:)
      ELSE 
        DBXPERPDY => CEMETERYP(0,:)
      END IF
      IF (LDBXPERPDZ) THEN
        DBXPERPDZ => PLSTLS(NADDP(72)+1,:)
      ELSE 
        DBXPERPDZ => CEMETERYP(0,:)
      END IF
      IF (LDBYPERPDX) THEN
        DBYPERPDX => PLSTLS(NADDP(73)+1,:)
      ELSE 
        DBYPERPDX => CEMETERYP(0,:)
      END IF
      IF (LDBYPERPDY) THEN
        DBYPERPDY => PLSTLS(NADDP(74)+1,:)
      ELSE 
        DBYPERPDY => CEMETERYP(0,:)
      END IF
      IF (LDBYPERPDZ) THEN
        DBYPERPDZ => PLSTLS(NADDP(75)+1,:)
      ELSE 
        DBYPERPDZ => CEMETERYP(0,:)
      END IF
      IF (LDEXDX) THEN
        DEXDX => PLSTLS(NADDP(76)+1,:)
      ELSE 
        DEXDX => CEMETERYP(0,:)
      END IF
      IF (LDEXDY) THEN
        DEXDY => PLSTLS(NADDP(77)+1,:)
      ELSE 
        DEXDY => CEMETERYP(0,:)
      END IF
      IF (LDEXDZ) THEN
        DEXDZ => PLSTLS(NADDP(78)+1,:)
      ELSE 
        DEXDZ => CEMETERYP(0,:)
      END IF
      IF (LDEYDX) THEN
        DEYDX => PLSTLS(NADDP(79)+1,:)
      ELSE 
        DEYDX => CEMETERYP(0,:)
      END IF
      IF (LDEYDY) THEN
        DEYDY => PLSTLS(NADDP(80)+1,:)
      ELSE 
        DEYDY => CEMETERYP(0,:)
      END IF
      IF (LDEYDZ) THEN
        DEYDZ => PLSTLS(NADDP(81)+1,:)
      ELSE 
        DEYDZ => CEMETERYP(0,:)
      END IF
      IF (LDEZDX) THEN
        DEZDX => PLSTLS(NADDP(82)+1,:)
      ELSE 
        DEZDX => CEMETERYP(0,:)
      END IF
      IF (LDEZDY) THEN
        DEZDY => PLSTLS(NADDP(83)+1,:)
      ELSE 
        DEZDY => CEMETERYP(0,:)
      END IF
      IF (LDEZDZ) THEN
        DEZDZ => PLSTLS(NADDP(84)+1,:)
      ELSE 
        DEZDZ => CEMETERYP(0,:)
      END IF
      IF (LDEFDX) THEN
        DEFDX => PLSTLS(NADDP(85)+1,:)
      ELSE 
        DEFDX => CEMETERYP(0,:)
      END IF
      IF (LDEFDY) THEN
        DEFDY => PLSTLS(NADDP(86)+1,:)
      ELSE 
        DEFDY => CEMETERYP(0,:)
      END IF
      IF (LDEFDZ) THEN
        DEFDZ => PLSTLS(NADDP(87)+1,:)
      ELSE 
        DEFDZ => CEMETERYP(0,:)
      END IF
      IF (LDPOTDX) THEN
        DPOTDX => PLSTLS(NADDP(88)+1,:)
      ELSE 
        DPOTDX => CEMETERYP(0,:)
      END IF
      IF (LDPOTDY) THEN
        DPOTDY => PLSTLS(NADDP(89)+1,:)
      ELSE 
        DPOTDY => CEMETERYP(0,:)
      END IF
      IF (LDPOTDZ) THEN
        DPOTDZ => PLSTLS(NADDP(90)+1,:)
      ELSE 
        DPOTDZ => CEMETERYP(0,:)
      END IF
      IF (LDBVINDX) THEN
        DBVINDX => PLSTLS(NADDP(91)+1:NADDP(92),:)
      ELSE 
        DBVINDX => CEMETERYP(0:0,:)
      END IF
      IF (LDBVINDY) THEN
        DBVINDY => PLSTLS(NADDP(92)+1:NADDP(93),:)
      ELSE 
        DBVINDY => CEMETERYP(0:0,:)
      END IF
      IF (LDBVINDZ) THEN
        DBVINDZ => PLSTLS(NADDP(93)+1:NADDP(94),:)
      ELSE 
        DBVINDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDPARMOMDX) THEN
        DPARMOMDX => PLSTLS(NADDP(94)+1:NADDP(95),:)
      ELSE 
        DPARMOMDX => CEMETERYP(0:0,:)
      END IF
      IF (LDPARMOMDY) THEN
        DPARMOMDY => PLSTLS(NADDP(95)+1:NADDP(96),:)
      ELSE 
        DPARMOMDY => CEMETERYP(0:0,:)
      END IF
      IF (LDPARMOMDZ) THEN
        DPARMOMDZ => PLSTLS(NADDP(96)+1:NINPTL,:)
      ELSE 
        DPARMOMDZ => CEMETERYP(0:0,:)
      END IF

      RETURN
      END SUBROUTINE EIRENE_ASSOCIATE_COMUSR


 
      SUBROUTINE EIRENE_ALLOC_CORNERS(IUNOUT)

      INTEGER, INTENT(IN) :: IUNOUT
      INTEGER :: N1DIM(12)
      INTEGER :: NTOT, ICO, NLST, I, J, NLSTTL, NTOT2
      

      IF (ALLOCATED(CORNER_PROFILES)) RETURN

      N1DIM = 0
cdr  are there any FEM interpolated background tallies in this run?

!  interpolation to vertices can only be done if input tally is available (active)
      LSMOPRO(1:NTALG) = LSMOPRO(1:NTALG) .AND. LIVTALI(1:NTALG)

c  NTOT2: total number of smoothed talles, counting also with species index 
      NTOT2 = 0
      DO I= 1, NTALG
        IF (LSMOPRO(I)) THEN
          NTOT2 = NTOT2 + NFRSTP(I)
        END IF
      END DO     

c  NADDCOR: cummulated index of position of smoothed tally J within all smoothed tallies
C  NLSTLL : highest tally index J amongst all smoothed tallies
      NADDCOR(1)=0
      DO 6 J=2,NTALG
        IF (LSMOPRO(J-1)) THEN
          NADDCOR(J)=NADDCOR(J-1)+NFRSTP(J-1)
          NLSTTL=J-1
        ELSE
          NADDCOR(J)=NADDCOR(J-1)
        END IF
 6    CONTINUE

      IF (LSMOPRO(NTALG)) NLSTTL = NTALG
C
c  NTOT: total number of smoothed talles, counting also with species index 
      NTOT = 0
      IF (ANY(LSMOPRO)) THEN
        NTOT = NADDCOR(NTALG)
        IF (NFRSTP(NLSTTL) > 1) NTOT = NTOT+NFRSTP(NLSTTL)
      END IF
      
      IF (NTOT > 0) THEN
cdr  allocate storage for background tallies on cell vertices
cdr  ncorner is set in GRID.f (levgeo=4,5) or in SNEIGH.f (levgeo=1,2,3)
        ALLOCATE (CORNER_PROFILES(NCORNER,NTOT))
      ELSE
        ALLOCATE (CORNER_PROFILES(1,1))
      END IF

       WRITE (55+IFOFF,'(A,T25,I15)')
     .        ' COMUSR(CORNERS) ',SIZE(CORNER_PROFILES)*8

      LDSMO = LDESMO .OR. LDISMO 
      LVSMO = LVXSMO .OR. LVYSMO .OR. LVZSMO .OR. LBVSMO
      LBSMO = LBXSMO .OR. LBYSMO .OR. LBZSMO .OR. LBFSMO
      LESMO = LEXSMO .OR. LEYSMO .OR. LEZSMO .OR. LEFSMO .OR.
     .        LPOTSMO

      IF (LTESMO) THEN
        TEINCORNER => CORNER_PROFILES(:,NADDCOR(1)+1)
      ELSE
        NULLIFY(TEINCORNER)
      END IF

      IF (LTISMO) THEN
        TIINCORNER => CORNER_PROFILES(:,NADDCOR(2)+1 : NADDCOR(3))
      ELSE
        NULLIFY(TIINCORNER)
      END IF

      IF (LDSMO) THEN
        IF (LDESMO) THEN
          DEINCORNER => CORNER_PROFILES(:,NADDCOR(3)+1)
        ELSE
          NULLIFY(DEINCORNER)
        END IF
        IF (LDISMO) THEN
          DIINCORNER => CORNER_PROFILES(:,NADDCOR(4)+1 : NADDCOR(5))
        ELSE
          NULLIFY(DIINCORNER)
        END IF
      ELSE
        NULLIFY(DEINCORNER)
        NULLIFY(DIINCORNER)
      END IF

      IF (LVSMO) THEN
        IF (LVXSMO) THEN
          VXINCORNER => CORNER_PROFILES(:,NADDCOR(5)+1 : NADDCOR(6))
        ELSE
          NULLIFY(VXINCORNER)
        END IF 
        IF (LVYSMO) THEN
          VYINCORNER => CORNER_PROFILES(:,NADDCOR(6)+1 : NADDCOR(7))
        ELSE
          NULLIFY(VYINCORNER)
        END IF 
        IF (LVZSMO) THEN
          VZINCORNER => CORNER_PROFILES(:,NADDCOR(7)+1 : NADDCOR(8))
        ELSE
          NULLIFY(VZINCORNER)
        END IF
        IF (LBVSMO) THEN
          BVINCORNER => CORNER_PROFILES(:,NADDCOR(22)+1 : NADDCOR(24))
        ELSE
          NULLIFY(BVINCORNER)
        END IF
      ELSE
        NULLIFY(VXINCORNER)
        NULLIFY(VYINCORNER)
        NULLIFY(VZINCORNER)
        NULLIFY(BVINCORNER)
      END IF

      IF (LBSMO) THEN
        IF (LBXSMO) THEN
          BXINCORNER => CORNER_PROFILES(:,NADDCOR(8)+1)
        ELSE
          NULLIFY(BXINCORNER)
        END IF
        IF (LBYSMO) THEN
          BYINCORNER => CORNER_PROFILES(:,NADDCOR(9)+1)
        ELSE
          NULLIFY(BYINCORNER)
        END IF
        IF (LBZSMO) THEN
          BZINCORNER => CORNER_PROFILES(:,NADDCOR(10)+1)
        ELSE
          NULLIFY(BZINCORNER)
        END IF
        IF (LBFSMO) THEN
          BFINCORNER => CORNER_PROFILES(:,NADDCOR(11)+1)
        ELSE
          NULLIFY(BFINCORNER)
        END IF
      ELSE
        NULLIFY(BXINCORNER)
        NULLIFY(BYINCORNER)
        NULLIFY(BZINCORNER)
        NULLIFY(BFINCORNER)
      END IF

      IF (LADSMO) THEN
        ADCORNER => CORNER_PROFILES(:,NADDCOR(12)+1 : NADDCOR(13))
      ELSE
        NULLIFY(ADCORNER)
      END IF

      IF (LEDRIFTSMO) THEN
        EDRIFTCORNER => CORNER_PROFILES(:,NADDCOR(13)+1 : NADDCOR(14))
      ELSE
        NULLIFY(EDRIFTCORNER)
      END IF

      IF (LVOLSMO) THEN
        VOLCORNER => CORNER_PROFILES(:,NADDCOR(14)+1 )
      ELSE
        NULLIFY(VOLCORNER)
      END IF

      IF (LWGHTSMO) THEN
        WGHTCORNER => CORNER_PROFILES(:,NADDCOR(15)+1 : NADDCOR(16))
      ELSE
        NULLIFY(WGHTCORNER)
      END IF

cdr these next two B field tallies should go into LBSMO
      IF (LBXPSMO) THEN
        BXPERPCORNER => CORNER_PROFILES(:,NADDCOR(16)+1 )
      ELSE
        NULLIFY(BXPERPCORNER)
      END IF

      IF (LBYPSMO) THEN
        BYPERPCORNER => CORNER_PROFILES(:,NADDCOR(17)+1 )
      ELSE
        NULLIFY(BYPERPCORNER)
      END IF

      IF (LESMO) THEN
        IF (LEXSMO) THEN
          EXCORNER => CORNER_PROFILES(:,NADDCOR(18)+1)
        ELSE
          NULLIFY(EXCORNER)
        END IF
        IF (LEYSMO) THEN
          EYCORNER => CORNER_PROFILES(:,NADDCOR(19)+1)
        ELSE
          NULLIFY(EYCORNER)
        END IF
        IF (LEZSMO) THEN
          EZCORNER => CORNER_PROFILES(:,NADDCOR(20)+1)
        ELSE
          NULLIFY(EZCORNER)
        END IF
        IF (LEFSMO) THEN
          EFCORNER => CORNER_PROFILES(:,NADDCOR(21)+1)
        ELSE
          NULLIFY(EFCORNER)
        END IF
        IF (LPOTSMO) THEN
          POTCORNER => CORNER_PROFILES(:,NADDCOR(22)+1)
        ELSE
          NULLIFY(POTCORNER)
        END IF
      ELSE
        NULLIFY(EXCORNER)
        NULLIFY(EYCORNER)
        NULLIFY(EZCORNER)
        NULLIFY(EFCORNER)
        NULLIFY(POTCORNER)
      END IF

      IF (LPARMOMSMO) THEN
        PARMOMCORNER => CORNER_PROFILES(:,NADDCOR(24+1) : NTOT )
      ELSE
        NULLIFY(PARMOMCORNER)
      END IF

      CORNER_PROFILES = 0._DP

      END SUBROUTINE EIRENE_ALLOC_CORNERS
 


      SUBROUTINE EIRENE_DEALLOC_COMUSR
C
      IF (.NOT.ALLOCATED(PLSTLS)) RETURN
 
      DEALLOCATE (PLSTLS)
      DEALLOCATE (CEMETERYP)
c
      DEALLOCATE (TEINL)
      DEALLOCATE (TIINL)
      DEALLOCATE (DEINL)
      DEALLOCATE (DIINL)

      DEALLOCATE (FLXOUT)
      DEALLOCATE (SAREA)

      
      DEALLOCATE (RMASSA)
      DEALLOCATE (RMASSM)
      DEALLOCATE (RMASSI)
      DEALLOCATE (RMASSPH)
      DEALLOCATE (RMASSP)

      DEALLOCATE (DIOD)
      DEALLOCATE (DATD)
      DEALLOCATE (DMLD)
      DEALLOCATE (DPLD)
      DEALLOCATE (DPHD)
      DEALLOCATE (DION)
      DEALLOCATE (DATM)
      DEALLOCATE (DMOL)
      DEALLOCATE (DPLS)
      DEALLOCATE (DPHOT)
 
      DEALLOCATE (TEXTS)
      DEALLOCATE (NMASSA)
      DEALLOCATE (NCHARA)
      DEALLOCATE (NFOLA)
      DEALLOCATE (NGENA)
      DEALLOCATE (NMASSM)
      DEALLOCATE (NCHARM)
      DEALLOCATE (NFOLM)
      DEALLOCATE (NGENM)
      DEALLOCATE (NMASSP)
      DEALLOCATE (NCHARP)
      DEALLOCATE (NCHRGP)
      DEALLOCATE (NMASSI)
      DEALLOCATE (NCHARI)
      DEALLOCATE (NCHRGI)
      DEALLOCATE (NFOLI)
      DEALLOCATE (NGENI)
      DEALLOCATE (NFOLPH)
      DEALLOCATE (NGENPH)
      DEALLOCATE (NPRT)
      DEALLOCATE (ISPEZ)
      DEALLOCATE (ISPEZI)
      DEALLOCATE (MPLSTI)
      DEALLOCATE (MPLSV)
      DEALLOCATE (ISPZ_BACK)
      DEALLOCATE (IADVE)
      DEALLOCATE (IADVS)
      DEALLOCATE (IADVT)
      DEALLOCATE (IADRC)
      DEALLOCATE (ICLVE)
      DEALLOCATE (ICLVS)
      DEALLOCATE (ICLVT)
      DEALLOCATE (ICLRC)
      DEALLOCATE (ISNVE)
      DEALLOCATE (ISNVS)
      DEALLOCATE (ISNVT)
      DEALLOCATE (ISNRC)
      DEALLOCATE (ICPVE)
      DEALLOCATE (ICPVS)
      DEALLOCATE (ICPVT)
      DEALLOCATE (ICPRC)
      DEALLOCATE (IBGVE)
      DEALLOCATE (IBGVS)
      DEALLOCATE (IBGVT)
      DEALLOCATE (IBGRC)
      DEALLOCATE (IADSE)
      DEALLOCATE (IADSS)
      DEALLOCATE (IADST)
      DEALLOCATE (IADSC)
      DEALLOCATE (NFRSTP)
      DEALLOCATE (NFSTPI)
      DEALLOCATE (NADDP)
      DEALLOCATE (NSPAN)
      DEALLOCATE (NSPEN)
      DEALLOCATE (NSPANW)
      DEALLOCATE (NSPENW)
      DEALLOCATE (INTLOPTS)
      DEALLOCATE (LGVAC)
      DEALLOCATE (LGDFT)
      DEALLOCATE (LSPCCLL)
      DEALLOCATE (LSMOPRO)
      DEALLOCATE (LIVTALI)
 
!pb      IF (NBACK_SPEC > 0) DEALLOCATE (BACK_SPEC)
      IF (ALLOCATED(BACK_SPEC)) DEALLOCATE (BACK_SPEC)

      IF (ALLOCATED(CORNER_PROFILES)) DEALLOCATE (CORNER_PROFILES)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_COMUSR
 
 
      SUBROUTINE EIRENE_INIT_COMUSR(ICAL)
 
      INTEGER, INTENT(IN) :: ICAL
      INTEGER, SAVE :: IFIRST=0

      IF (IFIRST == 0) THEN
        LSMOPRO = .FALSE.

        LTESMO      => LSMOPRO(1)
        LTISMO      => LSMOPRO(2)
        LDESMO      => LSMOPRO(3)
        LDISMO      => LSMOPRO(4)
c  plasma flow field
        LVXSMO      => LSMOPRO(5)
        LVYSMO      => LSMOPRO(6)
        LVZSMO      => LSMOPRO(7)
c  plasma flow parallel
        LBVSMO      => LSMOPRO(23)
        LPARMOMSMO  => LSMOPRO(24)

c  B field
        LBXSMO      => LSMOPRO(8)
        LBYSMO      => LSMOPRO(9)
        LBZSMO      => LSMOPRO(10)
        LBFSMO      => LSMOPRO(11)

c  B  perp
        LBXPSMO     => LSMOPRO(16)
        LBYPSMO     => LSMOPRO(17)

        LADSMO      => LSMOPRO(12)
        LEDRIFTSMO  => LSMOPRO(13)
        LVOLSMO     => LSMOPRO(14)
        LWGHTSMO    => LSMOPRO(15)

c  E field
        LEXSMO      => LSMOPRO(18)
        LEYSMO      => LSMOPRO(19)
        LEZSMO      => LSMOPRO(20)
        LEFSMO      => LSMOPRO(21)
        LPOTSMO     => LSMOPRO(22)


        IFIRST = 1
      ENDIF
 
      IF (ICAL == 1) THEN
cdr oct 18: initialization of input volumetric tallies moved to ICAL==2
        
        RMASSA = 0._DP
        RMASSM = 0._DP
        RMASSI = 0._DP
        RMASSPH = 0._DP
        RMASSP = 0._DP

        DIOD   = 0._DP
        DATD   = 0._DP
        DMLD   = 0._DP
        DPLD   = 0._DP
        DPHD   = 0._DP
        DION   = 0._DP
        DATM   = 0._DP
        DMOL   = 0._DP
        DPLS   = 0._DP
        DPHOT  = 0._DP
 
        TEXTS  = ' '
        NMASSA = 0
        NCHARA = 0
        NFOLA  = 0
        NGENA  = 0
        NMASSM = 0
        NCHARM = 0
        NFOLM  = 0
        NGENM  = 0
        NMASSP = 0
        NCHARP = 0
        NCHRGP = 0
        NMASSI = 0
        NCHARI = 0
        NCHRGI = 0
        NFOLI  = 0
        NGENI  = 0
        NFOLPH = 0
        NGENPH = 0
        NPRT   = 0
        ISPEZ  = 0
        ISPEZI = 0
        MPLSTI = 0
        MPLSV  = 0
        ISPZ_BACK = 0
        IADVE  = 0
        IADVS  = 0
        IADVT  = 0
        IADRC  = 0
        ICLVE  = 0
        ICLVS  = 0
        ICLVT  = 0
        ICLRC  = 0
        ISNVE  = 0
        ISNVS  = 0
        ISNVT  = 0
        ISNRC  = 0
        IADSE  = 0
        IADSS  = 0
        IADST  = 0
        IADSC  = 0
        NFRSTP = 0
        NFSTPI = 0
        NADDP  = 0
        NSPAN  = 0
        NSPEN  = 0
        NSPANW = 0
        NSPENW = 0
        INTLOPTS  = 0
        LGVAC  = .FALSE.
        LGDFT  = .FALSE.
        LSPCCLL = .FALSE.
        LIVTALI = .TRUE.

        
        LTEIN      => LIVTALI(1)
        LTIIN      => LIVTALI(2)
        LDEIN      => LIVTALI(3)
        LDIIN      => LIVTALI(4)
        LVXIN      => LIVTALI(5)
        LVYIN      => LIVTALI(6)
        LVZIN      => LIVTALI(7)
        LBXIN      => LIVTALI(8)
        LBYIN      => LIVTALI(9)
        LBZIN      => LIVTALI(10)
        LBFIN      => LIVTALI(11)
        LADIN      => LIVTALI(12)
        LEDRIFT    => LIVTALI(13)
        LVOL       => LIVTALI(14)
        LWGHT      => LIVTALI(15)
        LBXPERP    => LIVTALI(16)
        LBYPERP    => LIVTALI(17)
        LEXIN      => LIVTALI(18)
        LEYIN      => LIVTALI(19)
        LEZIN      => LIVTALI(20)
        LEFIN      => LIVTALI(21)
        LPOT       => LIVTALI(22)
        LBVIN      => LIVTALI(23)
        LPARMOM    => LIVTALI(24)

        LDTEDX     => LIVTALI(25)
        LDTEDY     => LIVTALI(26)
        LDTEDZ     => LIVTALI(27)
        LDTIDX     => LIVTALI(28)
        LDTIDY     => LIVTALI(29)
        LDTIDZ     => LIVTALI(30)
        LDDEDX     => LIVTALI(31)
        LDDEDY     => LIVTALI(32)
        LDDEDZ     => LIVTALI(33)
        LDDIDX     => LIVTALI(34)
        LDDIDY     => LIVTALI(35)
        LDDIDZ     => LIVTALI(36)
        LDVXDX     => LIVTALI(37)
        LDVXDY     => LIVTALI(38)
        LDVXDZ     => LIVTALI(39)
        LDVYDX     => LIVTALI(40)
        LDVYDY     => LIVTALI(41)
        LDVYDZ     => LIVTALI(42)
        LDVZDX     => LIVTALI(43)
        LDVZDY     => LIVTALI(44)
        LDVZDZ     => LIVTALI(45)
        LDBXDX     => LIVTALI(46)
        LDBXDY     => LIVTALI(47)
        LDBXDZ     => LIVTALI(48)
        LDBYDX     => LIVTALI(49)
        LDBYDY     => LIVTALI(50)
        LDBYDZ     => LIVTALI(51)
        LDBZDX     => LIVTALI(52)
        LDBZDY     => LIVTALI(53)
        LDBZDZ     => LIVTALI(54)
        LDBFDX     => LIVTALI(55)
        LDBFDY     => LIVTALI(56)
        LDBFDZ     => LIVTALI(57)
        LDADINDX   => LIVTALI(58)
        LDADINDY   => LIVTALI(59)
        LDADINDZ   => LIVTALI(60)
        LDEDRIFTDX => LIVTALI(61)
        LDEDRIFTDY => LIVTALI(62)
        LDEDRIFTDZ => LIVTALI(63)
        LDVOLDX    => LIVTALI(64)
        LDVOLDY    => LIVTALI(65)
        LDVOLDZ    => LIVTALI(66)
        LDWGHTDX   => LIVTALI(67)
        LDWGHTDY   => LIVTALI(68)
        LDWGHTDZ   => LIVTALI(69)
        LDBXPERPDX => LIVTALI(70)
        LDBXPERPDY => LIVTALI(71)
        LDBXPERPDZ => LIVTALI(72)
        LDBYPERPDX => LIVTALI(73)
        LDBYPERPDY => LIVTALI(74)
        LDBYPERPDZ => LIVTALI(75)
        LDEXDX     => LIVTALI(76)
        LDEXDY     => LIVTALI(77)
        LDEXDZ     => LIVTALI(78)
        LDEYDX     => LIVTALI(79)
        LDEYDY     => LIVTALI(80)
        LDEYDZ     => LIVTALI(81)
        LDEZDX     => LIVTALI(82)
        LDEZDY     => LIVTALI(83)
        LDEZDZ     => LIVTALI(84)
        LDEFDX     => LIVTALI(85)
        LDEFDY     => LIVTALI(86)
        LDEFDZ     => LIVTALI(87)
        LDPOTDX    => LIVTALI(88)
        LDPOTDY    => LIVTALI(89)
        LDPOTDZ    => LIVTALI(90)
        LDBVINDX   => LIVTALI(91)
        LDBVINDY   => LIVTALI(92)
        LDBVINDZ   => LIVTALI(93)
        LDPARMOMDX => LIVTALI(94)
        LDPARMOMDY => LIVTALI(95)
        LDPARMOMDZ => LIVTALI(96)
 
      ELSE IF (ICAL == 2) THEN
c  Active volumetric input tallies
        PLSTLS = 0._DP
c  Cemetery for inactive input tallies (no storage)
        CEMETERYP = 0._DP

        TEINL  = 0._DP
        TIINL  = 0._DP
        DEINL  = 0._DP
        DIINL  = 0._DP

        ICPVE  = 0
        ICPVS  = 0
        ICPVT  = 0
        ICPRC  = 0
        IBGVE  = 0
        IBGVS  = 0
        IBGVT  = 0
        IBGRC  = 0
 
      ELSE IF (ICAL == 3) THEN
 
        FLXOUT = 0._DP
        SAREA  = 666._DP
 
      END IF
 
      RETURN

      ENTRY EIRENE_COMUSR_REINIT
      IFIRST = 0
      RETURN

      END SUBROUTINE EIRENE_INIT_COMUSR
 
      END MODULE EIRMOD_COMUSR
