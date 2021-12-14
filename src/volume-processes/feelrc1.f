!pb  19.12.06: bremsstrahlung added
!dr  31.07.07: bug fix: tein(j) --> tein(k)
cdr  nov.14:  function brems, replaces Gaunt factor function,
cdr           reaction scaling factor factkk removed from bremsstrahlung
cdr  nov 17:  comments. bremsstrahlung handling is currently connected
cdr           to 2d tab format. This is not generally correct, only if
cdr           the 2d tab comes from ADAS.  Should be handled in read_tab2d?
cdr           or as iftflg(..4) option.
cdr sept 18:  try to revive storage save mode.
cdr           Rationalization with options in feelei1.f:
cdr           as for EI processes: KK < 0: default, minimal (here: -1,-2)
cdr                                 KK > 0: kk=kread, read from datadase
cdr                                 KK = 0: else, simple models (via jelrrc)

      FUNCTION EIRENE_FEELRC1 (IRRC,K)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRRC, K
      REAL(DP) :: PLS, DELE, EE, EIRENE_FEELRC1, ZX, corsum,
     .            EIRENE_FTABRC1,
     .            DEIMIN, ELRC, EIRENE_ENERGY_RATE_COEFF, BREMS, Z,
     .            eirene_brems
      INTEGER :: KK
      LOGICAL :: LADAS

      EIRENE_FEELRC1=0.D0
      KK=NELRRC(IRRC)

      IF (KK < 0) THEN  ! DEFAULT RECOMBINATION MODELS, -1, -2

        SELECT CASE (KK)
        CASE (-1)  ! E + H+  --> H + hv
            ZX=EIONH/MAX(1.E-5_DP,TEIN(K))
            corsum=(-0.5_dp*zx+0.59)/(zx+0.59)
            EIRENE_FEELRC1=-(1.5+corsum)*TEIN(K)*EIRENE_FTABRC1(IRRC,K)
        CASE (-2)  ! E + He+  --> He + hv
            ZX=EIONHE/MAX(1.E-5_DP,TEIN(K))
            corsum=(-0.5_dp*zx+0.35)/(zx+0.35)
            EIRENE_FEELRC1=-(1.5+corsum)*TEIN(K)*EIRENE_FTABRC1(IRRC,K)
        CASE DEFAULT
          GOTO 999
        END SELECT

      ELSE IF (KK > 0) THEN   ! KK = KREAD

        IF (JELRRC(IRRC) == 1) THEN
          ELRC = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),0._DP,.FALSE.,0)
          ELRC=ELRC+FACRRC(IRRC,2)
          ELRC=EXP(MAX(-100._DP,ELRC))
          EIRENE_FEELRC1=-ELRC*DEIN(K)
        ELSEIF (JELRRC(IRRC) == 9) THEN
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))
          ELRC = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),PLS,.FALSE.,1)
          EE=MAX(-100._DP,ELRC+DEINL(K)+FACRRC(IRRC,2))
          EIRENE_FEELRC1=-EXP(EE)
        ELSE
          GOTO 999
        END IF
c
c  in some case (e.g. ADAS electron cooling rate tables), bremsstrahlung is
c  added on top of free-bound radiation. Subtract this contribution here,
c  to avoid double-counting.
c
        LADAS = EIRENE_IS_RTCEW_TAB2D(KK)  ! ifit=3 <--> ladas.
        IF (LADAS.AND.(NCHRGP(IPLS) /= 0)) THEN
          Z = NCHRGP(IPLS)
          BREMS = EIRENE_BREMS(TEIN(K),DEIN(K),Z)/ELCHA   ! W per ion --> eV/s  per ion
          EIRENE_FEELRC1=EIRENE_FEELRC1 + BREMS
        END IF
c
c  turn radiation loss into an electron energy cooling/heating rate
c  by adding potential energy transfer contribution
        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          EIRENE_FEELRC1=EIRENE_FEELRC1+DELE*EIRENE_FTABRC1(IRRC,K)
        END IF

      ELSEIF (KK == 0) THEN
        if (jelrrc(irrc) == -1) then
          EIRENE_FEELRC1=EELRC1(IRRC,1)*EIRENE_FTABRC1(IRRC,K)
        elseif (jelrrc(irrc) == -2) then
          EIRENE_FEELRC1=-1.5*tein(k)*EIRENE_FTABRC1(IRRC,K)
        else
          goto 999
        endif
      END IF

      RETURN

  999 continue
      write (iunout,*) 'error in feelrc1, KK= 0 option'
      call eirene_exit_own(1)
      END FUNCTION EIRENE_FEELRC1
