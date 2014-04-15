cdr  feb 2014  :  iftflg :  select rate or rate coeficient 
cdr  march 2014:  density parameter added, H.4 option for rate coeff, ifit=2
 
      FUNCTION EIRENE_FTABPI3 (IRPI,K)
c  evaluate general heavy particle impact rate coefficient (or rate), 
c  for pi process no. IRPI,
c         in cell no. K
c  input via common: 
c         bulk collision partner: IPLS

c  FTABPI3 is only called in case MODCOL(4,2,IRPI)=1, i.e. rate depends only
c          on background parameters, not on test particle energy.
c 
c
c  hard wired: cut off (density) parameter for fits: 1e8
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IRPI, K
      REAL(DP) :: TBPIC(9), EIRENE_FTABPI3, DEIMIN, TII, PLS, TBPI, 
     .            EIRENE_RATE_COEFF, ERATE
      INTEGER :: J, I, II, KK, IPLSTI
 
      TBPI=0.D0
      KK = NREAPI(IRPI)


      DEIMIN=LOG(1.D8)
      PLS=MAX(DEIMIN,DEINL(K))
 
      IPLSTI=MPLSTI(IPLS)
      TII=TIINL(IPLSTI,K)+ADDPI(IRPI,IPLS)

c  erate: not needed (only intermediate for H-COL option. Remove ! 
c  input parameters for rate_coeff: ln(Ti), ln(ne)
c 
      TBPI = EIRENE_RATE_COEFF(KK,TII,PLS,.TRUE.,0,ERATE)*
     .       FACRPI(IRPI,1)
      IF (IFTFLG(KK,2) < 100) TBPI = TBPI*DIIN(IPLS,K)
 
      EIRENE_FTABPI3 = TBPI
 
      RETURN
      END
 
