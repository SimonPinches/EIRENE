c  introduced sept 2014. currenty largely identical to ftabpi3
c  already accommodates H.4 option (two parameter fits vs TII,PLS)

c  calls of ftabel3 in fpath..: not ready




      FUNCTION EIRENE_FTABEL3 (IREL,K)
c  evaluate elastic collision rate (1/s), 
c  for el process no. IREL,
c         in cell no. K
c  input via common: 
c         bulk collision partner: IPLS
c 
c  ftabel3 is currently not called. 

c  soon (first stage):
c  FTABEL3 is currently only called in case MODCOL(3,2,IREL)=1, i.e. rate depends only
c          on background parameters, not on test particle energy.
C  EVEN THAT MAY NOT BE TRUE:::
c 
c
c  hard wired: cut off (density) parameter for fits: 1e8
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IREL, K
      REAL(DP) :: EIRENE_FTABEL3, DEIMIN, TII, PLS, TBEL, 
     .            EIRENE_RATE_COEFF
      INTEGER :: KK, IPLSTI
 
      TBEL=0.D0
      KK = NREAEL(IREL)


      DEIMIN=LOG(1.D8)
      PLS=MAX(DEIMIN,DIINL(IPLS,K))
 
      IPLSTI=MPLSTI(IPLS)
      TII=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)

c  input parameters for rate_coeff: tii=ln(Ti), pls=ln(ni)
c 
      TBEL = EIRENE_RATE_COEFF(KK,K,TII,PLS,.TRUE.,0)*
     .       FACREL(IREL,1)
      IF (IFTFLG(KK,2) < 100) TBEL = TBEL*DIIN(IPLS,K)
 
      EIRENE_FTABEL3 = TBEL
 
      RETURN
      END
 
