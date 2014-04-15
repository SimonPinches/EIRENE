c  introduced march 2014. currenty identical to ftabpi3

      FUNCTION EIRENE_FTABCX3 (IRCX,K)
c  evaluate charge exchange rate coefficient (or rate), 
c  for cx process no. IRCX,
c         in cell no. K
c  input via common: 
c         bulk collision partner: IPLS

c  FTABCX3 is only called in case MODCOL(3,2,IRCX)=1, i.e. rate depends only
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
 
      INTEGER, INTENT(IN) :: IRCX, K
      REAL(DP) :: TBPIC(9), EIRENE_FTABCX3, DEIMIN, TII, PLS, TBCX, 
     .            EIRENE_RATE_COEFF, ERATE
      INTEGER :: J, I, II, KK, IPLSTI
 
      TBCX=0.D0
      KK = NREACX(IRCX)


      DEIMIN=LOG(1.D8)
      PLS=MAX(DEIMIN,DEINL(K))
 
      IPLSTI=MPLSTI(IPLS)
      TII=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)

c  erate: not needed (only intermediate for H-COL option. Remove ! 
c  input parameters for rate_coeff: ln(Ti), ln(ne)
c 
      TBCX = EIRENE_RATE_COEFF(KK,TII,PLS,.TRUE.,0,ERATE)*
     .       FACRCX(IRCX,1)
      IF (IFTFLG(KK,2) < 100) TBCX = TBCX*DIIN(IPLS,K)
 
      EIRENE_FTABCX3 = TBCX
 
      RETURN
      END
 
