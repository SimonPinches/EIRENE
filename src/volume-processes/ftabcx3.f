cdr  introduced march 2014. currenty largely identical to ftabpi3
c    already accommodates H.4 option (two parameter fits vs TII,PLS)

cdr  calls of ftabcx3 in fpath..: not ready, only for modcol=1 option




      FUNCTION EIRENE_FTABCX3 (IRCX,K)
c  evaluate charge exchange rate (1/s), 
c  for cx process no. IRCX,
c         in cell no. K
c  input via common: 
c         bulk collision partner: IPLS
c 

cdr(first stage):
c  FTABCX3 is currently only called in case MODCOL(3,2,IRCX)=1, i.e. rate depends only
c          on background parameters, not on test particle energy.
c 
c
cdr  hard wired: cut off (density) parameter for fits: 1e8. TO BE CHECKED WITH XSTCX.F
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IRCX, K
      REAL(DP) :: EIRENE_FTABCX3, DEIMIN, TII, PLS, TBCX, 
     .            EIRENE_RATE_COEFF
      INTEGER :: KK, IPLSTI
 
      TBCX=0.D0
      KK = NREACX(IRCX)


      DEIMIN=LOG(1.D8)
      PLS=MAX(DEIMIN,DIINL(IPLS,K))
 
      IPLSTI=MPLSTI(IPLS)
      TII=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)

cdr  :  check these next lines are only for modcol=1 !

c  erate: not needed (only intermediate for H-COL option. Remove ! 
c  input parameters for rate_coeff: ln(Ti), ln(ne)
c 
      TBCX = EIRENE_RATE_COEFF(KK,K,TII,PLS,.TRUE.,0)*
     .       FACRCX(IRCX,1)
      IF (IFTFLG(KK,2) < 100) TBCX = TBCX*DIIN(IPLS,K)
 
      EIRENE_FTABCX3 = TBCX
 
      RETURN
      END
 
