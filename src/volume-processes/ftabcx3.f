cdr  introduced march 2014. currenty largely identical to ftabpi3
c    already accommodates H.4 option (two parameter fits vs TII,PLS)

cdr  calls of ftabcx3 in fpath..: not ready, only for modcol=1 option




      FUNCTION EIRENE_FTABCX3 (IRCX,K)
c  evaluate charge exchange rate (1/s),
c  for CX process no. IRCX,
c         in cell no. K
c  input via common:
c         bulk collision partner: IPLS
c         process number in REACDAT structure:  KK=NREACX(IRCX)
c

cdr  (first stage):
c  FTABCX3 is currently only called in case MODC=eirene_idez(MODCLF(KK),3,5))=1,
c          i.e. rate depends only
c          on one single background parameter, Ti, and not on test particle energy.
c
cdr Sept 22: extended to MODC=3:  now two background parameters: ne and T (=Te=Ti)
cdr          tbd:  MODC=2: dependence on E0 and Ti
c
c
cdr  hard-wired: cut-off (density) parameter for fits: 1e8. TO BE CHECKED WITH XSTCX.F

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT  ! iunout
      USE EIRMOD_COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRCX, K
      REAL(DP) :: EIRENE_FTABCX3, DEIMIN, TII, PLS, TBCX,
     .            EIRENE_RATE_COEFF
      INTEGER :: KK, IPLSTI, MODC, IFLG
      INTEGER, EXTERNAL :: EIRENE_IDEZ

      TBCX=0.D0

      KK = NREACX(IRCX)
      MODC=EIRENE_IDEZ(MODCLF(KK),3,5)

cdr  for electron density dependence in condensed rates, e.g. MAR, MAD
cdr  option MODC=3
      DEIMIN=LOG(1.D8)
      PLS=MAX(DEIMIN,DEINL(K))

      IPLSTI=MPLSTI(IPLS)
      TII=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)

c
      if (modc.eq.1) then
c  input parameters for rate_coeff: ln(Ti) only
        IFLG=0
        TBCX = EIRENE_RATE_COEFF(KK,K,TII,0.0_DP,.TRUE.,IFLG)*
     .         FACRCX(IRCX,1)
        IF (IFTFLG(KK,2) < 100) TBCX = TBCX*DIIN(IPLS,K)

        EIRENE_FTABCX3 = TBCX

      elseif (modc.eq.3) then
c  input parameters for rate_coeff: ln(Ti), ln(ne)
c  isotopic shift included in TII  (hence: incorrectly also in Te)
        iflg=1
        TBCX = EIRENE_RATE_COEFF(KK,K,TII,PLS,.TRUE.,IFLG)*
     .       FACRCX(IRCX,1)
      IF (IFTFLG(KK,2) < 100) TBCX = TBCX*DIIN(IPLS,K)

      EIRENE_FTABCX3 = TBCX

      elseif (modc.eq.2) then
c  input parameters for rate_coeff: ln(Ti) and ln(E0)
        write (iunout,*) 'ftabcx3 called with modc=2 '
        write (iunout,*) 'option not ready, exit'
        call eirene_exit_own(1)
c  to be written
      else
        write (iunout,*) 'ftabcx3 called with invalid MODC, exit '
        write (iunout,*) 'kk, modc, ipls ',kk,modc,ipls
        call eirene_exit_own(1)
      endif

      RETURN
      END FUNCTION EIRENE_FTABCX3
