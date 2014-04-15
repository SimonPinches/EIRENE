      SUBROUTINE EIRENE_OUTLST

C  DIAGNOSTIC OUTPUT AFTER LAST HISTORY, BUT STILL INSIDE PARTICLE LOOP 

C  MOSTLY: CPU EFFICIANCY ISSUES



 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CLAST
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      REAL(DP) :: SMMEAN
      INTEGER :: IWR, IRCX, IREL, IRPI

C  1) EFFICIANCY OF RECECTION SAMPLING IN VELOCX, VELOEL, VELOPI 
      call EIRENE_leer(1)
      iwr=0
      do ircx=1,nrcxi
        if (iflrcx(ircx).gt.0) then
          if (iwr.eq.0) then
            WRITE (iunout,*) 'REJECTION SAMPLING EFFICIENCY IN VELOCX '
            write (iunout,*)
     .      'IRCX, TOTAL NO. OF CALLS TO VELOCX, MEAN NO. OF SAMPLING'
            iwr=1
          endif
          SMMEAN=xcmean(ircx)/(ncmean(ircx)+eps60)
          CALL EIRENE_MASJ2R('IRCX, NCMEAN, SMMEAN    ',
     .                        IRCX, NCMEAN(IRCX),SMMEAN) 
        endif
      enddo
      call EIRENE_leer(1)
      iwr=0
      do irel=1,nreli
        if (iflrel(irel).gt.0) then
          if (iwr.eq.0) then
            WRITE (iunout,*) 'REJECTION SAMPLING EFFICIENCY IN VELOEL '
            write (iunout,*)
     .      'IREL, TOTAL NO. OF CALLS TO VELOEL, MEAN NO. OF SAMPLING'
            iwr=1
          endif
          SMMEAN=xemean(irel)/(nemean(irel)+eps60)
          CALL EIRENE_MASJ2R('IREL, NEMEAN, SMMEAN    ',
     .                        IREL, NEMEAN(IREL),SMMEAN) 
        endif
      enddo
      call EIRENE_leer(1)
      iwr=0
      do irpi=1,nrpii
        if (iflrpi(irpi).gt.0) then
          if (iwr.eq.0) then
            WRITE (iunout,*) 'REJECTION SAMPLING EFFICIENCY IN VELOPI '
            write (iunout,*)
     .      'IRPI, TOTAL NO. OF CALLS TO VELOPI, MEAN NO. OF SAMPLING'
            iwr=1
          endif
          SMMEAN=xpmean(irpi)/(npmean(irpi)+eps60)
          CALL EIRENE_MASJ2R('IRPI, NPMEAN, SMMEAN    ',
     .                        IRPI, NPMEAN(IRPI),SMMEAN) 
        endif
      enddo
      call EIRENE_leer(1)
      RETURN
      END
