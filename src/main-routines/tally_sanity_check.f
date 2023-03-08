      subroutine eirene_tally_sanity_check(istra)

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_cestim
      use eirmod_coutau
      use eirmod_ccona
      USE EIRMOD_COMPRT
     , ,ONLY : IUNOUT
      USE EIRMOD_CLOGAU
     , ,ONLY : NLSPCSCL, NLSPCSCL_ATM, NLSPCSCL_MOL, NLSPCSCL_ION,
     ,         NLSPCSCL_PHOT
      use eirmod_comusr
     , ,ONLY : NATMI, NMOLI, NIONI, NPHOTI
      implicit none
      integer, intent(in) :: istra

      if (nlspcscl_atm) then
        if (lpaat) 
     .    call eirene_calc_tally_error(natmi, natm, natmi, 
     .                                 paati(:,istra), 'PAATI')
        if (lpaml) 
     .    call eirene_calc_tally_error(nmoli, nmol, natmi,
     .                                  pamli(:,istra), 'PAMLI')
        if (lpaio) 
     .    call eirene_calc_tally_error(nioni, nion, natmi,
     .                                  paioi(:,istra), 'PAIOI')
        if (lpapht) 
     .    call eirene_calc_tally_error(nphoti, nphot, natmi,
     .                                  paphti(:,istra),'PAPHTI')
        IF (lprfaat)
     .    call eirene_calc_tally_error(natmi, natm, natmi,
     .                                  prfaai(:,istra),'PRFAAI')
        IF (lprfaml)
     .    call eirene_calc_tally_error(nmoli, nmol, natmi,
     .                                  prfami(:,istra),'PRFAMI')
        IF (lprfaio)
     .    call eirene_calc_tally_error(nioni, nion, natmi,
     .                                  prfaii(:,istra),'PRFAII')
        IF (lprfapht)
     .    call eirene_calc_tally_error(nphoti, nphot, natmi,
     .                                  prfaphti(:,istra),'PRFAPHTI')
      end if

      if (nlspcscl_mol) then
        if (lpmat) 
     .    call eirene_calc_tally_error(natmi, natm, nmoli, 
     .                                 pmati(:,istra), 'PMATI')
        if (lpmml) 
     .    call eirene_calc_tally_error(nmoli, nmol, nmoli,
     .                                  pmmli(:,istra), 'PMMLI')
        if (lpmio) 
     .    call eirene_calc_tally_error(nioni, nion, nmoli,
     .                                  pmioi(:,istra), 'PMIOI')
        if (lpmpht) 
     .    call eirene_calc_tally_error(nphoti, nphot, nmoli,
     .                                  pmphti(:,istra),'PMPHTI')
        IF (lprfmat)
     .    call eirene_calc_tally_error(natmi, natm, nmoli,
     .                                  prfmai(:,istra),'PRFMAI')
        IF (lprfmml)
     .    call eirene_calc_tally_error(nmoli, nmol, nmoli,
     .                                  prfmmi(:,istra),'PRFMMI')
        IF (lprfmio)
     .    call eirene_calc_tally_error(nioni, nion, nmoli,
     .                                  prfmii(:,istra),'PRFMII')
        IF (lprfmpht)
     .    call eirene_calc_tally_error(nphoti, nphot, nmoli,
     .                                  prfmphti(:,istra),'PRFMPHTI')
      end if

      if (nlspcscl_ion) then
        if (lpiat) 
     .    call eirene_calc_tally_error(natmi, natm, nioni, 
     .                                 piati(:,istra), 'PIATI')
        if (lpiml) 
     .    call eirene_calc_tally_error(nmoli, nmol, nioni,
     .                                  pimli(:,istra), 'PIMLI')
        if (lpiio) 
     .    call eirene_calc_tally_error(nioni, nion, nioni,
     .                                  piioi(:,istra), 'PIIOI')
        if (lpipht) 
     .    call eirene_calc_tally_error(nphoti, nphot, nioni,
     .                                  piphti(:,istra),'PIPHTI')
        IF (lprfiat)
     .    call eirene_calc_tally_error(natmi, natm, nioni,
     .                                  prfiai(:,istra),'PRFIAI')
        IF (lprfiml)
     .    call eirene_calc_tally_error(nmoli, nmol, nioni,
     .                                  prfimi(:,istra),'PRFIMI')
        IF (lprfiio)
     .    call eirene_calc_tally_error(nioni, nion, nioni,
     .                                  prfiii(:,istra),'PRFIII')
        IF (lprfipht)
     .    call eirene_calc_tally_error(nphoti, nphot, nmoli,
     .                                  prfiphti(:,istra),'PRFIPHTI')
      end if


      if (nlspcscl_phot) then
        if (lpphat) 
     .    call eirene_calc_tally_error(natmi, natm, nphoti, 
     .                                 pphati(:,istra), 'PPHATI')
        if (lpphml) 
     .    call eirene_calc_tally_error(nmoli, nmol, nphoti,
     .                                  pphmli(:,istra), 'PPHMLI')
        if (lpphio) 
     .    call eirene_calc_tally_error(nioni, nion, nphoti,
     .                                  pphioi(:,istra), 'PPHIOI')
        if (lpphpht) 
     .    call eirene_calc_tally_error(nphoti, nphot, nphoti,
     .                                  pphphti(:,istra),'PPHPHTI')
        IF (lprfphat)
     .    call eirene_calc_tally_error(natmi, natm, nphoti,
     .                                  prfiai(:,istra),'PRFIAI')
        IF (lprfphml)
     .    call eirene_calc_tally_error(nmoli, nmol, nphoti,
     .                                  prfphmi(:,istra),'PRFPHMI')
        IF (lprfphio)
     .    call eirene_calc_tally_error(nioni, nion, nphoti,
     .                                  prfphii(:,istra),'PRFPHII')
        IF (lprfphpht)
     .    call eirene_calc_tally_error(nphoti, nphot, nphoti,
     .                                  prfphphti(:,istra),'PRFPHPHTI')
      end if

      return
      
      contains

      subroutine eirene_calc_tally_error(n1, n1d, n2, ar, arname)
      use eirmod_ccona
      implicit none
      integer, intent(in) :: n1, n1d, n2
      real(dp), target, intent(in) :: ar((n1d+1)*(n2+1))
      character(*), intent(in) :: arname
      real(dp), pointer :: p2(:,:)
      real(dp) :: error, sm, rel_err
      integer :: i

      p2(0:n1d,0:n2) => ar
      do i = 1, n1
        sm = sum(p2(i,1:n2))
        error = p2(i,0) - sm
        rel_err = error / max(abs(p2(i,0)),abs(sm),eps60)
        if (abs(error) > eps10) then
          write (iunout,*) 'PROBLEM FOUND IN TALLY_SANITY_CHECK'
          write (iunout,*) 'DISCREPANCY FOUND IN ARRAY ',arname
          write (iunout,*) 'INDEX = ',i
          write (iunout,*) 'ERROR = ',error
          write (iunout,*) 'REL. ERROR = ',rel_err
          write (iunout,'(a,a1,i2,a6,es25.16)') 
     .                     arname,'(',i,'0) = ',p2(i,0)
          write (iunout,'(a4,a,a1,i2,a2,i2,a4,es25.16)') 
     .                     'sum(',arname,'(',i,'1:',n2,') = ',sm
        end if       
      end do

      return
      end subroutine eirene_calc_tally_error

      end subroutine eirene_tally_sanity_check
