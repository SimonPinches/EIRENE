#ifndef NO_JSON
      subroutine eirene_wr0_json(json,this)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCOUPL
      USE EIRMOD_COMUSR
      USE EIRMOD_CLOGAU
      USE EIRMOD_CTEXT
      use json_module           !IGNORE

      implicit none

      class(json_core), intent(inout) :: json
      type(json_value), pointer, intent(inout) :: this
      type(json_value), pointer :: trgs, trg, flds, fld, conds, cond
      type(json_value), pointer :: intls, intl, ottls, ottl
      integer :: i

      call json%add(this,'LSYMET',lsymet)
      call json%add(this,'LBALAN',lbalan)
      call json%add(this,'LCOARSE',lcoarse)
      call json%add(this,'LSPRCL',lsprcl)

      call json%add(this,'NFLA',nfla)
      call json%add(this,'NCOND',ncond)
      call json%add(this,'NTRFRM',ntrfrm)
      call json%add(this,'NFULL',nfull)
      call json%add(this,'IBRAD',ibrad)
      call json%add(this,'IBPOL',ibpol)
      call json%add(this,'IBTOR',ibtor)

      call json%create_array(flds,'B2FLUIDS')
      do i = 1, nplsi
        call json%create_object(fld,'')
        call json%add(fld,'IPL',i)
        call json%add(fld,'IFLB',iflb(i))
        call json%add(fld,'FCTE',fcte(i))
        call json%add(fld,'BMASS',bmass(i))
        call json%add(fld,'LKINDP',lkindp(i))
        call json%add(flds,fld)
      end do
      call json%add(this,flds)

      IF (NCOND.GT.0) THEN
        call json%create_array(conds,'CONDENSE')
        DO I = 1, NCOND
          call json%create_object(cond,'')
          call json%add(cond,'CONDENSE_FROM',condense_from(i))
          call json%add(cond,'CONDENSE_TO',condense_to(i))
          call json%add(conds,cond)
        ENDDO
        call json%add(this,conds)
      ENDIF

      call json%add(this,'NCVP',ncvp)
      call json%add(this,'NCIP',ncip)
      call json%add(this,'NFCP',nfcp)
      call json%add(this,'NTARGI',ntargi)

      if (ntargi > 0) then

        call json%create_array(trgs,'TARGETS')
        do i = 1, ntargi
          call json%create_object(trg,'')
          call json%add(trg,'ITARG',i)
          call json%add(trg,'NDT',ndt(i))
          call json%add(trg,'NPTC',nptc(i))
          call json%add(trg,'NPTCM',nptcm(i))
          call json%add(trg,'NSPZI',nspzi(i))
          call json%add(trg,'NSPZE',nspze(i))
          call json%add(trg,'NEMOD',nemod(i))
          call json%add(trg,'ITYPE',itype(i))
          call json%add(trgs,trg)
        end do
        call json%add(this,trgs)
      end if

      call json%add(this,'CHGP',chgp)
      call json%add(this,'CHGEE',chgee)
      call json%add(this,'CHGEI',chgei)
      call json%add(this,'CHGMOM',chgmom)

      call json%add(this,'NAINB',nainb)
      if (nainb > 0) then
        call json%create_array(intls,'ADD_IN_TAL')
        do i = 1, nainb
          call json%create_object(intl,'')
          call json%add(intl,'IAIN',i)
          call json%add(intl,'NAINS',nains(i))
          call json%add(intl,'NAINT',naint(i))
          call json%add(intl,'TXTPLS',txtpls(i,12))
          call json%add(intl,'TXTPSP',txtpsp(i,12))
          call json%add(intl,'TXTPUN',txtpun(i,12))
          call json%add(intls,intl)
        end do
        call json%add(this,intls)
      end if

      call json%add(this,'NAOTB',naotb)
      if (naotb > 0) then
        call json%create_array(ottls,'ADD_OUT_TAL')
        do i = 1, naotb
          call json%create_object(ottl,'')
          call json%add(ottl,'IAOT',i)
          call json%add(ottl,'NAOTS',naots(i))
          call json%add(ottl,'NAOTT',naott(i))
          call json%add(ottls,ottl)
        end do
        call json%add(this,ottls)
      end if

      return
      end subroutine eirene_wr0_json
#else
      subroutine eirene_wr0_json
      return
      end subroutine eirene_wr0_json
#endif
