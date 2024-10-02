      subroutine eirene_wr0_json(json,this)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCOUPL
      USE EIRMOD_COMUSR
      USE EIRMOD_CTEXT
      use json_module           !IGNORE

      implicit none

      class(json_core),intent(inout) :: json
      type(json_value),pointer, intent(inout) :: this
      type(json_value),pointer :: trgs, trg, flds, fld, prts, prt
      type(json_value),pointer :: intls, intl, ottls, ottl
      integer :: i, j

      call json%add(this,'LSYMET',lsymet)
      call json%add(this,'LBALAN',lbalan)
      call json%add(this,'LCOARSE',lcoarse)

      call json%add(this,'NFLA',nfla)
      call json%add(this,'NCUTB',ncutb)
      call json%add(this,'NCUTL',ncutl)
      call json%add(this,'IMF',mshfrm)

      call json%create_array(flds,'B2FLUIDS')
      do i = 1, nplsi
        call json%create_object(fld,'')
        call json%add(fld,'IPL',i)
        call json%add(fld,'IFLB',iflb(i))
        call json%add(fld,'FCTE',fcte(i))
        call json%add(fld,'BMASS',bmass(i))
        call json%add(flds,fld)
      end do
      call json%add(this,flds)

      call json%add(this,'NDXA',ndxa)
      call json%add(this,'NDYA',ndya)
      call json%add(this,'NTARGI',ntargi)

      if (ntargi > 0) then
        call json%add(this,'NTGPRT',ntgprt(1:ntargi))

        call json%create_array(trgs,'TARGETS')
        do i = 1, ntargi
          call json%create_object(trg,'')
          call json%add(trg,'ITARG',i)
          call json%create_array(prts,'PARTS')
          do j = 1, ntgprt(i)
            call json%create_object(prt,'')
            call json%add(prt,'NDT',ndt(i,j))
            call json%add(prt,'NINCT',ninct(i,j))
            call json%add(prt,'NIXY',nixy(i,j))
            call json%add(prt,'NTIN',ntin(i,j))
            call json%add(prt,'NTEN',nten(i,j))
            call json%add(prt,'NIFLG',niflg(i,j))
            call json%add(prt,'NPTC',nptc(i,j))
            call json%add(prt,'NPTCM',nptcm(i,j))
            call json%add(prt,'NSPZI',nspzi(i,j))
            call json%add(prt,'NSPZE',nspze(i,j))
            call json%add(prt,'NEMOD',nemod(i,j))
            call json%add(prts,prt)
          end do
          call json%add(trg,prts)
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

      end subroutine eirene_wr0_json
