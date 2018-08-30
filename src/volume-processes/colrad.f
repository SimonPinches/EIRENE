cdr jan 18:  distinct from solps4.3 version: e_alpcr correct now.
cdr          (electron cooling/heating terms associated with recombination
cdr feb 18:  l_ext, q_ext, lopaque, pop_esc: must not change,
c            after first call, otherwise: reset LVIS 
c            so far: q_ext not connected (l_ext=.false.)
cdr may 18:  add population escape factors pop_esc(40,40), for hydrogen atom.
cdr          default: optically thin: pop_esc=1


      subroutine eirene_colrad (ir, icrm, ivar, 
     .                          icell, p1, p2, res)

!   driver routine for collisional-radiative models
!     calls internal CR code no. icrm, for cell no. ICELL
!     keeps all results from this call and 
!     marks the cells already visited (for icrm=1: lvis_h) to avoid douple calls
!     for one and the same cell, for two or more differenct CRM output quantities 

!   input:
!   ir:        reaction number, as stored in eirene input arrays.
!   icrm:      choice of internal CR model. Currently icrm=1: H-colrad
!              Soon: 
!              icrm=4: He-colrad, icrm=2:  H2-colrad
!   ivar:      this call to colrad pick one particular CR variable for cell icell.
!              Currently, for H_COLRAD, there are
!
!   icell:     cell for which collisional-radiative model should be calculated
!   p1:        first parameter (usually:  log_e temperature,...)
!   p2:        second parameter  (if any, e.g.  log_e (density),...,log_e(test particle energy),...)

!   output:
!   res:       result, for cell no. icell. 
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_comprt, only: iunout
 
      implicit none
 
      integer, intent(in) :: ir, icell, icrm, ivar
      real(dp), intent(in) :: p1, p2
      real(dp), intent(out) :: res

      real(dp) :: ALPCR, SCR, SCR_EXT, E_ALPCR, E_SCR, E_SCR_EXT
ctt  .           ,E_ALPCR_T, E_SCR_T, E_SCR_EXT_T   these arrays are for testing only
      integer :: i, irow_esc, icol_esc, irc

      real(dp), allocatable, save :: pop0(:), pop1(:), pop_ext(:), 
     .                               q_ext(:), 
     .                               pop_esc(:,:)
      real(dp), allocatable, save :: h_stor(:,:)
      logical, allocatable, save :: lvis_h(:)
      logical :: l_ext
 
c  try to avoid repeated calls to CR model in same plasma grid cell
c      for the current run/iteration/time-cycle     
      if (.not. allocated(lvis_h)) then
        allocate (lvis_h(nrad))
        allocate (h_stor(nhcol_store,nrad))
        lvis_h = .false.
      end if

      if (icrm == 1) then

! COLLISIONAL-RADIATIVE MODEL OF ATOMIC HYDROGEN
        if (.not.allocated(pop0)) then
          allocate(pop0(40))
          allocate(pop1(40))
          allocate(pop_ext(40))
          allocate(q_ext(40))      !   e.g. photo excitation rate for H*(n)
          allocate(pop_esc(40,40)) !   line population escape factor (default:==1)

          POP_ESC =1.0_DP          !   default: all transitions are optically thin

cdr  cummulate all population escape factors for internal CR model.
cdr  either read  
cdr              via reaction cards (block 4) 
cdr           or via line-emission cards (block 12) 
cdr  
          do irc = 1, nreac
cdr  scan over all reaction decks (from block 4 and/or block 12)
            if (reacdat(irc)%lrtc) then  ! rate coeff 
              if (reacdat(irc)%rtc%ifit == 5) then ! internal crm 
                if (reacdat(irc)%rtc%crm%iflav == 1) then ! h-col
                  irow_esc = reacdat(irc)%rtc%crm%irow_esc
                  icol_esc = reacdat(irc)%rtc%crm%icol_esc
                  if ((irow_esc > 0) .and. (icol_esc > 0)) then
                    pop_esc(irow_esc,icol_esc) = 
     .                      reacdat(irc)%rtc%crm%pop_esc
                    write (iunout,*) 'H-crm: pop_esc set for transition'
                    write (iunout,*) irc,irow_esc,'p-->',icol_esc, 
     .                               'to ',
     .                               pop_esc(irow_esc,icol_esc)          
                  end if
                end if
              end if
            end if
            if (reacdat(irc)%lrtcew) then ! energy rate coeff 
              if (reacdat(irc)%rtcew%ifit == 5) then ! crm 
                if (reacdat(irc)%rtcew%crm%iflav == 1) then ! h-col
                  irow_esc = reacdat(irc)%rtcew%crm%irow_esc
                  icol_esc = reacdat(irc)%rtcew%crm%icol_esc
                  if ((irow_esc > 0) .and. (icol_esc > 0)) then
                    pop_esc(irow_esc,icol_esc) = 
     .                      reacdat(irc)%rtcew%crm%pop_esc
                    write (iunout,*) 'H-crm: pop_esc set for transition'
                    write (iunout,*) irc,irow_esc,'e-->',icol_esc, 
     .                               'to ',
     .                               pop_esc(irow_esc,icol_esc)          
                  end if
                end if
              end if
            end if
            if (reacdat(irc)%loth) then ! other rate coeff, pop_coef 
              if (reacdat(irc)%oth%ifit == 5) then ! crm 
                if (reacdat(irc)%oth%crm%iflav == 1) then ! h-col
                  irow_esc = reacdat(irc)%oth%crm%irow_esc
                  icol_esc = reacdat(irc)%oth%crm%icol_esc
                  if ((irow_esc > 0) .and. (icol_esc > 0)) then
                    pop_esc(irow_esc,icol_esc) = 
     .                      reacdat(irc)%oth%crm%pop_esc
                    write (iunout,*) 'H-crm: pop_esc set for transition'
                    write (iunout,*) irc,irow_esc,'o-->',icol_esc, 
     .                               'to ', 
     .                               pop_esc(irow_esc,icol_esc)        
                  end if
                end if
              end if
            end if
         end do
        end if

        Q_EXT = 0._DP
        L_EXT = .FALSE.

        if (.not.lvis_h(icell))  then
! cell number ICELL has not yet been visited so far in this run
! cr-model needs to be calculated.
! In later calls, for this ICELL, 
!    we assume Q_EXT, L_EXT. LOPAQUE, POP_ESC to be unchanged !

          CALL EIRENE_H_COLRAD(P1, P2, Q_EXT, L_EXT,
     .                         POP0, POP1, POP_EXT,
     .                         ALPCR,    SCR,    SCR_EXT,
     .                         E_ALPCR,  E_SCR,  E_SCR_EXT,
ctt  .                        ,E_ALPCR_T,E_SCR_T,E_SCR_EXT_T
     .                         POP_ESC)
c
c  up to nhcol_store parameters from the cr-model are stored in cell ICELL
C  TBD: if .NOT.L_EXT: only case(1) to case(16) are available  
          do i = 1, nhcol_store
          
            select case(m_hcol(i))
c  effective ionisation rate
            case (1)                      ! H.4  2.1.5
              h_stor(i,icell) = scr
            case (2)                      ! H.10 2.1.5
c  electron cooling rate coeff. e_scr is negative from h-colrad
c  note: with delpot=-13.6 (input): this becomes the radiation loss rate coeff. alone
              h_stor(i,icell) = -e_scr
c  effective recombination rate
            case (3)                      ! H.4  2.1.8
              h_stor(i,icell) = alpcr
            case (4)                      ! H.10 2.1.8
c  electron cooling/heating rate coeff. (both signs possible. loss: negative e_alpcr))
c  note: with delpot=+13.6 (input): this becomes the radiation loss rate coeff. alone 
              h_stor(i,icell) = -e_alpcr
c  external source driven ionisation rate, e.g. photo-excitation driven ionisation
            case (5)                      ! H.4  2.1.5PH
              h_stor(i,icell) = scr_ext
            case (6)                      ! H.10 2.1.5PH
              h_stor(i,icell) = -e_scr_ext
c  population coefficients, coupling to ground state H(1) atom
            case (7)                      ! H.4  2.1.5a
              h_stor(i,icell) = pop1(3)
            case (8)                      ! H.4  2.1.5b
              h_stor(i,icell) = pop1(2)
            case (9)                      ! H.4  2.1.5c
              h_stor(i,icell) = pop1(4)
            case (10)                     ! H.4  2.1.5d
              h_stor(i,icell) = pop1(5)
            case (11)                     ! H.4  2.1.5e
              h_stor(i,icell) = pop1(6)
c  population coefficients, coupling to H+ ion
            case (12)                     ! H.4  2.1.8a
              h_stor(i,icell) = pop0(3)
            case (13)                     ! H.4  2.1.8b
              h_stor(i,icell) = pop0(2)
            case (14)                     ! H.4  2.1.8c
              h_stor(i,icell) = pop0(4)
            case (15)                     ! H.4  2.1.8d
              h_stor(i,icell) = pop0(5)
            case (16)                     ! H.4  2.1.8e
              h_stor(i,icell) = pop0(6)
c  population coefficients, coupling to external source of excitation (e.g. photons)
c  only availabel if L_EXT=.TRUE. in call to H_COLRAD
            case (17)                     ! H.4  2.1.5PHa
              h_stor(i,icell) = pop_ext(3)
            case (18)                     ! H.4  2.1.5PHb
              h_stor(i,icell) = pop_ext(2)
            case (19)                     ! H.4  2.1.5PHc
              h_stor(i,icell) = pop_ext(4)
            case (20)                     ! H.4  2.1.5PHd
              h_stor(i,icell) = pop_ext(5)
            case (21)                     ! H.4  2.1.5PHe
              h_stor(i,icell) = pop_ext(6)
            case default
              write (iunout,*) ' ERROR IN COLRAD, M_HCOL(I) '
              write (iunout,*) ' REQUESTED RATE FROM H_COLRAD ?? '
              call eirene_exit_own(1)
            end select

          end do

          lvis_h(icell) = .true.
        end if

cdr  Result RES from h_colrad is now calculated.
cdr  it may be a rate, an energly loss rate or a reduced population coefficient
        res = h_stor(ivar,icell)      
        return      
        

      else   ! icrm .ne.1
         write (iunout,*) ' REQUESTED COLLISIONAL-RADIATIVE MODEL' //
     .                    ' NOT AVAILABLE '
         WRITE (iunout,*) 'icrm ',icrm
         WRITE (iunout,*) 'ir   ',ir
         call eirene_exit_own
      end if

      entry eirene_colrad_reinit
cdr this must be done after each internal iteration or time cycle
      
      if (allocated(lvis_h)) then
         lvis_h = .false.
      end if
      
      if (allocated(h_stor)) then
         h_stor = 0._dp
      end if
      
      return

      entry eirene_dealloc_colrad

      if (allocated(lvis_h)) deallocate (lvis_h)
      if (allocated(h_stor)) deallocate (h_stor)
      if (allocated(pop0)) then
        deallocate (pop0)
        deallocate (pop1)
        deallocate (pop_ext)
        deallocate (q_ext)
        deallocate (pop_esc)
      end if

      return

      end subroutine eirene_colrad
