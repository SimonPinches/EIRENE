      subroutine eirene_colrad (ir, del, iflavor, ivar, 
     .                          icell, p1, p2, res)

!   driver routine for collisional-radiative models
!     calls internal CR code for cell no. ICELL
!     keeps all results from this call and 
!     marks the cells already visited (lvis_h) to avoid douple calls
!     for one and the same cell, for two or more differenct CRM output quantities 

!   input:
!   ir:        reaction number, as stored in eirene input arrays.
!   del:       energy shift (delpot parameter in input-card (eV).
!   iflavor:   choice of internal CR model. Currently iflavor=1: H-colrad
!              Soon: 
!              iflavor=4: He-colrad, iflavor=2:  H2-colrad
!   ivar:      
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
 
      integer, intent(in) :: ir, icell, iflavor, ivar
      real(dp), intent(in) :: p1, p2, del
      real(dp), intent(out) :: res

      real(dp) :: ALPCR, SCR, SCR_EXT, E_ALPCR, E_SCR, E_SCR_EXT,
     .            E_ALPCR_T, E_SCR_T, E_SCR_EXT_T
      integer :: i

      real(dp), allocatable, save :: pop0(:), pop1(:), pop_ext(:), 
     .                               q_ext(:)
      real(dp), allocatable, save :: h_stor(:,:)
      logical, allocatable, save :: lvis_h(:)
 
      
      if (.not. allocated(lvis_h)) then
        allocate (lvis_h(nrad))
        allocate (h_stor(nhcol_store,nrad))
        lvis_h = .false.
      end if

      if (iflavor == 1) then

! COLLISIONAL-RADIATIVE MODEL OF ATOMIC HYDROGEN
        if (.not.allocated(pop0)) then
          allocate(pop0(40))
          allocate(pop1(40))
          allocate(pop_ext(40))

          allocate(q_ext(40))    !   e.g. photo excitation rate for H*(n)
          Q_EXT = 0._DP
        end if

        if (.not.lvis_h(icell))  then
! cell number icell has not yet been visited so far in this run
! cr-model needs to be calculated

          CALL EIRENE_H_COLRAD(P1, P2, Q_EXT, POP0, POP1, POP_EXT,
     .                         ALPCR,    SCR,    SCR_EXT,
     .                         E_ALPCR,  E_SCR,  E_SCR_EXT,
     .                         E_ALPCR_T,E_SCR_T,E_SCR_EXT_T)
c
c  up to nhcol_store parameters from the cr-model are stored in cell ICELL 
          do i = 1, nhcol_store
          
            select case(m_hcol(i))
c  effective ionisation rate
            case (1)                      ! H.4  2.1.5
              h_stor(i,icell) = scr
            case (2)                      ! H.10 2.1.5
c  electron cooling rate coeff. 
c  note: with del=-13.6: this becomes the radiation loss rate coeff. alone
              h_stor(i,icell) = e_scr+del*scr
c  effective recombination rate
            case (3)                      ! H.4  2.1.8
              h_stor(i,icell) = alpcr
            case (4)                      ! H.10 2.1.8
c  electron cooling/heating rate coeff. (both signs possible)
c  note: with del=+13.6: this becomes the radiation loss rate coeff. alone 
              h_stor(i,icell) = e_alpcr+del*alpcr
c  external source driven ionisation rate, e.g. photo-excitation driven ionisation
            case (5)                      ! H.4  2.1.5PH
              h_stor(i,icell) = scr_ext
            case (6)                      ! H.10 2.1.5PH
              h_stor(i,icell) = e_scr_ext
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
c  population coefficients, coupling to H+ ion
              h_stor(i,icell) = pop1(6)
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
              write (iunout,*) ' ERROR IN COLRAD '
              write (iunout,*) ' REQUESTED RATE IS UNKNOWN '
              call eirene_exit_own(1)
            end select

          end do

          lvis_h(icell) = .true.
        end if

! rate is calculated
        res = h_stor(ivar,icell)      
        return      
        

      else   ! iflavor .ne.1
         write (iunout,*) ' REQUESTED COLLISIONAL-RADIATIVE MODEL' //
     .                    ' NOT AVAILABLE '
         call eirene_exit_own
      end if

      entry eirene_colrad_reinit
      
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
      end if

      return

      end subroutine eirene_colrad
