      subroutine eirene_colrad (ir, icell, iflavor, ivar, p1, p2, rate)

!   driver routine for collisional-radiative models

!   input:
!   ir:        reaction number, as stored in eirene arrays.
!   icell:     cell for which collisional-radiative model should be calulated
!   p1:        first parameter (usually:  log_e temperature,...)
!   p2:        second parameter  (if any, e.g.  log_e (density),...,log_e(test particle energy),...) 
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_comprt, only: iunout
 
      implicit none
 
      integer, intent(in) :: ir, icell, iflavor, ivar
      real(dp), intent(in) :: p1, p2
      real(dp), intent(out) :: rate

      real(dp) :: ALPCR, SCR, SCR_EXT, E_ALPCR, E_SCR, E_SCR_EXT,
     .            E_ALPCR_T, E_SCR_T, E_SCR_EXT_T
      integer :: i

      real(dp), allocatable, save :: pop0(:), pop1(:), pop2(:), q_ext(:)
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
          allocate(pop2(40))

          allocate(q_ext(40))    !   e.g. photo excitation rate for H*(n)
          Q_EXT = 0._DP
        end if

        if (.not.lvis_h(icell))  then

! rate needs to be calculated

          CALL EIRENE_H_COLRAD(P1, P2, Q_EXT, POP0, POP1, POP2,
     .                         ALPCR,    SCR,    SCR_EXT,
     .                         E_ALPCR,  E_SCR,  E_SCR_EXT,
     .                         E_ALPCR_T,E_SCR_T,E_SCR_EXT_T)

          do i = 1, nhcol_store
          
            select case(m_hcol(i))
            case (1)                      ! H.4  2.1.5
              h_stor(i,icell) = scr
            case (2)                      ! H.10 2.1.5
              h_stor(i,icell) = e_scr
            case (3)                      ! H.4  2.1.8
              h_stor(i,icell) = alpcr
            case (4)                      ! H.10 2.1.8
              h_stor(i,icell) = e_alpcr
            case (5)                      ! H.4  2.1.5PH
              h_stor(i,icell) = scr_ext
            case (6)                      ! H.10 2.1.5PH
              h_stor(i,icell) = e_scr_ext
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
            case (17)                     ! H.4  2.1.5PHa
              h_stor(i,icell) = pop2(3)
            case (18)                     ! H.4  2.1.5PHb
              h_stor(i,icell) = pop2(2)
            case (19)                     ! H.4  2.1.5PHc
              h_stor(i,icell) = pop2(4)
            case (20)                     ! H.4  2.1.5PHd
              h_stor(i,icell) = pop2(5)
            case (21)                     ! H.4  2.1.5PHe
              h_stor(i,icell) = pop2(6)
            case default
              write (iunout,*) ' ERROR IN COLRAD '
              write (iunout,*) ' REQUESTED RATE IS UNKNOWN '
              call eirene_exit_own(1)
            end select

          end do

          lvis_h(icell) = .true.
        end if

! rate is calculated
        rate = h_stor(ivar,icell)      
        return      
        

      else
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
        deallocate (pop2)
        deallocate (q_ext)
      end if

      return

      end subroutine eirene_colrad
