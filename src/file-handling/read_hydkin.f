cdr   This routine is only for internal use at FJZ
cdr   Purpose:  establish an interface to online A&M data repository and toolbox

cdr   called from:  setup_hydkin_reactions

cdr   feb 2014:  only started to add comments, then copied to read_table1_hydkin
cdr              for generalization

cdr   current implementation: started for rate coefficients (%RTC$) only.

cdr   july 16:   more error exits, to avoid code crashes when reading 1D tabulated data
cdr              currently this routine expects hard wired HYDKIN, CxHy format. 
cdr   to be done: distuingish between reading data, and automatted construction
cdr               of full blocks 4a,b,c,d,5 from a HYDKIN output

      subroutine EIRENE_read_hydkin (ir,filename,h123,reac,crc,
     .                        r1mn,r1mx,
     .                        e_el,e_k,lffl)

c  input:
c          ir:  reaction nunmber input block 4
c
c
c          e_el:  ?? available only for hydrocarbons?
c          e_k :  ?? available only for hydrocarbons?
c          lffl:  ??
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_comprt, only: iunout
 
      implicit none
 
      integer, intent(in) :: ir
      character(len=*), intent(in) :: reac, filename
      logical, intent(in) :: lffl
      character(4), intent(inout) :: h123
      character(3), intent(inout) :: crc
      real(dp) , intent(out) :: r1mn, r1mx, e_el, e_k
      character(132) :: zeile
      character(len=len(reac)+10) :: cpreac
      integer :: ianf, iend, ll, io, ie, iflg

      type(hydkin_data), pointer :: hp
 
      open (unit=28+ifoff,file=filename)

c  skip blank lines at top of file
  
      zeile = repeat(' ',len(zeile))
 

      do while (index(zeile,'Default energy mesh') == 0)
         read (28+ifoff,'(A132)',iostat=io,end=100) zeile
      end do
 
100   allocate (hp)
      hp%reacname = reac

! find number of temperatures 
      read (28+ifoff,'(A132)') zeile
      ianf = index(zeile,'=')
      read (zeile(ianf+1:),*) hp%ntemps
 
      do while (index(zeile,'EeVDef') == 0)
        read (28+ifoff,'(A132)',iostat=io,end=110) zeile
      end do
 
110   allocate (hp%temps(hp%ntemps))
      allocate (hp%rates(hp%ntemps))
      allocate (hp%ratio(hp%ntemps))
 
! read energies
      do ie=1, hp%ntemps
        read (28+ifoff,*) hp%temps(ie)
      end do
 
      r1mn = hp%temps(1)
      r1mx = hp%temps(hp%ntemps)
 
      e_el = 0._dp
      e_k = 0._dp
 
! find specified reaction
 
      cpreac = 'Eirname = '//adjustl(reac)
      ll = len_trim(cpreac)
 
      zeile = repeat(' ',len(zeile))
      do while (index(zeile,cpreac) == 0)
        read (28+ifoff,'(A132)',iostat=io,end=120) zeile
      end do
 
120   if (io .ne. 0) then
        write (iunout,*) ' ERROR READING REACTION FROM HYDKIN DATABASE '
        write (iunout,*) ' FILE ',filename, ' NOT FOUND'
        write (iunout,*) ' REACTION IS ',reac
        call EIRENE_exit_own(1)
      end if
 
! reaction found
! now read data
      zeile = repeat(' ',len(zeile))
      do while (index(zeile,'E_el') == 0)
        read (28+ifoff,'(A132)',iostat=io,end=130) zeile
      end do

130   if (io .ne. 0) then
        write (iunout,*) ' ERROR READING REACTION FROM HYDKIN DATABASE '
        write (iunout,*) ' FILE IS ',filename
        write (iunout,*) ' E_el NOT FOUND'
        call EIRENE_exit_own(1)
      end if

      ianf = index(zeile,'=')+1
      read (zeile(ianf:),*) e_el
 
      zeile = repeat(' ',len(zeile))
      do while (index(zeile,'E_K') == 0)
        read (28+ifoff,'(A132)',iostat=io,end=140) zeile
      end do

140   if (io .ne. 0) then
        write (iunout,*) ' ERROR READING REACTION FROM HYDKIN DATABASE '
        write (iunout,*) ' FILE IS ',filename
        write (iunout,*) ' E_K NOT FOUND'
        call EIRENE_exit_own(1)
      end if

      ianf = index(zeile,'=')+1
      read (zeile(ianf:),*) e_k

c  next: type of reaction:  EI, (=DS), CX, EL, RC, PI (=II) 
      zeile = repeat(' ',len(zeile))
      do while (index(zeile,'RPrT') == 0)
        read (28+ifoff,'(A132)',iostat=io,end=150) zeile
      end do

150   if (io .ne. 0) then
        write (iunout,*) ' ERROR READING REACTION FROM HYDKIN DATABASE '
        write (iunout,*) ' FILE IS ',filename
        write (iunout,*) ' RPrT NOT FOUND'
        call EIRENE_exit_own(1)
      end if
 
      ianf = index(zeile,'''')+1
      iend = ianf-1 + index(zeile(ianf:),'''') -1
      hp%rprt = zeile(ianf:iend)
 
      zeile = repeat(' ',len(zeile))
      do while (index(zeile,'RName') == 0)
        read (28+ifoff,'(A132)',iostat=io,end=160) zeile
      end do

160   if (io .ne. 0) then
        write (iunout,*) ' ERROR READING REACTION FROM HYDKIN DATABASE '
        write (iunout,*) ' FILE IS ',filename
        write (iunout,*) ' RName NOT FOUND'
        call EIRENE_exit_own(1)
      end if 

      ianf = index(zeile,'''')+1
      iend = ianf-1 + index(zeile(ianf:),'''') -1
      hp%reac_string = zeile(ianf:iend)
 
      zeile = repeat(' ',len(zeile))
      do while (index(zeile,'RData') == 0)
        read (28+ifoff,'(A132)',iostat=io,end=170) zeile
      end do

170   if (io .ne. 0) then
        write (iunout,*) ' ERROR READING REACTION FROM HYDKIN DATABASE '
        write (iunout,*) ' FILE IS ',filename
        write (iunout,*) ' RData NOT FOUND'
        call EIRENE_exit_own(1)
      end if 
 
! read rate coefficients:  cm**3/s
      do ie=1, hp%ntemps
        read (28+ifoff,*) hp%rates(ie)
      end do
 
      close (unit=28+ifoff)
 
      do ie=1, hp%ntemps-1
        hp%ratio(ie) = (hp%rates(ie+1) - hp%rates(ie)) /
     .                 (hp%temps(ie+1) - hp%temps(ie))
      end do
 
      if (reacdat(ir)%lrtc) then
        write (iunout,*) ' RATE COEFFICIENT ALREADY SPECIFIED',
     .                   ' FOR REACTION ',ir
        write (iunout,*) ' CHECK SPECIFICATION OF REACTIONS'
        deallocate (hp)
        call EIRENE_exit_own(1)
      end if
 
      reacdat(ir)%lrtc = .true.
      allocate(reacdat(ir)%rtc)
      nullify (reacdat(ir)%rtc%adas)
      nullify (reacdat(ir)%rtc%line)
      nullify (reacdat(ir)%rtc%poly)

      reacdat(ir)%rtc%hyd => hp
      reacdat(ir)%rtc%ifit = 4
      REACDAT(IR)%RTC%RC1MIN = hp%temps(1)
      REACDAT(IR)%RTC%RC1MAX = hp%temps(hp%ntemps)
      REACDAT(IR)%RTC%RC2MIN = 0._dp
      REACDAT(IR)%RTC%RC2MAX = huge(1._dp)
      REACDAT(IR)%RTC%FP1L = 0._DP
      REACDAT(IR)%RTC%FP1R = 0._DP
      REACDAT(IR)%RTC%FP2B = 0._DP
      REACDAT(IR)%RTC%FP2T = 0._DP
      REACDAT(IR)%RTC%JFEX1MN = 0
      REACDAT(IR)%RTC%JFEX1MX = 0
      REACDAT(IR)%RTC%JFEX2MN = 0
      REACDAT(IR)%RTC%JFEX2MX = 0
 
      if (lffl) then
        h123 = 'H.2 '
        MODCLF(IR)=MODCLF(IR)+100
        IFLG=2
C  DEFAULT RATE COEFFICIENT
        IFTFLG(IR,IFLG)=0
 
! find crc
 
        if (index(hp%rprt,'CX') > 0) then
           crc = 'CX '
           iswr(ir) = 3
        else if (index(hp%rprt,'_R') + index(hp%rprt,'R-DR') > 0) then
           crc = 'RC '
           iswr(ir) = 6
        else if (index(hp%rprt,'_DE') + index(hp%rprt,'_DI') +
     .          index(hp%rprt,'_I') + index(hp%rprt,'I-DI') +
     .          index(hp%rprt,'_CAD') > 0) then
           crc = 'EI '
           iswr(ir) = 1
        else
           write (iunout,*) ' UNKNOWN REACTION TYPE ',HP%RPRT
           write (iunout,*) ' USED IN REACTION ',reac
           write (iunout,*) ' CHECK SPECIFICATION OF REACTIONS'
           call EIRENE_exit_own(1)
        end if
      endif
 
      return
      end subroutine EIRENE_read_hydkin
