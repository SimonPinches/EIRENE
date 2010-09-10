 
 
      subroutine EIRENE_replace_string (inchar,rem,rep,iunout)

      use eirmod_precision
      use eirmod_parmmod

      implicit none

      character(len=*), intent(inout) :: inchar
      character(len=*), intent(in) :: rem, rep
      character(len=2*len(inchar)) :: outchar
      integer, intent(in) :: iunout
      integer :: lin, lout, lrem, lrep, i, io, ii, linc
 
      linc = len(inchar)
 
      lout = len(outchar)
      outchar = repeat(' ',lout)
 
      lin=len_trim(inchar)
      lrem=len_trim(rem)
      lrep=len_trim(rep)
 
      io = 0
      ii = 0
 
      i = index(inchar(ii+1:lin),rem(1:lrem)) - 1
 
      do while (i >= 0)
 
        if (io+i > lout) exit
        outchar(io+1:io+i) = inchar(ii+1:ii+i)
        io = io + i
        if (io+lrep > lout) exit
        outchar(io+1:io+lrep) = rep(1:lrep)
        io = io + lrep
 
        ii = ii + i + lrem
 
        i = index(inchar(ii+1:lin),rem(1:lrem)) - 1
 
      end do
 
      if (i > = 0) then
        write (iunout,*) ' ERROR IN REPLACE_STRING '
        write (iunout,*) ' STRING IS TOO SHORT TO HOLD ALL REPLACEMENTS'
        write (iunout,*) ' INCHAR = ',inchar
        write (iunout,*) ' REMCHAR = ',rem
        write (iunout,*) ' REPCHAR = ',rep
        write (iunout,*) ' STRING SHORTENED TO ',outchar(1:linc)
        inchar(1:linc) = outchar(1:linc)
        return
      end if
 
      outchar(io+1:io+lin-ii) = inchar(ii+1:lin)
      io = io + lin-ii
 
      if (io > linc) then
        write (iunout,*) ' ERROR IN REPLACE_STRING '
        write (iunout,*) ' STRING IS TOO SHORT TO HOLD ALL REPLACEMENTS'
        write (iunout,*) ' INCHAR = ',inchar
        write (iunout,*) ' REMCHAR = ',rem
        write (iunout,*) ' REPCHAR = ',rep
        write (iunout,*) ' STRING SHORTENED TO ',outchar(1:linc)
      end if
 
      inchar = outchar(1:linc)
 
      return
      end subroutine EIRENE_replace_string
 
