      SUBROUTINE EIRENE_GEOUSR_BIASED_GARCHING
c
c  version : 12.03.98 21:02
c
c======================================================================
C***  PREPARE DATA FOR ADDITIONAL SURFACES
c======================================================================
      USE EIRMOD_PRECISION
      use EIRMOD_PARMMOD
      use EIRMOD_CADGEO
      use EIRMOD_CTRCEI
      use EIRMOD_CCONA
      use EIRMOD_CGEOM
      use EIRMOD_CGRID
      use EIRMOD_CLGIN
      use EIRMOD_CPOLYG
      use EIRMOD_CINIT
      use EIRMOD_COMPRT
      use EIRMOD_CLOGAU

      IMPLICIT NONE

      logical :: first, normalcase, hlp_found
      integer :: onetwo(8),limpos(8),xpolpos(8),ypolpos(8)
      character(80) :: geometry_comment
      REAL(DP), PARAMETER :: hlp_tol=0.001

      INTEGER :: I, J, K, NBITS, M, N, L, NASMOD
      REAL(DP) :: HLP_P1, HLP_P2
      save first,onetwo,limpos,geometry_comment
      data first /.true./
c======================================================================
c*** At the first invocation, read the data from the Eirene input file
c*** and define the grid corners
c
      IF (NLPLG) THEN
      if(first) then
!       read(iunin,'(a80)') 'geometry_comment', already done in calling program
        write (iunout,*) 'GEOUSR: VERSION: BIASED-GARCHING'
        do i=1,8
          onetwo(i)=0
          limpos(i)=0
        end do
C INNER LEFT TARGET
        xpolpos(1)=1
        ypolpos(1)=npoint(1,1)
        xpolpos(2)=nr1st
        ypolpos(2)=npoint(1,1)
C OUTER RIGHT TARGET
        if(npplg.le.3) then
C  single null
          xpolpos(3)=1
          ypolpos(3)=npoint(2,npplg)
          xpolpos(4)=nr1st
          ypolpos(4)=npoint(2,npplg)
c  double null
        else if(npplg.eq.6) then
          xpolpos(3)=1
          ypolpos(3)=npoint(2,3)-1
          xpolpos(4)=nr1st
          ypolpos(4)=npoint(2,3)-1
        end if

c  double null only
        if(npplg.eq.6) then
C OUTER LEFT TARGET
          xpolpos(5)=1
          ypolpos(5)=npoint(1,4)+1
          xpolpos(6)=nr1st
          ypolpos(6)=npoint(1,4)+1
C INNER RIGHT TARGET
          xpolpos(7)=1
          ypolpos(7)=npoint(2,6)
          xpolpos(8)=nr1st
          ypolpos(8)=npoint(2,6)
        end if

        normalcase=.true.

        nasmod=max(npplg/3,1)*4


        do i=1,nasmod

cdr either read onetwo, limpos
cdr or     read limpos, onetwo


          read(iunin,*) onetwo(i),limpos(i)
          if(onetwo(i).lt.0) then
            normalcase=.false.
            onetwo(i)=-onetwo(i)
            read(iunin,*) xpolpos(i),ypolpos(i)
          end if
        end do
C
        first=.false.
!       write(iunout,*) 'GEOMETRY FOR'
!       write(iunout,'(a80)') geometry_comment

        if (npplg.le.3) then
          if(onetwo(1).eq.0) then
            write(iunout,*) 'Geometry fixup skipped'
            goto 1001
          end if

          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO INNER LEFT TARGET'')') onetwo(1),limpos(1)
          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO OUTER LEFT TARGET'')') onetwo(2),limpos(2)
          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO INNER RIGHT TARGET'')') onetwo(3),limpos(3)
          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO OUTER RIGHT TARGET'')') onetwo(4),limpos(4)
c  double null
        else if (npplg.eq.6) then
          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO INNER LEFT TARGET'')') onetwo(1),limpos(1)
          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO INNER LEFT TARGET'')') onetwo(2),limpos(2)
          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO OUTER RIGHT TARGET'')') onetwo(3),limpos(3)
          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO OUTER RIGHT TARGET'')') onetwo(4),limpos(4)
          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO OUTER LEFT TARGET'')') onetwo(5),limpos(5)
          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO OUTER LEFT TARGET'')') onetwo(6),limpos(6)
          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO INNER RIGHT TARGET'')') onetwo(7),limpos(7)
          write(iunout,'(''P'',i1,'' FOR SEGMENT '',i3,
     1     '' LINKED TO INNER RIGHT TARGET'')') onetwo(8),limpos(8)
        else
          write(iunout,*) 'Case NPPLG = ',NPPLG,' not coded. '
        end if
      end if
c
c*** Switch off the additional surfaces corresponding to the targets,
c*** that is, the surfaces between the ones to be linked to the grid
c*** corners.
c
      n=max(npplg/3,1)*4
      if (n.ne.nasmod) then
        write (iunout,*) 'ERROR IN GEOUSR_biased-garching'
        write (iunout,*) 'N, NASMOD =', N,NASMOD
        call eirene_exit_own(1)
      endif
c     print '(/(2i8))',(limpos(i),onetwo(i),i=1,n)
      do 990 i=1,n
        if(onetwo(i).eq.2) then
          m=limpos(i)
          hlp_p1=p2(1,m)
          hlp_p2=p2(2,m)
c          print *,'onetwo=2. igjum0= ',igjum0(j),',  i,hlp_p1,hlp_p2 =
c          print *,i,hlp_p1,hlp_p2
          do 980 l=1,nlimi
            do j=1,nlimi
              hlp_found=.false.
c              print *,'lgjum0,j,m = ',lgjum0(j),j,m
              if(igjum0(j)==0 .and. j.ne.m .and. iliin(j).eq.1) then
                if(abs(hlp_p1-p1(1,j)).le.hlp_tol .and.
     .                             abs(hlp_p2-p1(2,j)).le.hlp_tol) then
                  hlp_found=.true.
                  hlp_p1=p2(1,j)
                  hlp_p2=p2(2,j)
                else if(abs(hlp_p1-p2(1,j)).le.hlp_tol .and.
     .                             abs(hlp_p2-p2(2,j)).le.hlp_tol) then
                hlp_found=.true.
                hlp_p1=p1(1,j)
                hlp_p2=p1(2,j)
              end if
              if(hlp_found) then
c                print *,j,hlp_p1,hlp_p2
c*** Check whether this segment is marked as a target edge
                do k=1,n
                  if(j.eq.limpos(k)) then
                    if(onetwo(k).eq.1) then
                      go to 990
                    else
                      write(iunout,*) 'geousr_biased:',
     ,                             ' something is wrong with ',
     ,                             'the target chain definition.'
                      write(iunout,*)
     .                      'Check the data on the target edges ',
     ,                      'at the very end of the Eirene input file.'
                      call EIRENE_exit_own(1)
                      end if
                    end if
                  end do
c*** Switch off the segment
                  igjum0(j)=1
                  print *,'geousr_biased: segment ',j,'  is turned off'
                  go to 980
                end if
              end if
            end do
c*** The chain is broken
            print *,'geousr_biased: the chain is broken'
            go to 990
  980     continue
        end if
  990 continue
 1001 continue
C
C  ANFANG: MODIFY GEOMETRY
C
      n=max(npplg/3,1)*4
      if (n.ne.nasmod) then
        write (iunout,*) 'ERROR IN GEOUSR, Vs.: geousr_biased-garching'
        write (iunout,*) 'N, NASMOD =', N,NASMOD
        call eirene_exit_own(1)
      endif
c
      do i=1,max(npplg/3,1)*4
        if(onetwo(i).eq.1) then
          p1(1,limpos(i))=xpol(xpolpos(i),ypolpos(i))
          p1(2,limpos(i))=ypol(xpolpos(i),ypolpos(i))
        else if(onetwo(i).eq.2) then
          p2(1,limpos(i))=xpol(xpolpos(i),ypolpos(i))
          p2(2,limpos(i))=ypol(xpolpos(i),ypolpos(i))
        end if
      end do
c=====================================================
C  ABSCHALTEN NICHT ERREICHBARER ODER DOPPELT VORHANDENER FLAECHEN
C
C   HIERHER: LGJUM1, LGJUM2 SETZEN ZUR BESCHLEUNIGUNG (NICHT UNBEDINGT
C   NOETIG)
C
C   LGJUM1(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN AUF J SITZT
C
C   LGJUM2(J,I)=.TRUE. :
C   ABSCHALTEN DES ERSTEN SCHNITTPUNKTES MIT FLAECHE I, FALLS
C   TEILCHEN AUF J SITZT (FALLS I EINE FLAECHE ZWEITER ORDNUNG IST)
C
C   DEFAULTS: LGJUM1(J,J)=.TRUE. FUER EBENE FLAECHEN,
C             LGJUM2(J,J)=.TRUE. FUER FLAECHEN ZWEITER ORDNUNG
C
C
C   LGJUM3(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN IN ZELLE J SITZT
C
      IF (NOPTIM >= NSURF) THEN

      if(normalcase) then
        write(iunout,*) 'Setting IGJUM3 to 1 for',NSURF,NLIMI
        if (nlimpb.ge.nlimps) then
          do I=1,NLIMI
!pb            do J=1,NOPTIM
            do J=1,NSURF
              IGJUM3(J,I)=1
            end do
          end do
        else
          nbits=bit_size(1)
          do I=1,NLIMI
!pb            do J=1,NOPTIM
            do J=1,NSURF
              call EIRENE_bitset(igjum3,0,noptim,j,i,1,nbits)
            end do
          end do
        endif
      else
        write(iunout,*) 'Setting IGJUM3 to 0 for ',NSURF,NLIMI
        if (nlimpb.ge.nlimps) then
          do I=1,NLIMI
!pb            do J=1,NOPTIM
            do J=1,NSURF
              IGJUM3(J,I)=0
            end do
          end do
        else
          nbits=bit_size(1)
          do I=1,NLIMI
!pb            do J=1,NOPTIM
            do J=1,NSURF
              call EIRENE_bitset(igjum3,0,noptim,j,i,0,nbits)
            end do
          end do
        endif
      endif
      END IF

      END IF
C
C  SET SOME VOLUMES EXPLICITLY
C
C
C  MODIFY REFLECTION MODEL AT TARGET PLATES
C
C
      RETURN
      END
