C 27.6.05:  PHV_NROTA, PHV_NROTPH REMOVED
cdr  nov. 15:  comments,  irds --> irei
<<<<<<< HEAD
<<<<<<< HEAD
cdr  april 16:  added: fail safe (exit) step in case of more than one (distinct) bulk 
cdr             secondaries.
cdr             This is temporarily necessary, as a consequence of making the
cdr             (bulk) ion energy sources eapl, empl, eipl species dependent
cdr             We are not aware of any application of eirene, in which this new error exit
cdr             would be activated.  
=======
!pb  APR  16:  ipplds -> ipplei
<<<<<<< HEAD
>>>>>>> variable IPPLDS renamed to IPPLEI
=======
=======
!pb  APR  16:  ipplds -> ipplei, pplds -> pplei
<<<<<<< HEAD
>>>>>>> variable PPLDS renamed to PPLEI
!pb  APR  16:  ipatds -> ipatei
<<<<<<< HEAD
>>>>>>> variable IPATDS renamed to IPATEI
=======
=======
!pb  APR  16:  ipatds -> ipatei, patds -> patei
<<<<<<< HEAD
>>>>>>> variable PATDS renamed to PATEI
!pb  APR  16:  ipmlds -> ipmlei
<<<<<<< HEAD
>>>>>>> variable IPMLDS renamed to IPMLEI
=======
=======
!pb  APR  16:  ipmlds -> ipmlei, pmlds -> pmlei
<<<<<<< HEAD
>>>>>>> variable PMLDS renamed to PMLEI
!pb  APR  16:  ipiods -> ipioei
>>>>>>> variable IPIODS renamed to IPIOEI
=======
!pb  APR  16:  ipiods -> ipioei, piods -> pioei
<<<<<<< HEAD
>>>>>>> variable PIODS renamed to PIOEI
=======
!pb  MAY  16:  nrds   -> nrei
>>>>>>> variable NRDS renamed to NREI
C
      SUBROUTINE EIRENE_SETAMD(ICAL)
C
C  SET ATOMIC AND MOLECULAR DATA: DRIVER
C
CDR  CALLED IN INITIALIZATION PHASE OF RUN
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMXS
      USE EIRMOD_COMSOU
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY : IUNOUT
      USE EIRMOD_CZT1
      USE EIRMOD_PHOTON
 
      IMPLICIT NONE
 
      real(dp) :: tpb1, tpb2, second_own
      INTEGER, INTENT(IN) :: ICAL
      INTEGER :: I, IRPI, IREI, IERROR
 
!pb      tpb1 = second_own()
 
      IF (ICAL == 0) THEN
        NRCX=0
        NREL=0
        NRPI=0
        NREI=0
        NREC=0
        NBGV=0
        NROT=0
        CALL EIRENE_XSECTA_PARAM
        CALL EIRENE_XSECTM_PARAM
        CALL EIRENE_XSECTI_PARAM
        CALL EIRENE_XSECTP_PARAM
        CALL EIRENE_XSECTPH_PARAM
 
        NRCX=MAX(1,NRCX)
        NREL=MAX(1,NREL)
        NRPI=MAX(1,NRPI)
        NREI=MAX(1,NREI)
        NREC=MAX(1,NREC)
        NBGV=MAX(1,NBGV)
        NROT=MAX(1,NROT)
 
        CALL EIRENE_SET_PARMMOD(2)
        CALL EIRENE_ALLOC_COMXS(2)
        CALL EIRENE_ALLOC_COMSOU(2)
        CALL EIRENE_ALLOC_CZT1(2)

        MAXSPC(0:4) = (/ NPHOTI,NATMI,NMOLI,NIONI,NPLSI /)
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for setamd(0) ',tpb2-tpb1
!pb        tpb1 = tpb2
 
        RETURN
      ELSE
        CALL EIRENE_INIT_CMDTA(2)
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for init_cmdta ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      END IF
 
      NRCXI=0
      NRELI=0
      NRPII=0
      NREII=0
      NRRCI=0
      NRBGI=0
csw 27jul2011
      NPBGKP=0 !VK
      NPBGKA=0 !VK
      NPBGKM=0 !VK
      NPBGKI=0 !VK
csw

      CALL EIRENE_XSECTA
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsecta ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      CALL EIRENE_XSECTM
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsectm ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      CALL EIRENE_XSECTI
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsecti ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      CALL EIRENE_XSECTP
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsectp ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      CALL EIRENE_XSECTPH
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsectph ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      CALL EIRENE_CONDENSE
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for condense ',tpb2-tpb1
!pb        tpb1 = tpb2
 
c
cdr  set some further assistant arrays, for ei and pi processes:
cdr  accumulated information from A, M, I ,P and PH for particle processes 'ei' and 'pi'. 
cdr  These array are stored in comxs and are used for scoring
cdr  tallies in update.f (tracklength) and collide.f (coll. estim) exclusively
cdr  They are for indirect indexing, in loops over secondary species.
cdr  e.g. rather than 
cdr                   do iat=1,natmi
cdr         now:      
cdr                   do i   =1,ipatds(irei,0)   (<=natmi,  possibly much shorter loop) 
cdr                      iat = ipatds(irei,i)    (now we know: iat is a secondary indeed)
cdr                      inum= patds(irei,iat)   (there are inum secondaries of species iat)
cdr                      ...
cdr                   enddo
cdr 
      IERROR = 0

      IPATEI = 0
      IPMLEI = 0
      IPIOEI = 0
cdr   IPPHDS = 0   ARRAY IPPHDS IS STILL MISSING, NO PHOTON SECONDARIES IN EI REACTIONS.
      IPPLEI = 0
      DO IREI=1,NREI
        ipatei(IREI,0)=COUNT(PATEI(IREI,1:) > 0)  ! amongst all natm species there are ipatei(...,0) (<= natm)
cdr                                                 atomic species which appear as secondaries, 
cdr                                                 with one or more per atomic species iatm 
        IF (ipatei(IREI,0).GT.0) THEN
             IPATEI(IREI,1:ipatei(IREI,0))=PACK( (/ (i,i=1,natm) /),
     .                                     PATEI(IREI,1:) > 0)
cdr  IPATEI(IREI,...)=iatm means:  one or more secondaries of species iatm
c
cdr  the arrays patei,...,pplei, and p2nd, contain the further information:
cdr  "how many" of this secondary species iatm arise after process irei.
        END IF
        ipmlei(IREI,0)=COUNT(PMLEI(IREI,1:) > 0)
        IF (ipmlei(IREI,0).GT.0) THEN
             IPMLEI(IREI,1:ipmlei(IREI,0))=PACK( (/ (i,i=1,nmol) /),
     .                                     PMLEI(IREI,1:) > 0)
        END IF
        ipioei(IREI,0)=COUNT(PIOEI(IREI,1:) > 0)
        IF (ipioei(IREI,0).GT.0) THEN
             IPIOEI(IREI,1:ipioei(IREI,0))=PACK( (/ (i,i=1,nion) /),
     .                                     PIOEI(IREI,1:) > 0)
        END IF
        ipplei(IREI,0)=COUNT(PPLEI(IREI,1:) > 0)       
        IF (ipplei(IREI,0).GT.0) THEN
             IPPLEI(IREI,1:ipplei(IREI,0))=PACK( (/ (i,i=1,npls) /),
     .                                     PPLEI(IREI,1:) > 0)
        END IF
        if (ipplei(IREI,0) > 1) then
          IERROR = IERROR + 1
          write (iunout,*) 'MORE THAN ONE BULK ION SPECIES SPECIFIED ',
     .          'AS SECONDARY PARTICLE OF EI REACTION IREI = ',IREI
        end if
      END DO

cdr:  same as above, for PI processes 
      IPATPI = 0
      IPMLPI = 0
      IPIOPI = 0
cdr   IPPHPI = 0   ARRAY IPPHPI IS STILL MISSING, NO PHOTON SECONDARIES IN PI REACTIONS.
      IPPLPI = 0
      DO IRPI=1,NRPI
        ipatpi(IRPI,0)=COUNT(PATPI(IRPI,1:) > 0)
        IF (ipatpi(IRPI,0).GT.0) then         
          IPATPI(IRPI,1:ipatpi(IRPI,0))=PACK( (/ (i,i=1,natm) /),
     .                                  PATPI(IRPI,1:) > 0)
        endif
        ipmlpi(IRPI,0)=COUNT(PMLPI(IRPI,1:) > 0)
        IF (ipmlpi(IRPI,0).GT.0) then        
          IPMLPI(IRPI,1:ipmlpi(IRPI,0))=PACK( (/ (i,i=1,nmol) /),
     .                                  PMLPI(IRPI,1:) > 0)
        endif
        ipiopi(IRPI,0)=COUNT(PIOPI(IRPI,1:) > 0)
        IF (ipiopi(IRPI,0).GT.0) then        
          IPIOPI(IRPI,1:ipiopi(IRPI,0))=PACK( (/ (i,i=1,nion) /),
     .                                  PIOPI(IRPI,1:) > 0)
        endif
        ipplpi(IRPI,0)=COUNT(PPLPI(IRPI,1:) > 0)
        IF (ipplpi(IRPI,0).GT.0) then        
          IPPLPI(IRPI,1:ipplpi(IRPI,0))=PACK( (/ (i,i=1,npls) /),
     .                                  PPLPI(IRPI,1:) > 0)
        endif
        if (ipplpi(IRPI,0) > 1) then
          IERROR = IERROR + 1
          write (iunout,*) 'MORE THAN ONE BULK ION SPECIES SPECIFIED ',
     .          'AS SECONDARY PARTICLE OF PI REACTION IRPI = ',IRPI
        end if
      END DO
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for packs und counts ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      if (ierror > 0) then
         write (iunout,*) 'only a temporary fail safe step'
         write (iunout,*) 'contact eirene group at fzj, if this occurs'  
         write (iunout,*) 'CALCULATION ABANDONNED '
         CALL EIRENE_EXIT_OWN(1)
      end if

      RETURN
      END
