C 27.6.05:  PHV_NROTA, PHV_NROTPH REMOVED
cdr  nov. 15:  comments,  irds --> irei
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
      USE EIRMOD_CZT1
      USE EIRMOD_PHOTON
 
      IMPLICIT NONE
 
      real(dp) :: tpb1, tpb2, second_own
      INTEGER, INTENT(IN) :: ICAL
      INTEGER :: I, IRPI, IREI
 
!pb      tpb1 = second_own()
 
      IF (ICAL == 0) THEN
        NRCX=0
        NREL=0
        NRPI=0
        NRDS=0
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
        NRDS=MAX(1,NRDS)
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

cdr 
      IPATDS = 0
      IPMLDS = 0
      IPIODS = 0
      IPPLDS = 0
      DO IREI=1,NRDS
        ipatds(IREI,0)=COUNT(PATDS(IREI,1:) > 0)  ! amongst all natm species there are ipatds (<= natm) 
cdr                                                 atomic species which appear a secondaries, 
cdr                                                 with one or more per atomic species iatm 
        IF (ipatds(IREI,0).GT.0) THEN          ! inserted by Derek Harting 26.03
             IPATDS(IREI,1:ipatds(IREI,0))=PACK( (/ (i,i=1,natm) /),
     .                                     PATDS(IREI,1:) > 0)
cdr  IPATDS(IREI,...)=iatm means:  one or more secondaries of species iatm
cdr  the arrays patds,,...p2nd, contain furthermore the information: 
cdr  "how many" of this secondary species iatm arise after process irei.
        END IF
        ipmlds(IREI,0)=COUNT(PMLDS(IREI,1:) > 0)
        IF (ipmlds(IREI,0).GT.0) THEN          ! inserted by Derek Harting 26.03
             IPMLDS(IREI,1:ipmlds(IREI,0))=PACK( (/ (i,i=1,nmol) /),
     .                                     PMLDS(IREI,1:) > 0)
        END IF
        ipiods(IREI,0)=COUNT(PIODS(IREI,1:) > 0)
        IF (ipiods(IREI,0).GT.0) THEN         ! inserted by Derek Harting 26.03.
             IPIODS(IREI,1:ipiods(IREI,0))=PACK( (/ (i,i=1,nion) /),
     .                                     PIODS(IREI,1:) > 0)
        END IF
        ipplds(IREI,0)=COUNT(PPLDS(IREI,1:) > 0)
        IF (ipplds(IREI,0).GT.0) THEN         ! inserted by Derek Harting 26.03.
             IPPLDS(IREI,1:ipplds(IREI,0))=PACK( (/ (i,i=1,npls) /),
     .                                     PPLDS(IREI,1:) > 0)
        END IF
      END DO

cdr:  same as above, for PI processes 
      IPATPI = 0
      IPMLPI = 0
      IPIOPI = 0
      IPPLPI = 0
      DO IRPI=1,NRPI
        ipatpi(IRPI,0)=COUNT(PATPI(IRPI,1:) > 0)
        IF (ipatpi(IRPI,0).GT.0) then         ! inserted by Derek Harting 26.03.
          IPATPI(IRPI,1:ipatpi(IRPI,0))=PACK( (/ (i,i=1,natm) /),
     .                                  PATPI(IRPI,1:) > 0)
        endif
        ipmlpi(IRPI,0)=COUNT(PMLPI(IRPI,1:) > 0)
        IF (ipmlpi(IRPI,0).GT.0) then         ! inserted by Derek Harting 26.03.
          IPMLPI(IRPI,1:ipmlpi(IRPI,0))=PACK( (/ (i,i=1,nmol) /),
     .                                  PMLPI(IRPI,1:) > 0)
        endif
        ipiopi(IRPI,0)=COUNT(PIOPI(IRPI,1:) > 0)
        IF (ipiopi(IRPI,0).GT.0) then         ! inserted by Derek Harting 26.03.
          IPIOPI(IRPI,1:ipiopi(IRPI,0))=PACK( (/ (i,i=1,nion) /),
     .                                  PIOPI(IRPI,1:) > 0)
        endif
        ipplpi(IRPI,0)=COUNT(PPLPI(IRPI,1:) > 0)
        IF (ipplpi(IRPI,0).GT.0) then         ! inserted by Derek Harting 26.03.
          IPPLPI(IRPI,1:ipplpi(IRPI,0))=PACK( (/ (i,i=1,npls) /),
     .                                  PPLPI(IRPI,1:) > 0)
        endif
      END DO
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for packs und counts ',tpb2-tpb1
!pb        tpb1 = tpb2
 
 
      RETURN
      END
