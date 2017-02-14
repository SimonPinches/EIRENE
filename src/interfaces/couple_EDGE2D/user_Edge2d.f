C ===== buildscript: user_Edge2d
C ===== SOURCE: broad_usr.f


      SUBROUTINE EIRENE_BROAD_USR
      IMPLICIT NONE
      RETURN
      END
C ===== SOURCE: geomd.f
C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D
C=======================================================================
      SUBROUTINE EIRENE_GEOMD(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST

      character(110) :: zeile
      REAL(DP) :: br(0:ndxp,0:ndyp,4),bz(0:ndxp,0:ndyp,4)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) ::
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
      REAL(DP) :: DX, DY
      INTEGER :: I0, I0E, IX, IY, I1, I2, I3, I4, IPART, I, J

      ndxa=0
      ndya=0

      read (30+ifoff,*)
      read (30+ifoff,*)
      read (30+ifoff,*)
      read (30+ifoff,*)

1     continue
      read (30+ifoff,'(a110)',end=99) zeile
      i0=index(zeile,'(')
      i0e=index(zeile,')')
      read (zeile(i0+1:i0e-1),*) ix,iy
      ndxa=max(ndxa,ix)
      ndya=max(ndya,iy)
      i1=index(zeile,': (')
      i2=index(zeile(i1+3:),')')+i1+2
      read (zeile(i1+3:i2-1),*) br(ix,iy,4),bz(ix,iy,4)
      i3=index(zeile(i2+1:),'(')+i2
      i4=i3+index(zeile(i3+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,3),bz(ix,iy,3)

      read (30+ifoff,'(a110)') zeile

      read (30+ifoff,'(a110)') zeile
      i1=index(zeile,'(')
      i2=index(zeile,')')
      read (zeile(i1+1:i2-1),*) br(ix,iy,1),bz(ix,iy,1)
      i3=i2+index(zeile(i2+1:),'(')
      i4=i2+index(zeile(i2+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,2),bz(ix,iy,2)

      read (30+ifoff,*)
      goto 1


99    continue
      ndxa=ndxa-1
      ndya=ndya-1
C
      DO 1015 IY=1,NDYA
        DO 1014 IX=1,NDXA
          X1(IX)=br(ix,iy,1)
          Y1(IX)=bz(ix,iy,1)
          X2(IX)=br(ix,iy,2)
          Y2(IX)=bz(ix,iy,2)
          X3(IX)=br(ix,iy,4)
          Y3(IX)=bz(ix,iy,4)
          X4(IX)=br(ix,iy,3)
          Y4(IX)=bz(ix,iy,3)
1014    CONTINUE
        CALL EIRENE_MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,
     .       PUX,PUY,PVX,PVY,NDXA,
     .       NR1ST,IY)
1015  CONTINUE
C
C SEARCH FOR THE CUTS
C
      IPART=1
      NPOINT(1,IPART)=1
      IY=1
      DO IX=1,NDXA
        DX=BR(IX+1,IY,1)-BR(IX,IY,2)
        DY=BZ(IX+1,IY,1)-BZ(IX,IY,2)
        IF (DX*DX+DY*DY.GT.EPS10) THEN
C CUT GEFUNDEN
          NPOINT(2,IPART)=IX+1+IPART-1
          IPART=IPART+1
          NPOINT(1,IPART)=IX+1+IPART-1
        ENDIF
      ENDDO
      NPOINT(2,IPART)=NDXA+IPART
C
      NPLP=IPART
C
      DO IY=1,NDYA
        DO IPART=1,NPLP
          DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
            XPOL(IY,IX)=BR(IX-(IPART-1),IY,1)
            YPOL(IY,IX)=BZ(IX-(IPART-1),IY,1)
          ENDDO
          XPOL(IY,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,IY,2)
          YPOL(IY,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,IY,2)
        ENDDO
      ENDDO
C INTRODUCE OUTERMOST RADIAL POLYGON
      DO IPART=1,NPLP
        DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
          XPOL(NDYA+1,IX)=BR(IX-(IPART-1),NDYA,4)
          YPOL(NDYA+1,IX)=BZ(IX-(IPART-1),NDYA,4)
        ENDDO
        XPOL(NDYA+1,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,NDYA,3)
        YPOL(NDYA+1,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,NDYA,3)
      ENDDO
C
      DO J=1,NDYA+1
        DO I=1,NPOINT(2,NPLP)
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
        END DO
      END DO
C
      ndxa=npoint(2,nplp)-1
c
C     do j=1,ndya+1
C       write (iunout,*)
C       write (iunout,*) 'in geomd polygon ',j
C       write (iunout,'(1p,6e12.4)') 
C    .        (xpol(j,i),ypol(j,i),i=1,npoint(2,nplp))
C     enddo
C
      RETURN
*//END GEOMD//
      END
C ===== SOURCE: geousr.f
C
C
C
      SUBROUTINE EIRENE_GEOUSR
C
C   PREPARE DATA FOR LIMITER-SURFACES
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CADGEO
      USE EIRMOD_COMPRT
      USE EIRMOD_CTRCEI
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      USE EIRMOD_CGRID
      USE EIRMOD_CLGIN
      USE EIRMOD_CINIT
      USE EIRMOD_CPOLYG
      IMPLICIT NONE
      CHARACTER(80) :: ZEILE
      REAL(DP) :: XCOOR, YCOOR, ZCOOR, xan, xen, yan, yen, zan, zen
      INTEGER :: NADMOD, NASMOD, NORMOD, NRS, IPUNKT, I,NSSIR, NSSIP,
     .           IDIR, IR, IP, IT, IC, IN, NAS
      logical :: lx1,lx2,lx3,lx4,ly1,ly2,ly3,ly4
C
C MODIFY GEOMETRY
C
C
cswx 24sep07
      return
cswx
      READ (IUNIN,'(A80)') ZEILE
      READ (IUNIN,'(3I6)') NADMOD,NASMOD,NORMOD

      DO I=1,NADMOD
        READ (IUNIN,'(2I6,3E12.4)') NRS,IPUNKT,XCOOR,YCOOR,ZCOOR

        SELECT CASE(IPUNKT)

        CASE DEFAULT
           WRITE (iunout,*) 'WRONG POINTNUMBER IN ADDUSR '
           WRITE (iunout,*) 'INPUT LINE READING'
           WRITE (iunout,'(2I6,1P,3E12.4)') NRS,IPUNKT,XCOOR,YCOOR,ZCOOR
           WRITE (iunout,*) ' IS IGNORED '

        CASE (1)
           P1(1,NRS)=XCOOR
           P1(2,NRS)=YCOOR
           P1(3,NRS)=ZCOOR

        CASE (2)
           P2(1,NRS)=XCOOR
           P2(2,NRS)=YCOOR
           P2(3,NRS)=ZCOOR

        CASE (3)
           P3(1,NRS)=XCOOR
           P3(2,NRS)=YCOOR
           P3(3,NRS)=ZCOOR

        CASE (4)
           P4(1,NRS)=XCOOR
           P4(2,NRS)=YCOOR
           P4(3,NRS)=ZCOOR

        CASE (5)
           P5(1,NRS)=XCOOR
           P5(2,NRS)=YCOOR
           P5(3,NRS)=ZCOOR

        CASE (6)
           P6(1,NRS)=XCOOR
           P6(2,NRS)=YCOOR
           P6(3,NRS)=ZCOOR

        END SELECT
      ENDDO

      DO I=1,NASMOD
        READ (IUNIN,'(5I6)') NAS,IPUNKT,NSSIR,NSSIP
        IF (IPUNKT.EQ.1) THEN
          P1(1,NAS)=XPOL(NSSIR,NSSIP)
          P1(2,NAS)=YPOL(NSSIR,NSSIP)
        ELSEIF (IPUNKT.EQ.2) THEN
          P2(1,NAS)=XPOL(NSSIR,NSSIP)
          P2(2,NAS)=YPOL(NSSIR,NSSIP)
        ELSE
          WRITE (iunout,*) 'WRONG POINTNUMBER IN ADDUSR '
          WRITE (iunout,*) 'INPUT LINE READING'
          WRITE (iunout,'(5I6)') NAS,IPUNKT,NSSIR,NSSIP
          WRITE (iunout,*) ' IS IGNORED '
        ENDIF
      ENDDO

      DO I = 1, NORMOD
        READ (IUNIN,'(5I6)') IDIR,IR,IP
        IF (IDIR == 1) THEN
          PLNX(IR,IP) = -PLNX(IR,IP)
          PLNY(IR,IP) = -PLNY(IR,IP)
        ELSE IF (IDIR == 2) THEN
          PPLNX(IR,IP) = -PPLNX(IR,IP)
          PPLNY(IR,IP) = -PPLNY(IR,IP)
        ELSE
          WRITE (iunout,*) ' IDIR =',IDIR,' NOT FORESEEN IN GEOUSR '
          WRITE (iunout,*) IDIR,IR,IP
          WRITE (iunout,*) ' IS IGNORED '
        END IF
      END DO
C
C
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
C
C  SET SOME VOLUMES EXPLIZIT
C
C
C  MODIFY REFLECTION MODEL AT TARGET PLATES
C
C     do i=1,nlimps
C       do isp=1,natmi+nmoli+nioni
C         recyct(isp,i)=1.
C       enddo
C     enddo

      RETURN
      END
C ===== SOURCE: iniusr.f


      SUBROUTINE EIRENE_iniUSR
      IMPLICIT NONE
      RETURN
      END
C ===== SOURCE: leausr.f


      FUNCTION EIRENE_LEAUSR(A,B,C)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A, B, C
      INTEGER :: EIRENE_LEAUSR
      EIRENE_LEAUSR=1
      RETURN
      END
C ===== SOURCE: locstr_usr.f
csw
csw routine to find a specific string in unit fp
csw mar2006
csw s.wiesen@fz-juelich.de
csw
      subroutine EIRENE_locstr_usr(fp,sstr,ier)
      implicit none
      integer, intent(in) :: fp
      integer, intent(out) :: ier
      character(*), intent(in) :: sstr
      character(256) :: line
      integer :: ierror

      ier=1
      rewind(fp)
      do
         read(fp,'(a)',iostat=ierror) line
         if(ierror /= 0) return

         if( index(line,sstr) == 1) then
            ier = 0
            return
         endif
      enddo
      end

C ===== SOURCE: modusr.f
c
c
      subroutine EIRENE_modusr
      return
      end
C ===== SOURCE: mshadj.f
C
C
      SUBROUTINE EIRENE_MSHADJ (X1,Y1,X2,Y2,XPLG,YPLG,
     &     NPLP,NPOINT,M1,NDX,IR)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1(*),Y1(*),X2(*),Y2(*)
      INTEGER, INTENT(IN) ::M1, NDX, IR
      REAL(DP), INTENT(INOUT) :: XPLG(M1,*),YPLG(M1,*)
      INTEGER, INTENT(INOUT) :: NPOINT(2,*), NPLP
      REAL(DP) :: EPS, D, D1
      INTEGER :: NPRT, NP, I
      LOGICAL :: LDUMCL,LPRT
C  LPRT=.TRUE.  : VALID PART
C  LDUMCL=.TRUE.: DUMMY CELL
C
      EPS=1.E-20
      NPRT=0
      NP=0
      LDUMCL=.FALSE.
      LPRT=.FALSE.
C
      DO 1 I=1,NDX
        D=ABS((X1(I)-X2(I))**2+(Y1(I)-Y2(I))**2)
        IF (D.LE.EPS) THEN
          IF (LPRT) THEN
C  ENDE DER VALID ZELLEN
            NPOINT(2,NPRT)=NP
            LPRT=.FALSE.
            LDUMCL=.TRUE.
          ENDIF
        ELSEIF (D.GT.EPS) THEN
C  STARTPUNKT DIESES POLYGONS = ERSTER VALID PUNKT?
          IF (.NOT.LPRT.AND..NOT.LDUMCL) THEN
            LPRT=.TRUE.
            NPRT=NPRT+1
            NP=NP+1
            XPLG(IR,NP)=X1(I)
            YPLG(IR,NP)=Y1(I)
            NPOINT(1,NPRT)=NP
          ELSEIF (.NOT.LPRT.AND.LDUMCL) THEN
            D1=ABS((X1(I+1)-X2(I+1))**2+(Y1(I+1)-Y2(I+1))**2)
            IF (D1.GT.EPS) THEN
C  ENDE DER DUMMYZELLEN: SUCHE NAECHSTE ZELLE MIT D1.GT.0.
              LDUMCL=.FALSE.
              LPRT=.TRUE.
              NPRT=NPRT+1
              NPOINT(1,NPRT)=NP
            ENDIF
          ENDIF
        ENDIF
        NP=NP+1
        XPLG(IR,NP)=X2(I)
        YPLG(IR,NP)=Y2(I)
1     CONTINUE
      NPOINT(2,NPRT)=NP
      NPLP=NPRT
      RETURN
      END
C ===== SOURCE: outusr.f
C+---------------------------------------------------------------+
C| Purpose:                                                      |
C| --------                                                      |
C| Write the eirene.tranfer file to give the EIRENE results back |
C| to EDGE2D.                                                    |
C| Write also eirene.chemFluxDep file, to store the neutral      |
C| particle fluxes on to the walls for the next EIRENE call      |
C+---------------------------------------------------------------+
C| Modifications:                                                |
C| --------------                                                |
C| 16/07/2010   D.Harting    Added writing of neutral fluxes to  |
C|                           file eirene.chemFluxDep. Added also |
C|                           two variables to eirene_user        |
C|                           namelist.                           |
C| 23/11/2010   D.Harting    If EDGE2D is used with density      |
C|                           control by puff+recycling, the      |
C|                           number of puffeing surfaces and thus|
C|                           the number of additional surfaces   |
C|                           (NLIM) may vary. Before, the neutral|
C|                           flux file eirene.chemFluxDep from   |
C|                           previous run was checked to have the|
C|                           same number of add. surfaces as the |
C|                           actual run. This forced a stopping  |
C|                           of the code. Now the actual triangle|
C|                           number and its side is checked, to  |
C|                           asure that the right neutral flux   |
C|                           is used.                            |
C| 24/11/2010   D.Harting    Added for backward compatibility    |
C|                           a version number to the neutral flux|
C|                           file eirene.chemFluxDep. If the     |
C|                           Version number in the file and in   |
C|                           code are not matching, the file is  |
C|                           not read and zero neutral flux is   |
C|                           assumed in the actual run. At the   |
C|                           end of the eirene run, a new neutral|
C|                           flux file with the current version  |
C|                           number is generated.                |
C+---------------------------------------------------------------+
      subroutine EIRENE_outusr
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comusr
      use eirmod_ccona
      use eirmod_cestim
      use eirmod_ctrcei
      use eirmod_comxs
      use eirmod_cgeom
      use eirmod_cinit
      use eirmod_ctrig
      use eirmod_csdvi
      use eirmod_clogau
      use eirmod_comsou
      use eirmod_comprt
      use eirmod_czt1
      use eirmod_coutau
      implicit none
      integer :: fp,mtri,ir,is

      real(dp) :: pa,pi,pm,pph,vvol,val,sumvol,sumval,vdenpara,bb
      integer :: j,i
      real(dp) :: c1,c2,c3,c4

csw ratecoeff.dat
      integer :: iaei, iacx, imei, imcx, iiei, iicx, iirc
      integer :: irei, ircx, irrc
      real(dp) :: de,di
      integer :: idsc,ipl,kk,nrc,ireac,k,np,nr,msg
csw 25oct07
      real(dp) :: x1,x2,y1,y2,ar,xc
      real(dp),allocatable :: sumpotpl(:)
csw
      logical :: lcxsigma

      integer, external :: EIRENE_idez

      logical, save :: ldebug

      integer, save :: eirene_nbirth,eirene_njetto
      character(len=256), save :: eirene_fbirth,eirene_ftransfer,
     &     eirene_fstoreneutflux
      real(dp) :: eirene_phi_offsets(9)
      integer :: eirene_wallFluxModel ! calculation of wall fluxes for chemical sputtering
c                = 0: no wall fluxes are used (old edge2d model)
c                = 1: only ion fluxes are used
c                = 2: ion fluxes and neutral fluxes from last eirene iteration are used
c                = 3: ion and neutral fluxes are used, and EIRENE is iterated to give 
c                     converged neutral fluxes.
      logical :: eirene_use_elstepdat_bug
      logical  :: lfound
      real(dp) :: neutralFluxFileVersion
      namelist /eirene_user/eirene_nbirth,eirene_njetto,
     .                      eirene_fbirth,eirene_ftransfer,
     .                      eirene_phi_offsets,
     .                      eirene_fstoreneutflux,
     .                      eirene_wallFluxModel,
     .                      eirene_use_elstepdat_bug

      NeutralFluxFileVersion = 1.0
      ldebug=.false.
csw
      eirene_ftransfer = 'eirene.transfer'
      eirene_njetto=0
cdmh
      eirene_fstoreneutflux = 'eirene.chemFluxDep'
      eirene_wallFluxModel = 1
cdmh      
      open(unit=9998,file='eirene_user.namelist')
      read(9998,eirene_user)
      close(9998)
csw

csw 25oct07
      allocate(sumpotpl(npls))
csw
      do k=1,npls
         sumpotpl(k) = 0.
         do np=1,3
            do nr=1,ntrii
                if(  inmti(np,nr) == 2+1 .or.
     .               inmti(np,nr) == 3+1) then
                   msg = nlim+nsts+inspat(np,nr)
                   
                   x1 = xtrian(necke(np,nr))
                   y1 = ytrian(necke(np,nr))
                   if(np == 3) then
                      x2 = xtrian(necke(1,nr))
                      y2 = ytrian(necke(1,nr))
                   else
                      x2 = xtrian(necke(np+1,nr))
                      y2 = ytrian(necke(np+1,nr))
                   endif
                   
                   ar = sqrt( (x1-x2)**2 + (y1-y2)**2)
                   xc = (x1+x2)/2.d0
c     potpl in 1/s/cm
                   sumpotpl(k) = sumpotpl(k) +
     .                  estims(naddw(25)+k,msg)
     .                  /1.6022e-19/2.d0/pia/xc
                endif
            enddo
         enddo
      enddo
csw


      fp = 4999
      open(unit=fp,file=trim(eirene_ftransfer),access='sequential',
     .     status='replace')

      write(fp,'(a)') '* ntrii,nrad  :'
      write(fp,'(3(1x,i6))') ntrii,nrad

      write(fp,'(a)') '* nstordr:'
      write(fp,'(1x,i6)') nstordr

      write(fp,'(a)') '* natm,nmol,nion,nphot  :'
      write(fp,'(4(1x,i6))') natm,nmol,nion,nphot     

      write(fp,'(a)') '* npls:'
      write(fp,'(3(1x,i6))') npls

      write(fp,'(a)') '* nlimps,nlim,nsts:'
      write(fp,'(3(1x,i6))') nlimps,nlim,nsts

c---------------------------------------
      write(fp,'(a,i6)') '* BULK SPECIES NPLS = ',npls      
      do ipls=1,npls
         ityp=4
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IPLS = ',ipls
         do ir=1,ntrii
            c1 = mapl(ipls,ir)
            c2 = mmpl(ipls,ir)
            c3 = mipl(ipls,ir)
c            c4 = mphpl(ipls,ir)
            c4 = 0.
            write(fp,'(i6,30(1x,e14.6))') ir,
     .           papl(ipls,ir),
     .           pmpl(ipls,ir),
     .           pipl(ipls,ir),
c     .           pphpl(ipls,ir),
     .           0.,
c pspl:
     .           papl(ipls,ir)+pmpl(ipls,ir)+pipl(ipls,ir)
c     .                        +pphpl(ipls,ir),
     .                        +0.,
c mapl:
     .           c1,
c mmpl:
     .           c2,
c mipl:
     .           c3,
c mphpl:
     .           c4,
c mspl:
     .           c1+c2+c3+c4,

c not used anymore... (changed background profiles in extra file?)
     .           diin(ipls,ir),
     .           vxin(ipls,ir),
     .           vyin(ipls,ir),
     .           vzin(ipls,ir),
     .           bvin(ipls,ir),
     .           tiin(ipls,ir),
     .           edrift(ipls,ir),

c     .           eapl(ir),empl(ir),eipl(ir),ephpl(ir),
     .           eapl(ir),empl(ir),eipl(ir),0.,
c espl:
c     .           eapl(ir)+empl(ir)+eipl(ir)+ephpl(ir),
     .           eapl(ir)+empl(ir)+eipl(ir)+0.,

c     .           eael(ir),emel(ir),eiel(ir),ephel(ir),
     .           eael(ir),emel(ir),eiel(ir),0.,
c esel:
c     .           eael(ir)+emel(ir)+eiel(ir)+ephel(ir)
     .           eael(ir)+emel(ir)+eiel(ir)+0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IPLS = ',ipls
         do is=1,nlimps
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potpl(ipls,is),
     .           eotpl(ipls,is),
     .           sptpl(ipls,is),
     .           spump(ispz,is)
         enddo
         write(fp,'(20(1x,e14.6))') sumpotpl(ipls),srec(ipls,0)
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* ATOMIC SPECIES, NATM = ',natm
      do iatm=1,natm
         ityp=1
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IATM = ',iatm
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdena(iatm,ir)*bxin(ir) + vydena(iatm,ir)
     .              *byin(ir) + vzdena(iatm,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir, 
     .           pdena(iatm,ir),            
     .           vxdena(iatm,ir),
     .           vydena(iatm,ir),
     .           vzdena(iatm,ir),
     .           vdenpara,
     .           edena(iatm,ir),
     .           edena(iatm,ir)/(pdena(iatm,ir)+1.d-60),
     .           0.,
     .           sigma(1,ir),
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IATM = ',iatm
         do is=1,nlimps
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potat(iatm,is),
     .           prfaat(iatm,is),
     .           prfmat(iatm,is),
     .           prfiat(iatm,is),
c     .           prfphat(iatm,is),
     .           0.,
     .           prfpat(iatm,is),
     .           eotat(iatm,is),
     .           erfaat(iatm,is),
     .           erfmat(iatm,is),
     .           erfiat(iatm,is),
c     .           erfphat(iatm,is),
     .           0.,
     .           erfpat(iatm,is),
     .           sptat(iatm,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* MOLECULAR SPECIES, NMOL = ',nmol
      do imol=1,nmol
         ityp=2
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IMOL = ',imol
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdenm(imol,ir)*bxin(ir) + vydenm(imol,ir)
     .              *byin(ir) + vzdenm(imol,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir, 
     .           pdenm(imol,ir),
     .           vxdenm(imol,ir),
     .           vydenm(imol,ir),
     .           vzdenm(imol,ir),
     .           vdenpara,
     .           edenm(imol,ir),
     .           edenm(imol,ir)/(pdenm(imol,ir)+1.d-60),
     .           0.,
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IMOL = ',imol
         do is=1,nlimps
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potml(imol,is),
     .           prfaml(imol,is),
     .           prfmml(imol,is),
     .           prfiml(imol,is),
c     .           prfphml(imol,is),
     .           0.,
     .           prfpml(imol,is),
     .           eotml(imol,is),
     .           erfaml(imol,is),
     .           erfmml(imol,is),
     .           erfiml(imol,is),
c     .           erfphml(imol,is),
     .           0.,
     .           erfpml(imol,is),
     .           sptml(imol,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* IONIC SPECIES, NION = ',nion
      do iion=1,nion
         ityp=3
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IION = ',iion
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdeni(iion,ir)*bxin(ir) + vydeni(iion,ir)
     .              *byin(ir) + vzdeni(iion,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir, 
     .           pdeni(iion,ir),
     .           vxdeni(iion,ir),
     .           vydeni(iion,ir),
     .           vzdeni(iion,ir),
     .           vdenpara,
     .           edeni(iion,ir),
     .           edeni(iion,ir)/(pdeni(iion,ir)+1.d-60),
     .           0.,
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IION = ',iion
         do is=1,nlimps
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potio(iion,is),
     .           prfaio(iion,is),
     .           prfmio(iion,is),
     .           prfiio(iion,is),
c     .           prfphio(iion,is),
     .           0.,
     .           prfpio(iion,is),
     .           eotio(iion,is),
     .           erfaio(iion,is),
     .           erfmio(iion,is),
     .           erfiio(iion,is),
c     .           erfphio(iion,is),
     .           0.,
     .           erfpio(iion,is),
     .           sptio(iion,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* PHOTONIC SPECIES, NPHOT = ',nphot
      do iphot=1,nphot
         ityp=0
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IPHOT = ',iphot
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdenph(iphot,ir)*bxin(ir) + vydenph(iphot,ir)
     .              *byin(ir) + vzdenph(iphot,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir, 
     .           pdenph(iphot,ir),
     .           vxdenph(iphot,ir),
     .           vydenph(iphot,ir),
     .           vzdenph(iphot,ir),
     .           vdenpara,
     .           edenph(iphot,ir),
     .           edenph(iphot,ir)/(pdenph(iphot,ir)+1.d-60),
     .           0.,
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IPHOT = ',iphot
         do is=1,nlimps
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potpht(iphot,is),
     .           prfapht(iphot,is),
     .           prfmpht(iphot,is),
     .           prfipht(iphot,is),
     .           prfphpht(iphot,is),
     .           prfppht(iphot,is),
     .           eotpht(iphot,is),
     .           erfapht(iphot,is),
     .           erfmpht(iphot,is),
     .           erfipht(iphot,is),
     .           erfphpht(iphot,is),
     .           erfppht(iphot,is),
     .           sptpht(iphot,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* MISC DATA'
      do ir=1,ntrii
         write(fp,'(i6,20(1x,e14.6))') ir,
     .        dble(ncltal(ir)),0.,0.,
     .        vol(ir),voltal(ir),
     .        dein(ir),tein(ir),bxin(ir),byin(ir),bzin(ir),bfin(ir)
      enddo

csw 20dec07-----------------------------
      write(fp,'(a)') '* STRATUM DATA'
      write(fp,'(i6)') nstra
      do istra=1,nstra
        do iatm=1,natm
          write(fp,'(2i6,20(1x,e14.6))') istra,iatm,
     .              wtota(iatm,istra),0.
        enddo
        do imol=1,nmol
          write(fp,'(2i6,20(1x,e14.6))') istra,imol,
     .              wtotm(imol,istra),0.
        enddo
        do iion=1,nion
          write(fp,'(2i6,20(1x,e14.6))') istra,iion,
     .              wtoti(iion,istra),0.
        enddo
        do iphot=1,nphot
          write(fp,'(2i6,20(1x,e14.6))') istra,iphot,
     .              wtotph(iphot,istra),0.
        enddo
        do ipls=1,npls
          write(fp,'(2i6,20(1x,e14.6))') istra,ipls,
     .              wtotp(ipls,istra),0.
        enddo
      enddo

csw 24jan08
c---------------------------------------
      if(eirene_njetto .gt. 0) then
        write(fp,'(a)') '* ADDITIONAL TALLY DATA FOR JETTO NEUTRALS'
        write(fp,'(6i6)') nadv
        do ir=1,ntrii
          write(fp,'(i6,20(1x,e14.6))') ir,( addv(i,ir), i=1,nadv ),
     .                                     ( ppat(i,ir), i=1,natm )
        enddo
      endif

      if(eirene_nbirth .gt. 0) then
        write(fp,'(a)') '* ADDITIONAL TALLY DATA FOR JETTO NBI'
        write(fp,'(6i6)') nstra
        do is=1,nstra
          write(fp,'(i6,20(1x,e14.6))') is, etota(is)
        enddo
      endif
 
      close(fp)
csw 25oct07
      deallocate(sumpotpl)
csw

cdmh 15jun10
c----------------------------------------
c     store neutral particle fluxes [A] on wall
      fp = 4999
      open(unit=fp,file=trim(eirene_fstoreneutflux),access='sequential',
     .     status='replace')
      write(fp,'(a28,f14.6)') "* Neutral flux file version:",
     &     NeutralFluxFileVersion
      write(fp,'(a,a)') '*  NLIM,   NSTS,  NGITT, NGSTAL,',
     &     '   NATM, NLMPGS,  NTRII'
      write(fp,'(7i8)') NLIM, NSTS, NGITT, NGSTAL, NATM, NLMPGS, NTRII
      do iatm=1,NATM
         write(fp,'(a,i0)') '* neutral fluxes from atom species ',iatm
         write(fp,'(a)') "*  idx,  neutral_flux, ITRIA, ISIDE, ISURF"
         do is=1,NLMPGS
c           find corresponding triangle
            lfound = .false.
            nr = 0
            np = 0
            do i=1,ntrii
               do j=1,3
                  if ((INSPAT(j,i).eq. is -(NLIM+NSTS))
     &                 .and.(INSPAT(j,i).ne.0) )then
                     if (lfound) then                        
                        write(iunout,*)"* EIRENE_OUTUSR:"
                        write(iunout,*)"* Edge twice found"
                        call EIRENE_exit_own(1)
                     endif
                     lfound=.true.
                     nr = i
                     np = j
                  endif
               enddo
            enddo                 
            if ((nr.le.0).or.(np.le.0)) then
               write(fp,'(i6,1x,e14.6,1x,i6,1x,i6,1x,i6)') 
     &              is,POTAT(iatm,is),nr,np,0
            else
               write(fp,'(i6,1x,e14.6,1x,i6,1x,i6,1x,i6)') 
     &              is,POTAT(iatm,is),nr,np,INMTI(np,nr)
            endif
         enddo                  !is
      enddo                     !iatm
      close(fp)

      end
 
C ===== SOURCE: plausr.f
c
c modified: s.wiesen@fz-juelich.de
c
C+---------------------------------------------------------------+
C| Purpose:                                                      |
C| --------                                                      |
C| reads target plasma information from casename.zplasma         |
C| sets step-functions                                           |
C| reads also neutral particle fluxes from 'eirene.chemFluxDep'  |
C| and sets FLXOUT for fluxdependency of chemical sputtering.    |
C+---------------------------------------------------------------+
C| Modifications:                                                |
C| --------------                                                |
C| 16/07/2010   D.Harting    Added reading of neutral fluxes from|
C|                           file eirene.chemFluxDep. Added also |
C|                           two variables to eirene_user        |
C|                           namelist.                           |
C| 23/11/2010   D.Harting    If EDGE2D is used with density      |
C|                           control by puff+recycling, the      |
C|                           number of puffeing surfaces and thus|
C|                           the number of additional surfaces   |
C|                           (NLIM) may vary. Before, the neutral|
C|                           flux file eirene.chemFluxDep from   |
C|                           previous run was checked to have the|
C|                           same number of add. surfaces as the |
C|                           actual run. This forced a stopping  |
C|                           of the code. Now the actual triangle|
C|                           number and its side is checked, to  |
C|                           asure that the right neutral flux   |
C|                           is used.                            |
C| 24/11/2010   D.Harting    Added for backward compatibility    |
C|                           a version number to the neutral flux|
C|                           file eirene.chemFluxDep. If the     |
C|                           Version number in the file and in   |
C|                           code are not matching, the file is  |
C|                           not read and zero neutral flux is   |
C|                           assumed in the actual run. At the   |
C|                           end of the eirene run, a new neutral|
C|                           flux file with the current version  |
C|                           number is generated.                |
C| 22/03/2011   D.Harting    Do not exit anymore if one of the   |
C|                           geometrical parameters in           |
C|                           eirene.chemFluxDep are not matching.|
C|                           Just ignore the neutral flux from   |
C|                           the previous run and continue.      |
C+---------------------------------------------------------------+
      SUBROUTINE EIRENE_PLAUSR
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CSTEP
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CINIT
      USE EIRMOD_COMSOU
      USE EIRMOD_CTRIG
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      IMPLICIT NONE
      REAL(DP) :: FACTOR, EIRENE_STEP
      INTEGER :: NLINES, ITRI, ISIDE, I, NBIN, ISTRA, ISRFS, ISOR, JJJ,
     .           ITEC1, ITEC2, ITEC3, ISTEP, INDSRF, IS1, IERROR, IPLS
      INTEGER :: EIRENE_IDEZ
      INTEGER, ALLOCATABLE :: KSTEP(:), INOSRC(:), IPLAN(:), IPLEN(:)
      REAL(DP) :: FLX, TE, TI, DE, MC, FE, FI, FSH, VP, FEL, DUM
      REAL(DP) :: DELR, FL, MCC, FFEL, CS, vx,vy,vz,di,usrval
      real(dp) :: xref, yref, bzref, facbz, x, y, rad, bx, by, bz, bf,
     .            errbx, errby, errbz, errbf, bxmax, bymax, bzmax, 
     .            bfmax, dfdx, dfdy, dfdz, xref2, yref2, bzref2, facbz2,
     .            xref3, yref3, bzref3, facbz3
      integer :: nref, icell, EIRENE_learc1, nplcll, ipolg, nref2, nref3
      CHARACTER(256) :: line,sstr,filename
      character(2) :: cstr2
      integer, allocatable :: indextmp(:),itritmp(:),isidetmp(:)
      integer, allocatable :: isegtmp(:)
      real(dp), allocatable :: tetmp(:),fetmp(:),fshtmp(:)
      REAL(DP), ALLOCATABLE :: PSI(:), PSI_CORNER(:), copy(:)

      integer, parameter :: fp=31
      integer :: ll,ier,j,iseg,ind, fp2

      integer, save :: eirene_nbirth,eirene_njetto
      character(len=256), save :: eirene_fbirth,eirene_ftransfer,
     &     eirene_fstoreneutflux
      real(dp) :: eirene_phi_offsets(9)

      integer :: ntr, NLIM_tmp, NSTS_tmp, NGITT_tmp, NGSTAL_tmp, 
     &     NATM_tmp, IPLS_tmp, NTR_tmp, NLMPGS_tmp
      REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: hydIonFLX, EIRENE_wall_area
      REAL(DP),ALLOCATABLE,DIMENSION(:)   :: hydNeutFLX
      INTEGER, ALLOCATABLE,DIMENSION(:,:) :: hydNeutFLX_info
      logical :: lex,dbg_out
      real(dp):: leng, twopi, tmp
      integer :: eirene_wallFluxModel ! calculation of wall fluxes for chemical sputtering
c                = 0: no wall fluxes are used (old edge2d model)
c                = 1: only ion fluxes are used
c                = 2: ion fluxes and neutral fluxes from last eirene iteration are used
c                = 3: ion and neutral fluxes are used, and EIRENE is iterated to give 
c                     converged neutral fluxes.
      logical :: eirene_use_elstepdat_bug
      real(dp) :: neutralFluxFileVersion

      namelist /eirene_user/eirene_nbirth,eirene_njetto,
     .                      eirene_fbirth,eirene_ftransfer,
     .                      eirene_phi_offsets,
     .                      eirene_fstoreneutflux,
     .                      eirene_wallFluxModel,
     .                      eirene_use_elstepdat_bug



      interface
        subroutine EIRENE_cell_to_corner (f, fcorner)
          use eirmod_precision
          implicit none
          real(dp), intent(in) :: f(:)
          real(dp), intent(out) :: fcorner(:)
        end subroutine EIRENE_cell_to_corner

        subroutine EIRENE_df_dxyz (fcorner, icell, x, y, z, 
     .                      dfdx, dfdy, dfdz) 
          use eirmod_precision
          implicit none
          real(dp), intent(in) :: fcorner(:), x, y, z
          real(dp), intent(out) :: dfdx, dfdy, dfdz
          integer, intent(in) :: icell
        end subroutine EIRENE_df_dxyz 
      end interface

      CALL EIRENE_ALLOC_CSTEP
      ALLOCATE (KSTEP(NSTEP))
      ALLOCATE (INOSRC(NSTEP))
      ALLOCATE (IPLAN(NSTEP))
      ALLOCATE (IPLEN(NSTEP))
      KSTEP = 0
      INOSRC = 0
      IPLAN = 0
      IPLEN = 0

      dbg_out=.true.

c     set default namelist values
      eirene_ftransfer = 'eirene.transfer'
      eirene_njetto=0
cdmh
      eirene_fstoreneutflux = 'eirene.chemFluxDep'
      eirene_wallFluxModel = 1
      NeutralFluxFileVersion = 1.0
cdmh      

c     get namelist config
      open(unit=9998,file='eirene_user.namelist')
      read(9998,eirene_user)
      close(9998)

c begin dmh added 21.06.2010
      ll=len_trim(casename)
      filename=casename(1:ll) // '.zplasma'
      open (unit=fp+ifoff,file=filename,access='sequential',
     .     form='formatted')

c first get number of triangles from misc plasma data:
        write(sstr,'(a20)') 
     .          '*** MISC PLASMA DATA'
        CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
        if(ier /=0) then
           write(*,*) 'PLAUSR: ',sstr,' not found'
           close(fp+ifoff)
           call EIRENE_exit_own(1)
        endif
        do j=1,2
           read(fp+ifoff,'(a)') line
        enddo

        read(fp+ifoff,'(i7)') ntr
        if (ntr /= nr1st-1) then
           write (*,*) 'PLAUSR:', sstr
           write (*,*) ' wrong number of triangles in plasma file'
           write (*,*) ' check for correct number in file ',filename
           call EIRENE_exit_own(1)
        endif
c allocate temporary array to store plasma flux (hydrogen isotope)
        allocate(hydIonFLX(3,ntr))
        hydIonFLX(:,:) = 0.0
c allocate temporary array to store wall area of elements
        allocate(EIRENE_wall_area(3,ntr))
        EIRENE_wall_area(:,:) = 0.0
c end dmh added 21.06.2010

c firstly, read misc target data
      write(sstr,'(a20)') '*** MISC TARGET DATA'
      CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
      if(ier /=0) then
         write(*,*) 'PLAUSR: ',sstr,' not found'
         close(fp+ifoff)
         CALL EIRENE_exit_own(1)
      endif
      do j=1,2
         read(fp+ifoff,'(a)') line
      enddo
      READ (fp+ifoff,*) NLINES
cswx 24sep07
      if(nlines == 0) then
        close(fp+ifoff)
        return
      endif
cswx
      if(nlines <= 0) then
         write(*,*) ' PLAUSR: error, nlines=',nlines
         close(fp+ifoff)
         call EIRENE_exit_own(1)
      endif

      allocate(indextmp(nlines))
      allocate(itritmp(nlines))
      allocate(isidetmp(nlines))
      allocate(tetmp(nlines))
      allocate(fetmp(nlines))
      allocate(fshtmp(nlines))
      allocate(isegtmp(nlines))
      do j=1,nlines
         read(fp+ifoff,'(3i7,3(1x,e14.7))') 
     .        indextmp(j),
     .        itritmp(j),
     .        isidetmp(j),
     .        tetmp(j),
     .        fetmp(j),
     .        fshtmp(j)
         isegtmp(indextmp(j)) = j
      enddo


c now read species dependent target data
      DO ISTRA=1,NSTRAI
         IF (.NOT.NLSRF(ISTRA)) CYCLE
         DO ISRFS=1,NSRFSI(ISTRA)

c     get step function index ISTEP
            ISOR=SORLIM(ISRFS,ISTRA)
            ITEC1=EIRENE_IDEZ(ISOR,1,4)
            ITEC2=EIRENE_IDEZ(ISOR,2,4)
            ITEC3=EIRENE_IDEZ(ISOR,3,4)
            IF ((ITEC1 /= 4).AND.(ITEC2 /= 4).AND.(ITEC3 /= 4)) CYCLE
            ISTEP=SORIND(ISRFS,ISTRA)
            IF (ISTEP.EQ.0) THEN
               WRITE (*,*) 'ERROR IN PRIMARY SOURCE DATA '
               WRITE (*,*) 'STEPFUNCTION REQUESTED FOR SOURCE SURFACE '
               WRITE (*,*) 'NO. ',INSOR(ISRFS,ISTRA),' BUT SORIND.EQ.0.'
               CALL EIRENE_EXIT_own(1)
            ELSEIF (ISTEP.GT.NSTEP) THEN
               CALL EIRENE_MASPRM('NSTEP',5,NSTEP,'ISTEP',5,ISTEP,IERROR)
               CALL EIRENE_EXIT_own(1)
            ENDIF

c     get surface index INDSRF
            INDSRF=INSOR(ISRFS,ISTRA)
            IF (INDSRF < 0) INDSRF=NLIM+ABS(INDSRF)

c     get species index/indices IPLAN(ISTEP) --> IPLEN(ISTEP)
            IF (NSPEZ(ISTRA) <= 0) THEN
               IPLAN(ISTEP)=1
               IPLEN(ISTEP)=NPLSI
               ipls = 0
            ELSE
               IPLAN(ISTEP)=NSPEZ(ISTRA)
               IPLEN(ISTEP)=NSPEZ(ISTRA)
               ipls=nspez(istra)
            END IF
c     fudge species index for atomic impurity flux (get it from NEMODS index K)
            IPLS_tmp = EIRENE_IDEZ(NEMODS(ISTRA),4,4)
            IF (IPLS_tmp.gt.1) then 
               ipls=IPLS_tmp
            ENDIF
c     search target tag in .zplasma file
           IPLS_tmp = IPLS 
           write(cstr2,'(i2.2)') ipls
           write(sstr,'(a23)') 
     .          '*** ION #'//cstr2//' TARGET DATA'
           CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
           if(ier /=0) then
              write(*,*) 'PROUSR: tag ',sstr,'not found'
              close(fp+ifoff)
              call EIRENE_exit_own(1)
           endif
           do j=1,2
              read(fp+ifoff,'(a)') line
           enddo

c     read step functions
            READ (fp+ifoff,*) NLINES
            DO I=1, NLINES
               read(fp+ifoff,'(i7,11(1x,e14.7))')
     .              ind,
     .              flx,ti,di,
     .              vx,vy,vz,
     .              fi,fel,
     .              vp,mc,
     .              usrval

               iseg = isegtmp(ind)
               itri = itritmp(iseg)
               iside= isidetmp(iseg)
               te   = tetmp(iseg)
               fe   = fetmp(iseg)
               fsh  = fshtmp(iseg)

c begin added dmh 21.06.2010
c     store ion flux in [A] of main plasma
               if(ipls_tmp.eq.1) then
                  hydIonFLX(iside,itri) = flx
               endif
c end added dmh 21.06.2010

               IF (INMTI(ISIDE,ITRI) == INDSRF) THEN
                  IF (KSTEP(ISTEP) == 0) RRSTEP(ISTEP,1) = 0._DP
                  KSTEP(ISTEP) = KSTEP(ISTEP) + 1
                  INOSRC(ISTEP) = ISTRA
                  IS1 = ISIDE + 1
                  IF (IS1.GT.3) IS1=1
                  IRSTEP(ISTEP,KSTEP(ISTEP))=ITRI
                  IPSTEP(ISTEP,KSTEP(ISTEP))=ISIDE
                  ITSTEP(ISTEP,KSTEP(ISTEP))=1
                  IASTEP(ISTEP,KSTEP(ISTEP))=0
                  IBSTEP(ISTEP,KSTEP(ISTEP))=1
                  DELR =  SQRT(
     .                 (XTRIAN(NECKE(ISIDE,ITRI))
     .                 -XTRIAN(NECKE(IS1,ITRI)))**2+
     .                 (YTRIAN(NECKE(ISIDE,ITRI))
     .                 -YTRIAN(NECKE(IS1,ITRI)))**2)
                  RRSTEP(ISTEP,KSTEP(ISTEP)+1)=
     .                 RRSTEP(ISTEP,KSTEP(ISTEP)) + DELR         
                  TESTEP(ISTEP,KSTEP(ISTEP)) = TE
                  FESTEP(ISTEP,KSTEP(ISTEP)) = FE
C     IF NO SHEATH POTENTIAL SPECIFIED, DERIVE IT FROM ELECTRON ENERGY
C     FLUX BY SUBTRACTING THE KINETIC CONTRIBUTION 2.0*TE
                  IF (FSH.EQ.0..AND.FE.GE.2.0) FSH=FE-2.0
                  SHSTEP(ISTEP,KSTEP(ISTEP)) = FSH
                  DO IPLS=IPLAN(ISTEP), IPLEN(ISTEP)
                     TISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = TI ! eV
                     DISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = DI ! 1/cm**3
                     FISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = FI !   1
                     VPSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = abs(VP) ! cm/s 
C     VP OVERRULES MC, IF VP IS GIVEN and MC=0
                     MCC=0.0
                     IF (VP.NE.0.0) THEN
                        CS=SQRT((TI+TE)/RMASSP(IPLS))*CVEL2A
                        IF (CS.GT.0.) MCC=VP/CS
                     ENDIF
                     jjj=kstep(istep)
                     IF (MC.EQ.0.) MC=MCC
!pb                     MCSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = abs(MC) ! 1
C     THIS NEXT VECTOR IS V-PARALLEL, IN CARTESIAN COORDINATES 			   
                     VXSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VXIN(IPLS,ITRI)
                     VYSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VYIN(IPLS,ITRI)
                     VZSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VZIN(IPLS,ITRI)
c                     VXSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VX
c                     VYSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VY
c                     VZSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VZ
                     FLSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = ABS(FLX)/DELR
C     IF NO ION KINETIC ENERGY FLUX IS SPECIFIED, DERIVE IT FROM  FI,MC,VP
                     ffel=0.
                     IF (FI.gt.0..or.mc.gt.0.) then
                        ffel=(FI*TI+0.5*MC*MC*(TE+TI))*abs(FLX)
                     endif
                     IF (FEL.EQ.0.) FEL=FFEL
                     IF (eirene_use_elstepdat_bug) THEN
                        ELSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = FEL
                     ELSE
                        ELSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = FEL/DELR
                     ENDIF
csw
                     usrstep(ipls,istep,kstep(istep)) = usrval
csw
c                 enddo ipls,iplan
                  ENDDO

c              endif inmti
               ENDIF

c           enddo nlines
            ENDDO 

c        enddo isrfs
         ENDDO 

c     enddo istra
      ENDDO



      DO ISTEP = 1, NSTEP
         IF (KSTEP(ISTEP) > 0) THEN
            NBIN=KSTEP(ISTEP)+1
            FL=EIRENE_STEP(IPLAN(ISTEP),IPLEN(ISTEP),NBIN,ISTEP)
            FLUX(INOSRC(ISTEP))=FL
         END IF
      END DO
 
c begin added dmh 21.06.2010 for flux dependency of chemical sputtering
c     read neutral flux [A] to target and walls from last EIRENE run

      allocate(hydNeutFLX(NLMPGS))
      hydNeutFLX(:) = 0.0         
      allocate(hydNeutFLX_info(NLMPGS,3))
      hydNeutFLX_info = 0
      NLIM_tmp   = NLIM
      NSTS_tmp   = NSTS
      NGITT_tmp  = NGITT
      NGSTAL_tmp = NGSTAL
      NATM_tmp   = NATM
      NLMPGS_tmp = NLMPGS
      NTR_tmp    = NTR
      fp2 = 4999
      inquire(file=trim(eirene_fstoreneutflux),exist=lex)
      if (lex) then
          open(unit=fp2,file=trim(eirene_fstoreneutflux),
     &        access='sequential')
         read(fp2,'(a)') line
         read(line,'(a28)') sstr

c        check if netral flux file is compatible with actual code version
         if (index(sstr,"* Neutral flux file version:").ne.1) then
            WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning"
            write(IUNOUT,*) "Found obsolete neutral flux file: ",
     &           trim(eirene_fstoreneutflux)
            write(IUNOUT,*) "Ignoring neutral flux from previous run"
            WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning end"
            lex =.false.
         else
            read(line,'(a28,f14.6)') sstr,tmp
            if (tmp.ne.NeutralFluxFileVersion) then
               WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning"
               write(IUNOUT,'(a,a,a,f14.6)') 
     &              "Found obsolete neutral flux file: ",
     &              trim(eirene_fstoreneutflux),"; version:",tmp
               write(IUNOUT,'(a,f14.6)') 
     &              "but actual neutral flux file should be version:",
     &              NeutralFluxFileVersion
               write(IUNOUT,*) "Ignoring neutral flux from previous run"
               WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning end"
               lex =.false.
            else
c              neutral flux file is compatible with actual code version, so read it
               read(fp2,'(a)') line
               read(fp2,'(7i8)') NLIM_tmp, NSTS_tmp, NGITT_tmp,
     &              NGSTAL_tmp, NATM_tmp, NLMPGS_tmp, NTR_tmp
               if ((NLMPGS_tmp-NLIM_tmp-NSTS_tmp.ne.NLMPGS-NLIM-NSTS)
     &              .or.(NSTS_tmp.ne.NSTS).or.(NTR_tmp.ne.NTR).or.
     &              (NGITT_tmp.ne.NGITT).or.(NGSTAL_tmp.ne.NGSTAL)) then
                  WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning"
                  WRITE(IUNOUT,*) "Neutral flux file ",
     &                 trim(eirene_fstoreneutflux),
     &                 " does not fit to simulation!"
                  WRITE(IUNOUT,*)"NLIM = ",NLIM,"; NLIM_tmp = ",NLIM_tmp
                  WRITE(IUNOUT,*)"NSTS = ",NSTS,"; NSTS_tmp = ",NSTS_tmp
                  WRITE(IUNOUT,*)"NGITT = ",NGITT,
     &                 "; NGITT_tmp = ",NGITT_tmp
                  WRITE(IUNOUT,*)"NGSTAL = ",NGSTAL,
     &                 "; NGSTAL_tmp = ",NGSTAL_tmp
                  WRITE(IUNOUT,*)"NATM = ",NATM,"; NATM_tmp = ",NATM_tmp
                  WRITE(IUNOUT,*)"NTRII= ",NTR,"; NTRII_tmp = ",NTR_tmp
                  WRITE(IUNOUT,*)"NLMPGS_tmp-NLIM_tmp-NSTS_tmp = ",
     &                 NLMPGS_tmp-NLIM_tmp-NSTS_tmp,
     &                 "NLMPGS-NLIM-NSTS = ",NLMPGS-NLIM-NSTS
                  write(IUNOUT,*) 
     &                 "Ignoring neutral flux from previous run"
                  WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning end"
                  lex =.false.
               endif

               if (lex) then
                  deallocate(hydNeutFLX, hydNeutFLX_info)
                  allocate(hydNeutFLX(NLMPGS_tmp))
                  hydNeutFLX(:) = 0.0            
                  allocate(hydNeutFLX_info(NLMPGS_tmp,3))
                  hydNeutFLX_info = 0
                  
                  read(fp2,'(a)') line
                  read(fp2,'(a)') line
                  do i=1,NLMPGS_tmp
                     read(fp2,'(i6,1x,e14.6,1x,i6,1x,i6,1x,i6)') 
     &                    j, hydNeutFLX(i),
     &                    hydNeutFLX_info(i,1), hydNeutFLX_info(i,2),
     &                    hydNeutFLX_info(i,3)
                     if (i.ne.j) then
                        WRITE(IUNOUT,*) "* EIRENE_PLAUSR:"
                        write(IUNOUT,*) "Error reading from file: ",
     &                       trim(eirene_fstoreneutflux)
                        write(IUNOUT,*) "i = ",i,";  in file j = ",j
                        CALL EIRENE_EXIT_own(1)
                     endif
                  enddo         !i
               endif            !(lex)
            endif               !(tmp.ne.NeutralFluxFileVersion)
         endif                  !(index(sstr,"* Neutral flux file version:").ne.1)
         close(fp2)
      endif 
      
c     calculate area of wall surface elements
      twopi = 2.d0*dabs(dacos(-1.0))
      do itri=1,ntr
         do iside=1,3
            IS1 = ISIDE + 1
            IF (IS1.GT.3) IS1=1
            leng  = SQRT(
     .           (XTRIAN(NECKE(ISIDE,ITRI))
     .           -XTRIAN(NECKE(IS1,ITRI)))**2+
     .           (YTRIAN(NECKE(ISIDE,ITRI))
     .           -YTRIAN(NECKE(IS1,ITRI)))**2)
            EIRENE_wall_area(iside,itri) = twopi*leng
     &           *( XTRIAN(NECKE(ISIDE,ITRI)) 
     &           +  XTRIAN(NECKE(IS1,ITRI)) )/2.D0
         enddo                  ! iside
      enddo                     ! itri              

c     calculate target flux for chemical sputtering of Roth-formula
      FLXOUT(:) = 0.0

      if (eirene_wallFluxModel.ge.1) then
c     use hydrogen ion flux [A] to wall
         do itri=1,ntr
            do iside=1,3
               if ( (hydIonFLX(iside,itri).ne.0)
     &              .and.(INMTI(iside,itri).ne.0) )then
                  FLXOUT(NLIM+NSTS + INSPAT(iside,itri)) = 
     &                 FLXOUT(NLIM+NSTS + INSPAT(iside,itri))
     &                 + DABS(hydIonFLX(iside,itri))
               endif
            enddo               ! iside
         enddo                  ! itri
      endif

      if ((eirene_wallFluxModel.ge.2).and.lex) then
c     use hydrogen neutral flux [A] to wall
         do i=NLIM_tmp+NSTS_tmp+1,NLMPGS_tmp
            j = i-NLIM_tmp-NSTS_tmp+NLIM+NSTS            
            if (hydNeutFLX(i).ne.0) then
               if (hydNeutFLX_info(i,1).ne.0) then
                  lex=.false.
                  do itri=1,ntr
                     do iside=1,3
                        if ((INSPAT(iside,itri).eq. j -(NLIM+NSTS))
     &                       .and.(INSPAT(iside,itri).ne.0) )then
                           if (lex) then
                              WRITE(IUNOUT,*) "* EIRENE_PLAUSR:"
                              write(IUNOUT,*) "* Edge twice found"
                              CALL EIRENE_EXIT_own(1)
                           endif
                           lex=.true.
                           
                           if  ((itri .eq. hydNeutFLX_info(i,1)).or.
     &                          (iside .eq. hydNeutFLX_info(i,2)).or.
     &                          (INMTI(iside,itri)-NLIM-NSTS.eq.
     &                          hydNeutFLX_info(i,3)-NLIM_tmp-NSTS_tmp))
     &                          then
                              FLXOUT(j) = FLXOUT(j)
     &                             + DABS(hydNeutFLX(i))
                           else
                              WRITE(IUNOUT,*) "* EIRENE_PLAUSR:"
                              write(IUNOUT,*) 
     &                             "Neutral flux from previous run ",
     &                             "is not associated with the same ",
     &                             "triangle in this run"
                              write(IUNOUT,*) "itri_old = ",
     &                             hydNeutFLX_info(i,1),
     &                             "; itri_new = ",itri
                              write(IUNOUT,*) "iside_old = ",
     &                             hydNeutFLX_info(i,2),
     &                             "; iside_new = ",iside
                              write(IUNOUT,*) "isurf_old = ",
     &                             hydNeutFLX_info(i,3),
     &                             "; isurf_new = ",itri
                              CALL EIRENE_EXIT_own(1)
                           endif ! same triangle associated
                        endif   ! triangle found
                     enddo      ! iside
                  enddo         ! itri
               else
                  WRITE(IUNOUT,*) "* EIRENE_PLAUSR:"
                  write(IUNOUT,*) 
     &                 "* No triangle associated with neutral flux"
                  CALL EIRENE_EXIT_own(1)
               endif            ! hydNeutFLX_info(i).ne.0
            endif               ! hydNeutFLX(i).ne.0
         enddo                  ! i
      endif

c     Debug output of boundary
      if (dbg_out) then
      open(unit=fp2, file="eirene.chemSput_dbgOut",access='sequential')
      write(fp2,*) NLMPGS
      write(fp2,'(a,a,a)') 
     &     "* index, nsurf,           R1,           R2,",
     &     "           Z1,           Z2,       FLXOUT,",
     &     "      NeutFLX,       IonFLX,         area"
      do i=1,NLMPGS
         lex=.false.
         do itri=1,ntr
            do iside=1,3
               if ((INSPAT(iside,itri).eq. i -(NLIM+NSTS))
     &              .and.(INSPAT(iside,itri).ne.0) )then
                  if (lex) write(fp2,*)"*Edge twice found"
                  lex=.true.
                  IS1 = ISIDE + 1
                  IF (IS1.GT.3) IS1=1                 
                  write(fp2,'(I8,I7,8(x,e13.6))')
     &                 i,INMTI(iside,itri),
     &                 XTRIAN(NECKE(ISIDE,ITRI)),
     &                 XTRIAN(NECKE(IS1,ITRI)),
     &                 YTRIAN(NECKE(ISIDE,ITRI)),
     &                 YTRIAN(NECKE(IS1,ITRI)), 
     &                 FLXOUT(i)/
     &                 (1.6022D-19*EIRENE_wall_area(iside,itri)),
     &                 hydNeutFLX(i)
     &                 /(1.6022D-19*EIRENE_wall_area(iside,itri)),
     &                 hydIonFLX(iside,itri)
     &                 /(1.6022D-19*EIRENE_wall_area(iside,itri)),
     &                 EIRENE_wall_area(iside,itri)
               endif
            enddo
         enddo
         if (.not.lex) then
            write(fp2,'(I8,I7,8(x,e13.6))')
     &           i,0,0.D0,0.D0,0.D0,0.D0,FLXOUT(i),hydNeutFLX(i),0.D0
         endif
      enddo
      close(fp2)

      
      open(unit=fp2, file="eirene.chemSput_dbgOut2",access='sequential')
      write(fp2,*) ntr
      write(fp2,'(a,a)') 
     &     "*   itri, iside, nsurf,           R1,           R2,",
     &     "           Z1,           Z2"
      do itri=1,ntr
         do iside=1,3
            IS1 = ISIDE + 1
            IF (IS1.GT.3) IS1=1                 
            write(fp2,'(I8,x,I6,x,I6,4(x,e13.6))')
     &                 itri,iside,INMTI(iside,itri),
     &                 XTRIAN(NECKE(ISIDE,ITRI)),
     &                 XTRIAN(NECKE(IS1,ITRI)),
     &                 YTRIAN(NECKE(ISIDE,ITRI)),
     &                 YTRIAN(NECKE(IS1,ITRI))
         enddo
      enddo
      close(fp2)
      endif

c     FLXOUT is needed in #/(cm^2 s) (covert from A to #/(cm^2 s))
      do itri=1,ntr
         do iside=1,3
            if (INMTI(iside,itri).ne.0) then
               FLXOUT(NLIM+NSTS + INSPAT(iside,itri)) = 
     &              FLXOUT(NLIM+NSTS + INSPAT(iside,itri))
     &              /(1.6022D-19*EIRENE_wall_area(iside,itri))
            endif
         enddo                  ! iside
      enddo                     ! itri


c end added dmh 21.06.2010 for flux dependency of chemical sputtering

c dmh begin added output of step function data
      if (dbg_out) then
      open(unit=fp2, file="eirene.stepdat_dbg",access='sequential')
      write(fp2,"(A,A,A)") "ipls istep i  IRSTEP(ITRI)  IPSTEP(ISIDE)",
     &     "  X1   Y1  DR  FLSTEP  ELSTEP  SHSTEP  TISTEP",
     &     "  TESTEP"
      do ipls=1,npls
         do istep=1,nstep
            do i=1,KSTEP(ISTEP)
               ITRI = IRSTEP(ISTEP,I)
               ISIDE= IPSTEP(ISTEP,I)
               IS1 = ISIDE+1
               if(IS1.GT.3) IS1=1
               write(fp2,"(2(I3,x),I4,x,I6,x,I1,8(x,e13.6))")
     &              ipls, istep, i, ITRI, ISIDE, 
     &              XTRIAN(NECKE(ISIDE,ITRI)),YTRIAN(NECKE(ISIDE,ITRI)),
     &              (RRSTEP(istep,i+1)-RRSTEP(istep,i)),
     &              FLSTEP(IPLS,ISTEP,i), ELSTEP(IPLS,ISTEP,i), 
     &              SHSTEP(ISTEP,i), TISTEP(IPLS,ISTEP,i), 
     &              TESTEP(ISTEP,i)
               write(fp2,"(2(I3,x),I4,x,I6,x,I1,8(x,e13.6))")
     &              ipls, istep, i, ITRI, ISIDE, 
     &              XTRIAN(NECKE(IS1,ITRI)), YTRIAN(NECKE(IS1,ITRI)),
     &              (RRSTEP(istep,i+1)-RRSTEP(istep,i)),
     &              FLSTEP(IPLS,ISTEP,i), ELSTEP(IPLS,ISTEP,i), 
     &              SHSTEP(ISTEP,i), TISTEP(IPLS,ISTEP,i), 
     &              TESTEP(ISTEP,i)
               write(fp2,*)
            enddo
         enddo
      enddo
      close(fp2)
      endif
c dmh begin added output of step function data


c     cleanup
      DEALLOCATE (KSTEP)
      DEALLOCATE (INOSRC)
      DEALLOCATE (IPLAN)
      DEALLOCATE (IPLEN)
      deallocate(indextmp)
      deallocate(itritmp)
      deallocate(isidetmp)
      deallocate(tetmp)
      deallocate(fetmp)
      deallocate(fshtmp)
      deallocate(isegtmp)
      deallocate(hydIonFLX)
      deallocate(hydNeutFLX)
      deallocate(hydNeutFLX_info)
      deallocate(EIRENE_wall_area)

      close(fp+ifoff)


      return
!pb      return
 4711 continue

      if (nain < 1) then
         write (iunout,*) ' no additional input tally available '
         write (iunout,*) ' for PSI-function'
         write (iunout,*) ' calculation of B field from PSI function'
         write (iunout,*) ' abandonned '
         return
      end if

      ALLOCATE (PSI(NRAD))
      ALLOCATE (PSI_CORNER(NRAD))
      ALLOCATE (COPY(NRAD))
      PSI=0.
      PSI_CORNER=0.
      COPY=0._dp

      if (naini >= 8) then
        call EIRENE_prousr(copy,1+5*npls,0._dp,0._dp,
     .        0._dp,0._dp,0._dp,0._dp,0._dp,nsbox)
        adin(7,1:nsbox) = copy(1:nsbox)

        call EIRENE_prousr(copy,2+5*npls,0._dp,0._dp,
     .       0._dp,0._dp,0._dp,0._dp,0._dp,nsbox)
        adin(8,1:nsbox) = copy(1:nsbox)
      end if 

      call EIRENE_prousr(psi,5+5*npls,0._dp,0._dp,0._dp,
     .     0._dp,0._dp,0._dp,0._dp,nsbox)

      psi(1:nsbox) = psi(1:nsbox) * 1.e4_dp

      adin(1,1:nsbox) = psi(1:nsbox)

      call EIRENE_cell_to_corner (psi, psi_corner)

      xref = 200._dp
      yref = 0._dp
      NREF=EIRENE_LEARC1(XREF,YREF,0._DP,IPOLG,1,NR1STM,
     .     .FALSE.,.FALSE.,1,'PLAUSR      ')
      BZREF=BZIN(NREF)*BFIN(NREF)
      FACBZ=xcom(nref)*BZREF

      write (iunout,*) ' xref,  yref,  bzref,  facbz ',
     .                   xref,  yref,  bzref,  facbz

      xref2= 380._dp
      yref2 = 0._dp
      NREF2=EIRENE_LEARC1(XREF2,YREF2,0._DP,IPOLG,1,NR1STM,
     .     .FALSE.,.FALSE.,1,'PLAUSR      ')
      BZREF2=BZIN(NREF2)*BFIN(NREF2)
      FACBZ2=xcom(nref2)*BZREF2

      write (iunout,*) ' xref2, yref2, bzref2, facbz2 ',
     .                   xref2, yref2, bzref2, facbz2

      xref3= 250._dp
      yref3 = 150._dp
      NREF3=EIRENE_LEARC1(XREF3,YREF3,0._DP,IPOLG,1,NR1STM,
     .     .FALSE.,.FALSE.,1,'PLAUSR      ')
      BZREF3=BZIN(NREF3)*BFIN(NREF3)
      FACBZ3=xcom(nref3)*BZREF3

      write (iunout,*) ' xref3, yref3, bzref3, facbz3 ',
     .                   xref3, yref3, bzref3, facbz3

      errbx = 0._dp
      errby = 0._dp
      errbz = 0._dp
      errbf = 0._dp
      bxmax = 0._dp
      bymax = 0._dp
      bzmax = 0._dp
      bfmax = 0._dp
      nplcll = 0

      do icell = 1, ntrii
        x = xcom(icell)
        y = ycom(icell)
        rad = x
        call EIRENE_df_dxyz (psi_corner, icell, x, y, 
     &       0._dp, dfdx, dfdy, dfdz) 

        bx = -dfdy / rad
        by =  dfdx / rad
        bz = facbz / rad
        bf = sqrt(bx*bx + by*by + bz*bz)

        if (naini >= 6) then
        adin(3,icell) = bx
        adin(4,icell) = by
        adin(5,icell) = bz
        adin(6,icell) = bf

        else

        bx = bx / bf
        by = by / bf
        bz = bz / bf
        if (ixtri(icell)+iytri(icell) == 0) then
! outside plasma region: set b-field
          bxin(icell) = bx
          byin(icell) = by
          bzin(icell) = bz
          bfin(icell) = bf
        else
          nplcll = nplcll + 1
          bxmax = max(bxmax, abs(bxin(icell)))
          bymax = max(bymax, abs(byin(icell)))
          bzmax = max(bzmax, abs(bzin(icell)))
          bfmax = max(bfmax, abs(bfin(icell)))

          errbx = errbx + abs(bx-bxin(icell))
          errby = errby + abs(by-byin(icell))
          errbz = errbz + abs(bz-bzin(icell))
          errbf = errbf + abs(bf-bfin(icell))
          if (icell == nref) then
            write (iunout,*) ' reference cell ',nref
            write (iunout,*) ' bxin, byin, bzin, bfin ',
     .             bxin(icell), byin(icell), bzin(icell), bfin(icell)
            write (iunout,*) ' bx,   by,   bz,   bf   ',
     .             bx, by, bz, bf
          elseif (icell == nref2) then
            write (iunout,*) ' reference cell ',nref2
            write (iunout,*) ' bxin, byin, bzin, bfin ',
     .             bxin(icell), byin(icell), bzin(icell), bfin(icell)
            write (iunout,*) ' bx,   by,   bz,   bf   ',
     .             bx, by, bz, bf
          elseif (icell == nref3) then
            write (iunout,*) ' reference cell ',nref3
            write (iunout,*) ' bxin, byin, bzin, bfin ',
     .             bxin(icell), byin(icell), bzin(icell), bfin(icell)
            write (iunout,*) ' bx,   by,   bz,   bf   ',
     .             bx, by, bz, bf
          end if
        end if
        end if

      end do

      adin(2,1:nsbox) = bzin(1:nsbox) * bfin(1:nsbox)

      if (.false.) then
      errbx = errbx / nplcll
      errby = errby / nplcll
      errbz = errbz / nplcll
      errbf = errbf / nplcll

      write (iunout,*) ' no of cell in plasma region ', nplcll
      write (iunout,*) ' errbx, errby, errbz, errbf ',
     .                   errbx, errby, errbz, errbf

      write (iunout,*) ' bxmax, bymax, bzmax, bfmax ',
     .                   bxmax, bymax, bzmax, bfmax
      end if

 99   RETURN
      END
C ===== SOURCE: pltusr.f
C
C
      SUBROUTINE EIRENE_PLTUSR(PLABLE,J)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLGIN
      IMPLICIT NONE
      LOGICAL, INTENT(INOUT) :: PLABLE
      INTEGER, INTENT(IN) :: J
      INTEGER :: IB, MERK, IER
      REAL(DP) :: XS, YS, ZX0, ZY0, ZZ0, CX, CY, CZ, RZYLB, B0B, B1B,
     .          B2B, B3B, F0B, F1B, F2B, F3B, X0, Y0, B0, B1, B2, Z1, Z2
      INTEGER :: IN
      RETURN
      END
C ===== SOURCE: prousr.f
CDK USER
C
C   USER SUPPLIED SUBROUTINES
C
C           ************
C           *  EDGE2D  *  (fem-interface)
C           ************
C
C s.wiesen@fz-juelich.de (mar06)
C
C get plasma values for bulk species IPLS (module comprt, set by subr. plasma)
c
c reads casename.zplasma
c
      SUBROUTINE EIRENE_PROUSR (PRO,INDX,P0,P1,P2,P3,P4,P5,PROVAC,N)
C
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CINIT
      USE EIRMOD_CCONA
      USE EIRMOD_COMPRT

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: P0, P1, P2, P3, P4, P5, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      INTEGER, INTENT(IN) :: INDX, N
      CHARACTER(256) :: FILENAME, line,sstr
      INTEGER :: NTR, I, J, LL,ier
      REAL(DP), ALLOCATABLE, SAVE :: PLAS(:,:,:)
      REAL(DP) :: bnorm
      INTEGER, SAVE :: INDAR(11)=(/ (0, i=1,11) /)
      integer, parameter :: fp=31
      character(2) :: cstr2
      integer :: idum
      real(dp) :: rdum
csw
c     reset?
      if(indx < 0) then
         if(allocated(plas)) deallocate(plas)
         indar(1:11) = 0
c just in case...:
         close(fp+ifoff)
         return
      endif
csw      
      
c read in plasma data from fort.31 ?
      if (.not.allocated(plas)) then
        allocate(plas(11,nrad,0:npls))
        plas=0.

        ll=len_trim(casename)
        filename=casename(1:ll) // '.zplasma'
        open (unit=fp+ifoff,file=filename,access='sequential',
     .        form='formatted')


c misc plasma data:
        write(sstr,'(a20)') 
     .          '*** MISC PLASMA DATA'
        CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
        if(ier /=0) then
           write(*,*) 'PROUSR: ',sstr,' not found'
           close(fp+ifoff)
           call EIRENE_exit_own(1)
        endif
        do j=1,2
           read(fp+ifoff,'(a)') line
        enddo

        read(fp+ifoff,'(i7)') ntr
        if (ntr /= nr1st-1) then
           write (*,*) 'PROUSR:', sstr
           write (*,*) ' wrong number of triangles in plasma file'
           write (*,*) ' check for correct number in file ',filename
           call EIRENE_exit_own(1)
        endif

        do j=1,ntr
           read(fp+ifoff,'(i7,7(1x,e14.7))') idum,
c     .          te,ne,bx,by,bz,pot,psi
     .          plas(1,j,0),
     .          rdum,
     .          plas(7,j,0),
     .          plas(8,j,0),
     .          plas(9,j,0),
     .          rdum,
!pb     .          rdum
     .          plas(6,j,0)

           bnorm=sqrt(plas(7,j,0)**2 + plas(8,j,0)**2 + plas(9,j,0)**2)
           if (bnorm < eps12) then
              plas(7,j,0) = 0._dp
              plas(8,j,0) = 0._dp
              plas(9,j,0) = 1._dp
              plas(10,j,0) = 0._dp
           else
              plas(7:9,j,0) = plas(7:9,j,0)/bnorm
              plas(10,j,0) =  bnorm
           endif
        enddo
        plas(9,ntr+1:nrad,0) = 1._dp

c loop over species:
        do i=1,npls
           write(cstr2,'(i2.2)') i
           write(sstr,'(a23)') 
     .          '*** ION #'//cstr2//' PLASMA DATA'
           CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
           if(ier /=0) then
              write(*,*) 'PROUSR: tag ',sstr,'not found'
              close(fp+ifoff)
              call EIRENE_exit_own(1)
           endif

           do j=1,3
              read(fp+ifoff,'(a)') line
           enddo

           read(fp+ifoff,'(i7)') ntr
           if (ntr /= nr1st-1) then
              write (*,*) 'PROUSR:', sstr
              write (*,*) ' wrong number of triangles in plasma file'
              write (*,*) ' check for correct number in file ',filename
              call EIRENE_exit_own(1)
           endif

           do j=1,ntr
              read(fp+ifoff,'(i7,6(1x,e14.7))') idum,
c     .             ti,ni,vx,vy,vz,zi
     .             plas(2,j,i),
     .             plas(3,j,i),
     .             plas(4,j,i),
     .             plas(5,j,i),
     .             plas(6,j,i),
     .             plas(11,j,i)
           enddo
        enddo

        close(fp+ifoff)
      endif

c
c fort.31 read in now
c
c sort in values
c
      if(indx == 0) then
c te         
         pro(1:n) = plas(1,1:n,0)
         indar(1) = indar(1)+1

      elseif (indx == 1) then
c ti
         pro(1:n) = plas(2,1:n,ipls)
         indar(2) = indar(2)+1

      elseif (indx == 1+1*npls) then
c ni
         pro(1:n) = plas(3,1:n,ipls)
         indar(3) = indar(3) + 1

      elseif (indx == 1+2*npls) then
c vx
         pro(1:n) = plas(4,1:n,ipls)
         indar(4) = indar(4) + 1

      elseif (indx == 1+3*npls) then
c vy
         pro(1:n) = plas(5,1:n,ipls)     
         indar(5) = indar(5) + 1

      elseif (indx == 1+4*npls) then
c vz
         pro(1:n) = plas(6,1:n,ipls)
         indar(6) = indar(6) + 1

      elseif (indx == 1+5*npls) then
c bx
         pro(1:n) = plas(7,1:n,0)
         indar(7) = indar(7) + 1

      elseif (indx == 2+5*npls) then
c by
         pro(1:n) = plas(8,1:n,0)
         indar(8) = indar(8) + 1

      elseif (indx == 3+5*npls) then
c bz
         pro(1:n) = plas(9,1:n,0)
         indar(9) = indar(9) + 1

      elseif (indx == 4+5*npls) then
c bf
         pro(1:n) = plas(10,1:n,0)
         indar(10) = indar(10) + 1

      elseif (indx == 5+5*npls) then
! psi
         pro(1:n) = plas(6,1:n,0)

      elseif (indx == 7+5*npls) then
! zi
         pro(1:n) = plas(11,1:n,ipls)

      else
         write (iunout,*) ' prousr: no data provided for index ',indx
         pro(1:n) = 0._dp
      endif

      return
      end
C ===== SOURCE: refusr.f


      SUBROUTINE EIRENE_REFUSR
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XMW,XCW,XMP,XCP,ZCOS,ZSIN,EXPI,RPROB,
     .                        E0TERM
      INTEGER, INTENT(IN) :: IGASF,IGAST
      ENTRY EIRENE_RF0USR
      ENTRY EIRENE_SPTUSR
      ENTRY EIRENE_SP0USR
      ENTRY EIRENE_SP1USR
      ENTRY EIRENE_RF1USR (XMW,XCW,XMP,XCP,IGASF,IGAST,ZCOS,ZSIN,EXPI,
     .              RPROB,E0TERM,*,*,*,*)
      RETURN
      END
C ===== SOURCE: retusr.f
c
c
      subroutine EIRENE_retusr(sig)
      USE EIRMOD_PRECISION
      implicit none
      real(dp), intent(in) :: sig
      return
      end
C ===== SOURCE: samusr.f
C
C
      SUBROUTINE EIRENE_SAMUSR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)
C
C  SAMPLE INITAL COORDIANTES X,Y,Z ON ADDITIONAL SURFACE NLLI
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CADGEO
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6
      REAL(DP), INTENT(OUT) :: X0,Y0,Z0,TEWL,TIWL(*),DIWL(*),
     .                         VXWL(*),VYWL(*),VZWL(*),
     .                         EFWL(*), SHWL, WEISPZ(*)
      INTEGER, INTENT(IN) :: NLSF,is1, is2
      INTEGER, INTENT(OUT) :: IRUSR, IPUSR, ITUSR, IAUSR, IBUSR
      REAL(DP) :: X, Y, T, B0, B1, B2, Z1, Z2
      REAL(DP), EXTERNAL :: EIRENE_RANF_EIRENE
      INTEGER :: IER

      entry EIRENE_sm0usr (is1,is2,sorad1,sorad2,sorad3,
     &     sorad4,sorad5,sorad6)
      return

      entry EIRENE_SM1USR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)

      RETURN
      END
C ===== SOURCE: sigusr.f


      SUBROUTINE EIRENE_SIGUSR(IFIRST,JJJ,ZDS,DUMMY1,PSIG,DUMMY2,ARGST,
     .                  XD0,YD0,ZD0,XD1,YD1,ZD1)
C
C  INPUT:
C          IFIRST: FLAG FOR INITIALISATION
C          NCELL:   INDEX IN TALLY ARRAYS FOR CURRENT ZONE
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C  OUTPUT: CONTRIB. FROM CELL NCELL AND CHORD SEGMENT JJJ TO:
C          THE H ALPHA FLUX PSIG(I),I=0,4 CONTRIBUTIONS
C          FROM ATOMS, MOLECULES, TEST IONS AND BULK IONS
C          THE INTEGRANT ARGST SUCH THAT INTEGR.(ARGST*DL) = PSIG
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CSPEI
      USE EIRMOD_COMSOU
      USE EIRMOD_CLOGAU
      USE EIRMOD_COMXS
      USE EIRMOD_COMSIG
      USE EIRMOD_CUPD
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGRID
      USE EIRMOD_CCONA
      USE EIRMOD_CADGEO
      USE EIRMOD_CLGIN
      USE EIRMOD_CSDVI
      USE EIRMOD_COUTAU
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JJJ
      INTEGER, INTENT(INOUT) :: IFIRST
      REAL(DP), INTENT(INOUT) :: PSIG(0:NSPZ+10),ARGST(0:NSPZ+10,NRAD)
      REAL(DP), INTENT(IN) :: ZDS,DUMMY1,DUMMY2,XD0,YD0,ZD0,XD1,YD1,ZD1
      RETURN
      END
C ===== SOURCE: talusr.f
c
c
      subroutine EIRENE_talusr (ICOUNT,VECTOR,TALTOT,TALAV,
     .              TXTTL,TXTSP,TXTUN,ILAST,*)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      implicit NONE
      integer, intent(in) :: icount
      integer, intent(out) :: ilast
      real(dp), intent(in) :: vector(*), TALTOT, TALAV
      character(len=*) :: txttl,txtsp,txtun
      integer :: i

      ilast=1
      return 1
      end
C ===== SOURCE: timusr.f


      SUBROUTINE EIRENE_TIMUSR(N,X,Y,Z,VX,VY,VZ,N1,N2,T,IC,IE,NP,NL)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: X,Y,Z,VX,VY,VZ,T,cx,cy,cz,sc
      INTEGER, INTENT(IN) :: N, N1, N2, IC, IE, NP, IS, NRCELL
      LOGICAL :: NL

      ENTRY EIRENE_NORUSR(is,x,y,z,cx,cy,cz,sc,VX,VY,VZ,NRCELL)

      RETURN
      END
C ===== SOURCE: tmsusr.f


      SUBROUTINE EIRENE_TMSUSR (T0)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: T0
      RETURN
      END
C ===== SOURCE: upcusr.f
C
C
      SUBROUTINE EIRENE_UPCUSR(WS,IND)
C
C  USER SUPPLIED COLLISION ESTIMATOR, VOLUME AVERAGED
C
C+---------------------------------------------------------------+
C| Modifications:                                                |
C| --------------                                                |
C| 16/07/2010   D.Harting    Added two variables to eirene_user  |
C|                           namelist for use of fluxdependency  |
C|                           in chemical sputtering.             |
C+---------------------------------------------------------------+
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
      USE EIRMOD_CGRID
      USE EIRMOD_CCONA
      use EIRMOD_CUPD
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WS
      INTEGER, INTENT(IN) :: IND
cswx 24sep07
      logical,save :: lfirst=.false.
      integer, save :: num=0,eirene_nbirth,eirene_njetto
      real*8,save :: delang
      real*8 :: xx,yy,zz,d,pphi,alph,v0v
      character(len=256), save :: eirene_fbirth,eirene_ftransfer,
     &     eirene_fstoreneutflux
      real(dp) :: eirene_phi_offsets(9)
      integer :: eirene_wallFluxModel ! calculation of wall fluxes for chemical sputtering
c                = 0: no wall fluxes are used (old edge2d model)
c                = 1: only ion fluxes are used
c                = 2: ion fluxes and neutral fluxes from last eirene iteration are used
c                = 3: ion and neutral fluxes are used, and EIRENE is iterated to give 
c                     converged neutral fluxes.
      logical :: eirene_use_elstepdat_bug
      real*8, allocatable,save :: rdata(:,:)
      integer, allocatable,save :: idata(:,:)
      integer :: i,j
      real*8 :: cosrot, sinrot, my_velx, my_vely, my_velz
      
      namelist /eirene_user/eirene_nbirth,eirene_njetto,
     .                      eirene_fbirth,eirene_ftransfer,
     .                      eirene_phi_offsets,
     .                      eirene_fstoreneutflux,
     .                      eirene_wallFluxModel,
     .                      eirene_use_elstepdat_bug
cswx
C
C     WS=WEIGHT/SIGTOT=WEIGHT/(VEL*ZMFPI)=WEIGHT/(VEL*SIGMA,MACR.)
C
C  FOR PARTICLE DENSITY IN CELL NO. NCELL
C     COLV(1,NCELL)=COLV(1,NCELL)+WS
C
cswx 24sep07
      if(.not.lfirst) then
        lfirst=.true.
        num=0
        delang = pi2a/dble(nttram)

c default maximum number
        eirene_nbirth=2000
        eirene_fbirth = ' '
cdmh
      eirene_fstoreneutflux = 'eirene.chemFluxDep'
      eirene_wallFluxModel = 1
cdmh      

        open(unit=9998,file='eirene_user.namelist')
        read(9998,eirene_user)
        close(9998)
        if(eirene_nbirth .gt. 0) then
          allocate(idata(eirene_nbirth,2))
          allocate(rdata(eirene_nbirth,9))
        endif
      endif

      if(         eirene_nbirth .gt. 0
     .      .and. ind .eq. 1 
     .      .and. num .lt. eirene_nbirth 
     .      .and. weight .gt. 1.d-2
     .      .and. e0.gt.2.e4) then
        num=num+1

        d =sqrt(x0*x0 + z0*z0)
        alph = pi2a/dble(nt3rd-1)
        pphi = mod(phi+pi2a-(alph/2.d0),pi2a)
        !pphi = mod(phi+pi2a,pi2a)
        xx = cos(pphi) * d
        yy = y0
        zz = sin(pphi) * d

        alph=0.
        alph = alph - rotsav_torcol
        ! back-rotate from torcol, periodic boundary
        alph = alph +(360.d0/dble(nt3rd-1))/360d0  *pi2a
        ! add offsets from beamboxes
        alph = alph + eirene_phi_offsets(istra)/360.d0 * pi2a

        cosrot=cos(alph)
        sinrot=sin(alph)
        my_velx = velx*cosrot + velz*sinrot
        my_vely = vely
        my_velz =-velx*sinrot + velz*cosrot
c        write(598, '(a,20(1x,e13.6))') 'UPCUSR:',
c     .        rotsav_torcol, my_velx, my_vely, my_velz
        v0v=my_velx*bxin(ncell)+my_vely*byin(ncell)+my_velz*bzin(ncell)
        idata(num,1) = num
        idata(num,2) = 2
        rdata(num,1) = sqrt(xx*xx + zz*zz)/100.d0
        rdata(num,2) = pphi/degrad
        rdata(num,3) = yy/100.d0
        rdata(num,4) = -v0v ! not used in ascot anymore
        rdata(num,5) = e0
        rdata(num,6) = weight
        rdata(num,7) = my_velx*vel/100.d0
csw change ordering
        rdata(num,8) = -my_velz*vel/100.d0
        rdata(num,9) = my_vely*vel/100.d0
      endif
cswx
      RETURN

      entry EIRENE_upcusr_reinit
      if(lfirst) then
        lfirst=.false.

        if(eirene_nbirth .gt. 0) then
          open(unit=9998,file=trim(eirene_fbirth),
     .       form='formatted',access='sequential',status='replace')
          rewind(9998)
          write(9998,'(i7)') num
          write(9998,'(a1,a4,a4,9a14)') '#','NP','NZ','R[m]',
     .                                'PHI[deg]',
     .                                'Z[m]',
     .                                'V0/V','E[eV]','W',
     .                                'VX[m/s]','VY[m/s]','VZ[m/s]'
          do i=1,num
            write(9998,'(i5,i4,20e14.5)') 
     .         (idata(i,j), j=1,2),
     .         (rdata(i,j), j=1,9)
          enddo
          close(9998)
          deallocate(idata)
          deallocate(rdata)
        endif
      endif
      return
      END
C ===== SOURCE: upnusr.f
c
c
      subroutine EIRENE_upnusr
      return
      end
C ===== SOURCE: upsusr.f
C
C
      SUBROUTINE EIRENE_UPSUSR(WT,IND)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WT
      INTEGER, INTENT(IN) :: IND
      RETURN
      END
C ===== SOURCE: uptusr.f
C
C
      SUBROUTINE EIRENE_UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
C+---------------------------------------------------------------+
C| Modifications:                                                |
C| --------------                                                |
C| 23/08/2006                VPX, VPY, VRX, VRY changed to       |
C|                           ALLOCATABLE, SAVE to speed up       |
C|                           subroutine call (save time in       |
C|                           storage allocation)                 |
C| 16/07/2010   D.Harting    Added two variables to eirene_user  |
C|                           namelist for use of fluxdependency  |
C|                           in chemical sputtering.             |
C+---------------------------------------------------------------+
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CUPD
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEZ
      USE EIRMOD_CGRID
      USE EIRMOD_CLOGAU
      USE EIRMOD_CCONA
      USE EIRMOD_CPOLYG
      USE EIRMOD_CZT1
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                         XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP), ALLOCATABLE, SAVE :: CNDYNA(:),CNDYNP(:)
CDR
      REAL(DP), ALLOCATABLE, SAVE :: VPX(:),VPY(:),VRX(:),VRY(:)
CDR
      INTEGER :: IAT, IPL, I, IR, IP, IRD
      INTEGER, SAVE :: IFIRST, IA1, IA2, IA3, NA4, INDEXM, INDEXF
      DATA IFIRST/0/
csw
      real(dp) :: dist,wtr
      integer :: iaei,irei,iacx,ircx
      integer, save :: num=0,eirene_nbirth,eirene_njetto
      character(len=256), save :: eirene_fbirth,eirene_ftransfer,
     &     eirene_fstoreneutflux
      real(dp) :: eirene_phi_offsets(9)
      integer :: eirene_wallFluxModel ! calculation of wall fluxes for chemical sputtering
c                = 0: no wall fluxes are used (old edge2d model)
c                = 1: only ion fluxes are used
c                = 2: ion fluxes and neutral fluxes from last eirene iteration are used
c                = 3: ion and neutral fluxes are used, and EIRENE is iterated to give 
c                     converged neutral fluxes.
      logical :: eirene_use_elstepdat_bug
      namelist /eirene_user/eirene_nbirth,eirene_njetto,
     .                      eirene_fbirth,eirene_ftransfer,
     .                      eirene_phi_offsets,
     .                      eirene_fstoreneutflux,
     .                      eirene_wallFluxModel,
     .                      eirene_use_elstepdat_bug
csw

      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        ALLOCATE (CNDYNA(NATM))
        ALLOCATE (CNDYNP(NPLS))
        DO IAT=1,NATMI
          CNDYNA(IAT)=1.D3*AMUA*RMASSA(IAT)
        END DO
        DO IPL=1,NPLSI
          CNDYNP(IPL)=1.D3*AMUA*RMASSP(IPL)
        END DO
C
CDR
CDR  PROVIDE A RADIAL UNIT VECTOR PER CELL
CDR  VPX,VPY,  NEEDED FOR PROJECTING PARTICLE VELOCITIES
CDR  SAME FOR POLOIDAL UNIT VECTOR VRX,VRY
C
        ALLOCATE (VPX(NRAD))
        ALLOCATE (VPY(NRAD))
        ALLOCATE (VRX(NRAD))
        ALLOCATE (VRY(NRAD))
        DO I=1,NRAD
          VPX(I)=0.
          VPY(I)=0.
          VRX(I)=0.
          VRY(I)=0.
        END DO
        DO IR=1,NR1STM
          DO IP=1,NP2NDM
            IRD=IR+(IP-1)*NR1P2
            VPX(IRD)=PLNX(IR,IP)
            VPY(IRD)=PLNY(IR,IP)
            VRX(IRD)=PPLNX(IR,IP)
            VRY(IRD)=PPLNY(IR,IP)
          END DO
        END DO
        IA1=NATMI+NMOLI
        IA2=2*IA1
        IA3=3*IA1
        NA4=4*IA1
        INDEXM=NPLSI
        INDEXF=2*NPLSI
csw
        eirene_njetto=0
cdmh
        eirene_fstoreneutflux = 'eirene.chemFluxDep'
        eirene_wallFluxModel = 1
cdmh      
        open(unit=9998,file='eirene_user.namelist')
        read(9998,eirene_user)
        close(9998)
csw
      ENDIF
csw
csw 24jan08 additional tallies for JETTO cold neutral coupling
csw
      if(       eirene_njetto .gt. 0 .and.
     .          nadv .eq. 8 .and. ityp.eq.1 .and. iatm.eq.1) then
        do i=1,ncou
          dist = clpd(i)
          ird = nrcell+nupc(i)*NR1P2+NBLCKA
          wtr = wv*dist
c         sources/sinks due ionisation:
          iaei = 1
          irei=lgaei(iatm,iaei)
          addv(1,ird) = addv(1,ird)+wtr*sigvei(irei)
          addv(4,ird) = addv(4,ird)+wtr*sigvei(irei)*esigei(irei,5)
          addv(5,ird) = addv(5,ird)+wtr*sigvei(irei)*esigei(irei,4)

c         sources/sinks due recombination:
c           not included

c         net sources due CX:
          iacx = 1
          ircx=lgacx(iatm,iacx,0)
          addv(8,ird) = addv(8,ird)+wtr*sigvcx(ircx)*
     .               (E0 - esigcx(ircx,1))    
        enddo
      endif
      RETURN
      entry EIRENE_uptusr_reinit
      if(ifirst .ne. 0) then
        ifirst=0
        if(allocated(cndyna)) deallocate(cndyna)
        if(allocated(cndynp)) deallocate(cndynp)
        if(allocated(vpx)) deallocate(vpx)
        if(allocated(vpy)) deallocate(vpy)
        if(allocated(vrx)) deallocate(vrx)
        if(allocated(vry)) deallocate(vry)
      endif
      return
      END



C ===== SOURCE: vdion.f


      FUNCTION EIRENE_VDION (I)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      REAL(DP) :: EIRENE_VDION
      EIRENE_VDION=0.
      RETURN
      END
C ===== SOURCE: vecusr.f


      SUBROUTINE EIRENE_VECUSR (I,VX,VY,VZ,IPLS)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, IPLS
      REAL(DP), INTENT(IN) :: VX,VY,VZ
      RETURN
      END
C ===== SOURCE: volusr.f


      SUBROUTINE EIRENE_VOLUSR(N,A)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: A(*)
      INTEGER, INTENT(IN) :: N
      RETURN
      END
