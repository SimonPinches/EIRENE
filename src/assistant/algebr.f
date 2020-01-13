cdr  Aug 19  : Cleanup. The routine MECKER is
cdr            now contained here.
cdr            tbd:
cdr            Similarly: ZERLEG, RUKSUB, OPRAND,
cdr            SUBTIT, SCHRIT, SIGNOK can also
cdr            be made local to this present module.
cdr
Cdr  Sept. 16: Bug fix: added option: two or more constants next to each other
C
C-----------------------------------------------------------------------
      SUBROUTINE EIRENE_ALGEBR (TERM,OPER,IZIF,CONST,NOP)
C-----------------------------------------------------------------------
C
C     AUTOR:           ST. HUBER
C     IHK-KENNZIFFER:  121
C     DATUM:           25-NOV-1988
C
C     FUNKTION:
C
C     DAS PROGRAMM LIEST ARITHMETISCHE AUSDRUECKE EIN,
C     RUFT DAS UNTERPROGRAMM ZERLEG AUF
C     UND GIBT ENTWEDER DIE EINZELNEN ZERLEGUNGEN ODER DIE
C     REGELVERLETZUNG AUS
C
C-----------------------------------------------------------------------
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
C
C     ARGUMENT:
C
         CHARACTER(*), INTENT(INOUT) :: TERM
C           : EINZULESENDER AUSDRUCK

         CHARACTER(2), INTENT(OUT) :: OPER(*)
         INTEGER, INTENT(OUT) :: IZIF(4,*)
         REAL(DP), INTENT(INOUT) :: CONST(*)
         INTEGER, INTENT(OUT) :: NOP
C
C     KONSTANTENDEKLARATION :
C
         INTEGER, PARAMETER :: ZMAX = 20
C           : ANZAHL DER MAXIMALEN ZERLEGUNGEN

         INTEGER, PARAMETER :: MAXLEN = 72
C           : MAXIMALE STRINGLAENGE

C
C     LOKALE VARIABLEN :
C
         INTEGER ::  LAENGE
C           : AKTUELLE LAENGE VON TERM

         CHARACTER(MAXLEN) :: HLFTERM
C           : HILFSSTRING ZUM UMSPEICHERN

         CHARACTER(MAXLEN+2):: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

         INTEGER :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         INTEGER :: TEIL
C           : AKTUELLE ANZAHL DER ZERLEGUNGEN

         CHARACTER(MAXLEN) :: PART(ZMAX)
C           : FELD VON STRINGS, AUF DENEN DIE EINZELNEN
C             ELEMENTARZERLEGUNGEN FESTGEHALTEN WERDEN

         INTEGER :: IPART(ZMAX)
C           : AKTUELLE LAENGEN VON PART(ZMAX)

         CHARACTER(MAXLEN) :: ARITH(ZMAX)
C           : FELD VON STRINGS, AUF DENEN DIE TEIL-TE GENERATION
C             VON AUSDRU FESTGEHALTEN WIRD

         INTEGER :: IARITH(ZMAX)
C           : AKTUELLE LAENGEN VON ARITH(ZMAX)

         CHARACTER(MAXLEN) :: HILFE
C           : ARBEITSSPEICHER FUER UNTERPROGRAMM ZERLEG

         INTEGER :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

CHR
CHR      VARIABLEN ZUR MODIFIKATION DES PROGRAMMES
         INTEGER :: NR, ANFANG, ENDE, FELDIND, IK, IKM, IKP
         CHARACTER(4) :: FO
         CHARACTER(12) :: ERSETZ(10)
         character(10) :: buchst
chr
C
C     HILFSVARIABLEN :
C
         INTEGER :: I, IC
chr
chr   string, der die neuen variablennamen enthaelt
      buchst='ABCDEFGHIJ'
chr
!pb count number of constant terms
      ic = 0

!pb change TERM to uppercase
      call eirene_uppercase(term)
C
C        LESE TERM UND WERTE AUS
C
         IF (TERM .NE. ' ') THEN
C
C           VERARBEITUNG, FUER DEN FALL, DASS KEINE LEERZEILE
C           EINGELESEN WURDE
C
chr         vorbereiten des terms fuer die weitere verarbeitung, d.h.
chr         bringen der operanden in die vom programm verlangte
chr         zweistellige alphabetische form
            nr=0
  101       anfang=index(term,'<')
            if (anfang.ne.0) then
               nr=nr+1
               ende=index(term,'>')
chr            abspeichern der ersetzten operanden
chr            der i-te operand wird hierbei in der i-ten dimension
chr            des feldes ersetz abgelegt und durch den i-ten
chr            buchstaben im alphabet ersetzt. gleichheit von
chr            operanden wird hierbei nicht beruecksichtigt
               ersetz(nr)=term(anfang:ende)
               if (anfang.gt.1) then
               hlfterm=term(1:anfang-1)//buchst(nr:nr)//buchst(nr:nr)//
     .              term(ende+1:)
               else
               hlfterm=buchst(nr:nr)//buchst(nr:nr)//
     .              term(ende+1:)
               endif
               term=hlfterm
               goto 101
            endif
chr
C
C           ERMITTELN DER LAENGE VON TERM
C
            LAENGE=LEN_TRIM(TERM)

            AUSDRU=TERM
            AKTLEN=LAENGE

            CALL EIRENE_ZERLEG(AUSDRU, AKTLEN, IPART, PART, IARITH,
     .                         ARITH, TEIL, HILFE, ERROR)

            IF (ERROR .EQ. 0) THEN
C
C              AUSGABE DER ZERLEGUNG
C
               NOP=TEIL
               DO 30, I=1,TEIL
chr               ausgabe der zerlegung in der form:
chr                    operator ziffer1 ziffer2 ziffer3 ziffer4
chr               wobei jeweils die 1. und 2. sowie die 3. und 4.
chr               ziffer eine einheit bilden.
chr               entweder stellt eine solche einheit einen operanden
chr               dar oder ein zwischenergebnis mit der 1. ziffer als
chr               nummer und der 2. als 0 zur kennzeichnung des paares
chr               als zwischenergebnis
                  oper(i) = '  '
                  izif(1:4,i) = 0
                  if (part(i)(7:7).eq.'Q') then
                    oper(i) = part(i)(7:8)
                    feldind=index(buchst,part(i)(10:10))
                    if (feldind == 0) then
                      READ(PART(I)(12:12),'(I1)') IZIF(1,I)
                      IZIF(2,I)=0
                    else
                      IK=INDEX(ERSETZ(FELDIND),',')
                      IF (IK.EQ.0) THEN
                        IC = IC + 1
                        IZIF(1,I)=-I
                        IZIF(2,I)= IC
                        CALL EIRENE_RDCN (ERSETZ(FELDIND),CONST(IC))
                      ELSE
                        IKM=IK-1
                        IKP=IK+1
                        FO(1:4)='(I )'
                        WRITE (FO(3:3),'(I1)') IK-2
                        READ(ERSETZ(FELDIND)(2:IKM),FO) IZIF(1,I)
                        IF (ERSETZ(FELDIND)(IK+2:IK+2).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP),'(I1)')
     .                         IZIF(2,I)
                        ELSEIF (ERSETZ(FELDIND)(IK+3:IK+3).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP+1),'(I2)')
     .                         IZIF(2,I)
                        ELSEIF (ERSETZ(FELDIND)(IK+4:IK+4).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP+2),'(I3)')
     .                         IZIF(2,I)
                        ENDIF
                      ENDIF
                    ENDIF
                  elseif (part(i)(7:7).ne.'Z') then
                     OPER(I)=PART(I)(9:9)
                     feldind=index(buchst,part(i)(7:7))
                     IK=INDEX(ERSETZ(FELDIND),',')
                     IF (IK.EQ.0) THEN
                       IC = IC + 1
                       IZIF(1,I)=-I
                       IZIF(2,I)= IC
                       CALL EIRENE_RDCN (ERSETZ(FELDIND),CONST(IC))
                     ELSE
                       IKM=IK-1
                       IKP=IK+1
                       FO(1:4)='(I )'
                       WRITE (FO(3:3),'(I1)') IK-2
                       READ(ERSETZ(FELDIND)(2:IKM),FO) IZIF(1,I)
                       IF (ERSETZ(FELDIND)(IK+2:IK+2).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP),'(I1)')
     .                         IZIF(2,I)
                       ELSEIF (ERSETZ(FELDIND)(IK+3:IK+3).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP+1),'(I2)')
     .                         IZIF(2,I)
                       ELSEIF (ERSETZ(FELDIND)(IK+4:IK+4).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP+2),'(I3)')
     .                         IZIF(2,I)
                       ENDIF
                     ENDIF

                     IF (PART(I)(10:10).NE.'Z') THEN
                       FELDIND=INDEX(BUCHST,PART(I)(10:10))
                       IK=INDEX(ERSETZ(FELDIND),',')
                       IF (IK.EQ.0) THEN
                         IC = IC + 1
                         IZIF(3,I)=-I
                         IZIF(4,I)=IC
                         CALL EIRENE_RDCN (ERSETZ(FELDIND),CONST(IC))
                       ELSE
                         IKM=IK-1
                         IKP=IK+1
                         FO(1:4)='(I )'
                         WRITE (FO(3:3),'(I1)') IK-2
                         READ(ERSETZ(FELDIND)(2:IKM),FO) IZIF(3,I)
                         IF (ERSETZ(FELDIND)(IK+2:IK+2).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP),'(I1)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+3:IK+3).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+1),'(I2)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+4:IK+4).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+2),'(I3)')
     .                          IZIF(4,I)
                         ENDIF
                       ENDIF
                     ELSE
                        READ(PART(I)(12:12),'(I1)') IZIF(3,I)
                        IZIF(4,I)=0
                     ENDIF
                  ELSE
                     OPER(I)=PART(I)(10:10)
                     READ(PART(I)(9:9),'(I1)') IZIF(1,I)
                     IZIF(2,I)=0
                     if (part(i)(11:11).ne.'Z') then
                       FELDIND=INDEX(BUCHST,PART(I)(11:11))
                       IK=INDEX(ERSETZ(FELDIND),',')
                       IF (IK.EQ.0) THEN
                         IC = IC + 1
                         IZIF(3,I)=-I
                         IZIF(4,I)=IC
                         CALL EIRENE_RDCN (ERSETZ(FELDIND),CONST(IC))
                       ELSE
                         IKM=IK-1
                         IKP=IK+1
                         FO(1:4)='(I )'
                         WRITE (FO(3:3),'(I1)') IK-2
                         READ(ERSETZ(FELDIND)(2:IKM),FO) IZIF(3,I)
                         IF (ERSETZ(FELDIND)(IK+2:IK+2).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP),'(I1)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+3:IK+3).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+1),'(I2)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+4:IK+4).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+2),'(I3)')
     .                          IZIF(4,I)
                         ENDIF
                       ENDIF
                     else
                        READ(PART(I)(13:13),'(I1)') IZIF(3,I)
                        IZIF(4,I)=0
                     endif
                  endif
   30          CONTINUE
            ELSE
C
C              AUSGABE DER FEHLERMELDUNG
C
               WRITE(iunout,'(2A)') 'FOLGENDE REGELVERLETZUNG ',
     >                        'WURDE ERKANNT:'
               CALL EIRENE_MECKER(ERROR)
               NOP=0
            ENDIF
         ENDIF
C
C     ENDE VON ALGEBR
C
      RETURN

      CONTAINS

cdr This routine writes "complaint messages"
cdr from ALGEBR, i.e. from the routines
cdr that try to decipher the coded algebraic
cdr expressions for the algebraic tallies
cdr ALGV, ALGS specified in input block 10C and 10E, resp.
cdr Mecker.f should be moved into ALGEBR.f


C-----------------------------------------------------------------------
                SUBROUTINE EIRENE_MECKER(ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     AUSGABE DER REGELVERLETZUNGEN
C
C-----------------------------------------------------------------------
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN


      IF     (ERROR .EQ.  1) THEN
        WRITE(iunout,*) 'DER AUSDRUCK ENTHAELT EIN UNGUELTIGES ZEICHEN.'
      ELSEIF (ERROR .EQ.  2) THEN
        WRITE(iunout,*)
     >  'ZWISCHEN ZWEI KLAMMERAUSDRUECKEN BEFINDET SICH',
     >  ' KEIN OPERATOR.'
      ELSEIF (ERROR .EQ.  3) THEN
         WRITE(iunout,*) 'EIN KLAMMERAUSDRUCK IST LEER.'
      ELSEIF (ERROR .EQ.  4) THEN
         WRITE(iunout,*) 'DER AUSDRUCK IST FALSCH GEKLAMMERT.'
      ELSEIF (ERROR .EQ.  5) THEN
         WRITE(iunout,*) 'DER AUSDRUCK ENTHAELT MEHR ALS 3 INEINANDER',
     >              'GESCHACHTELTE KLAMMERN.'
      ELSEIF (ERROR .EQ.  6) THEN
        WRITE(iunout,*)
     >  'EIN OPERAND BESTEHT AUS MEHR ALS ZWEI BUCHSTABEN.'
      ELSEIF (ERROR .EQ.  7) THEN
         WRITE(iunout,*) 'ZWEI OPERATOREN STEHEN NEBENEINANDER.'
      ELSEIF (ERROR .EQ.  8) THEN
         WRITE(iunout,*) 'NACH EINEM OPERATOR FOLGT EINE SCHLIESSENDE',
     >              ' KLAMMER.'
      ELSEIF (ERROR .EQ.  9) THEN
        WRITE(iunout,*)
     >  'DAS ERSTES ZEICHEN DES AUSDRUCKS IST OPERATOR, ',
     >  'UND KEIN PRAEFIX.'
      ELSEIF (ERROR .EQ. 10) THEN
         WRITE(iunout,*) 'DER AUSDRUCK ENDET MIT EINEM OPERATOR.'
      ELSEIF (ERROR .EQ. 11) THEN
         WRITE(iunout,*) 'DAS ERSTE ZEICHEN IN EINEM KLAMMERAUSDRUCK',
     >              ' IST EIN OPERATOR.'
      ELSEIF (ERROR .EQ. 12) THEN
         WRITE(iunout,*) 'DER AUSDRUCK ENTHAELT MEHR ALS 15 OPERATOREN.'
      ELSEIF (ERROR .EQ. 13) THEN
         WRITE(iunout,*) 'DIE (ANZAHL DER OPERANDEN)-1 IST UNGLEICH ',
     >              ' DER (ANZAHL DER OPERATOREN).'
      ENDIF
C
C     ENDE VON MECKER
C
      END SUBROUTINE EIRENE_MECKER

      END
