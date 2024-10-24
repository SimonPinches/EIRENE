      module eirmod_algebra

      USE EIRMOD_PRECISION

      private

      public :: eirene_algebr

      integer, save :: iprout

      contains

!pb  Apr 20  : subroutine EIRENE_ALGEBRA is now contained in a module
!pb            all associated routines are placed in the module as well.
cdr  Aug 19  : Cleanup. The routine MECKER is
cdr            now contained here.
cdr            tbd:
cdr            Similarly: ZERLEG, RUKSUB, OPRAND,
cdr            SUBTIT, SCHRIT, SIGNOK can also
cdr            be made local to this present module.
cdr
Cdr  Sept. 16: Bug fix: added option: two or more constants next to each other
c  summer 2019: added function eirene_replace(...)
C
C-----------------------------------------------------------------------
      SUBROUTINE EIRENE_ALGEBR (TERM,OPER,IZIF,CONST,NOP,IOUT)
C-----------------------------------------------------------------------
C
C     AUTOR:           ST. HUBER
C     IHK-KENNZIFFER:  121
C     DATUM:           25-NOV-1988
C
C     FUNKTION:
C
C     DAS PROGRAMM LIEST ARITHMETISCHE AUSDRUECKE EIN,
C     RUFT DAS UNTERPROGRAMM 'ZERLEG' AUF
C     UND GIBT ENTWEDER DIE EINZELNEN ZERLEGUNGEN ODER DIE
C     REGELVERLETZUNG AUS
C
C-----------------------------------------------------------------------
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
         INTEGER, INTENT(IN) :: IOUT
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
         INTEGER :: LAENGE
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
         INTEGER :: NR, ANFANG, ENDE, FELDIND, IK
         CHARACTER(12) :: ERSETZ(10)
         character(10) :: buchst
chr
C
C     HILFSVARIABLEN :
C
         INTEGER :: I, IC, IDIG
         logical :: prefix
         EXTERNAL :: EIRENE_UPPERCASE

         IPROUT = IOUT
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
            anfang=index(term,'<')
            do while (anfang.ne.0)
              nr=nr+1
              ende=index(term,'>')
chr           abspeichern der ersetzten operanden
chr           der i-te operand wird hierbei im i-ten element
chr           des feldes ersetz abgelegt und durch den i-ten
chr           buchstaben im alphabet ersetzt. gleichheit von
chr           operanden wird hierbei nicht beruecksichtigt
              ersetz(nr)=term(anfang:ende)
              if (anfang.gt.1) then
                hlfterm=term(1:anfang-1)//buchst(nr:nr)//buchst(nr:nr)//
     .              term(ende+1:)
              else
                hlfterm=buchst(nr:nr)//buchst(nr:nr)//
     .              term(ende+1:)
              endif
              term=hlfterm
              anfang=index(term,'<')
            end do
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
               DO I=1,TEIL
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
                  prefix = .false.
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
                        CALL EIRENE_READ_TAL_SPECS
     .                       (ERSETZ(FELDIND), IZIF(1,I), IZIF(2,I))
                      ENDIF
                    ENDIF
                  elseif (part(i)(7:7).ne.'Z') then
                     feldind=index(buchst,part(i)(7:7))
                     if (feldind == 0) then
!  prefix gefunden
                       oper(i) = part(i)(7:7)
                       ik = 0
                       prefix = .true.
                     else
                       OPER(I)=PART(I)(9:9)
                       IK=INDEX(ERSETZ(FELDIND),',')
                     end if
                     IF (IK.EQ.0) THEN
                       IC = IC + 1
                       IZIF(1,I)=-I
                       IZIF(2,I)= IC
                       if (prefix) then
                         const(ic) = 0.d0
                       else
                         CALL EIRENE_RDCN (ERSETZ(FELDIND),CONST(IC))
                       end if
                     ELSE
                       CALL EIRENE_READ_TAL_SPECS
     .                      (ERSETZ(FELDIND), IZIF(1,I), IZIF(2,I))
                     ENDIF

                     idig = 10
                     if (prefix) idig = 8

                     IF (PART(I)(idig:idig).NE.'Z') THEN
                       FELDIND=INDEX(BUCHST,PART(I)(idig:idig))
                       IK=INDEX(ERSETZ(FELDIND),',')
                       IF (IK.EQ.0) THEN
                         IC = IC + 1
                         IZIF(3,I)=-I
                         IZIF(4,I)=IC
                         CALL EIRENE_RDCN (ERSETZ(FELDIND),CONST(IC))
                       ELSE
                         CALL EIRENE_READ_TAL_SPECS
     .                        (ERSETZ(FELDIND), IZIF(3,I), IZIF(4,I))
                       ENDIF
                     ELSE
                        READ(PART(I)(12:12),'(I1)') IZIF(3,I)
                        IZIF(4,I)=0
                     ENDIF
                   ELSE
                     idig = 10
                     OPER(I)=PART(I)(idig:idig)
                     READ(PART(I)(idig-1:idig-1),'(I1)') IZIF(1,I)
                     IZIF(2,I)=0
                     if (part(i)(idig+1:idig+1).ne.'Z') then
                       FELDIND=INDEX(BUCHST,PART(I)(idig+1:idig+1))
                       IK=INDEX(ERSETZ(FELDIND),',')
                       IF (IK.EQ.0) THEN
                         IC = IC + 1
                         IZIF(3,I)=-I
                         IZIF(4,I)=IC
                         CALL EIRENE_RDCN (ERSETZ(FELDIND),CONST(IC))
                       ELSE
                         CALL EIRENE_READ_TAL_SPECS
     .                        (ERSETZ(FELDIND), IZIF(3,I), IZIF(4,I))
                       ENDIF
                     else
                       READ(PART(I)(idig+3:idig+3),'(I1)') IZIF(3,I)
                       IZIF(4,I)=0
                     endif
                  endif
               END DO  !  I
            ELSE
C
C              AUSGABE DER FEHLERMELDUNG
C
               WRITE(iprout,'(2A)') 'FOLGENDE REGELVERLETZUNG ',
     >                        'WURDE ERKANNT:'
               CALL EIRENE_MECKER(ERROR)
               NOP=0
            ENDIF
         ENDIF
C
C     ENDE VON ALGEBR
C
      RETURN

      END SUBROUTINE EIRENE_ALGEBR

*************************************************************************

      SUBROUTINE EIRENE_READ_TAL_SPECS (TAL, ZIF1, ZIF2)
      IMPLICIT NONE
      CHARACTER(12), INTENT(IN) :: TAL
      INTEGER, INTENT(OUT) :: ZIF1, ZIF2
      INTEGER :: IA, IE

      IA = INDEX(TAL,'<')
      IE = INDEX(TAL,'>')

      ZIF1 = 0
      ZIF2 = 0
      IF (IA*IE > 0) READ (TAL(IA+1:IE-1),*) ZIF1, ZIF2

      RETURN
      END SUBROUTINE EIRENE_READ_TAL_SPECS

*************************************************************************

cdr This routine writes "complaint messages"
cdr from ALGEBR, i.e. from the routines
cdr that try to decipher the coded algebraic
cdr expressions for the algebraic tallies
cdr ALGV, ALGS specified in input block 10C and 10E, resp.


C-----------------------------------------------------------------------
                SUBROUTINE EIRENE_MECKER(ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     AUSGABE DER REGELVERLETZUNGEN
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN


      IF     (ERROR .EQ.  1) THEN
        WRITE(iprout,*) 'DER AUSDRUCK ENTHAELT EIN UNGUELTIGES ZEICHEN.'
      ELSEIF (ERROR .EQ.  2) THEN
        WRITE(iprout,*)
     >  'ZWISCHEN ZWEI KLAMMERAUSDRUECKEN BEFINDET SICH',
     >  ' KEIN OPERATOR.'
      ELSEIF (ERROR .EQ.  3) THEN
         WRITE(iprout,*) 'EIN KLAMMERAUSDRUCK IST LEER.'
      ELSEIF (ERROR .EQ.  4) THEN
         WRITE(iprout,*) 'DER AUSDRUCK IST FALSCH GEKLAMMERT.'
      ELSEIF (ERROR .EQ.  5) THEN
         WRITE(iprout,*) 'DER AUSDRUCK ENTHAELT MEHR ALS 3 INEINANDER',
     >              'GESCHACHTELTE KLAMMERN.'
      ELSEIF (ERROR .EQ.  6) THEN
        WRITE(iprout,*)
     >  'EIN OPERAND BESTEHT AUS MEHR ALS ZWEI BUCHSTABEN.'
      ELSEIF (ERROR .EQ.  7) THEN
         WRITE(iprout,*) 'ZWEI OPERATOREN STEHEN NEBENEINANDER.'
      ELSEIF (ERROR .EQ.  8) THEN
         WRITE(iprout,*) 'NACH EINEM OPERATOR FOLGT EINE SCHLIESSENDE',
     >              ' KLAMMER.'
      ELSEIF (ERROR .EQ.  9) THEN
        WRITE(iprout,*)
     >  'DAS ERSTES ZEICHEN DES AUSDRUCKS IST OPERATOR, ',
     >  'UND KEIN PRAEFIX.'
      ELSEIF (ERROR .EQ. 10) THEN
         WRITE(iprout,*) 'DER AUSDRUCK ENDET MIT EINEM OPERATOR.'
      ELSEIF (ERROR .EQ. 11) THEN
         WRITE(iprout,*) 'DAS ERSTE ZEICHEN IN EINEM KLAMMERAUSDRUCK',
     >              ' IST EIN OPERATOR.'
      ELSEIF (ERROR .EQ. 12) THEN
         WRITE(iprout,*) 'DER AUSDRUCK ENTHAELT MEHR ALS 15 OPERATOREN.'
      ELSEIF (ERROR .EQ. 13) THEN
         WRITE(iprout,*) 'DIE (ANZAHL DER OPERANDEN)-1 IST UNGLEICH ',
     >              ' DER (ANZAHL DER OPERATOREN).'
      ENDIF
C
C     ENDE VON MECKER
C
      END SUBROUTINE EIRENE_MECKER

*************************************************************************

C-----------------------------------------------------------------------
              SUBROUTINE EIRENE_FEHLER(AUSDRU,AKTLEN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     DIESES UNTERPROGRAMM UEBERPRUEFT MOEGLICHE REGELVERLETZUNGEN
C     VON AUSDRU
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         INTEGER, PARAMETER :: MAXOPT=15
C           : MAXIMAL ZULAESSIGE ANZAHL VON OPERATOREN IN AUSDRU

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

C
C     AUSGABEPARAMETER :
C
         INTEGER, INTENT(OUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     LOKALE VARIABLEN :
C
         INTEGER :: OANDEN
C           : ANZAHL DER OPERANDEN IN AUSDRU

         INTEGER :: OTOREN
C           : ANZAHL DER OPERATOREN IN AUSDRU


      ERROR=0
C
C     UEBERPRUEFUNG AUF GUELTIGE ZEICHEN
C
      CALL EIRENE_SIGNOK(AUSDRU,AKTLEN,ERROR)
      IF (ERROR .EQ. 0) THEN
C
C        UEBERPRUEFUNG AUF KORREKTE KLAMMERUNG
C
         CALL EIRENE_KLAMME(AUSDRU,AKTLEN,ERROR)
         IF (ERROR .EQ. 0) THEN
C
C           UEBERPRUEFUNG AUF KORREKTHEIT DER OPERANDEN
C
            CALL EIRENE_OPRAND(AUSDRU,AKTLEN,OANDEN,ERROR)
            IF (ERROR .EQ. 0) THEN
C
C           UEBERPRUEFUNG AUF KORREKTHEIT DER OPERATOREN
C
               CALL EIRENE_OPERAT(AUSDRU,AKTLEN,OTOREN,ERROR)
               IF (ERROR .EQ. 0) THEN
                  IF (OTOREN .GT. MAXOPT) THEN
C
C                    AUSDRU ENTHAELT MEHR ALS DIE MAXIMAL ZULAESSIGE
C                    OPERATORENANZAHL
C
                     ERROR=12
                  ELSEIF (OTOREN .NE. OANDEN-1) THEN
C
C                    DIE ANZAHL DER OPERATOREN IST VERSCHIEDEN VON
C                    DER ANZAHL DER OPERANDEN-1
C
                     ERROR=13
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C
C     ENDE VON FEHLER
C
      END SUBROUTINE EIRENE_FEHLER

*************************************************************************


C-----------------------------------------------------------------------
             SUBROUTINE EIRENE_KLAMME(AUSDRU,AKTLEN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     UEBERPRUEFUNG AUF KORREKTE KLAMMERUNG
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         INTEGER, PARAMETER :: KLAMAX=3
C           : MAXIMAL ZULAESSIGE ANZAHL VON KLAMMERSCHACHTELUNGEN

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(IN) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     LOKALE VARIABLEN :
C
         INTEGER :: ANZAHL
C           : VARIABLE, DIE JE NACH IHREM WERT ANGIBT
C             = 0  ANZAHL DER OEFFNENDEN KLAMMERN IST GLEICH DER
C                  ANZAHL DER SCHLIESSENDEN KLAMMERN
C             > 0  ANZAHL DER OEFFNENDEN KLAMMERN IST GROESSER DER
C                  ANZAHL DER SCHLIESSENDEN KLAMMERN
C             < 0  ANZAHL DER OEFFNENDEN KLAMMERN IST KLEINER DER
C                  ANZAHL DER SCHLIESSENDEN KLAMMERN

         INTEGER :: MINKLA
C           : MINIMUM VON ANZAHL

         INTEGER :: MAXKLA
C           : MAXIMUM VON ANZAHL

         INTEGER :: KLAMON
C           : POSITION DER ZULETZ GEFUNDENEN OEFFNENDEN KLAMMER

         INTEGER :: KLAMOF
C           : POSITION DER ZULETZ GEFUNDENEN SCHLIESSENDEN KLAMMER

C
C     HILFSVARIABLE :
C
         INTEGER :: I


      ANZAHL=0
      MAXKLA=0
      MINKLA=0
      KLAMON=0
      KLAMOF=-1

      DO 10, I=1,AKTLEN
         IF (AUSDRU(I:I) .EQ. '(' .AND. ERROR .EQ. 0 ) THEN
            ANZAHL=ANZAHL+1
            MAXKLA=MAX(MAXKLA,ANZAHL)
            KLAMON=I

            IF (KLAMON .EQ. KLAMOF+1) THEN
C
C              ZWISCHEN ZWEI KLAMMERAUSDRUECKEN STEHT KEIN OPERATOR
C
               ERROR=2
            ENDIF
         ELSEIF (AUSDRU(I:I) .EQ. ')' .AND.  ERROR .EQ. 0 ) THEN
            ANZAHL=ANZAHL-1
            MINKLA=MIN(MINKLA,ANZAHL)
            KLAMOF=I

            IF (KLAMOF .EQ. KLAMON+1) THEN
C
C              DER KLAMMERAUSDRUCK IST LEER
C
               ERROR=3
            ENDIF
         ENDIF
   10 CONTINUE

      IF ((ANZAHL .NE. 0  .OR. MINKLA .LT. 0) .AND. ERROR .EQ. 0) THEN
C
C        FALSCHE KLAMMERSETZUNG
C
         ERROR=4
      ELSEIF (MAXKLA .GT. KLAMAX .AND. ERROR .EQ. 0) THEN
C
C        ES WURDEN MEHR ALS KLAMAX KLAMMERN GESCHACHTELT
C
         ERROR=5
      ENDIF
C
C     ENDE VON KLAMME
C
      END SUBROUTINE EIRENE_KLAMME

*************************************************************************


C-----------------------------------------------------------------------
          SUBROUTINE EIRENE_OPERAT(AUSDRU,AKTLEN,OTOREN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     UEBERPRUEFUNG EINER REGELVERLETZUNG BEI OPERATOREN
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION
C
         CHARACTER(5), PARAMETER :: FAKTOR='^*/+-'

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(IN) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     AUSGABEPARAMETER :
C
         INTEGER, INTENT(OUT) :: OTOREN
C           : ANZAHL DER OPERATOREN IN AUSDRU

C
C     HILFSVARIABLEN :
C
         INTEGER :: I, POS

      OTOREN=0

      I=1
C
C     REPEAT
C
   10 CONTINUE
         POS= INDEX (FAKTOR, AUSDRU(I:I) )
         IF (POS .GT. 0  .AND.  ERROR .EQ. 0 ) THEN
            IF (INDEX(FAKTOR, AUSDRU(I+1:I+1)) .GT. 0) THEN
C
C              2 OPERATOREN STEHEN NEBENEINANDER
C
               ERROR=7
            ELSEIF (INDEX(')',AUSDRU(I+1:I+1)) .GT. 0) THEN
C
C              NACH EINEM OPERATOR FOLGT EINE SCHLIESSENDE KLAMMER
C
               ERROR=8
            ENDIF
            OTOREN=OTOREN+1
         ENDIF
         I=I+1
      IF (I .LT. AKTLEN  .AND.  ERROR .EQ. 0 ) GOTO 10
C
C     UNTIL : BIS AUSDRU VOLLSTAENDIG DURCHLAUFEN
C             ODER FEHLER AUFGETRETEN
C
C     UEBERPRUEFEN, OB DER 1. OPERAND IN AUSDRU EIN PRAEFIX BESITZT
C
      POS=INDEX(FAKTOR,AUSDRU(1:1))
      IF (POS .GT. 0 ) THEN
         IF ( POS .LE. 3) THEN
C
C           AUSDRU BEGINN MIT EINEM OPERATOR
C
            ERROR=9
         ELSE
            OTOREN=OTOREN-1
         ENDIF
      ENDIF

      IF (INDEX(FAKTOR,AUSDRU(AKTLEN:AKTLEN)).GT. 0  .AND.
     >     ERROR .EQ. 0) THEN
C
C        AUSDRU ENDET MIT EINEM OPERATOR
C
         ERROR=10
      ENDIF

      DO 20, I=1,AKTLEN-1
C
C        DIE ENDGUELTIGE ANZAHL DER OPERATOREN ERGIBT SICH AUS
C        DER BISHERIGEN ANZAHL DER OPERATOREN ABZUEGLICH DER
C        ANZAHL DER PRAEFIXE
C
         IF (AUSDRU(I:I) .EQ. '(' .AND.
     >        INDEX(FAKTOR,AUSDRU(I+1:I+1)) .GT. 0 ) THEN
            IF (INDEX(FAKTOR, AUSDRU(I+1:I+1) ) .GT. 3) THEN
               OTOREN=OTOREN-1
            ELSE
C
C              EINER GEOEFFNETEN KLAMMER FOLGT EIN OPERATOR
C
               ERROR=11
            ENDIF
         ENDIF
   20 CONTINUE
C
C     ENDE VON OPERAT
C
      END SUBROUTINE EIRENE_OPERAT

*************************************************************************


C-----------------------------------------------------------------------
           SUBROUTINE EIRENE_OPRAND(AUSDRU,AKTLEN,OANDEN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     UEBERPRUEFEN AUF ZULAESSIGE OPERANDEN
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         CHARACTER(26), PARAMETER ::
     P                  BUCHST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     AUSGABEPARAMETER :
C
         INTEGER, INTENT(OUT) :: OANDEN
C           : ANZAHL DER OPERANDEN IN AUSDRU

C
C     HILFSVARIABLEN :
C
         INTEGER :: I, POS


      OANDEN=0
      AUSDRU(AKTLEN+1:AKTLEN+2)='  '

      I=1
C
C     WHILE-1 : SOLANGE AUSDRU NOCH NICHT VOLLSTAENDIG DURCHLAUFEN
C               UND KEIN FEHLER AUFGETRETEN
C
   20 IF (I .LE. AKTLEN  .AND.  ERROR .EQ. 0 ) THEN
         POS= INDEX (BUCHST, AUSDRU(I:I) )
C
C        WHILE-2 : SOLANGE ANFANG VON EINEM OPERANDEN GEFUNDEN
C                  UND AUSDRU NOCH NICHT VOLLSTAENDIG DURCHLAUFEN
C                  UND KEIN FEHLER AUFGETRETEN
C
   10    IF (POS .GT. 0  .AND. I .LE. AKTLEN .AND. ERROR .EQ. 0 ) THEN
            IF (INDEX(BUCHST, AUSDRU(I+1:I+1)) .GT. 0  .AND.
     >          INDEX(BUCHST, AUSDRU(I+2:I+2)) .GT. 0 ) THEN
C
C              OPERAND BESTEHT AUS MEHR ALS 2 BUCHSTABEN
C               => ZU LANG
C
               ERROR=6
            ELSE
C
C              GEHE WEITER IM STRING
C
               if (AUSDRU(I:I) == AUSDRU(I+1:I+1)) OANDEN=OANDEN+1
               I=I+2
               POS=INDEX( BUCHST, AUSDRU(I:I) )
            ENDIF
            GOTO 10
         ENDIF
C
C        ENDWHILE-2
C
C        KEIN ANFANG VON EINEM OPERANDEN GEFUNDEN, GEHE WEITER IM STRING
C
         I=I+1
         GOTO 20
      ENDIF
C
C     ENDWHILE-1
C
C     ENDE VON OPRAND
C
      END SUBROUTINE EIRENE_OPRAND

*************************************************************************


C-----------------------------------------------------------------------
        SUBROUTINE EIRENE_RUKSUB(TEIL, IPART, PART, IARITH, ARITH,
     .  HILFE)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     EINFUEGEN VON 'Z' FUER ZERLEGUNG
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         CHARACTER(10), PARAMETER :: ZIFFER='0123456789'

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: TEIL
C           : AKTUELLE ANZAHL DER ZERLEGUNGEN

C
C     EIN/AUSGABEPARAMETER :
C
         CHARACTER(*), INTENT(INOUT) :: PART(*)
C           : FELD VON STRINGS, AUF DENEN DIE EINZELNEN
C             ELEMENTARZERLEGUNGEN FESTGEHALTEN WERDEN

         INTEGER, INTENT(INOUT) :: IPART(*)
C           : AKTUELLE LAENGEN VON PART(ZMAX)

         CHARACTER(*), INTENT(INOUT) :: ARITH(*)
C           : FELD VON STRINGS, AUF DENEN DIE TEIL-TE GENERATION
C             VON AUSDRU FESTGEHALTEN WIRD

         INTEGER, INTENT(INOUT) :: IARITH(*)
C           : AKTUELLE LAENGEN VON ARITH(ZMAX)

         CHARACTER(*), INTENT(INOUT) :: HILFE
C           : HILFSSTRING, DER FUER ZUWEISUNGEN BENOETIGT WIRD

C
C     HILFSVARIABLEN :
C
         INTEGER :: I, J

      DO J=1,TEIL
C
C        EINFUEGEN  VON 'Z' FUER ZERLEGUNG IN PART(J)
C
         I=1
         DO WHILE (I .LT. IPART(J))
            IF (INDEX(ZIFFER,PART(J)(I:I)) .GT. 0) THEN
               IF (I .EQ. 1) THEN
                  HILFE='Z'//PART(J)(I:IPART(J))//' '
               ELSE
                  HILFE=PART(J)(1:I-1)//'Z'//
     >                   PART(J)(I:IPART(J))//' '
               ENDIF

               I=I+2
               IPART(J)=IPART(J)+1
               PART(J)=HILFE
            ENDIF
            I=I+1
         END DO
C
C        EINFUEGEN VON 'Z' FUER ZERLEGUNG IN ARITH(J)
C
         I=1
         DO WHILE (I .LT. IARITH(J))
            IF (INDEX(ZIFFER,ARITH(J)(I:I)) .GT. 0) THEN
               IF (I .EQ. 1) THEN
                  HILFE='Z'//ARITH(J)(I:IARITH(J))//' '
               ELSE
                  HILFE=ARITH(J)(1:I-1)//'Z'//
     >                   ARITH(J)(I:IARITH(J))//' '
               ENDIF

               I=I+2
               IARITH(J)=IARITH(J)+1
               ARITH(J)=HILFE
            ENDIF
            I=I+1
         END DO

      END DO  !  J
C
C     ENDE VON RUKSUB
C
      END SUBROUTINE EIRENE_RUKSUB

*************************************************************************


C-----------------------------------------------------------------------
      SUBROUTINE EIRENE_SCHRIT(AUSDRU, AKTLEN, ALPHA, OMEGA, TEIL,
     >                  IPART, PART, IARITH, ARITH)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     ZERLEGUNG EINES TEILAUSDRUCKS,
C     DIESER AUSDRUCK ENTHAELT KEINE KLAMMERN
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
C     KONSTANTENDEKLARATION :
C
         CHARACTER(5), PARAMETER :: FAKTOR='^*/+-'
C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

         INTEGER, INTENT(INOUT) :: TEIL
C           : AKTUELLE ANZAHL DER ZERLEGUNGEN

         INTEGER, INTENT(INOUT) :: ALPHA
C           : POSITION DES ERSTEN ZEICHENS VOM TEILAUSDRUCK

         INTEGER, INTENT(INOUT) :: OMEGA
C           : POSITION DES LETZTEN ZEICHENS VOM TEILAUSDRUCK

         CHARACTER(*), INTENT(INOUT) :: PART(*)
C           : FELD VON STRINGS, AUF DENEN DIE EINZELNEN
C             ELEMENTARZERLEGUNGEN FESTGEHALTEN WERDEN

         INTEGER, INTENT(INOUT) :: IPART(*)
C           : AKTUELLE LAENGEN VON PART(ZMAX)

         CHARACTER(*), INTENT(INOUT) :: ARITH(*)
C           : FELD VON STRINGS, AUF DENEN DIE TEIL-TE GENERATION
C             VON AUSDRU FESTGEHALTEN WIRD

         INTEGER, INTENT(INOUT) :: IARITH(*)
C           : AKTUELLE LAENGEN VON ARITH(ZMAX)

C
C     LOKALE VARIABLEN :
C
         INTEGER :: BEGINN
C           : POSITION IN AUSDRU, BEI DER DIE ZERLEGUNG BEGINNT

         INTEGER :: ENDE
C           : POSITION IN AUSDRU, BEI DER DIE ZERLEGUNG BEENDET WIRD

         INTEGER :: OPER1
C           : POSITION DES 1. OPERATORS IN AUSDRU

         INTEGER :: OPER2
C           : POSITION DES 2. OPERATORS IN AUSDRU

         INTEGER :: RANK1
C           : RANG DES 1. OPERATORS

         INTEGER :: RANK2
C           : RANG DES 2. OPERATORS

         LOGICAL :: ERFOLG
C           : GIBT AN, OB EINE ZERLEGUNG ERFOLGEN DARF

         LOGICAL :: PREFIX
C           : GIBT AN, OB 1.OPERAND IN AUSDRUCK EIN PRAEFIX BESITZT

         CHARACTER(2) :: TEILCH
C           : INHALT VON TEIL

C
C     HILFSVARIABLEN :
C
         INTEGER :: I, POS, IS

      BEGINN=ALPHA
      ENDE=OMEGA
      ERFOLG=.FALSE.
      I=BEGINN-1
C
C     SUCHE NACH DEM 1. OPERATOR IM TEILSTRING
C
      IS = SCAN(AUSDRU(BEGINN:OMEGA),FAKTOR)
      I = BEGINN + IS - 1
      IF (IS == 0) I = OMEGA
      POS=INDEX(FAKTOR, ausdru(I:I))

      PREFIX = (I .EQ. BEGINN)
C
C     WHILE : SOLANGE TEILSTRING EINEN OPERATOR ENTHAELT
C
      DO WHILE (POS .GT. 0)
C
C       POSITION UND RANGBESTIMMUNG DES OPERATORS
C       RANK = 0 : ^
C            = 1 : * ODER /
C            = 2 : + ODER -
C
        OPER1=I
        RANK1=POS/2

C
C       REPEAT-1
C
        DO WHILE (.NOT.ERFOLG)
C
C         SUCHE NACH DEM 2. OPERATOR IM TEILSTRING
C         REPEAT-2
C
          IS = SCAN(AUSDRU(I+1:OMEGA),FAKTOR)
          I = I + IS
          IF (IS == 0) I = OMEGA
          POS=INDEX(FAKTOR, ausdru(I:I))
C
C           UNTIL-2 : BIS TEILSTRING VOLLSTAENDIG DURCHLAUFEN ODER
C                     WEITEREN OPERATOR IM TEILSTRING GEFUNDEN
C
          IF (POS .GT. 0) THEN
C
C              WEITEREN OPERATOR IN TEILSTRING GEFUNDEN,
C              POSITION UND RANGBESTIMMUNG DES OPERATORS
C
            OPER2=I
            RANK2=POS/2
C
C              UEBERPRUEFEN, OB NOCH EIN OPERATOR IM TEILSTRING
C              GESUCHT WERDEN MUSS, ODER EINE ZERLEGUNG ERFOLGEN KANN
C
            IF (RANK2 .LT. RANK1 .AND. .NOT. PREFIX) THEN
              BEGINN=OPER1+1
              OPER1=OPER2
              RANK1=RANK2
            ELSEIF (RANK2 .EQ. RANK1  .AND.  RANK1 .EQ. 0
     >              .AND. .NOT. PREFIX) THEN
              BEGINN=OPER1+1
              OPER1=OPER2
              RANK1=RANK2
            ELSE
              ENDE=OPER2-1
              ERFOLG=.TRUE.
C                 zerlegung von beginn bis ende
            ENDIF

          ELSE
C
C              ES IST KEIN WEITERER OPERATOR IN TEILSTRING VORHANDEN,
C              DIE POSITION DES LETZTEN OPERATORS IST OPER1
C
            ERFOLG=.TRUE.
          ENDIF
        END DO
C
C       UNTIL-1 : BIS EINE ZERLEGUNG ERFOLGEN DARF
C
        TEIL=TEIL+1
C
C       UMWANDELN DER ZAHL TEIL IN EINE ZEICHENKETTE TEILCH
C
        WRITE(TEILCH,'(I2.2)') TEIL
C
C       ABSPEICHERN DER TEIL-TEN ELEMENTARZERLEGUNG
C
        IPART(TEIL)= 5 + ENDE+1-BEGINN
        PART(TEIL)=TEILCH //' = '//AUSDRU(BEGINN:ENDE) // ' '
C
C       ABSPEICHERN DER TEIL-TEN GENERATION VON AUSDRU
C
        IARITH(TEIL)=AKTLEN - (ENDE+1-BEGINN) +2
        IF (BEGINN .EQ. 1) then
          if (ENDE .EQ. AKTLEN) THEN
            ARITH(TEIL)=TEILCH
          ELSE
            ARITH(TEIL)=TEILCH//AUSDRU(ENDE+1:AKTLEN)
          ENDIF
        ELSEIF (ENDE .EQ. AKTLEN) THEN
          ARITH(TEIL)=AUSDRU(1:BEGINN-1)//TEILCH
        ELSE
          ARITH(TEIL)=AUSDRU(1:BEGINN-1)//TEILCH
     >                   //AUSDRU(ENDE+1:AKTLEN)
        ENDIF
        AKTLEN=IARITH(TEIL)
        AUSDRU=ARITH(TEIL)
C
C       WERTE ZURUECKSETZEN
C
        OMEGA= OMEGA - (ENDE+1-BEGINN) +2
        ERFOLG=.FALSE.
        BEGINN=alpha
        ENDE=omega
        I=BEGINN-1
        PREFIX=.FALSE.
C
C       SUCHE NACH DEM 1. OPERATOR IM TEILSTRING
C

        IS = SCAN(AUSDRU(I+1:OMEGA),FAKTOR)
        I = I + IS
        IF (IS == 0) I = OMEGA
        POS=INDEX(FAKTOR, ausdru(I:I))

      END DO
C
C     ENDWHILE
C     ZERLEGUNG BEENDET, DA TEILSTRING KEIN OPERATOR MEHR ENTHAELT
C
C     ENDE VON SCHRIT
C
      END SUBROUTINE EIRENE_SCHRIT

*************************************************************************


C-----------------------------------------------------------------------
            SUBROUTINE EIRENE_SIGNOK(AUSDRU,AKTLEN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     DAS UNTERPROGRAMM UEBERPRUEFT AUSDRU AUF ZULAESSIGE ZEICHEN
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         CHARACTER(33), PARAMETER ::
     P                  GULTIG='ABCDEFGHIJKLMNOPQRSTUVWXYZ+-*/^()'

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(IN) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN


      IF (VERIFY(AUSDRU(1:AKTLEN),GULTIG) >= 1) THEN
C
C        UNGUELTIGES ZEICHEN IN AUSDRU GEFUNDEN
C
         ERROR = 1
       END IF
C
C     ENDE VON SIGNOK
C
      END SUBROUTINE EIRENE_SIGNOK

*************************************************************************


C-----------------------------------------------------------------------
               SUBROUTINE EIRENE_SUBTIT(AUSDRU, AKTLEN, HILFE)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     ELIMINATION VON BLANKS UND SUBTITUTION VON '**' DURCH '^'
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     EIN/AUSGABEPARAMETER :
C
      INTEGER, INTENT(INOUT) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

      CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

      CHARACTER(*), INTENT(INOUT) :: HILFE
C           : HILFSSTRING, DER FUER ZUWEISUNGEN BENOETIGT WIRD

C
C     HILFSVARIABLEN :
C
      INTEGER :: I

C
C     ELIMINATION VON BLANKS
C
      I=1
C
C     WHILE-1 : SOLANGE AUSDRU NOCH NICHT VOLLSTAENDIG DURCHLAUFEN
C
      DO WHILE (I .LT. AKTLEN)
C
C        ENTFERNEN VON BLANKS
C
         IF (AUSDRU(I:I) .EQ. ' ') THEN
            IF ( I .EQ. 1) THEN
               HILFE=AUSDRU(I+1:AKTLEN) // ' '
            ELSE
               HILFE=AUSDRU(1:I-1) // AUSDRU(I+1:AKTLEN) // ' '
            ENDIF
            AUSDRU=HILFE
            AKTLEN=AKTLEN-1
         ELSE
            I=I+1
         ENDIF
      END DO
C
C     ENDWHILE-1
C
C
C     SUBTITUTION VON '**' DURCH '^'
C
      CALL EIRENE_REPLACE(AUSDRU,'**', '^' ,AKTLEN)
      CALL EIRENE_REPLACE(AUSDRU,'DX', 'QA',AKTLEN)
      CALL EIRENE_REPLACE(AUSDRU,'DY', 'QB',AKTLEN)
      CALL EIRENE_REPLACE(AUSDRU,'DZ', 'QC',AKTLEN)
      CALL EIRENE_REPLACE(AUSDRU,'EXP','QD',AKTLEN)
      CALL EIRENE_REPLACE(AUSDRU,'LN', 'QE',AKTLEN)
      CALL EIRENE_REPLACE(AUSDRU,'LOG','QF',AKTLEN)
      CALL EIRENE_REPLACE(AUSDRU,'ABS','QG',AKTLEN)
C
C     ENDE VON SUBTIT
C
      CONTAINS

      SUBROUTINE EIRENE_REPLACE (STR, REP, WITH, LENGTH)
      IMPLICIT NONE
      CHARACTER(*), INTENT(INOUT) :: STR
      CHARACTER(*), INTENT(IN) :: REP, WITH
      INTEGER, INTENT(INOUT) :: LENGTH
      INTEGER :: LREP, LWITH, POS
#ifdef F2003
      CHARACTER(:), ALLOCATABLE :: CREP, CWITH, HILFE
#else
      INTEGER :: N
      INTEGER, PARAMETER :: STRMAX = 512
      CHARACTER(LEN=STRMAX) :: CREP, CWITH, HILFE
#endif

#ifndef F2003
      N = LEN(ADJUSTL(TRIM(REP)))
      IF (N.GT.STRMAX) THEN
        WRITE(IPROUT,*)
     .   'INCREASE SIZE OF STRMAX IN SUBTIT.EIRENE_REPLACE'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
      N = LEN(ADJUSTL(TRIM(WITH)))
      IF (N.GT.STRMAX) THEN
        WRITE(IPROUT,*)
     .   'INCREASE SIZE OF STRMAX IN SUBTIT.EIRENE_REPLACE'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
#endif
      CREP = ADJUSTL(TRIM(REP))
      CWITH = ADJUSTL(TRIM(WITH))
      LREP = LEN_TRIM(CREP)
      LWITH = LEN_TRIM(CWITH)

      POS=INDEX(STR,CREP(1:LREP))

      DO WHILE ( POS .GE. 1 .AND. POS .LT. LENGTH-1)

#ifndef F2003
         N = LWITH+LENGTH-LREP+1
         IF (N.GT.STRMAX) THEN
           WRITE(IPROUT,*)
     .      'INCREASE SIZE OF STRMAX IN SUBTIT.EIRENE_REPLACE'
           CALL EIRENE_EXIT_OWN(1)
         ENDIF
#endif
         IF (POS == 1) THEN
           HILFE = CWITH(1:LWITH) // STR(POS+LREP:LENGTH)
         ELSE
           HILFE=STR(1:POS-1) // CWITH(1:LWITH) // STR(POS+LREP:LENGTH)
         END IF
         STR=HILFE
         LENGTH=LENGTH - (LREP - LWITH)
         POS=INDEX(STR,CREP(1:LREP))
      END DO

      END SUBROUTINE EIRENE_REPLACE

      END SUBROUTINE EIRENE_SUBTIT

*************************************************************************


CDR  decompose (interpret) expression AUSDRU into
cdr  single algebraic operations on tallies.

C-----------------------------------------------------------------------
      SUBROUTINE EIRENE_ZERLEG(AUSDRU, AKTLEN,
     .                         IPART, PART, IARITH,
     .                         ARITH, TEIL, HILFE, ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     DAS UNTERPROGRAMM UBERPRUEFT DEN AUSDRUCK AUF REGELVERLETZUNG
C     UND ZERLEGT IHN IN EINZELOPERATIONEN, SO DASS DER AUSDRUCK
C     BEI SEQUENTIELLER ABARBEITUNG RICHTIG AUSGEWERTET WIRD
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     EIN/AUSGABEPARAMETER :
C
      CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

      INTEGER, INTENT(INOUT) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

C
C     AUSGABEPARAMETER :
C
      INTEGER, INTENT(OUT) :: TEIL
C           : AKTUELLE ANZAHL DER ZERLEGUNGEN

      CHARACTER(*), INTENT(OUT) :: PART(*)
C           : FELD VON STRINGS, AUF DENEN DIE EINZELNEN
C             ELEMENTARZERLEGUNGEN FESTGEHALTEN WERDEN

      INTEGER, INTENT(OUT) :: IPART(*)
C           : AKTUELLE LAENGEN VON PART(ZMAX)

      CHARACTER(*), INTENT(OUT) :: ARITH(*)
C           : FELD VON STRINGS, AUF DENEN DIE TEIL-TE GENERATION
C             VON AUSDRU FESTGEHALTEN WIRD

      INTEGER, INTENT(OUT) :: IARITH(*)
C           : AKTUELLE LAENGEN VON ARITH(ZMAX)

      CHARACTER(*), INTENT(OUT) :: HILFE
C           : HILFSSTRING, DER FUER ZUWEISUNGEN BENOETIGT WIRD

      INTEGER, INTENT(OUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     LOKALE VARIABLEN :
C
      INTEGER :: ALPHA
C           : POSITION DES ERSTEN ZEICHENS VOM TEILAUSDRUCK

      INTEGER :: OMEGA
C           : POSITION DES LEZTEN ZEICHENS VOM TEILAUSDRUCK

         INTEGER :: BEGINN
C           : POSITION IN AUSDRU, BEI DER DIE ZERLEGUNG BEGINNT

         INTEGER :: ENDE
C           : POSITION IN AUSDRU, BEI DER DIE ZERLEGUNG BEENDET WIRD

         CHARACTER(6), PARAMETER :: FUCHAR= 'ABCDEF'

         CHARACTER(2) :: TEILCH
C           : INHALT VON TEIL
C
C     HILFSVARIABLEN :
C
      INTEGER :: POS
      LOGICAL :: CHECK


C
C     ELIMINATION VON BLANKS UND ERSETZEN VON '**' DURCH '^'
C
      CALL EIRENE_SUBTIT(AUSDRU, AKTLEN, HILFE)
C
C     UEBERPRUEFUNG VERSCHIEDENER REGELVERLETZUNGEN
C
      CALL EIRENE_FEHLER(AUSDRU, AKTLEN, ERROR)
      IF (ERROR .EQ. 0) THEN
C
C        BEGINN DER ZERLEGUNG
C
         TEIL=0

         ALPHA=1
         OMEGA=AKTLEN
         POS=INDEX( AUSDRU(1:AKTLEN), ')' )
C
C        WHILE : SOLANGE NOCH KLAMMERN IN VARIABLE "AUSDRU" VORHANDEN
C
!   22    IF (POS .GT. 0) THEN
         DO WHILE (POS .GT. 0)
C
C           BESTIMME LETZTES ZEICHEN IM ERSTEN INNERSTEN KLAMMERAUSDRUCK
C
            OMEGA=POS-1
C
C           BESTIMME ERSTES ZEICHEN IM ERSTEN INNERSTEN KLAMMERAUSDRUCK
C
            POS=0
            POS=INDEX( AUSDRU(POS+1:OMEGA), '(', .true. ) +POS
            ALPHA=POS+1
C
C           ZERLEGUNG DES KLAMMERAUSDRUCKS
C
            CALL EIRENE_SCHRIT(AUSDRU, AKTLEN, ALPHA, OMEGA, TEIL,
     >                  IPART, PART, IARITH, ARITH)
C
C           UEBERPRUEFUNG, OB DIE KLAMMERN WEGFALLEN KOENNEN
C
            IF (OMEGA+1-ALPHA .LE. 2
     >          .AND. INDEX('+-',AUSDRU(ALPHA:ALPHA)) .EQ. 0) THEN

!              PRUEFE, OB EINE FUNKTION AUSGEWERTET WERDEN SOLL
               CHECK = .FALSE.
               IF (ALPHA > 3) THEN
                 CHECK = (AUSDRU(ALPHA-3:ALPHA-3) == 'Q') .AND.
     >             (INDEX(FUCHAR,AUSDRU(ALPHA-2:ALPHA-2)) > 0)
               END IF
               IF (CHECK) THEN

!                 FUNKTIONSAUFRUF
                  TEIL = TEIL + 1
                  WRITE(TEILCH,'(I2.2)') TEIL
                  BEGINN = ALPHA - 3
                  ENDE = OMEGA + 1
                  IPART(TEIL) = 5 + ENDE - BEGINN + 1
                  PART(TEIL) = TEILCH//' = '//AUSDRU(BEGINN:ENDE)//' '
                  IARITH(TEIL)=AKTLEN - (ENDE+1-BEGINN) +2
                  IF (BEGINN .EQ. 1) then
                     if (ENDE .EQ. AKTLEN) THEN
                        ARITH(TEIL)=TEILCH
                     ELSE
                        ARITH(TEIL)=TEILCH//AUSDRU(ENDE+1:AKTLEN)
                     ENDIF
                  ELSEIF (ENDE .EQ. AKTLEN) THEN
                     ARITH(TEIL)=AUSDRU(1:BEGINN-1)//TEILCH
                  ELSE
                     ARITH(TEIL)=AUSDRU(1:BEGINN-1)//TEILCH
     >                    //AUSDRU(ENDE+1:AKTLEN)
                  ENDIF
                  AKTLEN=IARITH(TEIL)
                  AUSDRU=ARITH(TEIL)

               ELSE
C
C                 ENTFERNUNG DER KLAMMERN
C
                  IF (ALPHA-1 .EQ. 1) THEN
                     IF (OMEGA+1 .EQ. AKTLEN) THEN
                       HILFE=AUSDRU(ALPHA:OMEGA)
                     ELSE
                       HILFE=AUSDRU(ALPHA:OMEGA)//AUSDRU(OMEGA+2:AKTLEN)
                     ENDIF
                  ELSEIF (OMEGA+1 .EQ. AKTLEN) THEN
                     HILFE=AUSDRU(1:ALPHA-2)//AUSDRU(ALPHA:OMEGA)
                  ELSE
                     HILFE=AUSDRU(1:ALPHA-2)//AUSDRU(ALPHA:OMEGA)
     >                    //AUSDRU(OMEGA+2:AKTLEN)
                  ENDIF
                  AKTLEN=AKTLEN-2
                  AUSDRU=HILFE
                  IARITH(TEIL)=AKTLEN
                  ARITH(TEIL)=AUSDRU
               END IF
            ENDIF

            ALPHA=1
            OMEGA=AKTLEN
            POS=INDEX( AUSDRU(1:AKTLEN), ')' )
!            GOTO 22
         END DO
C
C        ENDWHILE
C
C        AUSDRUCK ENTHAELT KEINE KLAMMERAUSDRUECKE MEHR
C          => ZERLEGUNG DES GESAMMTEN AUSDRUCKS
C
         CALL EIRENE_SCHRIT(AUSDRU, AKTLEN, ALPHA, OMEGA, TEIL,
     >               IPART, PART, IARITH, ARITH)

         IF (TEIL .EQ. 0) THEN
C
C           ES SIND KEINE ZERLEGUNGEN GEMACHT WORDEN, DA
C           DER AUSDRUCK NUR AUS EINEM OPERAND BESTEHT
C
            TEIL = 1
            PART(1) = 'Z01 = ' // AUSDRU
            IPART(1) = 6 + AKTLEN
            ARITH(1) = 'Z01'
            IARITH(1) = 3
         ELSE
C
C           ERSETZEN VON '^' DURCH '**' UND EINFUEGEN VON 'Z' FUER
C           ZERLEGUNG
C
            CALL EIRENE_RUKSUB(TEIL, IPART, PART, IARITH, ARITH, HILFE)
         ENDIF

      ENDIF
C
C     ENDE VON ZERLEG
C
      END SUBROUTINE EIRENE_ZERLEG

*************************************************************************

C
C
      SUBROUTINE EIRENE_RDCN (ERSETZ,CONST)
      IMPLICIT NONE
C
      CHARACTER(*), INTENT(IN) :: ERSETZ
      REAL(DP), INTENT(OUT) :: CONST
      INTEGER :: IEXPO, IW, IPUNKT, ILEN, ICON
      CHARACTER(10) :: FORM
C
      ILEN=INDEX(ERSETZ,'>')-2
      FORM='          '
C
C   INTEGER?
C
      IF (INDEX(ERSETZ,'.').EQ.0) THEN
        FORM(1:5)='(I  )'
        IF (ILEN.GT.9) WRITE (FORM(3:4),'(I2)') ILEN
        IF (ILEN.LE.9) WRITE (FORM(3:3),'(I1)') ILEN
        READ (ERSETZ(2:ILEN+1),FORM) ICON
        CONST=FLOAT(ICON)
C
      ELSE
C
C   REAL CONSTANT
C
        IPUNKT=INDEX(ERSETZ,'.')-1
        IEXPO=INDEX(ERSETZ,'E')
        IF (IEXPO.EQ.0) IEXPO=INDEX(ERSETZ,'D')
C
        IF (IEXPO.EQ.0) THEN
C   F-FORMAT
          FORM(1:2)='(F'
          IF (ILEN.GT.9) THEN
            WRITE (FORM(3:4),'(I2)') ILEN
            IW=4
          ELSE
            WRITE (FORM(3:3),'(I1)') ILEN
            IW=3
          ENDIF
          FORM(IW+1:IW+3)='. )'
          WRITE (FORM(IW+2:IW+2),'(I1)') ILEN-IPUNKT
          READ (ERSETZ(2:ILEN+1),FORM) CONST
C
        ELSE
C
C    E-FORMAT
          IEXPO=IEXPO-1
          FORM(1:2)='(E'
          IF (ILEN.GT.9) THEN
            WRITE (FORM(3:4),'(I2)') ILEN
            IW=4
          ELSE
            WRITE (FORM(3:3),'(I1)') ILEN
            IW=3
          ENDIF
          FORM(IW+1:IW+3)='. )'
          WRITE (FORM(IW+2:IW+2),'(I1)') IEXPO-(IPUNKT+1)
          READ (ERSETZ(2:ILEN+1),FORM) CONST
        ENDIF
      ENDIF
C
      RETURN
      END SUBROUTINE EIRENE_RDCN

*************************************************************************

      end module eirmod_algebra
