      SUBROUTINE TNOUVR( KNOM , KTYPE , NOTAMS , MCTAMS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : OUVRIR LE LEXIQUE OU TABLEAU MS OU TABLEAU DE TABLEAUX
C ----- DEFINI PAR SON NOM KNOM
C
C ENTREES :
C ---------
C KNOM   : LE NOM A OUVRIR.IL EST DEFINI PAR PLUSIEURS MOTS SEPARES
C          PAR LE CARACTERE > . PAR EXEMPLE  ADAM>LIGNE
C
C ATTENTION : SI KNOM(1:2) EST DIFFERENT DE ~> ALORS ~> EST AJOUTE
C
C SORTIES :
C ---------
C KTYPE  : TYPE DU TABLEAU TS 'TAMS' OU 'LEXI' OU 'TATA' OU 'RELA'
C NOTAMS : NUMERO DU TABLEAU MS DEFINI PAR KNOM
C          0 S'IL N'EXISTE PAS
C MCTAMS : SON ADRESSE DANS MCN OU MCK SELON SON TYPE CARACTERES OU NON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C.......................................................................
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / MSSFTA /  MSSF(28),NTADAM
      COMMON / UNITES /  LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*(*)      KNOM
      CHARACTER*100      NOMLX
      CHARACTER*80       KNOM1,KNOM2
C                        TABLEAUX NECESSAIRES SUITE A PB DE CONCATENATION
      CHARACTER*4        KTYPE
C
      CHARACTER*4        KAUX
      CHARACTER*1        SEPARE
      DATA               SEPARE/ '>' /
C
C     PROTECTION DE KNOM DANS NOMLX CHAINE LOCALE
      NOMLX = KNOM
C
      KAUX = '~>'
      N1 = 0
      L  = LEN( KNOM )
      IF( NOMLX(1:1) .EQ. '~' ) THEN
         N1 = 1
         IF( NOMLX(2:2) .EQ. SEPARE ) THEN
            N1 = 2
         ENDIF
      ELSE IF( NOMLX(1:1) .EQ. SEPARE ) THEN
         N1 = 1
      ENDIF
C
C     CALCUL INTERMEDIAIRE POUR EVITER LE PROBLEME DE LA CONCATENATION
      KAUX  = '~>' // NOMLX(N1+1:N1+2)
      KNOM1 = KAUX // NOMLX(N1+3:L)
      NOMLX = KNOM1
C
C     SI 'LINE' EST DANS CE TEXTE, IL EST REMPLACE PAR 'LIGNE' POUR ETRE TROUVE!
      N2 = INDEX( NOMLX, 'LINE' )
      IF( N2 .GT. 0 ) THEN
         print *
         print *,'tnouvr: n2=',n2,' NO 1-ER CARACTERE DE LINE in NOMLX'
         KNOM2 = NOMLX(1:N2-1) // 'LIGNE' // NOMLX(N2+4:L)
         print *,'tnouvr: knom2=',knom2
         NOMLX  = KNOM2
C        NOMLX A UN CARACTERE DE PLUS!
         L = L + 1
         print *,'tnouvr: l=',L,'  nomlx=',NOMLX
         print *,'tnouvr: On ne doit pas passer par la'
         CALL XVPAUSE
      ENDIF
C
C     LE LEXIQUE DE DEPART EST ADAM
      NT = NTADAM
      N1 = INDEX( NOMLX , SEPARE )
      IF( NOMLX(3:3) .EQ. ' ' ) THEN
C        PAS DE NOM DERRIERE > DONC LEXIQUE ADAM
         NOTAMS = NTADAM
         CALL TAMSOU( NTADAM , MCTAMS )
         KTYPE = 'LEXI'
         RETURN
      ENDIF
C
C     LA 1-ERE ET DERNIERE LETTRE APRES >
      N1 = N1 + 1
C
C     BOUCLE SUR LES DIFFERENTS MOTS DE NOMLX SEPARES PAR '>'
C     RECHERCHE DU 1-ER CARACTERE > RENCONTRE
 10   N2 = INDEX( NOMLX(N1:L) , SEPARE )
C
      IF( N2 .GT. 0 ) THEN
C        N1:N2-1 EST LE NOM DU LEXIQUE INTERMEDIAIRE
C        OUVERTURE DE CE LEXIQUE A PARTIR DU NOM
         IF( NT .LE. 0 ) GOTO 9000
         N2 = N1 + N2 - 1
         CALL LXLXOU( NT , NOMLX(N1:N2-1) , NOTAMS , MCTAMS )
         NT = NOTAMS
         N1 = N2 + 1
         GOTO 10
      ENDIF
C
C     NOMLX(N1:L) EST LE NOM TERMINAL
C     REPERAGE DE "...
C     OUVERTURE DE CE NOM
      IF( NT .LE. 0 ) GOTO 9000
      CALL LXNMOU( NT , NOMLX(N1:L),KTYPE,NOTAMS,MCTAMS )
      RETURN
C
C     ERREUR  TABLEAU OU LEXIQUE INEXISTANT
 9000 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TNOUVR: ABSENCE DE '//NOMLX(1:N1-2)
      ELSE
         KERR(1) = 'TNOUVR: UNKNOWN '//NOMLX(1:N1-2)
      ENDIF
      CALL LEREUR
      NOTAMS = 0
      MCTAMS = 0
C
      RETURN
      END
