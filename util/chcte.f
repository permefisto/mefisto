      SUBROUTINE CHCTE( NLD , NCD , NLF , NCF ,
     %                  NCODEV , DBLVAL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE A PARTIR DU CARACTERE (NLD,NCD) D'UNE CONSTANTE REELLE
C ----- CE PEUT NE PAS ETRE UNE TELLE CONSTANTE
C
C ENTREES :
C ---------
C NLD,NCD : POSITION DANS KLG DU PREMIER CARACTERE A TRAITER
C
C SORTIES :
C ---------
C NLF,NCF : SI NCODEV=1 POSITION DANS KLG DU DERNIER CARACTERE
C                       DE LA CONSTANTE REELLE
C           SINON NLF=0
C NCODEV : 0 DBLVAL N'EST PAS INITIALISEE
C          1 DBLVAL EST INITIALISEE
C DBLVAL : VALEUR REELLE DOUBLE PRECISION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  DBLVAL,D
      CHARACTER*1       CAR
C
C     LE PREMIER CARACTERE DOIT ETRE UN CHIFFRE
      CAR = KLG(NLD)(NCD:NCD)
      IF( CAR .LT. '0' .AND. CAR .GT. '9' ) THEN
C        ERREUR: LA CONSTANTE NE DEBUTE PAS PAR UN CHIFFRE
         GOTO 9000
      ENDIF
C
C     LECTURE DES CHIFFRES . e E d D + - CHIFFRES
      NL     = NLD
      NC     = NCD
      NBCD   = -1
      NBCE   = -1
      ISIGN  = 0
      D      = 0D0
      LEXPOS = 0
      LINITD = 0
C
 10   CAR    = KLG(NL)(NC:NC)
      LUNITE = LECHIF( CAR )
      IF( LUNITE .GE. 0 ) THEN
C        UN CHIFFRE DE PLUS
         D      = D * 10 + LUNITE
C        LA CONSTANTE PEUT EXISTER
         LINITD = 1
C        NBCD EST LE NOMBRE DE CHIFFRE APRES LE . ET AVANT E
         IF( NBCD .GE. 0 ) NBCD = NBCD + 1
C        PASSAGE AU CARACTERE SUIVANT
         NC = NC + 1
         IF( NC .GT. NBCALI ) GOTO 1000
         GOTO 10
      ELSE IF( CAR .EQ. '.' ) THEN
C        LE POINT DECIMAL EST TROUVE
         IF( NBCD .GE. 0 ) THEN
C           . RETROUVE 2 FOIS
            GOTO 9000
         ENDIF
         NBCD = 0
C        PASSAGE AU CARACTERE SUIVANT
         NC   = NC + 1
         IF( NC .GT. NBCALI ) GOTO 1000
         GOTO 10
      ELSE IF( (NC .GT. NCD) .AND.
     %         (CAR .EQ. 'E' .OR. CAR .EQ. 'e' .OR.
     %          CAR .EQ. 'D' .OR. CAR .EQ. 'd')    ) THEN
C        LE CARACTERE E D e d N'EST PAS LE PREMIER DE LA CONSTANTE
         NCAVEX = NC
         NBCE   = 0
C        RECHERCHE DE L'EXPOSANT
         NC     = NC + 1
 20      CAR    = KLG(NL)(NC:NC)
         LUNITE = LECHIF( CAR )
         IF( LUNITE .GE. 0 ) THEN
C           UN CHIFFRE DE PLUS
            LEXPOS = LEXPOS * 10 + LUNITE
            NBCE   = NBCE + 1
C           PASSAGE AU CARACTERE SUIVANT
            NC = NC + 1
            IF( NC .GT. NBCALI ) GOTO 1000
            GOTO 20
         ELSE IF( CAR .EQ. '-' ) THEN
            IF( ISIGN .EQ. 0 ) THEN
               ISIGN = -1
C              PASSAGE AU CARACTERE SUIVANT
               NC = NC + 1
               IF( NC .GT. NBCALI ) GOTO 1000
               GOTO 20
            ELSE
C              +- INTERDIT
               GOTO 9000
            ENDIF
         ELSE IF( CAR .EQ. '+' ) THEN
            IF( ISIGN .EQ. 0 ) THEN
               ISIGN = 1
C              PASSAGE AU CARACTERE SUIVANT
               NC = NC + 1
               IF( NC .GT. NBCALI ) GOTO 1000
               GOTO 20
            ELSE
C              +- INTERDIT
               GOTO 9000
            ENDIF
         ELSE
            IF( NC .EQ. NCAVEX+1 ) THEN
C              ERREUR E NON SUIVI D'UN CHIFFRE OU - OU +
               GOTO 9000
            ENDIF
         ENDIF
      ELSE
         IF( LINITD .GT. 0 ) THEN
C           FIN DE LA CONSTANTE REELLE
            GOTO 1000
         ELSE
C           CE N'EST PAS UNE CONSTANTE REELLE
            GOTO 9000
         ENDIF
      ENDIF
C
C     FIN DE LA CONSTANTE
 1000 NLF = NL
      NCF = NC-1
C
C     LE SIGNE NUL EST UNE ABSENCE DE SIGNE DONC +
      IF( ISIGN .EQ. 0 ) ISIGN = 1
C
C     S'IL EXISTE E SUIVI D'AUCUN CHIFFRE => ERREUR
      IF( NBCE .EQ. 0 ) GOTO 9000
C
C     LE NOMBRE DE CHIFFRES APRES . ET AVANT E EST NBCD
C     LA VALEUR DE LA CONSTANTE REELLE
      IF( NBCD .GT. 0 ) THEN
C        REEL AVEC UNE PARTIE DECIMALE
         IF( LEXPOS .GT. 0 ) THEN
C           REEL AVEC PARTIE DECIMALE ET EXPOSANT
            DBLVAL = D * ( 0.1D0 ** NBCD )*( 10D0 ** (ISIGN*LEXPOS) )
         ELSE
C           REEL AVEC PARTIE DECIMALE ET SANS EXPOSANT
            DBLVAL = D * ( 0.1D0 ** NBCD )
         ENDIF
      ELSE
C        REEL SANS PARTIE DECIMALE
         IF( LEXPOS .GT. 0 ) THEN
C           REEL SANS PARTIE DECIMALE ET AVEC EXPOSANT
            DBLVAL = D * ( 10D0 ** (ISIGN*LEXPOS) )
         ELSE
C           REEL SANS PARTIE DECIMALE ET SANS EXPOSANT
            DBLVAL = D
         ENDIF
      ENDIF
      NCODEV = 1
      RETURN
C
C     ERREUR : LA CHAINE N'EST PAS UNE CONSTANTE
 9000 NLF    = 0
      NCODEV = 0
      END
