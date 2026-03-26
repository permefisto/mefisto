      SUBROUTINE PTPLPR( NUTYOB, NUOBJ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TROUVER LE PLSVO CLIQUE SUR L'ECRAN GRAPHIQUE
C -----

C ENTREE :
C --------
C NUTYOB : NUMERO DU TYPE D'OBJET PLSVO 1:POINT, 2:LIGNE, ...

C SORTIE :
C --------
C NUOBJ  : >0 NUMERO DE L'OBJET DANS SON LEXIQUE
C          =0 SI LE TABLEAU DES ITEMS EST VIDE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C23456---------------------------------------------------------------012
      PARAMETER     (DIST0=1E28, RAYON=100)
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"

C     SAISIE D'UN POINT PAR CLIC DE LA SOURIS
C     ---------------------------------------
 10   CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
      IF( NOTYEV .LE. 0 ) GOTO 10

C     OBJET POSSIBLE DE TYPE D'OBJET
      DISTMI = DIST0
      NUOBJ  = 0

C     ADRESSE DE DEBUT DES ITEMS DE CES OBJETS PRESENTS SUR L'ECRAN
      MNIT = MNITEM( NUTYOB )
      IF( MNIT .LE. 0 ) RETURN
      IF( MCN(MNIT+2) .LE. 0 ) RETURN

C     LE TABLEAU D'ITEMS N'EST PAS VIDE
      MOTS = MCN( MNIT )

C     MNIT+2 POINTEUR SUR LE NOMBRE ACTUEL D'ITEMS
      DO 20 I=1,MCN(MNIT+2)

C        LES COORDONNEES ECRAN DE L'OBJET
         MNIT = MNIT + MOTS
         DIST = (NX-MCN(MNIT)) ** 2  + (NY-MCN(MNIT+1)) ** 2
         IF( DIST .GT. RAYON**2 ) GOTO 20
C        CLIC DANS LE CERCLE DE RAYON 100 PIXELS
         IF( DIST .LT. DISTMI ) THEN
C           OBJET PLUS PROCHE
            NUOBJ  = MCN(MNIT+2)
            DISTMI = DIST
         ENDIF

 20   ENDDO

      IF( NUOBJ .LE. 0 ) GOTO 10

      RETURN
      END
