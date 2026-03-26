      SUBROUTINE LIRLEX( NTLX , NCVALS , NOM , NVALEN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE UN NOM DANS UN LEXIQUE
C -----
C
C       LA DEMANDE EST TOUJOURS SATISFAITE CAR PAS DE RETOUR TANT QUE
C       LA VALEUR DU NOM EST INCORRECTE SAUF @ DEMANDANT L'ABANDON
C
C CF $MEFISTO/td/g/grammaire_lu DE DEFINITION DU LANGAGE UTILISATEUR
C
C ENTREE :
C --------
C NTLX   : NUMERO TMS DU LEXIQUE DANS LEQUEL LE NOM DOIT EXISTER
C
C ENTREE ET SORTIE :
C ------------------
C NCVALS : EN ENTREE : 11 LE NOM EST INITIALISE EN ENTREE
C                       0 SINON
C          EN SORTIE : 11 LE NOM EST RETROUVE DE NUMERO NVALEN
C                      -1 ABANDON DE LA LECTURE PAR ENTREE DE @
C
C SORTIE :
C --------
C NOM    : LE NOM LU
C NVALEN : LE NUMERO DU NOM DANS LE LEXIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1990
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      COMMON / LECLXQ / LECLEX
      DOUBLE PRECISION  D
      CHARACTER*(*)     NOM
      CHARACTER*14      NOMOBJ(1:7),NOMOBA(7)
      DATA              NOMOBJ/'POINT'    , 'LIGNE'  ,
     %                         'SURFACE'  , 'VOLUME' , 'OBJET' ,
     %                         'FONCTION' , 'TRANSFORMATION' /
      DATA              NOMOBA/'POINT'    , 'LINE'  ,
     %                         'SURFACE'  , 'VOLUME' , 'OBJECT' ,
     %                         'FUNCTION' , 'TRANSFORMATION' /
C
C     DEMANDE DE LECTURE D'UN NOM
      IF( NCVALS .NE. 0 ) THEN
         NOTYP = 2
      ELSE
         NOTYP = 0
      ENDIF
C
      IF( INTERA .GE. 3 .AND. LHLECT .EQ. 1 ) THEN
C
C        LECTURE A L'AIDE DE LA SOURIS
C        SOIT PAR CLIC SUR L'IMAGE DU CLAVIER
C        SOIT PAR CLIC SUR L'OBJET
C        ------------------------------------
C        L'INVITE EST EVENTUELLEMENT COMPLETEE
         IF( NBLGRC(NRINVI) .LT. MXLGIN ) THEN
            NBLGRC(NRINVI) = NBLGRC(NRINVI) + 1
         ENDIF
C
C        RECHERCHE DU TYPE DE L'OBJET DU LEXIQUE
         DO 20 I=1,7
            IF( NTLX .EQ. NTMN(I) ) GOTO 30
 20      CONTINUE
C        LEXIQUE DIFFERENT DE POINTS,...,OBJETS
         IF( LANGAG .EQ. 0 ) THEN
            KINVI(NBLGRC(NRINVI)) = 'LE NOM'
         ELSE
            KINVI(NBLGRC(NRINVI)) = 'The NAME'
         ENDIF
         GOTO 100
C
C        AJOUT DE L'IDENTIFICATEUR DE VALEUR DEMANDEE
 30      LECLEX = I
         IF( LANGAG .EQ. 0 ) THEN
            KINVI(NBLGRC(NRINVI)) = 'NOM de ' // NOMOBJ(LECLEX)
         ELSE
            KINVI(NBLGRC(NRINVI)) = 'NAME of ' // NOMOBA(LECLEX)
         ENDIF
C
 100     MDLGRC(NRINVI) = MXCANB( NBLGRC(NRINVI) , KINVI )
         MDRECT(NRINVI) = MDLGRC(NRINVI)
      ENDIF
C
C     LECTURE OU INTERPRETATION DIRECTE
 10   CALL DONNMF( NOTYP , 2 , NCVALS , D , NOM )
      IF( NCVALS .EQ. -1 ) GOTO 9000
      IF( NCVALS .NE.  2 ) GOTO 10
C
C     OUVERTURE VERIFICATION DE L'EXISTENCE DU NOM DANS LE LEXIQUE
      CALL LXNMNO( NTLX , NOM , NVALEN , MNLX )
      IF( NVALEN .LE. 0 ) THEN
C        NOM INCORRECT NON RETROUVE DANS LE LEXIQUE
         IF ( LHLECT.GT.1 .AND. LHLECT.LE.MXLECT ) THEN
C           LECTURE ACTUELLE SUR UN FICHIER AVEC UNE ERREUR
C           AFFICHAGE DU RETOUR AU CLAVIER
            CALL QUIFIC
         ELSE
            NBLGRC(NRERR) = 0
         ENDIF
         I=INDEX( NOM, ' ' )
         IF( I .LE. 0 ) I=1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(NBLGRC(NRERR)+1) = 'NOM INCONNU: ' // NOM(1:I)
            KERR(NBLGRC(NRERR)+2) = 'A choisir PARMI:'
         ELSE
            KERR(NBLGRC(NRERR)+1) = 'UNKNOWN NAME: ' // NOM(1:I)
            KERR(NBLGRC(NRERR)+2) = 'MUST be CHOSEN AMONG:'
         ENDIF
         NBLGRC(NRERR) = NBLGRC( NRERR ) + 2
         CALL LXIM0( MNLX )
         GOTO 10
      ENDIF
      NCVALS = 11
C
C     INVITE EFFACEE
 9000 NBLGRC( NRINVI ) = 0
      LECLEX = 0
C
CCC 8000  CALL CLAVIS
      RETURN
      END
