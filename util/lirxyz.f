      SUBROUTINE LIRXYZ( NCVALS , XYZ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE LES 3 COORDONNEES D'UN POINT
C -----
C
C       LA DEMANDE EST TOUJOURS SATISFAITE CAR PAS DE RETOUR TANT QUE
C       LA VALEUR DE XYZ EST INCORRECTE SAUF @ DEMANDANT L'ABANDON
C
C CF $MEFISTO/td/g/grammaire_lu DE DEFINITION DU LANGAGE UTILISATEUR
C
C ENTREE ET SORTIE :
C ------------------
C NCVALS : EN ENTREE : 12 LES XYZ SONT INITIALISES EN ENTREE
C                       0 SINON
C          EN SORTIE :  1 LES 3 VALEURS XYZ SONT INITIALISEES
C                      -1 ABANDON DE LA LECTURE PAR ENTREE DE @
C
C SORTIE :
C --------
C XYZ    : LA VALEUR DES 3 REELS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1990
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      DOUBLE PRECISION  D
      REAL              XYZ(3)
      INTRINSIC         REAL
C
CCC      CHARACTER*1       KCHAIN
CCC      CHARACTER*1       KXYZ(3)
CCC      DATA              KXYZ / 'X' , 'Y' , 'Z' /
C
C     DEMANDE DE LECTURE D'UN REEL DOUBLE PRECISION
      IF( NCVALS .NE. 0 ) THEN
         NOTYP = 5
      ELSE
         NOTYP = 0
      ENDIF

CCC      IF( INTERA .GE. 3 ) THEN
CCCC
CCCC        LE TEXTE EN BOUT DE LA DERNIERE LIGNE
CCC         IF( NBLGRC(NRINVI) .LT. MXLGIN ) THEN
CCC            NBLGRC(NRINVI) = NBLGRC(NRINVI) + 1
CCC         ENDIF
CCC         IF( NUIDEN .GT. 0 .AND. NUIDEN .LE. MXIDEN ) THEN
CCCC           AJOUT DE L'IDENTIFICATEUR DE VALEUR DEMANDEE
CCC            I      = NUDCNB( KIDENT(NUIDEN) )
CCC            KINVI(NBLGRC(NRINVI)) = 'X DE ' //
CCC     %                               KIDENT(NUIDEN)(1:I) // '=?'
CCC         ELSE
CCC            KINVI(NBLGRC(NRINVI)) = KXYZ(1) // '=?'
CCC         ENDIF
CCC         MDLGRC(NRINVI) = MXCANB( NBLGRC(NRINVI) , KINVI )
CCC         MDRECT(NRINVI) = MDLGRC(NRINVI)
CCC      ENDIF
C
      DO 100 NN=1,3
 10      IF( INTERA .GE. 3 ) THEN
C
C           TRAITEMENT DE L'INVITE
            IF( NN .EQ. 1 ) THEN
               CALL INVITE( 98 )
            ELSE IF( NN .EQ. 2 ) THEN
               CALL INVITE( 99 )
            ELSE
               CALL INVITE( 100 )
            ENDIF
CCC            CALL LEMENU
         ENDIF
C
         CALL DONNMF( NOTYP , 5 , NCVALS , D , KCHAIN )
         IF( NCVALS .EQ. -1 ) GOTO 9000
         IF( NCVALS .NE.  1 ) GOTO 10
         XYZ(NN) = REAL( D )

 100  CONTINUE
C
C     INVITE EFFACEE
 9000 NBLGRC( NRINVI ) = 0
C
CCC 8000  CALL CLAVIS
      RETURN
      END
