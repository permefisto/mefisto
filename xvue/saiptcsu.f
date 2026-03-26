      SUBROUTINE SAIPTCSU( NOTYEV, NOPXX, NOPXY, NOCHAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SI NOTYEV>0, RESTITUER LES COORDONNEES PIXELS DU CLIC SOURIS
C ----- SI NOTYEV<0, RESTITUER LE CARACTERE TAPE AU CLAVIER
C                    ET SI LE CARACTERE EST ECHAPPEMENT METTRE @
C                    POUR UNE EQUIVALENCE DE L'ABANDON
C       EN SURLIGNANT LA LIGNE DU MENU OU SE TROUVE LE MENU

C SORTIES :
C ---------
C NOTYEV  : NUMERO DU TYPE DE L'EVENEMENT
C           1 ou 2 ou 3 => NUMERO DU BOUTON DE LA SOURIS ET
C                     NOPXX, NOPXY SONT INITIALISES AU CLIC DE LA SOURIS
Cccc 17/12/2020  non
Cccc       LE BOUTON 2 DE LA SOURIS => -1 => FRAPPE CLAVIER DE @ ABANDON
Cccc 17/12/2020  non

C          -1 => CARACTERE TAPE AU CLAVIER AVEC NOCHAR
C           0 => ABANDON

C NOPXX, NOPXY : SI NOTYEV>0 LES COORDONNEES PIXELS DU POINT CLIQUE
C NOCHAR       : SI NOTYEV<0 NUMERO ASCII DU CARACTERE TAPE AU CLAVIER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN LJLL UPMC & St PIERRE DU PERRAY    AVRIL 2009
C MODIFS : PERRONNET ALAIN          Saint PIERRE DU PERRAY Decembre 2020
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      INTRINSIC      INT

      NUMERC0 = 0

C     ATTENDRE LA FRAPPE D'UNE TOUCHE DE LA SOURIS OU DU CLAVIER
C     ==========================================================
 10   CALL XVSOURIS( NOTYEV, NOCHAR, NOPXX, NOPXY )
C     RETOURNE NOTYEV
C        = 0 Si ABANDON demande par frappe de la touche Echappement ou @
C        = 1 Si CLIC ENFONCE et RELACHE D'UN BOUTON DE LA SOURIS
C               => NOPXX, NOPXY et NOCHAR=NUMERO DU BOUTON
C        =-1 Si CLIC SEULEMENT ENFONCE  D'UN BOUTON DE LA SOURIS
C               => NOPXX, NOPXY et NOCHAR=NUMERO DU BOUTON
C        =-2 Si le pointeur de la souris a bouge
C               => NOPXX, NOPXY
C        = 2 Si FRAPPE D'UN CARACTERE AU CLAVIER
C               => NOCHAR=NUMERO DU CARACTERE DANS LA TABLE ASCII

      IF( NOTYEV .EQ. 0 ) THEN

C        ABANDON DEMANDE
         RETURN

      ELSE IF( NOTYEV .LT. 0 ) THEN

C        SOURIS BOUGEE
         IF( NBLGRC(NRMENU) .GT. 0 ) THEN
C           LE POINTEUR EST IL SUR UNE DES LIGNES DU RECTANGLE MENU?
            IF( (NOPXX.GE.XRECT(NRMENU)) .AND.
     &          (NOPXX.LE.XRECT(NRMENU)+DXRECT(NRMENU)) .AND.
     &          (NOPXY.GE.YRECT(NRMENU)) .AND.
     &          (NOPXY.LE.YRECT(NRMENU)+DYRECT(NRMENU)) ) THEN
C              CLIC DANS LE MENU : NUMERO DE LA CASE ?
               NUMERC = INT( (NOPXY-YRECT(NRMENU))/DYLGRC(NRMENU)+ 1.0 )
C              AFFICHAGE DU MENU AVEC LA CASE NUMERC SURLIGNEE EN JAUNE
               IF( NUMERC .NE. NUMERC0 ) THEN
                  CALL RECTTXSU( NRMENU, KMENU, 2, NAMENU, NUMERC )
                  NUMERC0 = NUMERC
               ENDIF
            ENDIF
C           CLIC EN DEHORS DES LIMITES DU RECTANGLE DU MENU
         ENDIF

         GOTO 10

      ELSE IF( NOTYEV .EQ. 1 ) THEN

C        UN DES 3 BOUTONS DE LA SOURIS

ccc   17/12/2020   --------------------------------------------------------
ccc         IF( NOCHAR .EQ. 2 ) THEN
cccC           BOUTON 2 DE LA SOURIS <=> FRAPPE AU CLAVIER DU CARACTERE '@'
ccc            NOTYEV = -1
ccc            NOCHAR = 64
ccc         ELSE
cccC           BOUTON 1 OU 3 => PAS DE MODIFICATION
ccc            NOTYEV = NOCHAR
ccc         ENDIF
ccc   17/12/2020   --------------------------------------------------------

C        BOUTON 1 ou 2 ou 3 => PAS DE MODIFICATION
         NOTYEV = NOCHAR
CCC
CCC         PRINT *,'SAIPTCSU: EVENEMENT SOURIS : BOUTON ',NOTYEV,
CCC     %           ' POSITION DU CLIC X=',NOPXX,' Y=',NOPXY

      ELSE IF( NOTYEV .EQ. 2 ) THEN

C        FRAPPE D'UN CARACTERE AU CLAVIER
         NOTYEV = -1

CCC         PRINT *,'SAIPTCSU: EVENEMENT CLAVIER : NO ASCII CARACTERE ',
CCC     %            NOCHAR,' => CARACTERE=',CHAR(NOCHAR)

         IF( NOCHAR .EQ. 27 ) THEN

C           LE CARACTERE 'ECHAPPEMENT' DEVIENT LE CARACTERE '@'
            NOCHAR = 64

         ENDIF

      ELSE

         GOTO 10

      ENDIF
      END
