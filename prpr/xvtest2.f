      PROGRAM XVTEST2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TESTER LE GRAPHISME XVUE EN MODE OBJET 2D
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS  SEPTEMBRE 1995
C2345X7..............................................................012
C     DECLARATIONS NECESSAIRES POUR LES COULEURS ET FONTES
      include"./incl/trvari.inc"
      include"./incl/xvfontes.inc"
      COMMON /UNITES/  LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION DMCN
      COMMON           DMCN(100)
C
      PARAMETER       (NBSOM=6)
      REAL             XY(2,NBSOM)
      DATA             XY / -2.0,0.0,  -1.0,-1.0,   1.0,-1.0,
     %                       2.0,0.0,   1.0, 1.0,  -1.0, 1.0 /
C
      LECTEU = 5
      IMPRIM = 6
C
      PRINT *
      PRINT *,'=========================================='
      PRINT *,'Test du mode OBJET 2D de la librairie XVUE'
      PRINT *,'=========================================='
C
C     CREATION OUVERTURE DE LA FENETRE XVUE
C     =====================================
      CALL XVOUVRIR
C
C     DEMANDE DE TRACE POSTSCRIPT
C     ===========================
      CALL XVINITIERPS( 0 )
C
C     LA MISE A L'ECHELLE
C     ===================
      XMIN = -2.0
      XMAX =  2.0
      YMIN = -2.0
      YMAX =  1.0
      CALL ISOFENETRE( XMIN, XMAX, YMIN, YMAX )
C
C     TRACES DE 4 TRIANGLES
C     =====================
      CALL XVEPAISSEUR( 5 )
      CALL TRIANGLE( NCROUG, NCMAGE, XY, 1, 2, 6 )
      CALL TRIANGLE( NCBLAN, NCMAGE, XY, 2, 3, 5 )
      CALL TRIANGLE( NCBLAN, NCMAGE, XY, 2, 5, 6 )
      CALL TRIANGLE( NCBLEU, NCMAGE, XY, 3, 4, 5 )
C
C     TRACE DU SEGMENT XY(1) XY(4)
C     ============================
      CALL TRAIT2D( NCJAUN, XY(1,1), XY(2,1), XY(1,4), XY(2,4) )
C
C     TRACES DE TEXTES
C     ================
C     CHARGEMENT DE LA FONTE DE HAUTEUR LA PLUS PROCHE DE 15 PIXELS
C     ET DE PREFERENCE GRASSE ET 'courier'
      CALL CHOIXFONTE( 15 )
C     LE NUMERO DE LA FONTE CHARGEE EST NOFONT (cf ~/nef/incl/xvfontes.inc)
      CALL XVEPAISSEUR( 1 )
      X = -1.0
      Y = -1.2
      CALL TRAIT2D( NCBLAN, X-0.2, Y, X, Y )
      CALL TEXTE2D( NCBLAN, X, Y, 'TEXTE2D' )
      Y = -1.4
      CALL TRAIT2D(   NCVERT, X-0.2, Y, X, Y )
      CALL SYMBOLE2D( NCVERT, X, Y, 'SYMBOLE2D' )
      Y = -1.6
      CALL TRAIT2D(  NCCYAN, X-0.2, Y, X, Y )
      CALL ENTIER2D( NCCYAN, X, Y, 123 )
      Y = -1.8
      CALL TRAIT2D( NCROUG, X-0.2, Y, X, Y )
      CALL REEL2D(  NCROUG, X, Y, 3.14159, '(G13.5)' )
C
C     PASSAGE EN MODE PIXEL
C     =====================
      NPX = NUPXEX(  1.0 )
      NPY = NUPXEY( -1.2 )
      CALL XVTYPETRAIT( 1 )
      CALL XVTRAIT( NPX-100, NPY, NPX, NPY )
      CALL XVTEXTE( 'XVTEXTE', 8, NPX, NPY )
C
C     TRACE EFFECTIF
      CALL MEMPXFENETRE
      CALL XVVOIR
C
C     ARRET AVANT LA FIN
C     ==================
 100  PRINT *
      PRINT *,'ENTRER UN CARACTERE POUR FINIR'
      CALL XVSOURIS( NOTYEV, NC, NPX, NPY )
      IF( NOTYEV .NE. 2 ) GOTO 100
C
C     SAUVEGARDE DU TRACE DANS UN FICHIER POSTSCRIPT
C     ==============================================
      CALL XVSAUVERPS( 'XVUEtest2', 9 )
C
C     ENVOI SUR L'IMPRIMANTE DU FICHIER POSTSCRIPT
C     ============================================
      CALL XVIMPRIMERPS( 'XVUEtest2', 9 )
C
C     FERMETURE DE XVUE
C     =================
      CALL XVFERMER
      END
      SUBROUTINE TRIANGLE( NCFACE, NCARET, XY, NS1, NS2, NS3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE LE TRIANGLE DE SOMMETS (NS1, NS2, NS3) AVEC LA COULEUR
C -----  NCFACE ET LES ARETES SELON LA COULEUR NCARET
C
C ENTREES:
C --------
C NCFACE : NUMERO DE LA COULEUR DU TRIANGLE (SI <0 PAS DE TRACE)
C NCARET : NUMERO DE LA COULEUR DES 3 ARETES DU TRIANGLE
C          (SI <0 PAS DE TRACE)
C XY     :  LES 2 COORDONNEES DES SOMMETS DU TABLEAU XY
C NS1,NS2,NS3 : LES NUMEROS DANS XY DES 3 SOMMETS DU TRIANGLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS  SEPTEMBRE 1995
C2345X7..............................................................012
      REAL  XY(2,*), T(3,2)
C
C     REGROUPEMENT DES COORDONNEES DES 3 SOMMETS DU TRIANGLE
      DO 10 I=1,2
         T(1,I) = XY(I,NS1)
         T(2,I) = XY(I,NS2)
         T(3,I) = XY(I,NS3)
 10   CONTINUE
C
C     LE TRACE DU TRIANGLE ET DE SES ARETES EN MODE OBJET 2D
      CALL FACE2D( NCFACE, NCARET, 3, T(1,1), T(1,2) )
      END
