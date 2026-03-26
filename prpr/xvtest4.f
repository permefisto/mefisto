      PROGRAM XVTEST4
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TESTER LA FERMETURE DE LA FENETRE MEFISTO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS    JANVIER 1999
C2345X7..............................................................012
      COMMON /UNITES/  LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION DMCN
      COMMON           DMCN(1000)
C
      REAL             AXPTV(3), AXEIL(3)
      REAL             XYZR(3,4)
C
      PARAMETER       (NBSOM=8)
      REAL             XYZ(3,NBSOM)
      DATA             XYZ / 0,0,0, 1,0,0, 1,1,0, 0,1,0,
     %                       0,0,1, 1,0,1, 1,1,1, 0,1,1 /
C
      DATA             XYZR/ 0,0,0, 1.4,0,0, 0,1.4,0, 0,0,1.4 /
C
      LECTEU = 5
      IMPRIM = 6
C
      PRINT *
      PRINT *,'======================================================='
      PRINT *,'Test de la FERMETURE de la fenetre graphique de MEFISTO'
      PRINT *,'======================================================='
C
C     CREATION OUVERTURE DE LA FENETRE XVUE
C     =====================================
      CALL XVOUVRIR
C
C     CHARGEMENT D'UNE GROSSE FONTE (33 PIXELS DE HAUTEUR)
C     =============================
      CALL CHOIXFONTE( 33 )
C
C     LA DEFINITION DE L'AXONOMETRIE
C     ==============================
      AXPTV(1) = 0.0
      AXPTV(2) = 0.0
      AXPTV(3) = 0.0
      AXEIL(1) = 2.0
      AXEIL(2) = 2.0
      AXEIL(3) = 2.0
      AXLAR    = 1.4
      AXHAU    = 1.4
      AXARR    = 0.0
      AXAVA    = 0.0
      CALL AXONOMETRIE( AXPTV, AXEIL, AXLAR, AXHAU, AXARR, AXAVA )
C
C     TRACES DES TRIANGLES, FACES DE L'OBJET 3D
C     =========================================
      CALL XVTYPETRAIT( 1 )
      CALL XVEPAISSEUR( 3 )
      AXPTV(1) = 0.5
      AXPTV(2) = 0.5
      AXPTV(3) = 0.5
      CALL PTVLONLAT( AXPTV, 25.0, 35.0 )
C
      DO 10 I=1,320
         CALL OBJET3D( XYZR, XYZ )
 10   CONTINUE
C
C     TRACE EFFECTIF
      CALL MEMPXFENETRE
      CALL XVVOIR
 100  PRINT *
      PRINT *,'ENTRER UN CARACTERE POUR FINIR'
      CALL XVSOURIS( NOTYEV, NC, NPX, NPY )
      IF( NOTYEV .NE. 2 ) GOTO 100
C
C     FERMETURE DE XV
C     ================
      PRINT *,'APPEL XVFERMER'
      CALL XVFERMER
      PRINT *,'APRES XVFERMER'
      STOP 'FIN NORMALE'
      END
      SUBROUTINE TRIANGL3D( NCFACE, NCARET, XYZ, NS1, NS2, NS3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE LE TRIANGLE 3D DE SOMMETS (NS1, NS2, NS3) AVEC LA COULEUR
C -----  NCFACE ET LES ARETES SELON LA COULEUR NCARET
C
C ENTREES:
C --------
C NCFACE : NUMERO DE LA COULEUR DU TRIANGLE 3D (SI <0 PAS DE TRACE)
C NCARET : NUMERO DE LA COULEUR DES 3 ARETES DU TRIANGLE 3D
C          (SI <0 PAS DE TRACE)
C XYZ     :  LES 3 COORDONNEES DES SOMMETS DU TABLEAU XYZ
C NS1,NS2,NS3 : LES NUMEROS DANS XYZ DES 3 SOMMETS DU TRIANGLE 3D
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS  SEPTEMBRE 1995
C2345X7..............................................................012
      REAL  XYZ(3,*), T(3,3)
C
C     REGROUPEMENT DES COORDONNEES DES 3 SOMMETS DU TRIANGLE 3D
      DO 10 I=1,3
         T(I,1) = XYZ(I,NS1)
         T(I,2) = XYZ(I,NS2)
         T(I,3) = XYZ(I,NS3)
 10   CONTINUE
C
C     LE TRACE DU TRIANGLE 3D ET DE SES ARETES EN MODE OBJET 3D
      CALL FACE3D( NCFACE, NCARET, 3, T )
      END
      SUBROUTINE OBJET3D( XYZR, XYZ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE LE REPERE ET L'OBJET 3D A PARTIR DES 8 SOMMETS DU CUBE
C -----
C ENTREES:
C --------
C XYZ    :  LES 3 COORDONNEES DES 8 SOMMETS DU CUBE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS  SEPTEMBRE 1995
C2345X7..............................................................012
      include"./incl/trvari.inc"
      REAL  XYZR(3,4), XYZ(3,8)
C
C     TRACES DES FACES DE L'OBJET 3D
      CALL TRIANGL3D( NCROUG, NCMAGE,   XYZ, 4, 8, 5 )
      CALL TRIANGL3D( NCBLAN, NCMAGE,   XYZ, 2, 1, 6 )
      CALL TRIANGL3D( NCBLAN, NCMAGE,   XYZ, 1, 5, 6 )
      CALL TRIANGL3D( NCBLAN, NCMAGE,   XYZ, 1, 4, 5 )
      CALL TRIANGL3D( NCBLAN, -1,       XYZ, 3, 4, 5 )
      CALL TRIANGL3D( NCBLAN, -1,       XYZ, 3, 5, 6 )
      CALL TRIANGL3D( NCBLAN,  NCMAGE,  XYZ, 2, 3, 6 )
      CALL TRIANGL3D( NCBLEU,  NCMAGE,  XYZ, 3, 7, 6 )
      CALL TRIANGL3D( NCBLEU,  NCMAGE,  XYZ, 3, 8, 7 )
      CALL TRIANGL3D( NCBLEU,  NCMAGE,  XYZ, 6, 7, 8 )
      END
