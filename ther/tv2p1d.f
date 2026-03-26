      SUBROUTINE TV2P1D( NBJEUX, JEU, X, NODL,   NODLIB,
     &                   NOOBSF, NUMISU, NUMASU, LTDESU,
     &                   UG,     EXPOSANT,
     &                   VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL du VECTEUR ELEMENTAIRE de
C -----    INTEGRALE t[P(X)] COEF TEMP ([P(X)]{Ue})**EXPOSANT dX sur l'EF
C          LAGRANGE TRIANGLE 2P1D DE DEGRE 1 AVEC
C          INTEGRATION NUMERIQUE AUX 3 SOMMETS DU TRIANGLE
C
C ENTREES:
C --------
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C X      : LES 2 COORDONNEES DES 3 SOMMETS DE L'ELEMENT FINI
C NODL   : NUMERO DE DEGRE DE LIBERTE GLOBAL DES 3 DL LOCAUX
C NODLIB : NODLIB(I) = NUMERO DU DEGRE DE LIBERTE S'IL EST LIBRE
C                     -INDICE DANS LA LISTE DES DL BLOQUES S'IL EST BLOQUE
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CONDUCTIVITE DES SURFACES
C
C UG     : VECTEUR SOLUTION A ELEVER A LA PUISSANCE EXPOSANT
C EXPOSANT:ENTIER EXPOSANT DE Ue
C
C SORTIE :
C --------
C VE     : VECTEUR ELEMENTAIRE du TERME
C          INTEGRALE t[P(X)] COEF TEMP ([P(X)]{Ue})**EXPOSANT dX sur l'EF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET St PIERRE DU PERRAY & LJLL UPMC  FEVRIER 2010
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
C
      DOUBLE PRECISION UG(*), VE(3)
      REAL             X(3,2)
      INTEGER          LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
      INTEGER          EXPOSANT, NODL(3), NODLIB(*)
      DOUBLE PRECISION DELTA, XYZPI(3), COEFTE, UEL
      DOUBLE PRECISION X21, Y21, X31, Y31, X32, Y32
C
      IF( LTDESU(LPCOET,JEU,NOOBSF) .GT. 0 ) THEN
C
C        CONTRIBUTION DE LA SURFACE AU COEFFICIENT DE LA TEMPERATURE
C        ===========================================================
         X21 = X(2,1) - X(1,1)
         X31 = X(3,1) - X(1,1)
         X32 = X(3,1) - X(2,1)
C
         Y21 = X(2,2) - X(1,2)
         Y31 = X(3,2) - X(1,2)
         Y32 = X(3,2) - X(2,2)
C
C        Det DF * POIDS
         DELTA = ABS( X21 * Y31 - X31 * Y21 ) * 6D0
C
         DO L=1,3
C
C           CALCUL DE LA SOLUTION AU POINT D'INTEGRATION L = SOMMET L
C           NUMERO DU DL DANS LE MAILLAGE
            NDL = NODL(L)
C           NUMERO DU DL DANS LA NUMEROTATION DES DEGRES DE LIBERTE NON FIXES
            NDLL = NODLIB( NDL )
            IF( NDLL .GT. 0 ) THEN
C              DL NON FIXE
               UEL = UG( NDLL )
            ELSE
C              ATTENTION: LES AUTRES DL FIXES SONT SUPPOSES NULS
               UEL = 0D0
            ENDIF
            TEMPEL = UEL
C
C           RECHERCHE DU COEFFICIENT DE LA TEMPERATURE AU POINT D'INTEGRATION L
            XYZPI(1) = X(L,1)
            XYZPI(2) = X(L,2)
            XYZPI(3) = 0D0
            CALL RECOET( 3, NOOBSF, 3, XYZPI,
     %                   LTDESU(LPCOET,JEU,NOOBSF), COEFTE )
C
C           COEF TEMPERATURE = COEFTE * POIDS * DELTA
            COEFTE = COEFTE * DELTA
C
C           CONTRIBUTION A L'INTEGRALE DE
C           (  t[Pol(bl] COEFTEMP ([Pol(bl]{ue})**EXPOSANT  )
            VE(L) = COEFTE * ( UEL**EXPOSANT )
C
         ENDDO
C
      ELSE
C
C        PAS DE CONTRIBUTION DU VOLUME AU COEFFICIENT DE LA TEMPERATURE
C        ==============================================================
         DO L=1,3
            VE(L) = 0D0
         ENDDO
C
      ENDIF
C
      RETURN
      END
