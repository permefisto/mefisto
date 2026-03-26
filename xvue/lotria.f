      SUBROUTINE LOTRIA( NP,     AMPLI,
     %                   SOMMET, NS1, NS2, NS3, NCFACE, NCARET )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE LE TRIANGLE DE SOMMETS (NS1, NS2, NS3) AVEC LA COULEUR
C -----  NCFACE ET LES ARETES SELON LA COULEUR NCARET
C
C ENTREE :
C --------
C NP     : CENTRE EN PIXELS DU LOGO
C AMPLI  : FACTEUR D'AMPLIFICATION DES COORDONNEES DU LOGO
C SOMMET : XY EN INTEGER*2 DES 7 SOMMETS DU LOGO
C NS1,NS2,NS3 : LES NUMEROS DANS SOMMET DES 3 SOMMETS DU TRIANGLE
C NCFACE : NUMERO DE LA COULEUR DU TRIANGLE (SI <0 PAS DE TRACE)
C NCARET : NUMERO DE LA COULEUR DES 3 ARETES DU TRIANGLE
C          (SI <0 PAS DE TRACE)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS  SEPTEMBRE 1994
C2345X7..............................................................012
      INTEGER   NP(2)
      REAL      AMPLI
      INTEGER*2 SOMMET(2,7), T(2,3)
      INTRINSIC INT2
C
C     TRACE EFFECTIF DU REMPLISSAGE DU TRIANGLE
      T(1,1) = INT2( NP(1) + AMPLI * SOMMET(1,NS1) )
      T(1,2) = INT2( NP(1) + AMPLI * SOMMET(1,NS2) )
      T(1,3) = INT2( NP(1) + AMPLI * SOMMET(1,NS3) )
C
      T(2,1) = INT2( NP(2) - AMPLI * SOMMET(2,NS1) )
      T(2,2) = INT2( NP(2) - AMPLI * SOMMET(2,NS2) )
      T(2,3) = INT2( NP(2) - AMPLI * SOMMET(2,NS3) )
C
      IF( NCFACE .GE. 0 ) THEN
         CALL XVCOULEUR( NCFACE )
         CALL XVFACE( 3, T )
      ENDIF
C
C     TRACE DES ARETES DU TRIANGLE
      IF( NCARET .GE. 0 ) THEN
         CALL XVEPAISSEUR( 3 )
         CALL XVCOULEUR( NCARET )
         NX1 = T(1,3)
         NY1 = T(2,3)
         DO 20 I=1,3
            NX2 = T(1,I)
            NY2 = T(2,I)
            CALL XVTRAIT( NX1, NY1, NX2, NY2 )
            NX1 = NX2
            NY1 = NY2
20       CONTINUE
      ENDIF
C
      RETURN
      END
