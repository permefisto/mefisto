      SUBROUTINE CL3VED ( NBVAR, ALPHA1, B01, ALPHA2, B02, ALPHA3, B03,
     %                    B )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    B < =  ALPHA1 * B01 + ALPHA2 * B02 + ALPHA3 * B03
C -----    B, B01, B02, B03 VECTEURS DE NBVAR VARIABLES
C          ALPHA1, ALPHA2, ALPHA3 REELS DOUBLE PRECISION
C
C ENTREES:
C --------
C NBVAR  : NOMBRE DE VARIABLES DES VECTEURS B, B01, B02
C ALPHA1 : 1-ER  REEL         DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
C B01    : 1-ER  VECTEUR REEL DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
C ALPHA2 : 2-EME REEL         DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
C B02    : 2-EME VECTEUR REEL DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
C ALPHA3 : 3-EME REEL         DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
C B03    : 3-EME VECTEUR REEL DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
C
C SORTIE :
C --------
C B      : VECTEUR RESULTAT DE LA COMBINAISON LINEAIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     FEVRIER 1999
C2345X7..............................................................012
      DOUBLE PRECISION  B01(NBVAR),B02(NBVAR),B03(NBVAR),B(NBVAR),
     %                  ALPHA1,ALPHA2,ALPHA3
C
C     TESTS SUIVANT LES VALEURS DE ALPHA1 , ALPHA2, ALPHA3
C     POUR ELIMINER LES MULTIPLICATIONS PAR 1.D0
C     ====================================================
      IF( ALPHA1 .EQ. 1.D0 ) THEN
C
C        ALPHA1 = 1.D0
         IF( ALPHA2 .EQ. 1.D0 ) THEN
C
C           ALPHA1 = 1.D0 ET ALPHA2 = 1.D0
            IF( ALPHA3 .EQ. 1.D0 ) THEN
C
C              ALPHA1 = ALPHA2 = ALPHA3 = 1.D0
               DO 10 I = 1,NBVAR
                  B(I) = B01(I) + B02(I) + B03(I)
 10            CONTINUE
C
            ELSE
C
C              ALPHA1 = ALPHA2 = 1.D0 ET ALPHA3<>1.D0
               DO 11 I = 1,NBVAR
                  B(I) = B01(I) + B02(I) + ALPHA3*B03(I)
 11            CONTINUE
            ENDIF
C
         ELSE
C
C           ALPHA1 = 1.D0 ET ALPHA2<>1.D0
            IF( ALPHA3 .EQ. 1.D0 ) THEN
C
C              ALPHA1 = 1.D0 ET ALPHA2<>1.D0 ET ALPHA3 = 1.D0
               DO 20 I = 1,NBVAR
                  B(I) = B01(I) + ALPHA2*B02(I) + B03(I)
 20            CONTINUE
C
            ELSE
C
C              ALPHA1 = 1.D0 ET ALPHA2<>1.D0 ET ALPHA3<>1.D0
               DO 22 I = 1,NBVAR
                  B(I) = B01(I) + ALPHA2*B02(I) + ALPHA3*B03(I)
 22            CONTINUE
            ENDIF
C
         ENDIF
C
      ELSE
C
C        ALPHA1<>1D0
         IF( ALPHA2 .EQ. 1.D0 ) THEN
C
C           ALPHA1<>1.D0 ET ALPHA2 = 1.D0
            IF( ALPHA3 .EQ. 1.D0 ) THEN
C
C              ALPHA1<>1.D0 ET ALPHA2 = 1.D0 ET ALPHA3 = 1.D0
               DO 30 I = 1,NBVAR
                  B(I) = ALPHA1*B01(I) + B02(I) + B03(I)
 30            CONTINUE
C
            ELSE
C
C              ALPHA1<>1.D0 ET ALPHA2 = 1.D0 ET ALPHA3<>1.D0
               DO 33 I = 1,NBVAR
                  B(I) = ALPHA1*B01(I) + B02(I) + ALPHA3*B03(I)
 33            CONTINUE
            ENDIF
C
         ELSE
C
C           ALPHA1<>1.D0 ET ALPHA2<>1.D0
            IF( ALPHA3 .EQ. 1.D0 ) THEN
C
C              ALPHA1<>1.D0 ET ALPHA2<>1.D0 ET ALPHA3 = 1.D0
               DO 40 I = 1,NBVAR
                  B(I) = ALPHA1*B01(I) + ALPHA2*B02(I) + B03(I)
 40            CONTINUE
C
            ELSE
C
C              ALPHA1<>1.D0 ET ALPHA2<>1.D0 ET ALPHA3<>1.D0
               DO 44 I = 1,NBVAR
                  B(I) = ALPHA1*B01(I) + ALPHA2*B02(I) + ALPHA3*B03(I)
 44            CONTINUE
C
            ENDIF
         ENDIF
C
      ENDIF
      RETURN
      END
