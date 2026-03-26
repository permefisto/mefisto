      SUBROUTINE QUSTIN( NBS1, NBS2, STCARR, COSO  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER LES 3 COORDONNEES DES SOMMETS INTERNES DU
C -----     QUADRANGLE A MAILLER PAR LA METHODE DE GORDON-HALL
C           EN L'ABSENCE DE TANGENTES
C
C ENTREES:
C --------
C NBS1   : NOMBRE DE POINTS SUR LES LIGNES "HORIZONTALES"
C NBS2   : NOMBRE DE POINTS SUR LES LIGNES "VERTICALES"
C STCARR : COORDONNEES DES NBS1 x NBS2 SOMMETS DU CARRE UNITE
C COSO   : COORDONNEES DES SOMMETS DES 4 COTES (BORD) DU MAILLAGE
C
C SORTIES:
C --------
C COSO   : COORDONNEES DES NBS1 x NBS2 SOMMETS DU MAILLAGE DU QUADRANGLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NERIQUE UPMC PARIS     SEPTEMBRE 1996
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO,EPSXYZ
      REAL              STCARR(2,NBS1,NBS2)
      REAL              COSO(3,NBS1,NBS2)
      REAL              F(3,4)
C
      DO 100 J = 2, NBS2-1
C
         DO 90 I = 2, NBS1-1
C
C           LES 2 COORDONNEES DU SOMMET I,J DU CARRE UNITE
            XM = STCARR(1,I,J)
            YM = STCARR(2,I,J)
C
C           RECHERCHE DU POINT PROJETE SUR LA LIGNE 1
C           -----------------------------------------
            N2 = 1
            DO 10 K=2,NBS1
               IF( XM .LE. STCARR(1,K,1)+EPSXYZ ) THEN
                  N2 = K
                  GO TO 11
               END IF
 10         CONTINUE
C           LES NUMEROS DES SOMMETS DE L'ARETE CONTENANT LE POINT PROJETE
 11         N1 = N2 - 1
            S1 = (    XM           - STCARR(1,N1,1) )
     S          / ( STCARR(1,N2,1) - STCARR(1,N1,1) )
            S0 = 1.0 - S1
C           ARETE P1 => C0
            DO 14 K=1,3
               F(K,1) = COSO(K,N1,1) * S0 + COSO(K,N2,1) * S1
 14         CONTINUE
C
C           RECHERCHE DU POINT PROJETE SUR LA LIGNE 2
C           -----------------------------------------
            N2 = 1
            DO 20 K=2,NBS2
               IF( YM .LE. STCARR(2,NBS1,K)+EPSXYZ ) THEN
                  N2 = K
                  GO TO 21
               END IF
 20         CONTINUE
C           LES NUMEROS DES SOMMETS VOISINS
 21         N1 = N2 - 1
            S1 = (   YM               - STCARR(2,NBS1,N1) )
     S          / ( STCARR(2,NBS1,N2) - STCARR(2,NBS1,N1) )
            S0 = 1.0 - S1
C           ARETE P1 => C0
            DO 24 K=1,3
               F(K,2) = COSO(K,NBS1,N1) * S0 + COSO(K,NBS1,N2) * S1
 24         CONTINUE
C
C           RECHERCHE DU POINT PROJETE SUR LA LIGNE 3
C           -----------------------------------------
            N2 = 1
            DO 31 K=2,NBS1
               IF( XM .LE. STCARR(1,K,NBS2)+EPSXYZ ) THEN
                  N2 = K
                  GO TO 32
               END IF
 31         CONTINUE
C           LES NUMEROS DES SOMMETS VOISINS
 32         N1 = N2 - 1
            S1 = (   XM               - STCARR(1,N1,NBS2) )
     S          / ( STCARR(1,N2,NBS2) - STCARR(1,N1,NBS2) )
            S0 = 1.0 - S1
C           ARETE P1 => C0
            DO 34 K=1,3
               F(K,3) = COSO(K,N1,NBS2) * S0 + COSO(K,N2,NBS2) * S1
 34         CONTINUE
C
C           RECHERCHE DU POINT PROJETE SUR LA LIGNE 4
C           -----------------------------------------
            N2 = 1
            DO 41 K=2,NBS2
               IF( YM .LE. STCARR(2,1,K)+EPSXYZ ) THEN
                  N2 = K
                  GO TO 42
               END IF
 41         CONTINUE
C           LES NUMEROS DES SOMMETS
 42         N1 = N2 - 1
            S1 = (    YM            - STCARR(2,1,N1) )
     S          / (  STCARR(2,1,N2) - STCARR(2,1,N1) )
            S0 = 1.0 - S1
C           ARETE P1 => C0
            DO 44 K=1,3
               F(K,4) = COSO(K,1,N1) * S0 + COSO(K,1,N2) * S1
 44         CONTINUE
C
C           LE POINT (I,J) DU MAILLAGE D'APRES LA FORMULE DE GORDON ET HALL
C           ---------------------------------------------------------------
            XM1 = 1.0 - XM
            YM1 = 1.0 - YM
            DO 50 K=1,3
               COSO(K,I,J) = YM1 * ( F(K,1) - XM1 * COSO(K,   1,   1) )
     S                     + XM  * ( F(K,2) - YM1 * COSO(K,NBS1,   1) )
     S                     + YM  * ( F(K,3) - XM  * COSO(K,NBS1,NBS2) )
     S                     + XM1 * ( F(K,4) - YM  * COSO(K,   1,NBS2) )
 50         CONTINUE
C
 90      CONTINUE
C
 100  CONTINUE
      END
